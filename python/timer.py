import logging
import math
import time
from contextlib import ContextDecorator
from dataclasses import dataclass, field
from typing import Any, Callable, Dict, List, Optional, Union

logger = logging.getLogger(__name__)


class TimerError(Exception):
    """A custom exception used to report errors in use of Timer class"""


@dataclass
class Timer(ContextDecorator):
    """Time your code using a class, context manager, or decorator"""

    total_time: float = 0.0

    _start_time: Optional[float] = field(default=None, init=False, repr=False)
    name: Optional[str] = None
    text: str = "Elapsed time: {:0.4f} seconds"
    logger: Optional[Callable[[str], None]] = logger.info
    last: float = field(default=math.nan, init=False, repr=False)

    def start(self) -> None:
        """Start a new timer"""
        if self._start_time is not None:
            raise TimerError("Timer is running. Use .stop() to stop it")

        self._start_time = time.perf_counter()

    def stop(self) -> float:
        """Stop the timer, and report the elapsed time"""
        if self._start_time is None:
            raise TimerError("Timer is not running. Use .start() to start it")

        # Calculate elapsed time
        self.last = time.perf_counter() - self._start_time
        self._start_time = None

        # Report elapsed time
        if self.logger:
            self.logger(self.text.format(self.last))
        self.total_time += self.last

        return self.last

    def __enter__(self) -> "Timer":
        """Start a new timer as a context manager"""
        self.start()
        return self

    def __exit__(self, *exc_info: Any) -> None:
        """Stop the context manager timer"""
        self.stop()


class TimerCollection:
    def __init__(self, name: Optional[str] = None):
        self.name = name
        self.timers: Dict[str, Timer] = dict()
        self.lap_times: Dict[str, List[float]] = dict()
        self._current_timer = None

    def __repr__(self):
        s = f"{__class__.__name__}("
        args = ""
        if self.name is not None:
            args += self.name
        for name in self.timers.keys():
            args += f", {name}"
        s += args.strip(",").strip(" ") + ")"
        return s

    def default_name(self):
        try:
            name = next(iter(self.timers.keys()))
        except StopIteration:
            name = "0"
        finally:
            return name

    def start(self, name: Optional[Union[str, List[str]]] = None):
        if name is None:
            name = self.default_name()

        if isinstance(name, str):
            name = [name]
        for n in name:
            if n not in self.timers:
                self.timers[n] = Timer(name=n)
            if n not in self.lap_times:
                self.lap_times[n] = []
            self.timers[n].start()

    def stop(self, name: Optional[Union[str, List[str]]] = None):
        if name is None:
            name = self.default_name()

        if isinstance(name, str):
            name = [name]
        for n in name:
            if n not in self.timers:
                raise TimerError(
                    f"Timer {name} is not running. Use .start() to start it"
                )
            if n not in self.lap_times:
                raise TimerError(
                    f"Timer {name} is not running. Use .start() to start it"
                )
            self.lap_times[n].append(self.timers[n].stop())

    @property
    def total_time(self):
        return sum(self.total_times.values())

    @property
    def percentage_time(self):
        tot = self.total_time
        return {k: t.total_time / tot for k, t in self.timers.items()}

    @property
    def total_times(self):
        return {k: t.total_time for k, t in self.timers.items()}

    def report(self, name: Union[List[str], str] = None) -> str:

        return format_output("Total times", self.total_times)


def format_output(
    header: str,
    d: Dict[str, float],
    width: int = 50,
    precision: int = 5,
    key_value_ratio: float = 0.8,
):
    msg = "\n"
    msg += f"{header}:".center(width + 2, " ") + "\n"
    msg += "-" * (width + 2)

    key_width = int(width * key_value_ratio)
    value_width = width - key_width
    format_str = f"\n{{k:{key_width}}}: {{t:>{value_width}.{precision}f}} s"

    for k, t in d.items():
        msg += format_str.format(k=k, t=t)

    return msg
