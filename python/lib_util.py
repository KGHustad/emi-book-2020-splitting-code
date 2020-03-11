import ctypes
import ctypes.util
import os

import numpy as np


HERE = os.path.abspath(os.path.dirname(__file__))

class LibraryNotFoundException(Exception):
    pass

def load_library(filename: str) -> ctypes.CDLL:
    path = None
    for suffix in (".so", ".a", ".dylib"):
        if os.path.exists(filename + suffix):
            path = filename + suffix
            break
    #path = ctypes.util.find_library(filename)

    if path is None:
        raise LibraryNotFoundException(
            "Library file {0} does not exist. Try running rebuild_libs in the python/c directory to build the library files".format(
                filename
            )
        )

    lib = np.ctypeslib.load_library(path, HERE)
    return lib
