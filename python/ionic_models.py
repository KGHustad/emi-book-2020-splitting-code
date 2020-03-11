import os
from abc import ABC, abstractmethod

import numpy as np

import grandi_numpy
import util
try:
    import grandi_c
    has_grandi_c = True
except ImportError:
    has_grandi_c = False

grandi_num_states = 39
grandi_num_parameters = 112


class IonicModel(ABC):
    @abstractmethod
    def initialise(self, num_membrane_points, t_start=0, states=None, parameters=None):
        pass

    @abstractmethod
    def forward_n_steps(self, num_steps, V=None, scheme=None, dt=None):
        "Solves N time steps for the cell model, using the values in V as input to the first step, and later overwriting V with the output from the last step."
        pass

    @abstractmethod
    def state_index(self, state_name):
        pass

    @abstractmethod
    def parameter_index(self, parameter_name):
        pass

    @abstractmethod
    def init_parameter_values_single(self, **values):
        pass

    @abstractmethod
    def init_parameter_values_2d_from_array(self, parameters_single, num_membrane_points):
        pass

    @property
    @abstractmethod
    def t(self):
        pass

    @property
    @abstractmethod
    def V(self):
        pass

    @V.setter
    @abstractmethod
    def V(self, V):
        pass




class Grandi_NumPy(IonicModel):
    _valid_float_t = (np.float32, np.float64)
    _schemes = ["FE", "GRL1"]

    def __init__(self, float_t=np.float64):
        if float_t is not None and float_t not in self._valid_float_t:
            assert float_t
        self._float_t = float_t
        self._initialised = False
        self._states = None
        self._t = None
        self._V_idx = grandi_numpy.state_indices("V_m")

    def initialise(self, num_membrane_points, t_start=0, states=None, parameters=None):
        if states is None:
            # initialise states
            states_single = self.init_state_values_single()
            states = self.init_state_values_2d_from_array(states_single, num_membrane_points)
        else:
            assert states.shape == (grandi_num_states, num_membrane_points)

        if parameters is None:
            # initialise parameters
            parameters_single = self.init_parameter_values_single()
            parameters = self.init_parameter_values_2d_from_array(parameters_single, num_membrane_points)
        else:
            assert parameters.shape == (grandi_num_parameters, num_membrane_points)

        self._num_membrane_points = num_membrane_points
        self._states = states
        self._t = float(t_start)
        self._parameters = parameters

        self._initialised = True

    def _forward_step(self, scheme, dt):
        if scheme == "FE":
            grandi_numpy.FE(self._states, self._t, dt, self._parameters)
        elif scheme == "GRL1":
            grandi_numpy.GRL1(self._states, self._t, dt, self._parameters)

    def forward_n_steps(self, num_steps, V=None, scheme=None, dt=None):
        if dt is None:
            dt = 1e-3

        if scheme is None:
            scheme = self._schemes[0]
        assert scheme in self._schemes

        if V is not None:
            self.V = V
        for i in range(num_steps):
            self._forward_step(scheme, dt)
            self._t += dt
        if V is not None:
            V[:] = self.V

    @property
    def t(self):
        return self._t

    @property
    def V(self):
        V_idx = grandi_numpy.state_indices("V_m")
        return self._states[V_idx].copy()

    @V.setter
    def V(self, V):
        self._states[self._V_idx] = V[:]

    @property
    def num_membrane_points(self):
        return self._num_membrane_points

    @property
    def float_t(self):
        return self._float_t

    def state_index(self, state_name):
        return grandi_numpy.state_indices(state_name)

    def parameter_index(self, param_name):
        return grandi_numpy.parameter_indices(param_name)

    def init_state_values_single(self, **values):
        return grandi_numpy.init_state_values_single(*values)

    def init_state_values_2d_from_array(self, states_single, num_membrane_points):
        assert states_single.size == grandi_num_states
        states_single = states_single.reshape(grandi_num_states, 1)
        states = np.zeros((grandi_num_states, num_membrane_points), dtype=states_single.dtype)
        states += states_single
        return states

    def init_parameter_values_single(self, **values):
        return grandi_numpy.init_parameter_values_single(*values)

    def init_parameter_values_2d_from_array(self, parameters_single, num_membrane_points):
        assert parameters_single.size == grandi_num_parameters
        parameters_single = parameters_single.reshape(grandi_num_parameters, 1)
        parameters = np.zeros((grandi_num_parameters, num_membrane_points), dtype=parameters_single.dtype)
        parameters += parameters_single
        return parameters


class Grandi_Numba(Grandi_NumPy):
    def _forward_step(self, scheme, dt):
        if scheme == "FE":
            grandi_numba.FE(self._states, self._t, dt, self._parameters)
        elif scheme == "GRL1":
            grandi_numba.GRL1(self._states, self._t, dt, self._parameters)

if has_grandi_c:
    class Grandi_C(Grandi_NumPy):
        _valid_float_t = (np.float64,)
        _schemes = ["FE", "GRL1"]

        def _forward_step(self, scheme, dt):
            n_nodes = np.uint64(self.num_membrane_points)
            if scheme == "FE":
                grandi_c.FE(self._states, self._t, dt, self._parameters, n_nodes)
            elif scheme == "GRL1":
                grandi_c.GRL1(self._states, self._t, dt, self._parameters, n_nodes)

        def init_state_values_2d_from_array(self, states_single, num_membrane_points):
            assert states_single.size == grandi_num_states, f"got {states_single.size} states, expected {grandi_num_states}"
            # np.empty allocates memory without initialising it, which we need for performing first touch in C
            states = np.empty((grandi_num_states, num_membrane_points), dtype=states_single.dtype)
            n_nodes = np.uint64(num_membrane_points)
            grandi_c.init_state_values_2d_from_array(states, states_single, n_nodes)
            return states

        def init_parameter_values_2d_from_array(self, parameters_single, num_membrane_points):
            assert parameters_single.size == grandi_num_parameters, f"got {parameters_single.size} parameters, expected {grandi_num_parameters}"
            # np.empty allocates memory without initialising it, which we need for performing first touch in C
            parameters = np.empty((grandi_num_parameters, num_membrane_points), dtype=parameters_single.dtype)
            n_nodes = np.uint64(num_membrane_points)
            grandi_c.init_parameter_values_2d_from_array(parameters, parameters_single, n_nodes)
            return parameters

if __name__ == "__main__":
    import time

    N = 100_000
    num_steps = 10
    scheme = "FE"
    # scheme = "GRL1"

    grandi_implementations = [
        (Grandi_NumPy, "Grandi NumPy"),
    ]

    if has_grandi_c:
        pass
        grandi_implementations.append((Grandi_C, "Grandi C"))

    for G, desc in grandi_implementations:
        g = G()
        g.initialise(N)

        pre = time.perf_counter()
        g.forward_n_steps(num_steps, scheme=scheme)
        post = time.perf_counter()
        time_elapsed = post - pre
        cell_steps = N * num_steps
        print(
            "Cell steps per second ({1:s}): {0:g} (total time: {2:g} s)".format(
                cell_steps / time_elapsed, desc, time_elapsed
            )
        )
        V = g.V
        print(V.min(), V.mean(), V.max())
        print()
