import ctypes
import os
from ctypes import _CFuncPtr # type: ignore
from ctypes import c_int, c_int32, c_int64, c_uint, c_uint32, c_uint64, c_ulong, c_float, c_double, c_void_p, Structure

import numpy as np

from lib_util import load_library


float64_array = np.ctypeslib.ndpointer(dtype=c_double, ndim=1, flags="contiguous")
float64_array_2d = np.ctypeslib.ndpointer(dtype=c_double, ndim=2, flags="contiguous")


def _setup_init_state_values_2d_from_array(c_func: _CFuncPtr):
    "void init_state_values_2d_from_array(double* states, double *values, uint64_t n_nodes)"
    c_func.restype = None
    c_func.argtypes = [
        float64_array_2d, # states
        float64_array, # values
        c_uint64, # n_nodes
    ]

def _setup_init_parameter_values_2d_from_array(c_func: _CFuncPtr):
    "void init_parameter_values_2d_from_array(double* parameters, double *values, uint64_t n_nodes)"
    c_func.restype = None
    c_func.argtypes = [
        float64_array_2d, # parameters
        float64_array, # values
        c_uint64, # n_nodes
    ]

def _setup_forward_step(c_func: _CFuncPtr):
    c_func.restype = None
    c_func.argtypes = [
        float64_array_2d, # states
        c_double, # t
        c_double, # dt
        float64_array_2d, # parameters
        c_uint64, # n_nodes
    ]

def setup_lib(lib: ctypes.CDLL):
    _setup_init_state_values_2d_from_array(lib.init_state_values_2d_from_array)
    _setup_init_parameter_values_2d_from_array(lib.init_parameter_values_2d_from_array)
    _setup_forward_step(lib.FE)
    _setup_forward_step(lib.GRL1)


HERE = os.path.abspath(os.path.dirname(__file__))

_lib = load_library(
    os.path.abspath(os.path.join(HERE, "c/build/lib/libionic_model_grandi_cpu"))
)
setup_lib(_lib)

def FE(states: np.ndarray, t: float, dt: float, parameters: np.ndarray, n_nodes: np.uint64):
    _lib.FE(states, t, dt, parameters, n_nodes)

def GRL1(states: np.ndarray, t: float, dt: float, parameters: np.ndarray, n_nodes: np.uint64):
    _lib.GRL1(states, t, dt, parameters, n_nodes)

def init_state_values_2d_from_array(states: np.ndarray, values: np.ndarray, n_nodes: np.uint64):
    _lib.init_state_values_2d_from_array(states, values, n_nodes)

def init_parameter_values_2d_from_array(parameters: np.ndarray, values: np.ndarray, n_nodes: np.uint64):
    _lib.init_parameter_values_2d_from_array(parameters, values, n_nodes)
