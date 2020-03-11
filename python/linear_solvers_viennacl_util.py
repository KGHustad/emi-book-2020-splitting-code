import ctypes
import os
from ctypes import _CFuncPtr # type: ignore
from ctypes import c_int, c_int32, c_int64, c_uint32, c_uint64, c_ulong, c_float, c_double, c_void_p, Structure

import numpy as np
import scipy.sparse


c_idx_t = c_int32
c_idx_t_array = np.ctypeslib.ndpointer(dtype=c_idx_t, ndim=1, flags="contiguous")
float64_array = np.ctypeslib.ndpointer(dtype=c_double, ndim=1, flags="contiguous")
c_idx_t_ptr = ctypes.POINTER(c_idx_t)
float64_ptr = ctypes.POINTER(c_double)




class C_Matrix_CSR(Structure):
    _fields_ = [
        ("num_rows", c_uint64),
        ("num_cols", c_uint64),
        ("num_nonzeros", c_uint64),
        ("row_ptr", c_idx_t_ptr),
        ("col", c_idx_t_ptr),
        ("val", float64_ptr),
    ]


def errcheck(return_code, func, params):
    return_code_desc = [
        'success',
    ]
    if return_code != 0:
        msg = "ERROR: {0}".format(init_return_code_desc[return_code])
        raise RuntimeError(msg)

def _setup_solver_viennacl_init_cpu(c_func: _CFuncPtr):
    """
    int solver_viennacl_init(
        void **viennacl_solver_ctx_ptr_ptr,
        struct matrix_csr A_csr,
        int use_amg
    )
    """
    c_func.restype = c_int
    c_func.argtypes = [
        c_void_p, # viennacl_solver_ctx_ptr_ptr
        C_Matrix_CSR, # A_csr
        c_int, # use_amg
    ]
    c_func.errcheck = errcheck

def _setup_solver_viennacl_solve(c_func: _CFuncPtr):
    """
    int solver_viennacl_solve(
        viennacl_solver_ctx_t *ctx,
        ScalarType *b,
        ScalarType *x,
        idx_t N
    )
    """
    c_func.restype = c_int
    c_func.argtypes = [
        c_ulong, # viennacl_solver_ctx_ptr
        float64_array, # b
        float64_array, # x
        c_idx_t, # N
    ]
    c_func.errcheck = errcheck

def _setup_solver_viennacl_free(c_func: _CFuncPtr):
    """
    int solver_viennacl_free(viennacl_solver_ctx_t *ctx)
    """
    c_func.restype = c_int
    c_func.argtypes = [
        c_ulong, # viennacl_solver_ctx_ptr
    ]
    c_func.errcheck = errcheck


def setup_cpu_lib(lib: ctypes.CDLL):
    _setup_solver_viennacl_init_cpu(lib.solver_viennacl_init)
    _setup_solver_viennacl_solve(lib.solver_viennacl_solve)
    _setup_solver_viennacl_free(lib.solver_viennacl_free)
