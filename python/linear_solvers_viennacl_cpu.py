import ctypes
import ctypes.util
import os
from ctypes import _CFuncPtr # type: ignore
from ctypes import c_int, c_int32, c_int64, c_uint32, c_uint64, c_ulong, c_float, c_double, c_void_p, Structure

import numpy as np
import scipy.sparse

from lib_util import load_library
from linear_solvers import LinearSolver
from linear_solvers_viennacl_util import C_Matrix_CSR, setup_cpu_lib, c_idx_t, c_idx_t_ptr, float64_ptr


HERE = os.path.abspath(os.path.dirname(__file__))


_lib_cpu = load_library(os.path.abspath(os.path.join(HERE, "c/build/lib/liblinear_solver_viennacl_cpu")))
setup_cpu_lib(_lib_cpu)


class ViennaCL_CG_with_AMG_Solver_CPU(LinearSolver):
    _use_amg = True
    def __init__(self, A: scipy.sparse.spmatrix):
        if not isinstance(A, scipy.sparse.csr_matrix):
            A = A.tocsr()

        assert A.shape[0] == A.shape[1], "matrix must be square"
        self._N = A.shape[0]

        self._solver_viennacl_init(A, self._use_amg)


    def __del__(self):
        self._solver_viennacl_free()

    def solve(self, b: np.ndarray, x : np.ndarray = None):
        assert b.size == self._N
        if x is not None:
            assert isinstance(x, np.ndarray)
            assert x.size == b.size
        else:
            x = np.zeros_like(b)

        self._solver_viennacl_solve(b, x)
        return x

    def _solver_viennacl_init(self,
        A: scipy.sparse.csr_matrix,
        use_amg: bool
    ):
        solver_ctx_ptr = c_ulong(0)
        solver_ctx_ptr_ptr = ctypes.byref(solver_ctx_ptr)

        A_csr = C_Matrix_CSR(
            num_rows = A.shape[0],
            num_cols = A.shape[1],
            num_nonzeros = A.nnz,
            row_ptr = A.indptr.ctypes.data_as(c_idx_t_ptr),
            col = A.indices.ctypes.data_as(c_idx_t_ptr),
            val = A.data.ctypes.data_as(float64_ptr),
        )

        _lib_cpu.solver_viennacl_init(
            solver_ctx_ptr_ptr,
            A_csr,
            use_amg,
        )

        self._solver_ctx_ptr = solver_ctx_ptr

    def _solver_viennacl_solve(self, b, x):
        _lib_cpu.solver_viennacl_solve(self._solver_ctx_ptr, b, x, self._N)

    def _solver_viennacl_free(self):
        _lib_cpu.solver_viennacl_free(self._solver_ctx_ptr)


class ViennaCL_CG_Solver_CPU(ViennaCL_CG_with_AMG_Solver_CPU):
    _use_amg = False



if __name__ == '__main__':
    from binary_io import load_vec, load_csr_scipy

    synthetic_input = False
    if synthetic_input:
        N = 500
        A = scipy.sparse.spdiags(2.5*np.ones(N), 0, N, N)
        A = A.tocsr()

        b = np.full(N, 10)
    else:
        A = load_csr_scipy("A_csr.bin")
        b = load_vec("b.bin")

    print(A.data.shape)


    #s = ViennaCL_CG_Solver_CPU(A)
    s = ViennaCL_CG_with_AMG_Solver_CPU(A)

    x = s.solve(b)
    r = A@x - b
    print(np.linalg.norm(r) / np.linalg.norm(b))

    #x_ref = load_vec("../fdm-matlab/systems_extracellular_sym/8x8/x.bin")

    # from IPython import embed
    # embed()
    # exit()
