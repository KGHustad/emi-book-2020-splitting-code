from abc import ABC, abstractmethod

import numpy as np
import scipy.sparse.linalg


class LinearSolver(ABC):
    @abstractmethod
    def __init__(self, A: scipy.sparse.spmatrix):
        pass

    @abstractmethod
    def solve(self, b: np.ndarray):
        pass


class Scipy_Direct_Solver(LinearSolver):
    def __init__(self, A: scipy.sparse.spmatrix):
        self._A = A

    def solve(self, b: np.ndarray):
        return scipy.sparse.linalg.spsolve(self._A, b)


class Scipy_CG_Solver(LinearSolver):
    def __init__(self, A: scipy.sparse.spmatrix):
        self._A = A
        self._x_prev = None

    def solve(self, b: np.ndarray):
        x_prev = self._x_prev

        tol = 1e-5
        if x_prev is None:
            x, info = scipy.sparse.linalg.cg(self._A, b, tol=tol)
        else:
            x, info = scipy.sparse.linalg.cg(self._A, b, tol=tol, x0=x_prev)

        if info != 0:
            raise ValueError("Solution with CG did not succeed. Code", info)

        self._x_prev = x
        # return copy of x so that x_prev does not reference the same array as the one returned to the user
        return x.copy()


try:
    import scikits.umfpack
    scipy.sparse.linalg.use_solver(useUmfpack=True, assumeSortedIndices=True)
except ImportError:
    # umfpack is not installed, so we do not need to call scipy.sparse.linalg.use_solver
    pass
