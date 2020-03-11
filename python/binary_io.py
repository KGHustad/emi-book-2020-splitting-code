import numpy as np
import scipy.sparse


def load_vec(filename):
        with open(filename, 'r') as f:
            header = np.fromfile(f, dtype=np.uint64, count=1)
            N = header[0]
            v = np.fromfile(f, dtype=np.float64, count=N)
            return v

def load_csr(filename):
    with open(filename, 'rb') as f:
        header = np.fromfile(f, count = 3, dtype=np.uint64)
        N = header[0]
        NNZ = header[2]
        rowptr = np.fromfile(f, count=int(N+1), dtype=np.uint32)
        col = np.fromfile(f, count=NNZ, dtype=np.uint32)
        val = np.fromfile(f, count=NNZ, dtype=np.float64)

    return N, NNZ, rowptr, col, val

def load_csr_scipy(filename):
    N, NNZ, rowptr, col, val = load_csr(filename)
    return scipy.sparse.csr_matrix((val, col, rowptr), shape=(N,N))


def save_matrix_to_binary_csr(A, outfilename):
    with open(outfilename, 'wb') as f:
        header = np.zeros(3, dtype=np.uint64)
        header[0] = A.shape[0]
        header[1] = A.shape[1]
        header[2] = A.nnz

        header.tofile(f)

        A.indptr.astype(np.uint32).tofile(f)
        A.indices.astype(np.uint32).tofile(f)
        A.data.astype(np.float64).tofile(f)
