cimport cython
import numpy as np
cimport numpy as np

matrix = np.ndarray[np.float32_t, ndim=2]
ctypedef np.ndarray[np.float32_t, ndim=2] matrix_t

@cython.boundscheck(False)
@cython.wraparound(False)
cdef matrix_t times_(matrix_t a, matrix_t b):
    cdef matrix_t c
    cdef int n, p, m, i, j, k
    cdef np.float32_t s
    n, p, m = a.shape[0], a.shape[1], b.shape[1]
    c = np.zeros((n, m), dtype=np.float32)
    for i in xrange(n):
        for j in xrange(m):
            s = 0
            for k in xrange(p):
                s += a[i, k] * b[k, j]
            c[i, j] = s
    return c


def times(a, b):
    return times_(a, b)