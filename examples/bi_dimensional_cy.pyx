import cython
cimport numpy as cnp

@cython.wraparound(False)
@cython.boundscheck(False)
def f(cnp.ndarray[cnp.float64_t, ndim=1, mode='c'] x,
      cnp.ndarray[cnp.float64_t, ndim=1, mode='c'] fout):
    fout[0] = x[0] + (x[0] - x[1])**3/2 - 1
    fout[1] = x[1] + (x[1] - x[0])**3/2


@cython.wraparound(False)
@cython.boundscheck(False)
def j(cnp.ndarray[cnp.float64_t, ndim=1, mode='c'] x,
      cnp.ndarray[cnp.float64_t, ndim=2, mode='fortran'] Jout,
      cnp.ndarray[cnp.float64_t, ndim=1, mode='c'] fx):
    Jout[0, 0] = 1 + 3./2 * (x[0] - x[1])**2
    Jout[0, 1] = -3./2 * (x[0] - x[1])**2
    Jout[1, 0] = -3./2 * (x[1] - x[0])**2
    Jout[1, 1] = 1 + 3./2 * (x[1] - x[0])**2
