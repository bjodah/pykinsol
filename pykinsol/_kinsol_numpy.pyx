# -*- coding: utf-8; mode: cython -*-
# distutils: language = c++
# cython: language_level=3str

from cpython.object cimport PyObject
from libcpp cimport bool
cimport numpy as cnp
import numpy as np

from kinsol_numpy cimport PyKinsol

cnp.import_array()  # Numpy C-API initialization  (see include/kinsol_numpy.hpp)


cdef class KinsolSolver:

    cdef PyKinsol *thisptr

    def __cinit__(self, object f, object j, size_t nx, int ml=-1, int mu=-1):
        self.thisptr = new PyKinsol(<PyObject *>f, <PyObject *>j, nx, ml, mu)

    def __dealloc__(self):
        del self.thisptr

    def solve(self, cnp.ndarray[cnp.float64_t, ndim=1] x0, double fnormtol, double scsteptol,
              cnp.ndarray[cnp.float64_t, ndim=1] x_scale,
              cnp.ndarray[cnp.float64_t, ndim=1] f_scale,
              cnp.ndarray[cnp.float64_t, ndim=1] constraints,
              long int mxiter=200):
        if (x_scale.size != x0.size or
            f_scale.size != x0.size or
            constraints.size != x0.size):
            raise ValueError("Incompatible lengths")
        return <object>self.thisptr.solve(
            <cnp.PyArrayObject*>x0, fnormtol, scsteptol, mxiter,
            <cnp.PyArrayObject*>x_scale, <cnp.PyArrayObject*>f_scale, <cnp.PyArrayObject*>constraints)
