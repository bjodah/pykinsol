# -*- coding: utf-8; mode: cython -*-

from cpython.object cimport PyObject
from libcpp.vector cimport vector
from numpy cimport PyArrayObject

cdef extern from "kinsol_numpy.hpp" namespace "kinsol_numpy":
    cdef cppclass PyKinsol:
        const size_t nx
        const int mlower, mupper

        PyKinsol(PyObject*, PyObject*, size_t, int, int)
        PyObject* solve(PyArrayObject*, double, double, long int, PyArrayObject*, PyArrayObject*, PyArrayObject*) except +
