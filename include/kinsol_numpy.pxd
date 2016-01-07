# -*- coding: utf-8; mode: cython -*-

from cpython.object cimport PyObject
from libcpp.vector cimport vector

cdef extern from "kinsol_numpy.hpp" namespace "kinsol_numpy":
    cdef cppclass PyKinsol:
        const size_t nx
        const int mlower, mupper

        PyKinsol(PyObject*, PyObject*, size_t, int, int)
        PyObject* solve(PyObject*, double, double, long int, PyObject*, PyObject*, PyObject*) except +
