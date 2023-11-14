# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)

import numpy as np

from ._kinsol_numpy import KinsolSolver as _KinsolSolver


def _a(arr):
    return np.asarray(arr, dtype=np.float64)


def solve(f_cb, j_cb, x0, fnormtol=1e-6, scsteptol=1e-12, x_scale=None,
          f_scale=None, constraints=None, lband=-1, uband=-1, mxiter=200):
    """
    Solves a system of (non-linear) equations.

    Parameters
    ----------
    f_cb: callable
        Function with signature f(x, fout) which modifies fout *inplace*.
    j_cb: callable
        Function with signature j(x, jmat_out, fx) which modifies
        jmat_out *inplace*. fx contains current fout of f_cb(x).
    x0: array_like
        initial guess
    \\*\\*kwargs: dict

    Returns
    -------
    dict:
        'x': solution vector
        'success': bool
        'nfev': int
        'njev': int
    """
    x0 = np.asarray(x0, dtype=np.float64)
    if f_scale is None:
        f_scale = np.ones(x0.size)
    if x_scale is None:
        x_scale = np.ones(x0.size)
    if constraints is None:
        constraints = np.zeros(x0.size)
    x = np.array(x0, dtype=np.float64)
    solver = _KinsolSolver(f_cb, j_cb, x.size, lband, uband)
    return solver.solve(x, fnormtol, scsteptol, _a(x_scale), _a(f_scale), _a(constraints), mxiter)
