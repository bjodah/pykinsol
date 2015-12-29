# -*- coding: utf-8 -*-
"""
pykinsol provides a Python interface to KINSOL from the SUNDIALS package suite.
"""
from __future__ import division, absolute_import

from ._kinsol_numpy import solve as _solve
from ._release import __version__


def solve(f_cb, j_cb, nf, nx, x0, **kwargs):
    """
    Solves a system of (non-linear) equations.

    Parameters
    ----------
    f_cb: callable
        Function with signature f(x, params, fout) which modifies fout *inplace*.
    j_cb: callable
        Function with signature j(x, params, jmat_out) which modifies
        jmat_out *inplace*.
    x0: array_like
        initial guess
    \*\*kwargs: dict

    Returns
    -------
    dict:
        'x': solution vector
        'success': bool
        'nfev':
        'njev':
    """
    return _solve(f_cb, j_cb, nf, nx, x0, **kwargs)
