# -*- coding: utf-8 -*-

from pykinsol import solve


def f_cb(x, fout):
    # expression from example in SciPy documentation:
    # docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.root.html
    fout[0] = x[0] + (x[0] - x[1])**3/2 - 1
    fout[1] = x[1] + (x[1] - x[0])**3/2


def j_cb(x, Jout, fx):
    Jout[0, 0] = 1 + 3/2 * (x[0] - x[1])**2
    Jout[0, 1] = -3/2 * (x[0] - x[1])**2
    Jout[1, 0] = -3/2 * (x[1] - x[0])**2
    Jout[1, 1] = 1 + 3/2 * (x[1] - x[0])**2


def test_solve():
    result = solve(f_cb, j_cb, [0, 0])
    assert result['success']
    assert result['message'] == 'KIN_SUCCESS'
    assert result['status'] == 0
    assert result['nit'] > 1
    assert result['nfev'] > 1
    assert result['njev'] > 0
    assert result['time_cpu'] > 1e-15
    assert abs(result['x'][0] - 0.8411639) < 2e-7
    assert abs(result['x'][1] - 0.1588361) < 2e-7


def test_solve__failure():
    mxiter = 6
    result = solve(f_cb, j_cb, [.5, .5], mxiter=mxiter)
    assert not result['success']
    assert result['message'] == 'KIN_MAXITER_REACHED'
    assert result['nit'] == mxiter
