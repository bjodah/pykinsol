#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

from pykinsol import solve


def f_cb(x, fout):
    fout[0] = x[0] + (x[0] - x[1])**3/2 - 1
    fout[1] = x[1] + (x[1] - x[0])**3/2


def j_cb(x, Jout, fx):
    Jout[0, 0] = 1 + 3/2 * (x[0] - x[1])**2
    Jout[0, 1] = -3/2 * (x[0] - x[1])**2
    Jout[1, 0] = -3/2 * (x[1] - x[0])**2
    Jout[1, 1] = 1 + 3/2 * (x[1] - x[0])**2


def main(x=0., y=0., savetxt='None', verbose=False, cython=False, repeat=1):
    """
    Demonstrate how to solve a system of non-linear eqautions
    defined as SymPy expressions.
    """
    if cython:
        import pyximport
        pyximport.install()
        from bi_dimensional_cy import f, j
    else:
        f, j = f_cb, j_cb
    for i in range(repeat):
        result = solve(f, j, [x, y])
    assert result['success']
    if savetxt != 'None':
        np.savetxt(x, savetxt)
    else:
        if verbose:
            print(result)
        else:
            print(result['x'])


if __name__ == '__main__':
    try:
        import argh
        argh.dispatch_command(main, output_file=None)
    except ImportError:
        import sys
        if len(sys.argv) > 1:
            import warnings
            warnings.warn("Ignoring parameters, install argh by running "
                          "'python -m pip install --user argh'.")
        main()
