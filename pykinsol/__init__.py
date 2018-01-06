# -*- coding: utf-8 -*-
"""
pykinsol provides a Python interface to KINSOL from the SUNDIALS package suite.
"""
from __future__ import division, absolute_import

from ._release import __version__
from .core import solve


def get_include():
    from pkg_resources import resource_filename, Requirement
    return resource_filename(Requirement.parse(__name__),
                             '%s/include' % __name__)
