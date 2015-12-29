#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Tested with Sundials 2.6.2

import os
import sys
from setuptools import setup, Extension


pkg_name = 'pykinsol'

# Cythonize .pyx file if it exists (not in source distribution)
ext_modules = []

LLAPACK = os.environ.get('LLAPACK', 'lapack')

if len(sys.argv) > 1 and '--help' not in sys.argv[1:] and sys.argv[1] not in (
        '--help-commands', 'egg_info', 'clean', '--version'):
    import numpy as np
    USE_CYTHON = os.path.exists('pykinsol/_kinsol_numpy.pyx')  # not in sdist
    ext = '.pyx' if USE_CYTHON else '.cpp'
    ext_modules = [
        Extension('pykinsol._kinsol_numpy',
                  ['pykinsol/_kinsol_numpy'+ext],
                  language='c++', extra_compile_args=['-std=c++11'],
                  libraries=['sundials_kinsol', LLAPACK,
                             'sundials_nvecserial'],
                  include_dirs=[np.get_include(), './include'])
    ]
    if USE_CYTHON:
        from Cython.Build import cythonize
        ext_modules = cythonize(ext_modules, include_path=['./include'])

PYKINSOL_RELEASE_VERSION = os.environ.get('PYKINSOL_RELEASE_VERSION', '')

# http://conda.pydata.org/docs/build.html#environment-variables-set-during-the-build-process
CONDA_BUILD = os.environ.get('CONDA_BUILD', '0') == '1'
if CONDA_BUILD:
    try:
        PYKINSOL_RELEASE_VERSION = 'v' + open(
            '__conda_version__.txt', 'rt').readline().rstrip()
    except IOError:
        pass

release_py_path = os.path.join(pkg_name, '_release.py')

if (len(PYKINSOL_RELEASE_VERSION) > 1 and
   PYKINSOL_RELEASE_VERSION[0] == 'v'):
    TAGGED_RELEASE = True
    __version__ = PYKINSOL_RELEASE_VERSION[1:]
else:
    TAGGED_RELEASE = False
    # read __version__ attribute from _release.py:
    exec(open(release_py_path).read())

classifiers = [
    "Development Status :: 3 - Alpha",
    'License :: OSI Approved :: BSD License',
    'Operating System :: OS Independent',
    'Topic :: Scientific/Engineering',
    'Topic :: Scientific/Engineering :: Mathematics',
]

tests = [
    'pykinsol.tests',
]

descr = 'Python binding for kinsol from the sundials library.'
with open(pkg_name + '/__init__.py') as f:
    long_description = f.read().split('"""')[1]

setup_kwargs = dict(
    name=pkg_name,
    version=__version__,
    description=descr,
    long_description=long_description,
    classifiers=classifiers,
    author='Björn Dahlgren',
    author_email='bjodah@DELETEMEgmail.com',
    url='https://github.com/bjodah/' + pkg_name,
    license='BSD',
    packages=[pkg_name] + tests,
    ext_modules=ext_modules,
    install_requires=['numpy']
)

if __name__ == '__main__':
    import shutil
    try:
        if TAGGED_RELEASE:
            # Same commit should generate different sdist
            # depending on tagged version (set PYKINSOL_RELEASE_VERSION)
            # this will ensure source distributions contain the correct version
            shutil.move(release_py_path, release_py_path+'__temp__')
            open(release_py_path, 'wt').write(
                "__version__ = '{}'\n".format(__version__))
        setup(**setup_kwargs)
    finally:
        if TAGGED_RELEASE:
            shutil.move(release_py_path+'__temp__', release_py_path)
