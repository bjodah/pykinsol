#!/bin/bash
LLAPACK=openblas CPLUS_INCLUDE_PATH=${PREFIX}/include ${PYTHON} -m pip install --no-deps --ignore-installed .
