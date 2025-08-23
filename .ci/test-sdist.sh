#!/bin/bash
set -euxo pipefail
source $(compgen -G /opt-3/cpython-v3.*-apt-deb/bin/activate)
git clean -xfd
SUNDBASE=$(compgen -G /opt-3/sundials-6.*-release)
( export \
    CFLAGS="-isystem $SUNDBASE/include -isystem /usr/include/suitesparse" \
    LDFLAGS="-Wl,--disable-new-dtags -Wl,-rpath,$SUNDBASE/lib -L$SUNDBASE/lib" \
    CC=gcc CXX=g++ \
    ; python3 setup.py sdist \
    ; cd ./dist/ \
    ; pip install ${CI_REPO_NAME}-*.tar.gz )

( cd /; python3 -m pytest --pyargs pykinsol )

env \
    CFLAGS="-isystem $SUNDBASE/include -isystem /usr/include/suitesparse" \
    LDFLAGS="-Wl,--disable-new-dtags -Wl,-rpath,$SUNDBASE/lib -L$SUNDBASE/lib" \
    CC=gcc \
    CXX=g++ \
    python3 setup.py build_ext -i

pip install pytest-flakes pytest-cov matplotlib sphinx numpydoc sphinx-rtd-theme
./scripts/run_tests.sh --cov ${CI_REPO_NAME} --cov-report html
./scripts/coverage_badge.py htmlcov/ htmlcov/coverage.svg
./scripts/generate_docs.sh
