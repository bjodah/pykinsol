#!/bin/bash -x
PKG_NAME=${1:-${DRONE_REPO##*/}}
rm -r /usr/local/lib/python*/dist-packages/${PKG_NAME}*  # pip uninstall is useless
( cd /tmp; if ! ${PYTHON:-python3} -c "import $PKG_NAME"; then >&2 echo "Couldn't uninstall"; exit 1; fi )
set -e

for p in "${@:2}"
do
    export CPATH=$p/include:$CPATH LIBRARY_PATH=$p/lib:$LIBRARY_PATH LD_LIBRARY_PATH=$p/lib:$LD_LIBRARY_PATH
done

git clean -xfd

${PYTHON:-python3} setup.py sdist
if [ $NO_PIP == 1 ]; then
    ${PYTHON:-python3} setup.py build_ext -i
    export PYTHONPATH=$(pwd)
else
    PKG_VERSION=$(unset DISTUTILS_DEBUG; ${PYTHON:-python3} setup.py --version)
    (cd dist/; ${PYTHON:-python3} -m pip install $PKG_NAME-$PKG_VERSION.tar.gz)
fi
export LD_PRELOAD=$PY_LD_PRELOAD:$LD_PRELOAD
(cd /; ${PYTHON:-python3} -m pytest --pyargs $PKG_NAME)
(cd /; ${PYTHON:-python3} -c "from pykinsol import get_include as gi; import os; assert 'kinsol_numpy.pxd' in os.listdir(gi())")

#export PYTHONPATH=$(pwd)
./scripts/run_tests.sh
(cd examples/; jupyter nbconvert --to=html --ExecutePreprocessor.enabled=True --ExecutePreprocessor.timeout=300 *.ipynb)
(cd examples/; ../scripts/render_index.sh *.html)
