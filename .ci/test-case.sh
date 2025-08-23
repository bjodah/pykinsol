#!/bin/bash
set -euxo pipefail
show_help() {
    echo "Usage: ./$(basename $0) [opts] /path/to/sundials/prefix/e.g./usr/local"
    echo "--native"
    echo "--python <python-exe-path>"
}

NATIVE=0
TEST_ASAN=0
MAKE_TMP_DIR=0
RENDER=0
PYTHON=python3
while [ $# -gt 0 ]; do
    case $1 in
        --native)
            NATIVE=1
            shift
            ;;
        --python)
            shift
            PYTHON=$1
            shift
            ;;
        --asan)
            TEST_ASAN=1
            shift
            ;;
        --tmp)
            MAKE_TMP_DIR=1
            shift
            ;;
        --render)
            RENDER=1
            shift
            ;;
        --sundials-dir)
            shift
            SUNDBASE=$1
            shift
            ;;
        *)
            >&2 echo "Unrecognized parameter: $1"
            exit 1
    esac
done
if [ "$MAKE_TMP_DIR" = 1 ]; then
    REPO_TMP_DIR=$(mktemp -d)
    trap 'rm -rf -- "$REPO_TMP_DIR"' EXIT
    cp -ra . "$REPO_TMP_DIR/."
    cd "$REPO_TMP_DIR"
fi
if [ ! -e "$SUNDBASE/include/sundials/sundials_config.h" ]; then >&2 echo "No such directory: $SUNDBASE"; exit 1; fi
if [[ $SUNDBASE =~ *-extended || $SUNDBASE =~ *-single ]]; then
    >&2 echo "pykinsol currently only supports double precision"
    exit 1
fi
export CPATH=/usr/include/suitesparse  # include <klu.h>
export CXXFLAGS="${CXXFLAGS:-} -isystem $SUNDBASE/include"
export LDFLAGS="-Wl,--disable-new-dtags -Wl,-rpath,$SUNDBASE/lib -L$SUNDBASE/lib -lopenblas"
export LD_LIBRARY_PATH=$(compgen -G "/opt-2/llvm-*/lib")

if [ $TEST_ASAN -eq 1 ]; then
    export CC=clang
    export CXX=clang++
    LIBCXX_ASAN_ROOT=$(compgen -G "/opt-2/libcxx*-asan")
    LIBCXX_ASAN_INCLUDE=${LIBCXX_ASAN_ROOT}/include/c++/v1
    if [ ! -d "$LIBCXX_ASAN_INCLUDE" ]; then
        >&2 echo "Not a directory: $LIBCXX_ASAN_INCLUDE"
        exit 1
    fi
    export CXXFLAGS="$CXXFLAGS -fsanitize=address -stdlib++-isystem ${LIBCXX_ASAN_INCLUDE} -ferror-limit=5"
    export LDFLAGS="${LDFLAGS:-} -fsanitize=address -Wl,-rpath,${LIBCXX_ASAN_ROOT}/lib -L${LIBCXX_ASAN_ROOT}/lib -lc++ -lc++abi -stdlib=libc++"
    LLVM_ROOT=$(compgen -G "/opt-2/llvm-*")
    export LIBRARY_PATH="$LLVM_ROOT/lib:${LIBCXX_ASAN_ROOT}/lib:${LIBRARY_PATH:-}"
    #export PY_LD_PRELOAD="$(clang++ --print-file-name=libclang_rt.asan.so):$(clang++ --print-file-name=libstdc++.so)"
    #export PYTHON="env ASAN_OPTIONS=abort_on_error=1,detect_leaks=0 ${PYTHON:-python3}"
    export PYTHON_ENV="env ASAN_OPTIONS=abort_on_error=1,detect_leaks=0 LD_PRELOAD=$(clang++ --print-file-name=libclang_rt.asan.so):${LIBCXX_ASAN_ROOT}/lib/libc++.so.1.0:${LIBCXX_ASAN_ROOT}/lib/libc++abi.so.1.0:${LIBCXX_ASAN_ROOT}/lib/libunwind.so.1.0"  # Or this failure appears:
else
    export CC=gcc
    export CXX=g++
    #if $CXX --version | head -n 1 | grep -E 'g++\.*13\.[0-9]\.[0-9]$'; then exit 1; fi
    export CXXFLAGS="$CXXFLAGS -D_GLIBCXX_DEBUG -D_GLIBCXX_PEDANTIC"
    export CONTEXT="echo ''; echo ''; valgrind --error-exitcode=1"
    export PYTHON_ENV=""
fi


if [ -d ./build ]; then
    rm -r ./build
fi
if [ -d ./dist ]; then
    rm -r ./dist
fi

$PYTHON -m pip install build #--upgrade --upgrade-strategy=eager build setuptools==72.1.0 wheel
$PYTHON -m build . --sdist
$PYTHON -m pip uninstall -y pykinsol
cd dist/
CC=$CXX CFLAGS=$CXXFLAGS $PYTHON -m pip install *.tar.gz


# env \
#     LD_PRELOAD=${PY_LD_PRELOAD:-} \
   
$PYTHON_ENV $PYTHON -m pytest \
    -v \
    ${EXTRA_PYTEST_FLAGS:-} \
    --doctest-modules --pyargs pykinsol
cd -

if [[ $RENDER -eq 1 ]]; then
    $PYTHON -m doctest README.rst
    ( cd examples/; jupyter nbconvert \
                           --to=html \
                           --ExecutePreprocessor.enabled=True \
                           --ExecutePreprocessor.timeout=300 \
                           *.ipynb )
    ( cd examples/; ../scripts/render_index.sh *.html )
fi
