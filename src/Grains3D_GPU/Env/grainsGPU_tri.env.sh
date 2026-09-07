#!/bin/bash

# Login-node environment for this Grains3DGPU checkout.
# Uses the local nvhpc/25.1 module and targets H100 (sm_90).
if ! command -v module >/dev/null 2>&1; then
    source /etc/profile.d/modules.sh 2>/dev/null || true
    source /etc/profile.d/zz-60-lmod.sh 2>/dev/null || true
fi
module load nvhpc/25.1 2>/dev/null || true

_GCC12_BIN=/cvmfs/soft.computecanada.ca/gentoo/2023/x86-64-v3/usr/x86_64-pc-linux-gnu/gcc-bin/12
if [ -d "${_GCC12_BIN}" ]; then
    export PATH="${_GCC12_BIN}:${PATH}"
fi
unset _GCC12_BIN

export GRAINS_CPP_COMPILER=$(command -v g++ 2>/dev/null || echo g++)
export GRAINS_CPP_COMPILER_DIST="GNU"
export GRAINS_CPP_COMPILER_VERSION=$(g++ -dumpfullversion 2>/dev/null || echo 12.3.1)

export GRAINS_GPU_COMPILER=nvcc
export GRAINS_GPU_COMPILER_DIST="CUDA"
export GRAINS_GPU_COMPILER_VERSION="12.6"
export GRAINS_GPU_COMPILER_ROOT=${NVHPC}/Linux_x86_64/25.1/cuda/12.6
export GRAINS_GPU_COMPILER_INCDIR="${GRAINS_GPU_COMPILER_ROOT}/include"
export GRAINS_GPU_COMPILER_BINDIR="${GRAINS_GPU_COMPILER_ROOT}/bin"
export GRAINS_GPU_COMPILER_LIBDIR="${GRAINS_GPU_COMPILER_ROOT}/lib64"
export GRAINS_GPU_MATH_INCDIR=${NVHPC}/Linux_x86_64/25.1/math_libs/12.6/targets/x86_64-linux/include
export GRAINS_GPU_MATH_LIBDIR=${NVHPC}/Linux_x86_64/25.1/math_libs/12.6/targets/x86_64-linux/lib

export GRAINS_FULL_EXT=${GRAINS_CPP_COMPILER_DIST}-${GRAINS_CPP_COMPILER_VERSION}-${GRAINS_GPU_COMPILER_DIST}-${GRAINS_GPU_COMPILER_VERSION}
export GRAINS_HOME=/home/aliry95/Grains3DGPU
export GRAINS_ROOT=${GRAINS_HOME}/Grains
export GRAINS_INCDIR=${GRAINS_ROOT}/include
export GRAINS_OBJDIR=${GRAINS_ROOT}/obj${GRAINS_FULL_EXT}
export GRAINS_LIBDIR=${GRAINS_ROOT}/lib${GRAINS_FULL_EXT}

export GRAINS_XERCES_ROOT=${GRAINS_HOME}/XERCES-2.8.0
export GRAINS_XERCES_INCDIR="${GRAINS_XERCES_ROOT}/include"
export GRAINS_XERCES_BINDIR="${GRAINS_XERCES_ROOT}/bin"
export GRAINS_XERCES_LIBDIR="${GRAINS_XERCES_ROOT}/lib64-${GRAINS_CPP_COMPILER_DIST}-${GRAINS_CPP_COMPILER_VERSION}"

export GTEST_ROOT=/usr
export GTEST_INCLUDE_DIR=/usr/include
export GTEST_LIBRARY_DIR=/usr/lib64
export GRAINS_TEST_TIMEOUT=300
export GRAINS_TEST_PARALLEL_JOBS=8

echo -e '\033[31mGRAINS_HOME\033[0m =' $GRAINS_HOME
echo -e '\033[31mGRAINS_CPP_COMPILER\033[0m =' $GRAINS_CPP_COMPILER
echo -e '\033[31mGRAINS_CPP_COMPILER_DIST\033[0m =' $GRAINS_CPP_COMPILER_DIST
echo -e '\033[31mGRAINS_CPP_COMPILER_VERSION\033[0m =' $GRAINS_CPP_COMPILER_VERSION
echo -e '\033[31mGRAINS_GPU_COMPILER\033[0m =' $GRAINS_GPU_COMPILER
echo -e '\033[31mGRAINS_GPU_COMPILER_DIST\033[0m =' $GRAINS_GPU_COMPILER_DIST
echo -e '\033[31mGRAINS_GPU_COMPILER_VERSION\033[0m =' $GRAINS_GPU_COMPILER_VERSION
echo -e '\033[31mGRAINS_GPU_COMPILER_ROOT\033[0m =' $GRAINS_GPU_COMPILER_ROOT
echo -e '\033[31mGRAINS_FULL_EXT\033[0m =' $GRAINS_FULL_EXT
echo -e '\033[31mXERCES_ROOT\033[0m =' $GRAINS_XERCES_ROOT

export GRAINS_CPP_COMPILER_FLAGS="-m64 -O3 -fPIC -std=c++20 \
    -Wno-ctor-dtor-privacy \
    -Wall -Wextra -Wconversion -Wshadow -Wpedantic -Wwrite-strings \
    -fmax-errors=8 \
    -g"
export GRAINS_CPP_LINKER_FLAGS="${GRAINS_CPP_COMPILER_FLAGS} -shared"
export GRAINS_GPU_COMPILER="${GRAINS_GPU_COMPILER_BINDIR}/${GRAINS_GPU_COMPILER}"
export GRAINS_GPU_LINKER="${GRAINS_GPU_COMPILER_BINDIR}/${GRAINS_GPU_COMPILER}"
export GRAINS_NVCC_HOST_COMPILER=$(command -v g++ 2>/dev/null || echo g++)
export NVCC_PREPEND_FLAGS="-ccbin=${GRAINS_NVCC_HOST_COMPILER}"
export GRAINS_GPU_COMPILER_FLAGS_RELEASE="-x cu -m64 \
    -O3 -dlto -dc \
    -std=c++20 -arch=sm_90 \
    -cudart static -cudadevrt static \
    -use_fast_math -extra-device-vectorization -restrict \
    --extended-lambda --expt-relaxed-constexpr \
    -Xcompiler \"-rdynamic,-fPIC,-fopenmp\" \
    -g"
export GRAINS_GPU_LINKER_FLAGS_RELEASE="-O3 -dlto \
    -arch=sm_90 -lcudart \
    -use_fast_math -extra-device-vectorization -restrict \
    -L${GRAINS_GPU_MATH_LIBDIR} -lcudart -lcudadevrt -lcurand \
    -lgomp -lc \
    -g"

export GRAINS_GPU_COMPILER_FLAGS_DEBUG="-x cu -m64 \
    -O0 -G -dc \
    -std=c++20 -arch=sm_90 \
    -cudart static -cudadevrt static \
    --extended-lambda --expt-relaxed-constexpr \
    -Xcompiler \"-rdynamic,-fPIC,-fopenmp\" \
    -g -DDEBUG"
export GRAINS_GPU_LINKER_FLAGS_DEBUG="-O0 -G \
    -arch=sm_90 -lcudart \
    -L${GRAINS_GPU_MATH_LIBDIR} -lcudart -cudadevrt -lcurand \
    -lgomp -lc \
    -g"

export MODE=${MODE:-release}
if [ "$MODE" = "debug" ]; then
    export GRAINS_GPU_COMPILER_FLAGS="$GRAINS_GPU_COMPILER_FLAGS_DEBUG"
    export GRAINS_GPU_LINKER_FLAGS="$GRAINS_GPU_LINKER_FLAGS_DEBUG"
else
    export GRAINS_GPU_COMPILER_FLAGS="$GRAINS_GPU_COMPILER_FLAGS_RELEASE"
    export GRAINS_GPU_LINKER_FLAGS="$GRAINS_GPU_LINKER_FLAGS_RELEASE"
fi
export GRAINS_XERCES_FLAGS="-L${GRAINS_XERCES_LIBDIR} -lxerces-c -lxerces-depdom"
export GRAINS_Z_FLAGS="-lz"

export CMAKE_CXX_STANDARD=20
export CMAKE_BUILD_TYPE=Release
export CMAKE_CUDA_ARCHITECTURES=90
export CMAKE_CUDA_STANDARD=20
export CMAKE_PREFIX_PATH="${GRAINS_XERCES_ROOT}:${GRAINS_GPU_COMPILER_ROOT}:${CMAKE_PREFIX_PATH}"
export PKG_CONFIG_PATH="${GRAINS_XERCES_LIBDIR}/pkgconfig:${PKG_CONFIG_PATH}"

export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${GRAINS_XERCES_LIBDIR}
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${GRAINS_ROOT}/lib${GRAINS_FULL_EXT}
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${GRAINS_GPU_COMPILER_ROOT}/targets/x86_64-linux/lib
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${GRAINS_GPU_COMPILER_ROOT}/lib64
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${GRAINS_GPU_MATH_LIBDIR}

source $GRAINS_HOME/Env/grains_xerces.env.sh

if [ -d "${GRAINS_XERCES_ROOT}/lib64-GNU-12.3.1" ]; then
    export GRAINS_XERCES_LIBDIR="${GRAINS_XERCES_ROOT}/lib64-GNU-12.3.1"
elif [ -d "${GRAINS_XERCES_ROOT}/lib64-GNU-13.3.0" ]; then
    export GRAINS_XERCES_LIBDIR="${GRAINS_XERCES_ROOT}/lib64-GNU-13.3.0"
fi
export GRAINS_XERCES_FLAGS="-L${GRAINS_XERCES_LIBDIR} -lxerces-c -lxerces-depdom"
export PKG_CONFIG_PATH="${GRAINS_XERCES_LIBDIR}/pkgconfig:${PKG_CONFIG_PATH}"
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${GRAINS_XERCES_LIBDIR}