#!/bin/bash

# Sets up a debug and a release build outside of the source code with the paths
# that are needed on QBiG. The debug build will automatically set the `-g` flag
# for debugging, the release version will be performance and has `-O2 -DNDEBUG`
# set.

set -e
set -u
set -x

module load Architecture/KNL
module load CMake

# Load Intel compiler
module load Intel

module load IntelMPI
module load HDF5
module load Boost
module load Eigen

sourcedir=$(pwd)

for buildtype in release; do
    builddir=../sLapH-contractions-$buildtype
    rm -rf "$builddir"
    mkdir "$builddir"
    pushd "$builddir"

    cmake \
        "$sourcedir" \
        -DEIGEN3_INCLUDE_DIR='${EBROOTEIGEN}/include/' \
        -DCMAKE_C_COMPILER=icc \
        -DCMAKE_CXX_COMPILER=icpc \
        -DCMAKE_CXX_FLAGS_RELEASE='-xMIC-AVX512 -fma -qopenmp' \
        -DLIME_INCLUDE_DIRS='/homec/hbn28/hbn287/.local/include' \
        -DLIME_LIBRARIES='-L/homec/hbn28/hbn287/.local/lib -llime' \
        -DCMAKE_BUILD_TYPE=$buildtype

    make -j $(nproc)
    popd
done
