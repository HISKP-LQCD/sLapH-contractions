#!/bin/bash

# Sets up a debug and a release build outside of the source code with the paths
# that are needed on QBiG. The debug build will automatically set the `-g` flag
# for debugging, the release version will be performance and has `-O2 -DNDEBUG`
# set.

set -e
set -u
set -x

sourcedir=$(pwd)

for buildtype in release debug; do
    builddir=../sLapH-contractions-$buildtype
    rm -rf "$builddir"
    mkdir "$builddir"
    pushd "$builddir"

    cmake \
        "$sourcedir" \
        -DCMAKE_MODULE_PATH=/usr/share/cmake-3.0/Modules/ \
        -DCMAKE_BUILD_TYPE=$buildtype \
        -DLIME_INCLUDE_DIRS='/hadron/helmes/libraries/lime-1.3.2/include' \
        -DLIME_LIBRARIES='-L /hadron/helmes/libraries/lime-1.3.2/lib -llime' 

    make -j $(nproc)
    popd
done
