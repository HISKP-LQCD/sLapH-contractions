#!/bin/bash
# Copyright © 2017 Martin Ueding <dev@martin-ueding.de>
# Licensed under The MIT/Expat License

builddir=contraction-build
rm -rf "$builddir"
mkdir "$builddir"
cd "$builddir"

cmake \
    ../cntr.v0.1 \
    -DEIGEN3_INCLUDE_DIRS='/hiskp2/werner/libraries/eigen-3.2.7' \
    -DHDF5_INCLUDE_DIRS=/hiskp2/knippsch/hdf5-1.8.17/include \
    -DHDF5_LIBRARIES='-L/hiskp2/knippsch/hdf5-1.8.17/lib -lhdf5_cpp -lhdf5 -lsz -lz' \
    -DCMAKE_CXX_COMPILER=g++-4.7

make -j $(nproc)
