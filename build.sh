#!/usr/bin/env bash

export STRUMPACK_DIR=/global/homes/p/pghysels/cori/STRUMPACK/install

rm -rf build
rm -rf install
mkdir build
mkdir install
cd build

cmake ../ \
      -DCMAKE_BUILD_TYPE=Debug \
      -DCMAKE_CXX_COMPILER=CC \
      -DCMAKE_C_COMPILER=cc \
      -DCMAKE_Fortran_COMPILER=ftn

make
