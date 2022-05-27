#!/usr/bin/env bash

# export STRUMPACK_DIR=/global/homes/p/pghysels/cori/STRUMPACK/install
export STRUMPACK_DIR=${HOME}/STRUMPACK/install

rm -rf build
rm -rf install
mkdir build
mkdir install
cd build

cmake ../ \
      -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_CXX_COMPILER=mpiCC \
      -DCMAKE_C_COMPILER=mpicc \
      -DCMAKE_Fortran_COMPILER=mpif90 \
      -DCMAKE_CUDA_COMPILER=nvcc

make
