Solve large scale sparse linear systems from an rf-scidac applications
using STRUMPACK.


# Build instructions for SUMMIT

## Modules:
```
module unload cmake
module swap xl gcc/11.2.0
module load cuda/11.5.2
module load essl
module load spectrum-mpi
module load netlib-lapack
module load netlib-scalapack
module load cmake
module unload darshan-runtime
```

## Install metis
```
wget http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/metis-5.1.0.tar.gz
tar -xvzf metis-5.1.0.tar.gz
cd metis-5.1.0
mkdir install
make config shared=1 cc=gcc prefix=`pwd`/install       # specify your C compiler as cc=..
make
make install
export METIS_DIR=`pwd`/install    # this will help STRUMPACK find metis
```

## Build SLATE for GPU acceleration of ScaLAPACK functionality

```
git clone --recursive https://bitbucket.org/icl/slate

cd slate

export CXX=mpiCC      # or your preferred C++ compiler
export FC=mpif90      # or your preferred Fortran compiler
rm -rf build
rm -rf install
mkdir build
mkdir install
cd build
export SLATE_DIR=`pwd`/install

cmake .. \
      -DCMAKE_BUILD_TYPE=Release \
      -Dblas=essl \
      -Dgpu_backend=cuda \
      -DBUILD_SHARED_LIBS=OFF \
      -DCMAKE_INSTALL_PREFIX=../install \
      -DBLAS_LIBRARIES="${OLCF_ESSL_ROOT}/lib64/libessl.so;${OLCF_NETLIB_LAPACK_ROOT}/lib64/libblas.so" \
      -DLAPACK_LIBRARIES="${OLCF_ESSL_ROOT}/lib64/libessl.so;${OLCF_NETLIB_LAPACK_ROOT}/lib64/liblapack.so" \
      -DSCALAPACK_LIBRARIES="${OLCF_NETLIB_SCALAPACK_ROOT}/lib/libscalapack.so"

make -j8
make install
```


## Build STRUMPACK
```
git clone git@github.com:pghysels/STRUMPACK.git
cd STRUMPACK
export blaspp_DIR=$SLATE_DIR
export lapackpp_DIR=$SLATE_DIR
export slate_DIR=$SLATE_DIR

cmake ../ \
      -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_INSTALL_PREFIX=../install \
      -DCMAKE_CXX_COMPILER=mpiCC \
      -DCMAKE_C_COMPILER=mpicc \
      -DCMAKE_Fortran_COMPILER=mpif90 \
      -DCMAKE_CUDA_COMPILER=nvcc \
      -DCMAKE_CUDA_ARCHITECTURES="70" \
      -DSTRUMPACK_USE_CUDA=ON \
      -DTPL_BLAS_LIBRARIES="${OLCF_ESSL_ROOT}/lib64/libessl.so;${OLCF_NETLIB_LAPACK_ROOT}/lib64/libblas.so" \
      -DTPL_LAPACK_LIBRARIES="${OLCF_ESSL_ROOT}/lib64/libessl.so;${OLCF_NETLIB_LAPACK_ROOT}/lib64/liblapack.so" \
      -DTPL_SCALAPACK_LIBRARIES="${OLCF_NETLIB_SCALAPACK_ROOT}/lib/libscalapack.so" \
      -DSTRUMPACK_COUNT_FLOPS=ON \
      -DTPL_ENABLE_BPACK=OFF \
      -DTPL_ENABLE_ZFP=OFF \
      -DTPL_ENABLE_SLATE=ON
make -j8
make install
```

Finally run the driver using (adjust STRUMPACK path in `build.sh`)
```
bash build.sh
```

See `run_summit.sh` for a batch script.


