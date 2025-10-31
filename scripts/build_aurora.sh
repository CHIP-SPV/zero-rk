#!/bin/bash

module use /soft/modulefiles
module load chipStar/llvm17/20251010-17-5f25c7/release
module load cmake

export CHIP_JIT_FLAGS_OVERRIDE="-ze-opt-enable-auto-large-GRF-mode"
export CHIP_LOGLEVEL=off

export ZERORK_REACTOR_USE_HIPBLAS=1
export ZERORK_REACTOR_USE_LU=1

module list

export MPICH_CXX=hipcc
export MPICH_CC=hipcc

export CC=$(which mpicc)
export CXX=$(which mpicxx)

cmake -DCMAKE_EXE_LINKER_FLAGS="-lintlc" -DENABLE_GPU=ON -DENABLE_MAGMA=OFF -DCMAKE_CXX_COMPILER=hipcc -DCMAKE_C_COMPILER=hipcc -DZERORK_TESTS=OFF ../

make -j 8
ctest
make install
