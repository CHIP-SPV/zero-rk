#!/bin/bash

# loads modules by default but a flag that will avoid the module loads.
# this is useful for CI
LOAD_MODULES=true

  for arg in "$@"; do
      if [[ "$arg" == "--no-module-loads" ]]; then
          LOAD_MODULES=false
      fi
  done

  if $LOAD_MODULES; then
      module restore
      module use /soft/modulefiles
      module load chipStar/llvm19/20251107-19/release
      module load cmake
fi

export CHIP_JIT_FLAGS_OVERRIDE="-ze-opt-enable-auto-large-GRF-mode"
export CHIP_LOGLEVEL=off

module list

export MPICH_CXX=hipcc
export MPICH_CC=hipcc

export CC=$(which mpicc)
export CXX=$(which mpicxx)

cmake -DCMAKE_EXE_LINKER_FLAGS="-lintlc" -DENABLE_GPU=ON -DENABLE_MAGMA=OFF -DCMAKE_CXX_COMPILER=hipcc -DCMAKE_C_COMPILER=hipcc -DZERORK_TESTS=OFF ../

make -j
ctest
make install
