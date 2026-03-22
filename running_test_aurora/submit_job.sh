#!/bin/bash -l

echo Jobid: $PBS_JOBID
echo Running on host `hostname`
echo Running on nodes `cat $PBS_NODEFILE`

NNODES=`wc -l < $PBS_NODEFILE`
NRANKS=12       # Number of MPI ranks per node
NDEPTH=1        # Number of hardware threads per rank, spacing between MPI ranks on a node
NTHREADS=1      # Number of OMP threads per rank, given to OMP_NUM_THREADS
NTOTRANKS=$((NNODES*NRANKS))

export INST_DIR=$PWD/../build_aurora/inst_dir/
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${INST_DIR}/lib
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${INST_DIR}/lib64

ZERORK_TEST=${INST_DIR}/bin/zerork_random_reactors_test_gpu.x

module list

export CHIP_JIT_FLAGS_OVERRIDE="-ze-opt-enable-auto-large-GRF-mode"
#export ZERORK_GPUS_PER_NODE=12

echo "NUM_OF_NODES=${NNODES}  TOTAL_NUM_RANKS=${NTOTRANKS}  RANKS_PER_NODE=${NRANKS} THREADS_PER_RANK=${NTHREADS}"

echo "Binary path is $ZERORK_TEST"

ldd $ZERORK_TEST

echo "------ JOB STARTING ------"
date
START_TIME=$(date +"%Y-%m-%d %H:%M:%S")
mpiexec --np ${NTOTRANKS} -ppn ${NRANKS} -d ${NDEPTH} -env OMP_NUM_THREADS=${NTHREADS}  gpu_tile_compact.sh $ZERORK_TEST
error_code=$?
date
END_TIME=$(date +"%Y-%m-%d %H:%M:%S")
echo "------ JOB ENDED ------"

# Convert dates to seconds since epoch
DATE1_SECONDS=$(date -d "$START_TIME" +"%s")
DATE2_SECONDS=$(date -d "$END_TIME" +"%s")

# Calculate the difference in seconds
DURATION_SECONDS=$((DATE2_SECONDS - DATE1_SECONDS))

# Calculate hours, minutes, and seconds from the duration
HOURS=$((DURATION_SECONDS / 3600))
MINUTES=$(( (DURATION_SECONDS % 3600) / 60 ))
SECONDS=$((DURATION_SECONDS % 60))

echo "Job time in seconds: $DURATION_SECONDS"
exit $?
echo "Duration: $HOURS hours, $MINUTES minutes, and $SECONDS seconds."

