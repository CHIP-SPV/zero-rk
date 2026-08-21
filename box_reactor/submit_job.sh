#!/bin/bash -xe
# #PBS -A CombustTurbine
# #PBS -N zrk_h2_box_lu1
# #PBS -l walltime=01:00:00
# #PBS -l filesystems=flare
# #PBS -k doe
# #PBS -l place=scatter
# #PBS -q debug
# #PBS -l select=1
# #PBS -m be

# echo Working directory is $PBS_O_WORKDIR
# cd $PBS_O_WORKDIR

echo Jobid: $PBS_JOBID
echo Running on host `hostname`
echo Running on nodes `cat $PBS_NODEFILE`

NNODES=`wc -l < $PBS_NODEFILE`
NRANKS=12      # Number of MPI ranks per node
NDEPTH=1        # Number of hardware threads per rank, spacing between MPI ranks on a node
NTHREADS=1      # Number of OMP threads per rank, given to OMP_NUM_THREADS
NTOTRANKS=$((NNODES*NRANKS))

export INST_DIR=$PWD/../build_aurora/inst_dir/
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${INST_DIR}/lib
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${INST_DIR}/lib64

ZERORK_TEST=${INST_DIR}/bin/zerork_box_reactor_test_gpu.x

module list

export CHIP_JIT_FLAGS_OVERRIDE="-ze-opt-enable-auto-large-GRF-mode"
export CHIP_LOGLEVEL=off
#export NEOReadDebugKeys=1
#export EnableDeviceUsmAllocationPool=1
export ZERORK_REACTOR_USE_LU=1

echo "NUM_OF_NODES=${NNODES}  TOTAL_NUM_RANKS=${NTOTRANKS}  RANKS_PER_NODE=${NRANKS} THREADS_PER_RANK=${NTHREADS}"
echo "Binary path is $ZERORK_TEST"

ldd $ZERORK_TEST

echo "------ JOB STARTING ------"
date
START_TIME=$(date +"%Y-%m-%d %H:%M:%S")
#module load thapi
#export THAPI_SYNC_DAEMON=fs
#module load pti-gpu
#mpiexec --np ${NTOTRANKS} -ppn ${NRANKS} -d ${NDEPTH} -env OMP_NUM_THREADS=${NTHREADS} /home/applenco/thapi_devel_clean/build/ici/bin/iprof -- gpu_tile_compact.sh $ZERORK_TEST
#mpiexec --np ${NTOTRANKS} -ppn ${NRANKS} -d ${NDEPTH} -env OMP_NUM_THREADS=${NTHREADS} unitrace --device-timings -- gpu_tile_compact.sh $ZERORK_TEST
export ZERORK_TEST_N_STEPS=750
mpiexec --np ${NTOTRANKS} -ppn ${NRANKS} -d ${NDEPTH} -env OMP_NUM_THREADS=${NTHREADS}  gpu_tile_compact.sh $ZERORK_TEST

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
echo "Duration: $HOURS hours, $MINUTES minutes, and $SECONDS seconds."
