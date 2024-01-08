#!/bin/bash -l

#SBATCH -p course -t 15          # Specific course queue and max wallclock time
#SBATCH -D ./                    # Current working directory
#SBATCH --export=ALL

# Load Intel compiler module
module load compilers/intel/2019u5

echo "Node list: $SLURM_JOB_NODELIST"
echo "Number of nodes allocated: $SLURM_JOB_NUM_NODES"
echo "Number of threads or processes: $SLURM_NTASKS"
echo "Number of processes per node: $SLURM_TASKS_PER_NODE"
echo "Requested tasks per node: $SLURM_NTASKS_PER_NODE"
echo "Requested CPUs per task: $SLURM_CPUS_PER_TASK"
echo "Scheduling priority: $SLURM_PRIO_PROCESS"

SRC=$1
EXE=${SRC%%.c}.exe
rm -f ${EXE}

echo compiling $SRC to $EXE
icc -qopenmp -O2 -std=c99 $SRC coordReader.c ompcInsertion.c ompfInsertion.c ompnAddition.c -o $EXE

# Execute gomp-only with different coord files



./gnu-omp_ne.sh; ./gomp-only 9_coords.coord coStudent1_9.dat foStudent1_9.dat noStudent1_9.dat;
./gnu-omp_ne.sh; ./gomp-only 16_coords.coord coStudent1_16.dat foStudent1_16.dat noStudent1_16.dat;
./gnu-omp_ne.sh; ./gomp-only 512_coords.coord coStudent1_512.dat foStudent1_512.dat noStudent1_512.dat;

# Check if executable is present and run it
if test -x $EXE; then
    export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
    echo using ${OMP_NUM_THREADS} OpenMP threads
    for i in {1..10}; do ./${EXE}; done
else
    echo $SRC did not build to $EXE
fi
