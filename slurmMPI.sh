#!/bin/bash -l

#SBATCH -N 16                # Request 16 nodes
#SBATCH -n 32                # 32 MPI processes
#SBATCH -c 16                # 16 OpenMP threads per process
#SBATCH -p course            # Partition
#SBATCH --export=ALL
#SBATCH -D ./                # Working directory

# Load necessary modules
module load mpi/intel-mpi/2019u5/bin
module load compilers/intel/2019u5

# Set the number of OpenMP threads
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

# Compile the program
mpiicc -qopenmp coordReader.c main-mpi.c ompcInsertion.c ompfInsertion.c ompnAddition.c -std=c99 -lm

# Run the program
mpirun -np $procs ./a.out
