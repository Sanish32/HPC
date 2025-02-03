#!/bin/bash -l
#### Job submission script example to a CPU queue using hybrid programming model
#SBATCH --account courses
#SBATCH --partition courses
#SBATCH --time=00:05:00
#SBATCH --mem-per-cpu=4000
#### Number of nodes, number of MPI processes is nodes x ntasks
#SBATCH --nodes=1
#### Number of MPI processes per node
#SBATCH --ntasks-per-node=4
#### Number of openMP threads, number of total CPUs is nodes x ntasks x cpus
#SBATCH --cpus-per-task=4
##SBATCH --output=prog.out

export OMP_PROC_BIND=TRUE
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
module load gcc openmpi

mpicc -fopenmp -o OMP_MPI_2 OMP_MPI_2.c

time srun OMP_MPI_2
