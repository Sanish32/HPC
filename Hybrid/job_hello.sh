#!/bin/bash -l
#### (Remember for Slurm #SBATCH, one # is active, two or
#### more ## is commented out)
#SBATCH --account courses
#SBATCH --partition courses
#### Standard parameters
#SBATCH --time=00:05:00
#SBATCH --mem-per-cpu=4000
#### Number of nodes
#SBATCH --nodes=2
#### Number of MPI processes per node
#SBATCH --ntasks-per-node=2
#### Number of openMP threads per MPI process
#SBATCH --cpus-per-task=4
##SBATCH --output=prog.out

export OMP_PROC_BIND=TRUE
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
module purge   # unload all current modules
module load gcc openmpi

mpicc -fopenmp -o hello_class hello_class.c

time srun hello_class
