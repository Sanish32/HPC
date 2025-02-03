#!/bin/bash -l
#### Job submission script example to a CPU queue.
#### (Remember for Slurm #SBATCH, one # is active, two or
#### more ## is commented out)
#SBATCH --account=courses
#SBATCH --partition=courses
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem-per-cpu=4000
#SBATCH --time=00:05:00
##SBATCH --output=prog.out #You can optionally name output, otherwise slurm.jobid

module load gcc openmpi

mpicc -o SR4 MPI_SR_4.c

time srun SR4
