#!/bin/bash -l
#### Job submission script example to a CPU queue.
#### (Remember for Slurm #SBATCH, one # is active, two or
#### more ## is commented out)
####
#SBATCH --account=courses
#SBATCH --partition=courses
#SBATCH --mem-per-cpu=4000
#SBATCH --time=00:05:00
##SBATCH --output=prog.out #You can optionally name output, otherwise slurm.jobid

#SBATCH --nodes=2
#SBATCH --ntasks-per-node=2
module load gcc openmpi
mpicc -o SR1 MPI_SR_1.c
time srun SR1

#SBATCH --nodes=8
#SBATCH --ntasks-per-node=8
time srun --nodes=8 --ntasks-per-node=8 SR2

