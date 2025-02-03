#!/bin/bash -l
#### Job submission script example.
#### (Remember for Slurm #SBATCH, one # is active, two or
#### more ## is commented out)
#SBATCH --account courses
#SBATCH --partition courses
#### Standard parameters
#SBATCH --time=00:05:00
#SBATCH --mem-per-cpu=4000
#### For a small MPI job:
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
##SBATCH --output=prog.out

module purge
module load gcc openmpi

srun ../build/quicksort-serial


