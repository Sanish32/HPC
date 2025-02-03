#!/bin/bash -l
#SBATCH --time=00:05:00
#SBATCH --mem-per-cpu=4000
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
##SBATCH --output=vecAdd_CPU.out
#SBATCH --account courses
#SBATCH --partition courses

module load gcc

gcc -o vecAdd_CPU vecAdd_CPU.c

time srun vecAdd_CPU
