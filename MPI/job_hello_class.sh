#!/bin/bash -l
#SBATCH --account=courses
#SBATCH --partition=courses
#SBATCH --nodes=2
###What is the difference if you use here --ntasks=8?
#SBATCH --ntasks-per-node=8
#SBATCH --mem-per-cpu=4000
#SBATCH --time=00:05:00
##SBATCH --output=prog.out #You can optionally name output, otherwise slurm.jobid


module load gcc openmpi

mpicc -o hello_class hello_class.c

time srun hello_class
