#!/bin/bash -l
#SBATCH --gres=gpu:v100:1
#SBATCH --time=00:05:00
##SBATCH --output=BA.out
#SBATCH --account=courses
#SBATCH --partition=courses

module load gcc cuda 

nvcc -o BA Blocking_Async.cu

time srun BA
