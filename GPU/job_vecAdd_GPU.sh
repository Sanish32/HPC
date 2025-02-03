#!/bin/bash -l
#SBATCH --gres=gpu:v100:1
#SBATCH --time=00:05:00
##SBATCH --output=vecAdd_GPU.out
#SBATCH --account=courses
#SBATCH --partition=courses

module load gcc cuda 

nvcc -o vecAdd_GPU vecAdd_GPU.cu

time srun vecAdd_GPU
