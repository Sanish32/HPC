#!/bin/bash -l
#SBATCH --nodes=2  #Number of nodes, maximally 2
#SBATCH --ntasks-per-node=4 #Number of MPI ranks
#SBATCH --gres=gpu:v100:4 #Number of GPUs
#SBATCH -t 00:00:59
#SBATCH --account courses
#SBATCH --partition courses

module load gcc cuda openmpi cmake

nvcc -o HC hello_class.cu -lmpi

time srun HC
