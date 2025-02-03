#!/bin/bash -l
#SBATCH --account courses
#SBATCH --partition courses
#SBATCH --time=00:05:00
#SBATCH --mem-per-cpu=4000
#### Number of nodes, number of MPI processes is nodes x ntasks
#SBATCH --nodes=2
#### Number of MPI processes per node
#SBATCH --ntasks-per-node=2
##SBATCH --output=prog.out

module load gcc openmpi

mpicc -o MPIs_MPI_2 MPIs_MPI_2.c

time srun MPIs_MPI_2
