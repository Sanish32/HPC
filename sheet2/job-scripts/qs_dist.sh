#!/bin/bash -l
#SBATCH --account courses
#SBATCH --partition courses
#### Standard parameters
#SBATCH --time=00:05:00
#SBATCH --mem-per-cpu=4000
#### For a large MPI job:
#SBATCH --nodes=2         #Use two nodes
#SBATCH --ntasks-per-node=4  #Altogether eight tasks, evenly distributed over nodes
####End of large MPI job.
##SBATCH --output=prog.out

module purge
module load gcc openmpi cmake

srun ../build/quicksort-distributed


