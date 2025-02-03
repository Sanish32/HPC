## Exercise 2 - Distributed Quicksort
# Introduction
In this exercise, you will implement a stable version of parallel distributed quicksort across MPI processes to sort data common to all processes. 
In exercise 4 you will extend your implementation to work on multiple GPUs.
The exercise consists of two subtasks:

  1. Implement the function `quicksort` in `src/quicksort.cpp`. Your task here is to implement a  single threaded stable version of quicksort.

  2. Extend your implementation to be able to use multiple MPI processes. on a single node in `src/quicksort-distributed.cpp`. You need to implement `quicksort_distributed`  (you can are encouraged to re-use your solution from subtask 1. ) Now multiple processes hold a copy of the array and each process should end with the same sorted array.. Remember to allocate multiple MPI tasks  (-n should be larger than 1).

For each solution you have to implement there will be more detailed instructions in the corresponding file about what you need to do and what is expected from your implementations.

## Returnables
The solutions should be returned in a single zip file named <your student number>.zip, for example 12345.zip. The archive should contain at least `src/quicksort.cpp`, `src/quicksort_distributed.cpp`. You create the archive with the command zip <your student number>.zip -r src/ (note the -r flag: the archive must contain the src directory).

## Getting started
Run the following commands to get started:

```Bash
module load gcc cmake openmpi
cd pps-2024-exercises/exercise_sheets/sheet2
mkdir build && cd build
cmake .. && make -j
cd .. && mkdir yourrundir && cd yourrundir
sbatch ../job_scripts/qs_XX.sh
```
Two binaries should appear corresponding to tasks 1, 2: `quicksort-serial`, `quicksort-distributed`.
There are some batch scripts in `sheet2/job-scripts` to get you started.

## Hints

Only the functions `quicksort` and `quicksort_distributed` are used for grading. You can add additional helper functions as needed.



