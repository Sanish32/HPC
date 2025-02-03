# Exercise 4 - Distributed Computing on GPUs

## Introduction

In this exercise, you will implement functions to find the maximum integer in an array distributed across multiple GPUs (devices).
Additionally you will extend the quicksort you have implemented in exercise 2. to work on multiple GPUs


The exercise consists of three subtasks related to reductions:

  1. Implement the functions `reduce_kernel` and `reduce` in `src/reduce-single.cu`. Your task here is to find the maximum integer in an array allocated on a single device and return the result with `reduce`.
  
  2. Extend your implementation to work on multiple devices on a single node in `src/reduce-multi.cu`. The functions to implement are `reduce_kernel` (you can re-use your solution from subtask (1) and `reduce`. The array is now distributed across multiple devices and each device should compute a single result. You can perform the final reduction step on the host by transferring the data with the device-to-host variant of `cudaMemcpy`. In your slurm command, you should allocate a single task (`--ntasks=1`) but multiple devices (`--gres=gpu:v100:n`).

  3. Extend your implementation to work on multiple devices on multiple nodes using MPI in `src/reduce-mpi.cu`. The instructions are otherwise the same as in subtask (2), however, in this exercise you should allocate one process per device (for example `--ntasks-per-node=4 --nodes=2 --gres=gpu:v100:4`) and ensure your implementation works on multiple nodes.


And one subtask related to quicksort:

  4. Extend your quicksort implementation to work on GPUs. Start first by implementing `quicksort` in `src/quicksort_gpu.cu`. Then extend the GPU implementation to multiple GPU by implementing `quicksort_distributed` in `quicksort_distributed_gpu.cu` (you are allowed to and highly encouraged to re-use your previous solutions in exercise 2. for this).   In this exercise you should allocate one process per device (for example `--ntasks-per-node=4 --nodes=2 --gres=gpu:teslap100:4`).


> Note: The performance of your implementation is not graded and to get full points, you only need to ensure your implementations give the correct results. However, you should use the hardware relatively efficiently: an implementation that uses a single CUDA thread in task 1, or a single device in tasks 2 or 3, will receive 0 points. Additionally try to use CUDA-aware MPI instead of first moving data to the CPU and then communicating in task 4.
