// In this exercise, you will implement functions to find the maximum integer in an array distributed across multiple GPUs (devices).
// Additionally you will extend the quicksort you have implemented in exercise 2. to work on multiple GPUs
// The exercise consists of three subtasks related to reductions:


// Implement the functions reduce_kernel and reduce in src/reduce-single.cu. 
// Your task here is to find the maximum integer in an array allocated on a single device and return the result with reduce.

#include <stdio.h>
#include <stdlib.h>

#include "errchk.h"
#include "reduce.cuh"

// The device kernel for finding the maximum integer in `arr`
__global__ void
reduce_kernel(const int* in, const size_t count, int* out)
{
  // EXERCISE 4: Your code here
}

// The host function for calling the reduction kernel
int
reduce(const int* arr, const size_t initial_count)
{
  // EXERCISE 4: Your code here
  // Input:
  //  arr           - An array of integers
  //  initial_count - The number of integers in `arr`
  //
  // Return: the maximum integer in `arr`
  //
  // Task: allocate memory on the GPU, transfer `arr` into
  // the allocated space, and apply `reduce_kernel` iteratively on it
  // to find the maximum integer. Finally, move the result back to host
  // memory and return it.

  int num_devices, size = initial_count * sizeof(float);
  cudaGetDeviceCount(&num_devices);

  float 

  cudaMalloc((void**)&d_arr, size);
  cudaMemcpy(d_arr, arr, size, cudaMemcpyHostToDevice);

  reduce_kernel<<<1,256,size>>>(

  return -1;
}
