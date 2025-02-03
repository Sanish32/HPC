#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>


int
main(void)
{
// Check how many MPI tasks we have
  MPI_Init(NULL, NULL);
  int nprocs, pid;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &pid);
  if (pid == 0) printf("Num MPI processes: %d\n", nprocs);


// Check how many GPUs we have
  int num_devices;
  cudaGetDeviceCount(&num_devices);
  if (pid == 0) printf("Num devices: %d\n", num_devices);

// Where are we assigned to run?
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  int name_len;
  MPI_Get_processor_name(processor_name, &name_len);

// Make the processes to greet the class 
  printf("Hello class from node %s, rank %d out of %d MPI processes\n",
         processor_name, pid, nprocs);

// Make the GPUs greet the class
  
  for (int i = 0; i < num_devices; ++i) {
     cudaSetDevice(i);
     printf("Hello class, GPU %d active on node %s\n",i,processor_name);
  }

  MPI_Finalize();

}

