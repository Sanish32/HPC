#include <mpi.h>
#include <stdio.h>

int main(int argc, char** argv) {

  /* Blocking send and receive in between *two* processes */
  /* First send and receive */
  /* Start with bufsize 1, then keep on increasing it */
  int my_id, your_id, bufsize=10000;    // First try commenting out this
  //int my_id, your_id, bufsize=100000;    // Secondly comment out this line. What is happening and why?
  int sendbuf[bufsize],recvbuf[bufsize];
  MPI_Status status;

  printf("Initializing MPI...\n");
  printf("size of int= %d \n",sizeof(my_id));

  /* Initialize the MPI execution environment */
  MPI_Init(&argc,&argv);

  /* Ranks */
  MPI_Comm_rank(MPI_COMM_WORLD, &my_id);

  your_id = 1 - my_id;
  printf("Rank %d here, hello!\n",my_id);

  printf("Initiating send to rank %d\n",your_id);
  /* First send...*/
  MPI_Send(&sendbuf,bufsize,MPI_INT,your_id,0,MPI_COMM_WORLD); 

  printf("Waiting to receive from rank %d\n",your_id);  
  /* Then receive */
  MPI_Recv(&recvbuf,bufsize,MPI_INT,your_id,0,MPI_COMM_WORLD,&status);

  printf("Receive finished from %d\n",your_id);

  /* Finalize the MPI environment. */
  MPI_Finalize();
}

