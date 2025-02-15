#include <mpi.h>
#include <stdio.h>

int main(int argc, char** argv) {

  /* Blocking send and receive in between *two* processes */
  /* First receive and then send */
  /* What happens and why? */

  int my_id, your_id, bufsize=1;
  int sendbuf[bufsize],recvbuf[bufsize];
  MPI_Status status;

  /* Initialize the MPI execution environment */
  MPI_Init(&argc,&argv);

  /* Ranks */
  MPI_Comm_rank(MPI_COMM_WORLD, &my_id);

  your_id = 1 - my_id;
  printf("Rank %d here, hello!\n",my_id);

  /* Receive */
  printf("Waiting to receive from rank %d\n",your_id);  
  MPI_Recv(&recvbuf,bufsize,MPI_INT,your_id,0,MPI_COMM_WORLD,&status);

  printf("Receive finished from %d\n",your_id);

  printf("Initiating send to rank %d\n",your_id);
  /* Send...*/
  MPI_Send(&sendbuf,bufsize,MPI_INT,your_id,0,MPI_COMM_WORLD);

  printf("Send finished to %d\n",your_id);

  /* Finalize the MPI environment. */
  MPI_Finalize();
}

