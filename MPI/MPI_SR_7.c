#include <mpi.h>
#include <stdio.h>

int main(int argc, char** argv) {

  // Example code of using non-blocking Isend and Irecv
  // Similar to MPI_SR_1.c 
  // The new thing here is that in between send and wait you could do something else, as the communication is not blocking.
  // Try placing the Waits in different places. Does it matter where the Waits are?
  int my_id, your_id, bufsize=10000;
  int sendbuf[bufsize],recvbuf[bufsize];
  MPI_Status status;
  MPI_Request srequest,rrequest;

  /* Initialize the MPI execution environment */
  MPI_Init(&argc,&argv);

  /* Ranks */
  MPI_Comm_rank(MPI_COMM_WORLD, &my_id);

  your_id = 1 - my_id;

  sendbuf[5]=10*my_id-100;

  printf("Rank %d Initiating send to rank %d ...\n",my_id,your_id);
  /* First send...*/
  MPI_Isend(&sendbuf,bufsize,MPI_INT,your_id,0,MPI_COMM_WORLD,&srequest);
  printf("... rank %d returned from isend ...\n",my_id); 

  // 1st possible placement for Wait for sends
  MPI_Wait(&srequest,MPI_STATUS_IGNORE);
  printf("... rank %d finalized waiting of send requests!\n",my_id); 

  printf("Rank %d waiting to receive from rank %d ...\n",my_id,your_id);  
  /* Then receive */
  MPI_Irecv(&recvbuf,bufsize,MPI_INT,your_id,0,MPI_COMM_WORLD,&rrequest);
  printf("... rank %d returned from irecv ...\n",my_id); 

  // 2nd possible placement for Wait for sends
  //MPI_Wait(&srequest,MPI_STATUS_IGNORE);
  //printf("... rank %d finalized waiting of send requests!\n",my_id); 
  
  MPI_Wait(&rrequest,MPI_STATUS_IGNORE);
  printf("... rank %d finalized waiting of send requests and received this %d from %d!\n",my_id,recvbuf[5],your_id); 

  /* Finalize the MPI environment. */
  MPI_Finalize();
}

