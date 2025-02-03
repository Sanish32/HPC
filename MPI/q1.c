#include <mpi.h>
#include <stdio.h>
#include <limits.h>
#include <stdlib.h>
#include <unistd.h>  // For gethostname

int main(int argc, char** argv) {
    int my_id, your_id;
    int bufsize=2;
    int *sendbuf= (int*) malloc(bufsize * sizeof(int)), *recvbuf=(int*) malloc(bufsize * sizeof(int)); 
    MPI_Status status;
    MPI_Request srequest,rrequest;
    double start, end;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
    your_id = 1 - my_id;

    // Open the CSV file for writing
    const char *filenames[3] = {"send_recv_round.csv","isend_irecv_round.csv","bsend_recv_round.csv"};
    FILE *files[3];

    // Just want to make sure the length of filenames array
    int length_filenames = sizeof(filenames) / sizeof(filenames[0]);
    // printf("The length of the filenames array is %d\n",length_filenames);

    if (my_id == 0) {
        for (int i=0;i<length_filenames;i++)
        {
            files[i] = fopen(filenames[i],"w");

            if (files[i] == NULL) 
            {
                fprintf(stderr,"Unable to open '%s' file!\n", filenames[i]);

                // perror is bit useless because it can only take one string literal or agrument
                // perror("Unable to open file!\n");
                MPI_Abort(MPI_COMM_WORLD,1);
            }
            else if (files[i] != NULL)
            {
                // Write the header (first line) in the different csv files
                fprintf(files[i],"BufferSize,TimeTaken\n");
            }
        }
    }


    //Just trying out how to know which rank process is running on which node
    // if (my_id == 0) 
    // {   
    //     char buffer_hostname[256];
    //     // Get the hostname of the node
    //     gethostname(buffer_hostname, sizeof(buffer_hostname));
    //     printf("The node allocated for rank %d is %s\n",my_id,buffer_hostname);
    // }

    // else if(my_id==1) 
    // {       
    //     char buffer_hostname[256];
    //     // Get the hostname of the node
    //     gethostname(buffer_hostname, sizeof(buffer_hostname));
    //     printf("The node allocated for rank %d is %s\n",my_id,buffer_hostname);
    // }


    // First combination using classic MPI_Send and MPI_Recv
    for (bufsize = 1; bufsize <= 100000000; bufsize *=2) {

        sendbuf = (int*)realloc(sendbuf, bufsize * sizeof(int));
        recvbuf = (int*)realloc(recvbuf, bufsize * sizeof(int));

        // Synchronization to ensure both processes are ready before timing
        MPI_Barrier(MPI_COMM_WORLD);

        if (my_id == 0) {
            start = MPI_Wtime(); // Start timing
            MPI_Send(sendbuf, bufsize, MPI_INT, your_id, 0, MPI_COMM_WORLD);
            MPI_Recv(recvbuf, bufsize, MPI_INT, your_id, 0, MPI_COMM_WORLD, &status);
            end = MPI_Wtime(); // End timing

            // printf("Time taken for rank %d sending and receiving %d bytes: %f seconds\n", my_id, bufsize * sizeof(int), end - start);
            fprintf(files[0], "%d,%f\n", bufsize * sizeof(int), end - start); // Write packet size and time taken
        } else {
            MPI_Recv(recvbuf, bufsize, MPI_INT, your_id, 0, MPI_COMM_WORLD, &status);
            MPI_Send(recvbuf, bufsize, MPI_INT, your_id, 0, MPI_COMM_WORLD);   
        }
    }

    // Second combination using MPI_Isend and MPI_Irecv
    for (bufsize = 1; bufsize <= 100000000; bufsize *=2) {

        sendbuf = (int*)realloc(sendbuf, bufsize * sizeof(int));
        recvbuf = (int*)realloc(recvbuf, bufsize * sizeof(int));

        // Synchronization to ensure both processes are ready before timing
        MPI_Barrier(MPI_COMM_WORLD);

        if (my_id == 0) {
            start = MPI_Wtime(); // Start timing
            MPI_Isend(sendbuf, bufsize, MPI_INT, your_id, 0, MPI_COMM_WORLD,&srequest);

            // Wait for both operations to complete
            MPI_Wait(&srequest, MPI_STATUS_IGNORE);
            MPI_Irecv(recvbuf, bufsize, MPI_INT, your_id, 0, MPI_COMM_WORLD, &rrequest);

            // Wait for both operations to complete
            MPI_Wait(&rrequest, MPI_STATUS_IGNORE);

            end = MPI_Wtime(); // End timing
            // printf("Time taken for rank %d sending and receiving %d bytes: %f seconds\n", my_id, bufsize * sizeof(int), end - start);
            fprintf(files[1], "%d,%f\n", bufsize * sizeof(int), end - start); // Write packet size and time taken
        } else {
            MPI_Irecv(recvbuf, bufsize, MPI_INT, your_id, 0, MPI_COMM_WORLD, &rrequest);

            // Wait for both operations to complete
            MPI_Wait(&rrequest, MPI_STATUS_IGNORE);
            MPI_Isend(recvbuf, bufsize, MPI_INT, your_id, 0, MPI_COMM_WORLD,&srequest);

            // Wait for both operations to complete
            MPI_Wait(&srequest, MPI_STATUS_IGNORE);
        }
    }

    // Third combination using MPI_Bsend and MPI_Recv
    for (bufsize = 1; bufsize <= 100000000; bufsize *= 2) {
        sendbuf = (int*)realloc(sendbuf, bufsize * sizeof(int));
        recvbuf = (int*)realloc(recvbuf, bufsize * sizeof(int));

        // Synchronization to ensure both processes are ready before timing
        MPI_Barrier(MPI_COMM_WORLD);

        int buffer_size = MPI_BSEND_OVERHEAD + bufsize * sizeof(int);
        int *buffer = (int*) malloc(buffer_size);
        MPI_Buffer_attach(buffer, buffer_size);

        if (my_id == 0) {
            start = MPI_Wtime(); // Start timing

            MPI_Bsend(sendbuf, bufsize, MPI_INT, your_id, 0, MPI_COMM_WORLD);
            MPI_Recv(recvbuf, bufsize, MPI_INT, your_id, 0, MPI_COMM_WORLD, &status);
            end = MPI_Wtime(); // End timing
            // printf("Time taken for rank %d sending and receiving %d bytes: %f seconds\n", my_id, bufsize * sizeof(int), end - start);
            fprintf(files[2], "%d,%f\n", bufsize * sizeof(int), end - start); // Write packet size and time taken
        } else {
            MPI_Recv(recvbuf, bufsize, MPI_INT, your_id, 0, MPI_COMM_WORLD, &status);
            MPI_Bsend(recvbuf, bufsize, MPI_INT, your_id, 0, MPI_COMM_WORLD);
        }

        // Detach the buffer after use
        MPI_Buffer_detach(&buffer, &buffer_size);
        free(buffer);
    }


    // Close the CSV file after writing results
    for (int i = 0; i < length_filenames; i++) {
        if (files[i] != NULL) {
            fclose(files[i]); // Close only if successfully opened
        }
    }

    // Free allocated memory
    free(sendbuf);
    free(recvbuf);

    MPI_Finalize();
    return 0;
}
