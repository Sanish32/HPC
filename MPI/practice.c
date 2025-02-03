#include <mpi.h>
#include <stdio.h>

int main(int argc,char** argv)
{
    MPI_Init(&argc,&argv);

    int data;
    MPI_Win window;

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD,&world_rank);

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD,&world_size);

    MPI_Win_create(&data,sizeof(int),sizeof(int),MPI_INFO_NULL,MPI_COMM_WORLD,&window);

    MPI_Win_fence(0,window);

    if (world_rank<world_rank-1)
    {
        printf("Rank %d here and sending %d to %d\n",world_rank,world_rank,world_rank+1);
        MPI_Put(&world_rank,1,MPI_INT,world_rank+1)
    }
}