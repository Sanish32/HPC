#include <stdlib.h>
#include <stdio.h>
#include <algorithm>
#if DISTRIBUTED
    #include "quicksort_distributed.h"
    #include <mpi.h>
#else
    #include "quicksort.h"
#endif
bool is_sorted(float* data, int n){
    for(int i=1;i<n;i++){
        if(data[i]<data[i-1]){
            return false;
        }
    }
    return true;
}

int main(int argc, char* argv[]){
    int rank;
    #if DISTRIBUTED
        MPI_Init(&argc, &argv);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank); 
    #else
        rank = 0;
    #endif
    constexpr int size = 20000;
    // constexpr int size = 5;
    float data[size];
    float* c_data = (float*)malloc(size*sizeof(float));
    float* g_data;
    float* result;
    float data_gpu[size];

    srand(12345678);
    for(int i=0;i<size;i++){
        data[i] = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
        // data[i] = i%8; 
        // if (i%2==0)
        // {
        //     data[i] = static_cast <float> (2);
        // }
        // else
        // {
        //     if (i==2001)
        //     {
        //         data[i] = static_cast<float> (2);
        //     }
        //     else
        //     {
        //         data[i] = static_cast <float> (1);
        //     }   
        // }
    }
    for(int i=0;i<size;i++) c_data[i] = data[i];
    #if DISTRIBUTED
        quicksort_distributed(data[0],0,size,c_data,MPI_COMM_WORLD);
    #else
        quicksort(data[0],0,size,c_data);
    #endif
    result = c_data;

    // printf("The contents of data before are:\n");
    // for (int i=0;i<size;i++)
    // {
    //     printf("%d\n",data[i]);
    // }
    // printf("\n");
    
    std::stable_sort(data,data+size);
    bool sorted = is_sorted(result,size);

    // int temp;
    // temp = data[0];
    // data[0] = data[1];
    // data[1] = temp;

    // printf("The contents of data after are:\n");
    // for (int i=0;i<size;i++)
    // {
    //     printf("%d\n",data[i]);
    // }
    // printf("\n");

    // printf("The contents of sorted is %d:\n",sorted);

    if(sorted){
        printf("Is sorted at rank %d\n", rank);
    }else{
        printf("Not sorted at rank %d!!\n", rank);
    }
    bool is_correct = true;
    for(int i=0;i<size;i++) 
    {
        // printf("result[%d] is %f and data[%d] is %f\n",i,result[i],i,data[i]);
        is_correct &= (result[i] == data[i]);
    }

    if(is_correct){
        printf("Correct at rank: %d\n",rank);
    }else{
        printf("Incorrect at rank %d!!\n",rank);
    }
    #if DISTRIBUTED
	MPI_Barrier(MPI_COMM_WORLD);
        MPI_Finalize();
    #endif
    return 0;

}

