
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#include <iostream>
int partition(float pivot,int start, int end,float* &data)
{
	int i=start;
	float temp;

	for (int j=start+1;j<end;j++)
	{	
		if (pivot>data[j])
		{
			i++;
			if (i!=j)
			{
				temp = data[i];
				data[i] = data[j];
				data[j] = temp;
				// We dont care about the relative order for values which are not same as pivot value, because my quicksort/partitioning does it automatically
				// Smaller values are pushed front of the array and bigger values are pushed back linearly or sequentially

			}
		}
	}

	// Now, we know that value of i (which symbolizes the pivotal point
	// for next iteration of quicksort), anything before i is our focus/interest
	// because changing pivot may change the relative order for the same values

	// While on the other hand, anything after i can be ignored because there is no 
	// way that relative order can be disturbed. Swapping pivot point with i is when
	// relative order can be disturbed

	// Value of track should always be start and not 0 because second half of the quicksort works from i+1 position from previous array
	// or start of new array and not the 0th position (which is the only case for first half of the array)

	int track=start;
	for (int k=start;k<i;k++)
	{
		// The value of pivot is always same because we are swapping the same values, maintaining the relative order
		if (pivot==data[k])
		{
			temp = data[track];
			data[track] = data[k];
			data[k] = temp;

			// We are switching the relative order OR in simple words we are pushing 
			// same values ahead but ith position is never reached in this loop
			track = k;
		}
	}

	// It is possible that i and k position can have same value of pivot
	// If we swap those places, we will run into unrelative order
	
	if (data[track]!=data[i])
	{
		temp = data[track];
		data[track] = data[i];
		data[i] = temp;
	}
	// There is no else because if the values are same, they are already in relative order

	return i;
}

void quicksort(float pivot, int start, int end, float* &data)
{
/**
	Exercise 1: Your code here
	Input:
		pivot: a pivot value based on which to split the array in to less and greater elems
		start: starting index of the range to be sorted
		end: exclusive ending index of the range to be sorted
		data: array of floats to sort in range start till end
	Return:
		upon return the array range should be sorted
	Task: 	
		to sort the array using the idea of quicksort in a stable manner
		a sort is stable if it maintains the relative order of elements with equal values
**/

	// base case
	if (start>=end-1) return; 

	int i = partition(pivot,start,end,data);

	quicksort(data[start],start,i,data);
	quicksort(data[i+1],i+1,end,data);
}



void quicksort_distributed(float pivot, int start, int end, float* &data,  MPI_Comm comm)
{
/**
        Exercise 2: Your code here
        Input:
                pivot: a pivot value based on which to split the array in to less and greater elems
                start: starting index of the range to be sorted
                end: exclusive ending index of the range to be sorted
                data: array of floats to sort in range start till end
		      the array is the same for each MPI process
		comm: the communicator of the MPI processes
        Return:
                upon return each MPI process should have their array sorted
        Task:   
                To sort the array using the idea of quicksort in a stable manner using each MPI process.
                A sort is stable if it maintains the relative order of elements with equal values.

		Truly split the work among the process i.e. don't have each process independently sort the array.
		You should stick to the idea of quicksort also for the splitting of the array between processes.
		As an example you should not split the array into N equal sized parts (where N is the number of processors),
		sort them and then combine then (this would be following mergesort). Instead the splitting should be based on the
		pivot values
	
	       	 	

	Hints:
		why is the communicator given to you?
		maybe it could be a good idea to use MPI_Comm_split?
		could you reuse your quicksort implementation for the base case of a single process?

**/

	int i = partition(pivot,start,end,data);

	// Splitting section
	int nprocs;
	MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
	// printf("The total number of %d processors used.\n",nprocs);

	int comm_size;
	MPI_Comm_size(comm,&comm_size);
	
	int rank;
	MPI_Comm_rank(comm,&rank);

	if (nprocs==1 || comm_size==1)
	{
		quicksort(data[start],start,i,data);
		quicksort(data[i+1],i+1,end,data);
		return;
	}


	// float determine = static_cast<float> (i-start) / static_cast<float> (end-start);
	// int num_of_proc_on_left = static_cast<int> (determine) * comm_size;
	// int num_of_proc_on_right = comm_size - num_of_proc_on_left;

	int color = (rank<comm_size/2) ? 1:0;
	// if (rank<num_of_proc_on_left)
	// {
	// 	color=0;
	// }
	// else color=1;

	MPI_Comm new_comm;
	MPI_Comm_split(comm,color,rank,&new_comm);

	int root_rank = 0;
	int last_rank = comm_size-1;

	

	if (color==0)
	{  
		quicksort_distributed(data[start],start,i,data,new_comm);
		MPI_Bcast(&data[start],i-start,MPI_FLOAT,rank,comm);
	}
	else
	{  
		quicksort_distributed(data[i+1],i+1,end,data,new_comm);
		MPI_Bcast(&data[i+1],end-i-1,MPI_FLOAT,rank,comm);
	}

	MPI_Comm_free(&new_comm);
	return;
}
