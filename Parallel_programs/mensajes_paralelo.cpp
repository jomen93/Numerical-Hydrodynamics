#include <stdio.h>
#include "mpi.h"

int main( int argc, char *argv[])
{
	int rank;
	int size;

	// Initialize MPI
	MPI_Init(0, 0);
	// get the ID of each kernel (rank)
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	// get the number of run processes
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	printf("Hello world from processor %d of %d \n", rank, size);

	// send message to processor on the right
	// only if existe process on the right
	if (rank < (size-1))
	{
		// the message to send is the rank of face process
		float mE  = rank;    // the message to send is the ran of each process
		int count = 1;       // only one date is going to be sent 
		int to    = rank+1;  // processor on th right
		int tag   = 0;       // label, this must gonna be equal to recv    
		MPI_Send(&mE, count, MPI_FLOAT, to, tag, MPI_COMM_WORLD);
		printf("Process %d, sends --> %f to process %d\n", rank, mE, to);

	}

	if (rank > 0)
	{
		// message reception for processor on the left
		int count = 1;         // only one date is going to be receive 
		int from    = rank-1;  // rank of the left processor 
		int tag   = 0;         // label, this must gonna be equal to send   
		float mR;   		   // here saves the receive message
		MPI_Status stat;       // status variable
		MPI_Recv(&mR, count, MPI_FLOAT, from, tag, MPI_COMM_WORLD, &stat);
		printf("Process %d, recieves <-- %f from process %d\n", rank, mR, from);
	}

	// Finalize MPI
	MPI_Finalize();

	return 0; 

}