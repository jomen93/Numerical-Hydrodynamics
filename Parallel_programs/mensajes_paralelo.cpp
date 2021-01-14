#include <stdio.h> 
#include "mpi.h"

int main( int argc, char *argv[] )
{
    int rank;
    int size;
    //  initialization of MPI enviroment
    MPI_Init( 0, 0 );
    
    // we get the ID of each rank
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    // get the number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    

    // send messages to the processor on the right

    // only if exist processor on the right
    if (rank <= (size-1) )    
    {
    	float mE;
    	int count;
    	int to;
    	int tag;

    	//  el mensaje a enviar es el rank de cada proceso
    	//  the message to send is the rank of each process
    	mE = rank;       
    	// send only one single data
    	count = 1;         
    	// information arrival processor
    	to    = rank+1;   
    	// if the processor is the last, redefine to the first
    	if (to == size)
    	{
    		to = 0;
    	}
    	// label must be match to recv label
    	tag   = 0;         
    	printf("Process %d, sends -> %f to process %d\n", rank, mE, to);
    	MPI_Send(&mE, count, MPI_FLOAT, to, tag, MPI_COMM_WORLD ) ;
    }

    // receipt of message from left processor, only if the left processor exist
    if (rank >= 0 )           
    {
    	int count;
    	int from;
    	int tag;
    	// in this variable save the receive data
    	float mR;
    	// a single data will be received
    	count = 1;         
    	// rank of the left processor
    	from  = rank-1;    
    	// if data comes from the last processor, redefine from
    	if (from == -1)
    	{
    		from = 5;
    	}
    	// label, must be match with the send label
    	tag   = 0;         
    	// status variable
    	MPI_Status stat;   
    	MPI_Recv(&mR, count, MPI_FLOAT, from, tag, MPI_COMM_WORLD, &stat ) ;
    	printf("Process %d, recieved <- %f from process %d\n", rank, mR, from);
    }

    //  finaliza MPI
    MPI_Finalize();
} 
