//  Calculus in parallel using monte Carlo Method
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <time.h>

// Number of points
const int N_points = 1000000000;

// Generation the random numbers in the domain [0,1]X[0,1] 

int main(int argc, char const *argv[])
{
	int rank;
	int size;
	MPI_Status status;

	int count;
	double x;
	double y; 
	double pi;
	double pi_estim;
	double final_pi_estim;
	

	// Initialization of MPI
	MPI_Init(0, 0);

	// Get the rank of the parallel process
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	// Get the total number of process
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	// Obtain the number of points for this process
	const int p_points = N_points/size;

	// begin the random number generator
	srand(time(NULL)*(rank+1));

	// Monte Carlo
	printf("[%i] Generating %i points in [0,1]x[0,1]\n", rank, p_points);

	pi_estim = 0.0;
	count    = 0;

	for (int i = 1; i <= p_points; ++i)
	{
		x = (double)rand()/RAND_MAX;
		y = (double)rand()/RAND_MAX;
		if (pow(x, 2) + pow(y, 2) <= 1.0)
		{
			count += 1;
		}
	}
	pi_estim = 4*(double)count/p_points;

	printf("[%i] Estimate: %.8f\n", rank, pi_estim);
	
	// Processes that are not 0, send their partial estimate to rank 0
	if (rank != 0)
	{
		MPI_Send(&pi_estim, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	}

	// process 0 accumulates the partial estimates of the other
	// processes and produces the final estimate

	else if (rank == 0)
	{
		final_pi_estim = pi_estim;
		for (int r = 1; r < size; ++r)
		{
			MPI_Recv(&pi_estim, 1, MPI_DOUBLE, r, 0, MPI_COMM_WORLD, &status);
			final_pi_estim += pi_estim;
		}
	}
	final_pi_estim /= size;
	// barrier to wait every process to finish before printing result 
	MPI_Barrier(MPI_COMM_WORLD);
	
	// pi estimation
	if (rank == 0)
	{
		printf("\nTrue value: %.15f\n", acos(-1));
		printf("\nResult: %.15f\n", final_pi_estim);
	}

	// finalize MPI
	MPI_Finalize();
}