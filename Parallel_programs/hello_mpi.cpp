#include <stdio.h>
#include "mpi.h"

int main( int argc, char * argv[])
{
	int rank;
	int size;
	// inicializacion de MPI
	MPI_Init(0, 0);

	// Se obtiene el ID de cada núcleo (rank)
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	// se obtiene el número de procesos que etán coriendo
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	printf("Hola mundo desde el procesador %d de %d \n",rank ,size);

	// finaliza MPI
	MPI_Finalize();
}