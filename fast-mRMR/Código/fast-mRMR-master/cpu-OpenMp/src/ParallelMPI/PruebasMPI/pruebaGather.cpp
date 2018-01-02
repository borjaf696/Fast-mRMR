#include "mpi.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>


int main (int argc, char **argv)
{
	MPI_Init(&argc, &argv);

	int size, rank, cont = 0;
	int *vectorFinal, *vectoresParciales, * recvcounts, * displ;

	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	printf("Hola soy %d y empiezo\n",rank);

	vectoresParciales = (int*)malloc((2+1)*sizeof(int));
	vectorFinal = (int*)malloc((2*size+1)*sizeof(int));

	recvcounts = (int*) malloc(sizeof(int)*size);
	displ = (int*) malloc(sizeof(int)*size);

	for (int i = rank; i < rank+2; i++) {
		vectoresParciales[cont] = i;
		cont++;
	}
	if (rank == size-1)
		vectoresParciales[2] = 45;
	for (int i = 0; i < 2; i++)
		printf("Soy %d y mi valor es %d\n",rank,vectoresParciales[i]);

	for (int i = 0; i < size; ++i){
		recvcounts[i] = 2;
		if (i == size-1)
			recvcounts[i] = 3;
		displ[i] = i*2;
	}
	int send = 2;
	if (rank == size-1)
		send = 3;
	MPI_Allgatherv(vectoresParciales, send, MPI_INT,
        vectorFinal, recvcounts, displ, MPI_INT, MPI_COMM_WORLD);

	if (rank == 0)
		for (int i = 0; i< size*2+1;i++)
			printf("Valores finales %d\n",vectorFinal[i]);

	free(recvcounts);
	free(displ);
	free(vectoresParciales);
	free(vectorFinal);
	MPI_Finalize();
}
