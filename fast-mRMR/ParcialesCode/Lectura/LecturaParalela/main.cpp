#include "utils.h"
#include "RawData.h"

#include <mpi.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

int main (int argc, char * argv[]){
	char * file;
	int rank, nprocs;

	MPI_Init(&argc,&argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);


	if (argc >= 2){
		file = argv[1];
	}
	printf("Empiezo %s\n",file);

	RawData rw = RawData(file);

	printf("Termino");

	MPI_Finalize();
}
