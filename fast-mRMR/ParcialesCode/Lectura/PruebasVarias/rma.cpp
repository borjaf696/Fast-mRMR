#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <assert.h>

int main(int argc, char *argv[])
{
	int rank, n=2;
	MPI_Comm comm;
	int *inbuf,*outbuf;
	MPI_Win win;
	int i;

	if (argc != 2) {
    		fprintf(stderr, "Usage: simple_put_test num_test\n");
    		exit(1);
  	}

	int num_test = atoi(argv[1]);

	MPI_Init(NULL, NULL);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_split(MPI_COMM_WORLD, rank <=1, rank, &comm);

	/* Solo trabajan 0 y 1 */
	if (rank > 1){return 0;}
	/*Ventanas*/
	if(rank == 0){
		outbuf = (int *)malloc(sizeof(int) * n);
  		assert(outbuf != NULL);
		MPI_Win_create(MPI_BOTTOM, 0, sizeof(int),MPI_INFO_NULL, comm, &win);
		outbuf[0]=5;
		outbuf[1]=6;
		printf("Soy el proceso %d y tengo un %d y un %d\n",rank,outbuf[0],outbuf[1]);
	}
	else if (rank ==1){
	  	inbuf = (int *)malloc(sizeof(int) * n);
  		assert(inbuf != NULL);
	  	MPI_Win_create(inbuf, n * sizeof(int), sizeof(int),MPI_INFO_NULL, comm, &win);
		inbuf[0]=7;
		printf("Soy el proceso %d y tengo un %d\n",rank,inbuf[0]);
	}

	printf("Soy el proceso %d y hago %d iteraciones\n",rank,num_test);

	int results = 0;
	for(i=0;i<num_test;i++){
		/*Solo 0 escribe en proceso 1*/
		MPI_Win_fence(0,win);
		if(rank==0){
			MPI_Put(&outbuf[1], 1, MPI_INT, 1, 0, 1, MPI_INT, win);
			MPI_Put(&outbuf[1], 1, MPI_INT, 1, 1, 1, MPI_INT, win);
		}
		/*Completar la operación*/
		MPI_Win_fence(0,win);
		if(rank==1){
			if(inbuf[0]==6) results = results + 1;
			if (inbuf[1]==6) results = results + 1;
			printf("%d %d %d \n",inbuf[0],inbuf[1],inbuf[2]);
		}
	}
	if(rank==1) printf("Soy el proceso %d y encontre un 6 %d veces de %d(%f por ciento)\n",rank,results,num_test,results*100.0/num_test);


	/*Liberar el comunicador y ventana*/
	MPI_Barrier(comm);
	MPI_Comm_free(&comm);
	MPI_Win_free(&win);
	if(rank==1) free(inbuf);
	if(rank==0) free(outbuf);
	MPI_Barrier(MPI_COMM_WORLD);
  	MPI_Finalize();
}
