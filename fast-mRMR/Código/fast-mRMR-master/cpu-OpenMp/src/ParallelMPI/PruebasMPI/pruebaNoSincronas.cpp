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

    int send = 2;
    if (rank == size-1)
        send = 3;
    MPI_Request *requests;
    requests = (MPI_Request*)malloc((send)*sizeof(MPI_Request));
    for (int i= 0; i< 2; i++)
        MPI_Igather(&vectoresParciales[i],1,MPI_INT,
            &vectorFinal[size*i], 1, MPI_INT, 0, MPI_COMM_WORLD, &requests[i]);
    if (rank == size-1)
        MPI_Isend(&vectoresParciales[send-1], 1, MPI_INT,
            0,0, MPI_COMM_WORLD,&requests[send-1]);
    if (!rank)
        MPI_Recv(&vectorFinal[2*size],1,MPI_INT,
            MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);

    int flag = 0;
    while(!flag)
        MPI_Testall(send,requests, &flag, MPI_STATUSES_IGNORE);

    if (rank == 0)
        for (int i = 0; i< size*2+1;i++)
            printf("Valores finales %d\n",vectorFinal[i]);

    free(recvcounts);
    free(displ);
    free(vectoresParciales);
    free(vectorFinal);
    MPI_Finalize();
}
