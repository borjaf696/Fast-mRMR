#include "mpi.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>


int main (int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    int size, rank, cont = 0;
    int *vectorFinal, *vectoresParciales;
    MPI_Comm comunicadorExtra;

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    printf("Hola soy %d y empiezo\n",rank);

    vectoresParciales = (int*)malloc((2+1)*sizeof(int));
    vectorFinal = (int*)malloc((2*size+1)*sizeof(int));


    for (int i = rank; i < rank+2; i++) {
        vectoresParciales[cont] = i;
        cont++;
    }
    vectoresParciales[2] = 45+rank;

    for (int i = 0; i < 2; i++)
        printf("Soy %d y mi valor es %d\n",rank,vectoresParciales[i]);

    int send = (2*size+1) / size, sendGlobal = (2*size+1) / size;
    if (rank < (2*size+1)%size)
        for (int i = 0; i < ((2*size+1)%size); ++i)
            send++;
    MPI_Comm_split(MPI_COMM_WORLD, rank < ((2*size+1)%size), rank, &comunicadorExtra);
    printf("Soy %d y tengo que mandar %d\n",rank, send);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Request *requests;
    requests = (MPI_Request*)malloc((send)*sizeof(MPI_Request));
    for (int i= 0; i< send; i++)
        if (i != sendGlobal)
            MPI_Igather(&vectoresParciales[i],1,MPI_INT,
                    &vectorFinal[size*i], 1, MPI_INT, 0, MPI_COMM_WORLD, &requests[i]);
        else
            MPI_Igather(&vectoresParciales[i],1,MPI_INT,
                    &vectorFinal[size*i], 1, MPI_INT, 0, comunicadorExtra, &requests[i]);


    int flag = 0;
    while(!flag)
        MPI_Testall(send,requests, &flag, MPI_STATUSES_IGNORE);

    if (rank == 0)
        for (int i = 0; i< size*2+1;i++)
            printf("Valores finales %d\n",vectorFinal[i]);

    free(vectoresParciales);
    free(vectorFinal);
    MPI_Finalize();
}
