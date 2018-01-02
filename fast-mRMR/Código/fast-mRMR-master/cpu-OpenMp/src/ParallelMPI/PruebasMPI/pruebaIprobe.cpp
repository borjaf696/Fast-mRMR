#include "mpi.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>


int main (int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    int size, rank, cont = 0, mensajesFinalizacion = 0, msg, msg2 = -1, flag = 0,
            * bufferDst = (int*) malloc(4*sizeof(int)), count,buffer, msgReceived = 0;
    MPI_Request *request = (MPI_Request*)malloc(2*sizeof(MPI_Request));
    MPI_Status status;

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    msg = rank*1029;

    printf("Hola soy %d y empiezo\n",rank);

    MPI_Isend(&msg, 1, MPI_INT,
            0, 0, MPI_COMM_WORLD, &request[0]);

    if (!rank)
        while (mensajesFinalizacion != size-1){
            while (!flag)
                MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD,
                           &flag, &status);
            MPI_Get_count(&status, MPI_INT, &count);
            MPI_Recv(&buffer, count, MPI_INT, status.MPI_SOURCE,
                    status.MPI_TAG, MPI_COMM_WORLD, &status);
            if (buffer == -1)
                mensajesFinalizacion++;
            else
                bufferDst[status.MPI_SOURCE] = buffer;
            msgReceived++;
            flag = 0;
        }
    MPI_Isend(&msg2, 1, MPI_INT,
        0, 0, MPI_COMM_WORLD, &request[1]);
    while (!flag)
        MPI_Testall(2,request, &flag, MPI_STATUSES_IGNORE);
    if (!rank)
        for (int i = 0;i < mensajesFinalizacion; ++i)
            printf("Mensaje %d-%d\n",i,bufferDst[i]);
    printf("Mensajes recibidos %d\n",msgReceived);
    MPI_Finalize();
}
