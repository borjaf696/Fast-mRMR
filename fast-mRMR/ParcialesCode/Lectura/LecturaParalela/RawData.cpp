/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 *    http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
#include <mpi.h>
#include "RawData.h"
#include <string.h>


/**
 * Constructor that creates a rawData object.
 *
 * @param data_table this is a matrix of bytes containing the translated csv data.
 * @param ds the number of data samples.
 * @param fs the number of features that each sample has.
 */
RawData::RawData(char * filename) {
	dataFile = fopen(filename, "rb");
	MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, &dataFile2);
	calculateDSandFS();
	loadData();
	printf("Por probar\n");
	loadData2();
	//loadData3();
	calculateVR();
	MPI_File_close(&dataFile2);
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank == 0)
        comparacion();
}

RawData::~RawData() {
}

/**
 *
 */
void RawData::destroy() {
	free(valuesRange);
	free(data);
	free(data2);
}
/*Comparamos data2 y data*/
void RawData::comparacion(){
    printf("Comparamos\n");
    uint i,j;;
    for (i = 0; i < datasize; i++){
        for (j = 0; j < featuresSize; j++)
            if (data[j*datasize+i] > data2[j*datasize+i] || data[j*datasize+i] < data2[j*datasize+i])
                break;
        if (j != featuresSize){
            printf("Algo fue mal salimos j = %d i = %d data = %c data3=%c \n",j,i,data[j*datasize+i],data2[j*datasize+i]);
            break;
        }
    }
}
/**Calculates
 *	DataSize: Number of patterns or samples
 *  FeaturesSize: Number of features
 */
void RawData::calculateDSandFS() {
	uint featuresSizeBuffer[1];
	uint datasizeBuffer[1];
	fread(datasizeBuffer, sizeof(uint), 1, dataFile);
	fread(featuresSizeBuffer, sizeof(uint), 1, dataFile);
	datasize = datasizeBuffer[0];
	featuresSize = featuresSizeBuffer[0];
}

void RawData::loadData2(){
    //Uso de ventanas -> RMA
	uint i, j;
	int rank,size,count, FsPerProcess, DsPerProcess, rest;
	MPI_Win win;
	MPI_Offset sizeFile;
	MPI_Comm comm;
	MPI_Status status;
	//LLamadas tipicas
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_split(MPI_COMM_WORLD, rank <= size, rank, &comm);
	//Creamos una ventana para la variable data2 de manera que sea accesible por todos los procesos.
	/*Cada proceso puede exponer memoria sobre la misma MPI_Win lo que cambia (recordemos) es la
	posicion u offset dentro de esta. Esta funcion devuelve un objeto que representa a los
	procesos del comunicador que pueden acceder a la ventana.*/
	t_data buffer[1], *dataProcess;
	MPI_File_seek(dataFile2, 8, MPI_SEEK_SET);
	MPI_File_get_size(dataFile2, &sizeFile);
	MPI_File_seek(dataFile2, rank*sizeof(t_data), MPI_SEEK_CUR);
	if (rank == 0){
        data2 = (t_data*) calloc(featuresSize, sizeof(t_data) * datasize);
        //MPI_Win_allocate(datasize*featuresSize*sizeof(t_data), 0, MPI_INFO_NULL, MPI_COMM_WORLD, &data2, &win);
        MPI_Win_create(data2,sizeof(t_data)*featuresSize*datasize,sizeof(t_data),MPI_INFO_NULL,comm, &win);
        printf("Tengo en data2 %c\n",data2[10]);
    }else
		MPI_Win_create(MPI_BOTTOM, 0, sizeof(t_data),MPI_INFO_NULL, comm, &win);
	for (i = rank; i < datasize; i+=size)
		for (j = 0; j < featuresSize; j++){
			MPI_File_read(dataFile2, &buffer, 1, MPI_BYTE, &status);
			MPI_Win_fence(0,win);
			//El 4 parametro es el proceso, el quinto el desplazamiento dentro de este (si se ha puesto sizeof en la ventana no hace falta ponerlo aqui abajo de nuevo.
			MPI_Put(&buffer[0], 1, MPI_BYTE, 0, j*datasize+i, 1,MPI_BYTE, win);
			MPI_Win_fence(0,win);
			MPI_File_seek(dataFile2, (size)*sizeof(t_data), MPI_SEEK_CUR);
		}
    if (!rank)
        for (i = 0; i < datasize;i+=size)
            for (j = 0; j < featuresSize; j++)
                printf("Tengo en data2 %c\n",data2[j*datasize+i]);
}

void RawData::loadData3(){
    //Uso de igather
    int rank_mpi,size_mpi, flag = 0, *recvcount,*recvdisp;
    MPI_Offset sizeFile;
    uint i, j;
    t_data buffer[1],*individual_buffer;
    //MPI Declarations
    MPI_Request *request;
    MPI_Status status;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank_mpi);
    MPI_Comm_size(MPI_COMM_WORLD, &size_mpi);

    //Operaciones I/O: MPI_File_seek -> actualiza los punteros INDIVIDUALES de cada proceso
	MPI_File_get_size(dataFile2, &sizeFile);
    //MPI_File_seek(dataFile2, rank_mpi-1,MPI_SEEK_CUR);

    int modulo = featuresSize % size_mpi, featuresSizeInt = featuresSize-modulo;

    //Rank0(master) issues
     if (rank_mpi == 0){
        data3 = (t_data*) calloc(featuresSize, sizeof(t_data)*datasize);
        recvcount = (int*) malloc(size_mpi*sizeof(int));
        recvdisp = (int*) malloc(size_mpi*sizeof(int));
        for (int z = 0; z < size_mpi; ++z){
            recvcount[z] = 1;
            recvdisp[z] = datasize;
        }
    }
    individual_buffer = (t_data*)calloc(featuresSize,datasize*sizeof(t_data));
    request = (MPI_Request*)malloc(datasize*featuresSizeInt/size_mpi*sizeof(MPI_Request));
    int cont = 0;
    for (i = 0; i < datasize; i++)
        for (j = rank_mpi; j < featuresSizeInt; j+=size_mpi){
            MPI_File_read_at(dataFile2, (i*featuresSize+j)*sizeof(t_data)+8,&individual_buffer[j*datasize+i], 1, MPI_BYTE, &status);
            MPI_Igatherv(&individual_buffer[j*datasize+i],1,MPI_BYTE, &data3[i+datasize*j],recvcount, recvdisp, MPI_BYTE,0,MPI_COMM_WORLD, &request[cont]);
            //MPI_File_seek(dataFile2,size_mpi,MPI_SEEK_CUR);
            cont++;
        }
    MPI_Waitall(datasize*featuresSizeInt/size_mpi,request,MPI_STATUSES_IGNORE);
    printf("data3 %c\n",data3[1]);
    printf("Fin\n");
   /* if (!rank_mpi){
        for (i = 0; i < datasize; i++)
            for (j = featuresSizeInt; j < featuresSize; j++)
                MPI_File_read_at(dataFile2,i*featuresSize+j+8,&data3[i+datasize*j],1,MPI_BYTE,&status);
    }*/
    //Liberacion
    /*free(request);free(individual_buffer);*/
    if (!rank_mpi){free(recvcount);free(recvdisp);}
    printf("Esto si termina\n");
}

void RawData::loadData() {
	uint i, j;
	t_data buffer[1];
	data = (t_data*) calloc(featuresSize, sizeof(t_data) * datasize);
	fseek(dataFile, 8, 0);
	for (i = 0; i < datasize; i++) {
		for (j = 0; j < featuresSize; j++) {
			fread(buffer, sizeof(t_data), 1, dataFile);
			data[j * datasize + i] = buffer[0];
		}
	}
}

/**
 * Calculates how many different values each feature has.
 */
void RawData::calculateVR() {
	uint i, j;
	t_data dataReaded;
	uint vr;
	valuesRange = (uint*) calloc(featuresSize, sizeof(uint));
	for (i = 0; i < featuresSize; i++) {
		vr = 0;
		for (j = 0; j < datasize; j++) {
			dataReaded = data[i * datasize + j];
			if (dataReaded > vr) {
				vr++;
			}
		}
		valuesRange[i] = vr + 1;
	}
}

uint RawData::getDataSize() {
	return datasize;
}

uint RawData::getFeaturesSize() {
	return featuresSize;
}

/**
 * Returns how much values has a features FROM 1 to VALUES;
 */
uint RawData::getValuesRange(uint index) {
	return valuesRange[index];
}

/**
 * Returns an array with the number of possible values for each feature.
 */
uint * RawData::getValuesRangeArray() {
	return this->valuesRange;
}

/**
 * Returns a vector containing a feature.
 */
t_feature RawData::getFeature(int index) {
	return data + index * datasize;
}
