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

#include "RawData.h"
#include "omp.h"
#include <string.h>
#include <cmath>
#include "mpi.h"


/**
 * Constructor that creates a rawData object.
 *
 * @param data_table this is a matrix of bytes containing the translated csv data.
 * @param ds the number of data samples.
 * @param fs the number of features that each sample has.
 */
RawData::RawData(char * filename, int mpisize, int mpinum) {
	dataFile = fopen(filename, "rb"); //fopen devuelve un puntero al fichero filename
	MPI_File_open(MPI_COMM_WORLD, filename,MPI_MODE_RDONLY,MPI_INFO_NULL, &fh);
	mpiSize = mpisize;
	mpiNum = mpinum;
	calculateDSandFS();
	loadData();
    featureStudio = (t_feature) malloc(datasize*sizeof(t_data));
	calculateVR();
}

RawData::~RawData() {

}

size_t result;
/**
 *
 */
void RawData::destroy() {
	//free(valuesRange);
	//free(data);
	free(data2);
	free(valuesRange2);
}

/**Calculates
 *	DataSize: Number of patterns or samples
 *  FeaturesSize: Number of features
 */
/*Metodo que calcula el numero de Features y el numero de samples. (Para esto lee los 2 primeros valores que se encuentran
en el fichero de datos, .mrmr que se han colocado al principio FS y el DS -> luego hay que saltarlos de ahi el avance de
8 posiciones)*/
void RawData::calculateDSandFS() {
	uint featuresSizeBuffer[1];
	uint datasizeBuffer[1];
	//fread avanza el puntero en el fichero.
	result = fread(datasizeBuffer, sizeof(uint), 1, dataFile);
	if (result != 1)
		{fputs ("Reading error bufferDatafile\n",stderr); exit (3);}
	result = fread(featuresSizeBuffer, sizeof(uint), 1, dataFile);
	if (result != 1)
		{fputs ("Reading error featuresDatafile\n",stderr); exit (3);}
	datasize = datasizeBuffer[0];
	featuresSize = featuresSizeBuffer[0];
	featuresSize2 = (int) ceil(featuresSizeBuffer[0]/ mpiSize);
	uint sizePerProcess2 = featuresSize2, suma = 0, diff = featuresSize-featuresSize2 * mpiSize ;
	featuresSize2 = sizePerProcess2;
	inicio = featuresSize2*mpiNum;
	if ( diff > mpiNum) {
		featuresSize2++;
		suma = (mpiNum > 0) ? (mpiNum):0;
	}else
		suma = diff;
	inicio = inicio+suma;
	/*if (mpiNum == (mpiSize - 1))
		featuresSize2 = featuresSize-mpiNum*featuresSize2;*/

}

/*Rellena el atributo data con la informacion del fichero incluido en dataFile. Usa el tipo t_data -> que no es mas que
que un caracter. Se rellena data la info con ese valor para no limitarse a valores numericos lo que esta realmente bien.
dataFile parece que guarda la informacion contigua relativa a las features -> es decir el valor 0 leido contiene el valor
de la feature 0 para el sample 0, el valor 1 el valor 0 para la feature 1 para el sample 0. Es decir, parece que esta
dispuesto de modo FeaturesXData -> la traspuesta de la que se pide.*/
void RawData::loadData() {
	uint i, j = 0;
	int rank;
	long int index = 0;
	MPI_Status status;
	MPI_Request request;
	t_data buffer[1];
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	data2 = (t_data*) malloc(featuresSize2*sizeof(t_data) * datasize);
	if (data2 == NULL){
		printf("No such allocating allowed\n");
		exit(3);
	}
	//Avanza 8 bits en el fichero de datos, para saltar metaData
	fseek(dataFile,inicio+8,0);
	for (i = 0; i < datasize; i++) {
		for (j = 0; (j < featuresSize2); j++) {
			//Lee un bit de informacion
			result = fread(buffer, sizeof(t_data), 1, dataFile);
                        if (result != 1)
                                {fputs ("Reading error loadFile\n",stderr); exit (3);}
                        //Se va rellenando data -> para 1 mismo sample los distintos valores que toma para cada feature
				index =(long int) j*datasize+i;
                                data2[index] = buffer[0];
		}
		fseek(dataFile,featuresSize-featuresSize2,SEEK_CUR);
	}
	/*t_data buffer2[featuresSize2];
	MPI_File_seek( fh, inicio+8, MPI_SEEK_SET);
	printf("Inicio en %d leo %d\n",inicio,featuresSize2);
	for (i = 0; i < datasize; i++) {
		/*for (j = 0; (j < featuresSize2); j++) {
			//Lee un bit de informacion
			MPI_File_read( fh, buffer, 1, MPI_CHAR, &status );
			if (result != 1)
			{fputs ("Reading error loadFile\n",stderr); printf("Proceso %d\n",mpiNum);exit (3);}
			//Se va rellenando data -> para 1 mismo sample los distintos valores que toma para cada feature
			data2[j * datasize + i] = buffer[0];
		}
		MPI_File_iread(fh,buffer2,featuresSize2,MPI_CHAR,&request);
		for (j = 0; j < featuresSize2;j++)
			data2[j*datasize+i] = buffer2[j];
		MPI_File_seek(fh,featuresSize-featuresSize2,MPI_SEEK_CUR);
	}
	MPI_File_close(&fh);*/
}

/**
 * Calculates how many different values each feature has.
 */
/*Si que cuenta los valores diferentes -> caso extremo 1 solo valor 16, 100 veces. Vr llegara como maximo
a 16 -> pues vr == valor. Pero, y si lees primero 100 y luego veinte veces el valor 16 -> tienes 2 valores
y vas a crear un histograma con 16 valores. Hacer mi propia version*/
void RawData::calculateVR() {
	uint i, j;
	t_data dataReaded;
	long int index;
	uint vr;
	valuesRange2 = (uint*) calloc(featuresSize2, sizeof(uint));
	#pragma omp parallel for
    	for (i = 0; i< featuresSize2; i++){
		vr = 0;
		for (j = 0; j < datasize; j++){
			index = (long int)i*datasize+j;
			dataReaded = data2[index];
			if (dataReaded > vr)
				vr++;
		}
		valuesRange2[i] = vr + 1;
	}
}

uint RawData::getDataSize() {
	return datasize;
}

uint RawData::getFeaturesSize() {
	return featuresSize;
}

uint RawData::getFeaturesSizeReal() {
	return featuresSize2;
}

/**
 * Returns how much values has a features FROM 1 to VALUES;
 */
uint RawData::getValuesRange(uint index) {
	//return valuesRange[index];
	return valuesRange2[index];
}

/**
 * Returns an array with the number of possible values for each feature.
 */
uint * RawData::getValuesRangeArray() {
	return this->valuesRange2;
}

/**
 * Returns a vector containing a feature. Realmente devuelve un puntero
 */
t_feature RawData::getFeature(int index) {
	//return data + index * datasize;
	long int index2 = (long int) index*datasize;
	return data2+ index2;
	/*Como se apuntaba anteriormente cada feature contiene su informacion de manera contigua en memoria.*/
}

t_feature RawData::getFeatureStudio(t_feature data) {
    this->featureStudio = data;
}

t_feature RawData::getFeatureStudio() {
    return this->featureStudio;
}

void RawData::setValuesFeatureStudio(uint value){
    this->valuesRangeFeatureStudio = value;
}

uint RawData::getValuesFeatureStudio() {
    return this->valuesRangeFeatureStudio;
}
