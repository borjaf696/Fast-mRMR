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

/*Versión que usará la lectura única de un proceso o por parte de todos
 * como parametro de entrada.*/

#include "../../plibOp/JointProb.h"
#include "../../plibOp/MutualInfo.h"
#include "../../plibOp/utils.h"

#include <omp.h>
#include <vector>
#include <algorithm>
#include <limits>
#include <boost/bind.hpp>
#include <stdio.h>
#include <string.h>
#include <cmath>
#include <stdlib.h>
#include "mpi.h"


//Returns the index of the higher value in the classRelevances Vector different from classIndex
/*uint getMaxRelevance(vector<double> classRelevances, uint classIndex) {
	uint i = 0;
	uint newFeature = -1;
	double relevance = 0;
	for (i = 0; i < classRelevances.size(); ++i) {
		if (classRelevances[i] > relevance && i != classIndex) {
			relevance = classRelevances[i];
			newFeature = i;
		}
	}
	return newFeature;
}*/
uint getMaxRelevance(double* classRelevances, uint classIndex, uint rawDataSize) {
    uint i = 0;
    uint newFeature = -1;
    double relevance = 0;
    for (i = 0; i < rawDataSize ; ++i) {
        if (classRelevances[i] > relevance && i != classIndex) {
            relevance = classRelevances[i];
            newFeature = i;
        }
    }
    return newFeature;
}
void getMaxRelevance2(double* classRelevances, uint classIndex, uint rawDataSize, uint globalInteresValue, int mpi_rank, double *result) {
    uint i = 0;
    uint newFeature = -1;
    double relevance = 0;
    for (i = 0; i < rawDataSize ; ++i) {
        if (classRelevances[i] > relevance && i != classIndex) {
            relevance = classRelevances[i];
            newFeature = globalInteresValue*mpi_rank+i;
        }
    }
    result[0] = (double)newFeature;result[1] = relevance;
}
options parseOptions(int argc, char*argv[]) {
    options opts;
    opts.classIndex = 0;
    opts.selectedFeatures = 10;
    opts.file = (char *) malloc(20);
    strcpy(opts.file,"../data.mrmr");
    //opts.file = "../data.mrmr";

    if (argc > 1) {
        for (int i = 0; i < argc; ++i) {
            if (strcmp(argv[i], "-f") == 0) {
                opts.file = argv[i + 1];
            }
            if (strcmp(argv[i], "-a") == 0) {
                opts.selectedFeatures = atoi(argv[i + 1]) - 1;
            }
            if (strcmp(argv[i], "-c") == 0) {
                opts.classIndex = atoi(argv[i + 1]) - 1;
            }
            if (strcmp(argv[i], "-h") == 0) {
                printf(
                        "fast-mrmr:\nOptions:\n -f <inputfile>\t\tMRMR file generated using mrmrReader (default: data.mrmr).\n-c <classindex>\t\tIndicates the class index in the dataset (default: 0).\n-a <nfeatures>\t Indicates the number of features to select (default: 10).\n-h Prints this message");
                exit(0);
            }
        }
    }
    return opts;
}

int main(int argc, char* argv[]) {
    MPI_Init(&argc,&argv);
    int mpi_rank,mpi_size;

    //Inicializacion y obtencion de valores
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    //Presentaciones

    options opts;
    uint j = 0;
    uint newFeatureIndex = 0;
    uint lastFeatureIndex = 0;
    uint sizePerProcess = 0, sizePerProcessGlobal = 0;
    double mrmr = 0;
    double acum = 0;

    double * relevancesParciales;
    double * redundancesParciales;
    std::vector<int> selectedFeatures;
    int* visited;

    Timer tm;
    opts = parseOptions(argc, argv);
    tm.start();
    RawData rawData = RawData(opts.file);
    ProbTable prob = ProbTable(rawData);
    MutualInfo mutualInfo = MutualInfo(rawData, prob);
    visited = (int*) malloc(rawData.getFeaturesSize()*sizeof(int));

    for (uint z = 0;z < rawData.getFeaturesSize();z++)
        visited[z] = 0;


    //Variables de interes PorProceso
    sizePerProcessGlobal = (int) ceil(rawData.getFeaturesSize()/mpi_size);
    sizePerProcess = (int) ceil(rawData.getFeaturesSize()/mpi_size);
    if (mpi_rank == (mpi_size-1))
        sizePerProcess = rawData.getFeaturesSize()-mpi_rank*sizePerProcess;

    printf("Soy el proceso %d; tengo que ejecutar %d de %d\n",mpi_rank, sizePerProcess, rawData.getFeaturesSize());


    //Reserva y calculos parciales
    relevancesParciales = (double*) malloc(sizePerProcess*sizeof(double));
    redundancesParciales = (double*) malloc(sizePerProcess*sizeof(double));

    //Calculo de relevancias parciales a cada proceso.
    int iterator = 0;
    Timer timeRel;
    timeRel.start();

    int startPoint = mpi_rank*sizePerProcessGlobal;
    #pragma omp parallel for schedule(dynamic,10)
        for (uint l = startPoint; l < startPoint + sizePerProcess; ++l) {
            //printf("Thread %d\n",omp_get_thread_num());
            relevancesParciales[l-startPoint] = mutualInfo.get(opts.classIndex, l);
            redundancesParciales[l-startPoint] = 0;
        }
    /*if (!mpi_rank)
        printf("Valor del iterador %d\n",iterator);
    for (int i= 0; i< sizePerProcess; i++)
        printf("Relevancias %d Feature %d %f\n",mpi_rank,mpi_rank*sizePerProcess+i,relevancesParciales[i]);*/

    double timeMeasureRel = timeRel.stop();

    /*Para cada vector parcial se calcula su maximo, se envia al 0 y este discrimina cual es de interes.*/
    double * parInteres = (double*)malloc(2*sizeof(double)), * vectorInteres;
    getMaxRelevance2(relevancesParciales, opts.classIndex, sizePerProcess, sizePerProcessGlobal, mpi_rank, parInteres);

    if (!mpi_rank)
        vectorInteres = (double*)malloc(mpi_size*2* sizeof(double));

    MPI_Gather(parInteres, 2, MPI_DOUBLE,
               vectorInteres, 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if (!mpi_rank) {
        double max = vectorInteres[1];
        newFeatureIndex = vectorInteres[0];
        for (int i = 0; i < mpi_size * 2; i += 2)
            if (vectorInteres[i + 1] > max) {
                max = vectorInteres[i + 1];
                newFeatureIndex = (int) vectorInteres[i];
            }
    }
    MPI_Bcast(&newFeatureIndex,1, MPI_INT,
                0, MPI_COMM_WORLD);

    if (!mpi_rank)
        selectedFeatures.push_back(newFeatureIndex);
    lastFeatureIndex = newFeatureIndex;
    visited[lastFeatureIndex] = 1;

    //Liberamos la memoria reservada
    free(parInteres);
    if (!mpi_rank)
        free(vectorInteres);

    //Pruebas de tiempo sobre el original:
    Timer timeMrmr;
    timeMrmr.start();
    printf("Empieza lo bueno %d\n",mpi_rank);
    double * MrmrFeature, * MrmrFeaturesProcess;
    uint  etapas = 1;
    MrmrFeaturesProcess = (double*) malloc(2*sizeof(double));
    MrmrFeature = (double *) malloc(2 * mpi_size * sizeof(double));

    while (etapas < rawData.getFeaturesSize() - 1 /*-1 because class is discarded*/
           and etapas < opts.selectedFeatures) {
        //acum = infinito inicializacion simple y logica para meter el primer valor directamente.
        acum = -std::numeric_limits<double>::infinity();
        //Se divide en bloques de tamaño (en la medida de lo posible) homogeneos
        #pragma omp parallel for schedule(dynamic,10)
            for (j = startPoint; j < startPoint+sizePerProcess; ++j) {
                //printf("Thread %d\n",omp_get_thread_num());
                if (visited[j] == 0 && j != opts.classIndex) {
                    redundancesParciales[j-startPoint] += mutualInfo.get(lastFeatureIndex, j);
                    #pragma omp critical
                    {
                        mrmr = relevancesParciales[j-startPoint] - (redundancesParciales[j-startPoint] / (double) etapas);
                        if (mrmr > acum) {
                            acum = mrmr;
                            newFeatureIndex = j;
                        }
                    };
                }
            }
        MrmrFeaturesProcess[0] = acum;
        MrmrFeaturesProcess[1] = (double)newFeatureIndex;
        //Envio de la informacion
        MPI_Gather(MrmrFeaturesProcess, 2, MPI_DOUBLE,
            MrmrFeature, 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        if (!mpi_rank) {
            double max = -std::numeric_limits<double>::infinity();
            for (int i = 0; i < 2*mpi_size; i+=2) {
                if (MrmrFeature[i] > max) {
                    max = MrmrFeature[i];
                    newFeatureIndex = MrmrFeature[i + 1];
                }
            }
            selectedFeatures.push_back(newFeatureIndex);
            visited[newFeatureIndex] = 1;
            lastFeatureIndex = newFeatureIndex;
        }
        int featurePass = (int) lastFeatureIndex;
        MPI_Bcast(&featurePass, 1, MPI_INT, 0, MPI_COMM_WORLD);
        lastFeatureIndex = featurePass;
        visited[lastFeatureIndex] = 1;
        etapas++;
    }
    //Obtencion de tiempos estandar para Mrmr;Calculo de Relevancia;ProgramaTotal.
    double timeMeasureMrmr = timeMrmr.stop();
    printf("Tiempos medidos: tiempoRelevance: %f, tiempoMrmr %f\n",timeMeasureRel, timeMeasureMrmr);

    double timeTotal = tm.stop();
    printf("Tiempos totales medidos: %f\n",timeTotal);

    if (!mpi_rank)
        for (uint i=0; i < selectedFeatures.size(); i++)
            printf("%d,",selectedFeatures[i]);

    //Liberacion de memoria

    free(MrmrFeaturesProcess);
    free(MrmrFeature);

    free(visited);
    free(relevancesParciales);free(redundancesParciales);

    rawData.destroy();prob.destroy();
    printf("\n");
    MPI_Finalize();
    return (0);
}
