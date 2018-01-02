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

#include "../../plibOp/JointProb.h"
#include "../../plibOp/MutualInfo.h"
#include "../../plibOp/utils.h"

#include "omp.h"
#include <vector>
#include <algorithm>
#include <limits>
#include <boost/bind.hpp>
#include <stdio.h>
#include <string.h>
#include "omp.h"
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
    /*vector<double> relevances;
    vector<double> redundances;*/
    double * relevances, * relevancesParciales;
    double * redundances, * redundancesParciales;
    std::vector<int> selectedFeatures;
    int* visited, * recvcounts, * displ;

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

    //Reserva de memoria global por proceso
    relevances = (double*) malloc (rawData.getFeaturesSize()*sizeof(double));
    redundances = (double*) malloc (rawData.getFeaturesSize()*sizeof(double));

    //Reserva y calculos parciales
    relevancesParciales = (double*) malloc(sizePerProcess*sizeof(double));
    redundancesParciales = (double*) malloc(sizePerProcess*sizeof(double));

    recvcounts = (int*) malloc (mpi_size*sizeof(int));
    displ = (int*) malloc(mpi_size*sizeof(int));

    //Calculo de desplazamientos sobre cada proceso.
    for (int l = 0; l < mpi_size; ++l) {
        recvcounts[l] = (int) ceil(rawData.getFeaturesSize() / mpi_size);
        displ[l] = l * recvcounts[l];
        if (l == (mpi_size - 1))
            recvcounts[l] = rawData.getFeaturesSize() - l * recvcounts[l];
    }

    //Calculo de relevancias parciales a cada proceso.
    int iterator = 0;
    Timer timeRel;
    timeRel.start();

    #pragma omp parallel
    for (uint l = mpi_rank*sizePerProcessGlobal;l < displ[mpi_rank]+sizePerProcess; ++l){
        relevancesParciales[iterator] = mutualInfo.get(opts.classIndex, l);
        redundancesParciales[iterator] = 0;
        iterator++;
    }
    //Prueba del gather -> revisar todos los valores de relevance igual va a ser lo mejor.
    MPI_Allgatherv(relevancesParciales, sizePerProcess, MPI_DOUBLE,
                   relevances, recvcounts, displ, MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Allgatherv(redundancesParciales, sizePerProcess, MPI_DOUBLE,
                    redundances, recvcounts, displ, MPI_DOUBLE, MPI_COMM_WORLD);

    double timeMeasureRel = timeRel.stop();

    // Max relevance feature is added because no redundancy is possible.
    //newFeatureIndex = getMaxRelevance(relevances, opts.classIndex);
    /*Por lo de ahora vamos a dejar que todos lleven la cuenta de todo*/
    newFeatureIndex = getMaxRelevance(relevances, opts.classIndex, rawData.getFeaturesSize());
    if (!mpi_rank)
        selectedFeatures.push_back(newFeatureIndex);
    lastFeatureIndex = newFeatureIndex;
    visited[lastFeatureIndex] = 1;

    //Pruebas de tiempo sobre el original:
    Timer timeMrmr;
    timeMrmr.start();
    printf("Empieza lo bueno %d\n",mpi_rank);
    double * MrmrFeature, * MrmrFeaturesProcess;
    uint  etapas = 1;
    MrmrFeaturesProcess = (double*) malloc(2*sizeof(double));
    MrmrFeature = (double *) malloc(2 * mpi_size * sizeof(double));
    //MRMR Original
    while (etapas < rawData.getFeaturesSize() - 1 /*-1 because class is discarded*/
           and etapas < opts.selectedFeatures) {
        //acum = infinito inicializacion simple y logica para meter el primer valor directamente.
        acum = -std::numeric_limits<double>::infinity();
        //Se divide en bloques de tamaño (en la medida de lo posible) homogeneos
        for (j = mpi_rank*sizePerProcessGlobal;j < displ[mpi_rank]+sizePerProcess; ++j){
            if (visited[j] == 0 && j != opts.classIndex) {
                redundances[j] += mutualInfo.get(lastFeatureIndex, j);
                mrmr = relevances[j] - (redundances[j] /(double) etapas);
                if (mrmr > acum) {
                    acum = mrmr;
                    newFeatureIndex = j;
                }
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
    free(displ);free(recvcounts);
    free(relevancesParciales);free(redundancesParciales);
    free(relevances);free(redundances);

    rawData.destroy();prob.destroy();
    printf("\n");
    MPI_Finalize();
    return (0);
}
