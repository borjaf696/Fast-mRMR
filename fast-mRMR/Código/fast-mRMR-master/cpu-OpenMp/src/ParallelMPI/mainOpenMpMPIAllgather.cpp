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

    #pragma omp parallel for schedule(dynamic,10)
    for (uint l = mpi_rank*sizePerProcessGlobal;l < displ[mpi_rank]+sizePerProcess; ++l){
        relevancesParciales[iterator] = mutualInfo.get(opts.classIndex, l);
        redundancesParciales[iterator] = 0;
        iterator++;
    }
    //Prueba del gather
    MPI_Allgatherv(relevancesParciales, sizePerProcess, MPI_DOUBLE,
                   relevances, recvcounts, displ, MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Allgatherv(redundancesParciales, sizePerProcess, MPI_DOUBLE,
                    redundances, recvcounts, displ, MPI_DOUBLE, MPI_COMM_WORLD);

    double timeMeasureRel = timeRel.stop();

    // Max relevance feature is added because no redundancy is possible.
    //newFeatureIndex = getMaxRelevance(relevances, opts.classIndex);
    newFeatureIndex = getMaxRelevance(relevances, opts.classIndex, rawData.getFeaturesSize());
    selectedFeatures.push_back(newFeatureIndex);
    visited[newFeatureIndex] = 1;
    lastFeatureIndex = newFeatureIndex;

    // Modo de salvar cada caracteristica seleccionada.
    std::cout << newFeatureIndex << ",";

    //Pruebas de tiempo sobre el original:
    Timer timeMrmr;
    timeMrmr.start();
    //MRMR Original
    while (selectedFeatures.size() < rawData.getFeaturesSize() - 1 /*-1 because class is discarded*/ and selectedFeatures.size() < opts.selectedFeatures) {
        //acum = infinito inicializacion simple y logica para meter el primer valor directamente.
        acum = -std::numeric_limits<double>::infinity();
        #pragma omp parallel for schedule(dynamic,10)
        for (j = 0; j < rawData.getFeaturesSize(); ++j) {
            //If feature not in selected selectedFeatures -> se podrian usar conjuntos o eliminar de verdad la caracteristica -> ¿supone un ahorro sustancial?
            if (visited[j] == 0 && j != opts.classIndex) {
                redundances[j] += mutualInfo.get(lastFeatureIndex, j);
                mrmr = relevances[j] - (redundances[j] / selectedFeatures.size());
                if (mrmr > acum) {
                    acum = mrmr;
                    newFeatureIndex = j;
                }
            }
        }
        //Last feature doesn't prints comma.
        if ( (selectedFeatures.size() == (opts.selectedFeatures - 1)) or (selectedFeatures.size() == (rawData.getFeaturesSize() -2)) ){
            std::cout << newFeatureIndex << "\n";
        }else{
            std::cout << newFeatureIndex << "," ;
        }
        selectedFeatures.push_back(newFeatureIndex);
        visited[newFeatureIndex] = 1;
        lastFeatureIndex = newFeatureIndex;
    }
    //Obtencion de tiempos estandar para Mrmr;Calculo de Relevancia;ProgramaTotal.
    double timeMeasureMrmr = timeMrmr.stop();
    printf("Tiempos medidos: tiempoRelevance: %f, tiempoMrmr %f\n",timeMeasureRel, timeMeasureMrmr);

    double timeTotal = tm.stop();
    printf("Tiempos totales medidos: %f\n",timeTotal);

    //Liberacion de memoria
    free(visited);
    free(displ);free(recvcounts);
    free(relevancesParciales);free(redundancesParciales);
    free(relevances);free(redundances);

    rawData.destroy();prob.destroy();
    printf("\n");
    MPI_Finalize();
    return (0);
}
