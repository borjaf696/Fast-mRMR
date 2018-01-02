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

#include "../../plibOp2/JointProb.h"
#include "../../plibOp2/MutualInfo.h"
#include "../../plibOp2/utils.h"

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

    double * relevancesParcialesProcess;
    double * redundancesParcialesProcess;
    std::vector<int> selectedFeatures;
    int* visited;

    Timer tm;
    opts = parseOptions(argc, argv);
    tm.start();
    //Esto solo el 0
    RawData rawData;
    if (!mpi_rank)
        rawData = RawData(opts.file);
    //Parte nuevisima
    RawData rawData2;
    Envio envios;
    int numElEnv = 3;
    int *send = (int*)malloc(numElEnv*sizeof(int));
    int *bigSend, *displ, *scount, *displFeatures, *scountFeatures;
    uint dataSizes = 0;
    rawData2 = RawData();
    if (!mpi_rank) {
        bigSend = (int *) malloc(mpi_size * numElEnv * sizeof(int));
        int sizeProcess = (int) ceil(rawData.getFeaturesSize()/mpi_size);
        for (int i = 0; i < (numElEnv * mpi_size); i += numElEnv) {
            if (i == ((mpi_size-1)*numElEnv))
                bigSend[i] = rawData.getFeaturesSize()-(i/numElEnv)*sizeProcess;
            else
                bigSend[i] = sizeProcess;
            bigSend[i+1] = rawData.getDataSize();
            bigSend[i+2] = sizeProcess;
        }
    }
    MPI_Scatter(bigSend, numElEnv, MPI_INT, send, numElEnv, MPI_INT, 0, MPI_COMM_WORLD);
    envios.featuresSize = send[0];
    envios.dataSize = send[1];
    sizePerProcessGlobal = send[2];
    rawData2.setValues(envios);
    dataSizes = rawData2.getDataSize();
    if (!mpi_rank) {
        dataSizes = rawData.getDataSize();
        displ = (int*)malloc(mpi_size*sizeof(int));displ[0] = 0;
        scount = (int*)malloc(mpi_size*sizeof(int));scount[0] = send[0];
        displFeatures = (int*)malloc(mpi_size*sizeof(int));displFeatures[0] = displ[0]*dataSizes;
        scountFeatures = (int*)malloc(mpi_size*sizeof(int));scountFeatures[0] = scount[0]*dataSizes;
        for (int i = 1; i < mpi_size; i++) {
            displ[i] = displ[i - 1] + scount[i - 1];
            scount[i] = bigSend[numElEnv * i];
            displFeatures[i] = displ[i]*dataSizes;
            scountFeatures[i] = scount[i]*dataSizes;
        }
 	printf("Get values %d\n",rawData.getDataArray()[100]);
    }
    uint * ranges = (uint*) malloc(send[0]*sizeof(uint));
    t_data * datos = (t_data*)malloc(send[0]*sizeof(t_data)*dataSizes);
    MPI_Scatterv(rawData.getValuesRangeArray(), scount, displ, MPI_INT, ranges, send[0], MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Scatterv(rawData.getDataArray(),scountFeatures,displFeatures,MPI_CHAR,rawData2.getDataArray(),send[0]*dataSizes, MPI_CHAR, 0, MPI_COMM_WORLD);
    rawData2.setValuesRangeArray(ranges);
//    rawData2.setValuesDataArray(datos);
    printf("Process %d Ranges %d Datos %d\n", mpi_rank, rawData2.getValuesRange(0),rawData2.getFeature(0)[0]);

    if (!mpi_rank) {
        for (int i = 0; i < mpi_size; i++)
            printf("Process %d RangeValue %d Datos %d\n", i,rawData.getValuesRange(displ[i]),rawData.getFeature(displ[i])[0]);
        free(bigSend);
        free(displ);free(displFeatures);
        free(scount);free(scountFeatures);
    }
    free(ranges);free(datos);
    free(send);
    if (!mpi_rank)
        rawData.destroy();

    //Datos movidos: 3 envios. -> Continuar mañana con un nuevo try y corroborar resultados

    ProbTable prob2 = ProbTable(rawData2);

    MutualInfo mutualInfoProcess = MutualInfo(rawData2, prob2);
    visited = (int*) malloc(rawData2.getFeaturesSize()*sizeof(int));

    for (uint z = 0;z < rawData2.getFeaturesSize();z++)
        visited[z] = 0;

    //Reserva y calculos parciales
    int startPoint = mpi_rank*sizePerProcessGlobal;
    relevancesParcialesProcess = (double*) malloc(rawData2.getFeaturesSize()*sizeof(double));
    redundancesParcialesProcess = (double*) malloc(rawData2.getFeaturesSize()*sizeof(double));


    //Calculo de relevancias parciales a cada proceso -> envio del vector de valores de la feature de interes.
    t_feature vectInteres = (t_feature)malloc((dataSizes+1)*sizeof(t_data));
    double * probVector;
    int claseStudy = opts.classIndex;
    int proceso = (int)floor((float)claseStudy/sizePerProcessGlobal);
    if (proceso == mpi_size)
        proceso--;
    claseStudy -=startPoint;
    if (mpi_rank == proceso) {
        probVector = (double*) malloc(rawData2.getValuesRange(claseStudy)*sizeof(double));
        for (uint i = 0; i < dataSizes; i++)
            vectInteres[i] = rawData2.getFeature(claseStudy)[i];
        vectInteres[dataSizes] = rawData2.getValuesRange(claseStudy);
        for (uint i = 0; i < rawData2.getValuesRange(claseStudy); i++)
            probVector[i] = prob2.getProbability(claseStudy,i);
    }
    MPI_Bcast(vectInteres, dataSizes+1, MPI_CHAR, proceso, MPI_COMM_WORLD);
    if (mpi_rank != proceso)
        probVector = (double *) malloc(vectInteres[dataSizes] * sizeof(double));
    MPI_Bcast(probVector, vectInteres[dataSizes], MPI_DOUBLE, proceso, MPI_COMM_WORLD);
    //Envio del vectorProbabilidades, rangesFeature, ValuesFeature

    Timer timeRel;
    timeRel.start();

    #pragma omp parallel for
    	for (uint l = startPoint; l < startPoint + rawData2.getFeaturesSize(); ++l) {
        	relevancesParcialesProcess[l-startPoint] = mutualInfoProcess.get(vectInteres,l-startPoint,probVector);
	        redundancesParcialesProcess[l-startPoint] = 0;
    	}
    printf("Termine %d\n",mpi_rank);
    //Relevancias calculadas -> easy peasy

    double timeMeasureRel = timeRel.stop();

    /*Para cada vector parcial se calcula su maximo, se envia al 0 y este discrimina cual es de interes.*/
    double * parInteres = (double*)malloc(2*sizeof(double)), * vectorInteres;
    getMaxRelevance2(relevancesParcialesProcess, opts.classIndex, rawData2.getFeaturesSize(), sizePerProcessGlobal, mpi_rank, parInteres);

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
    if ((lastFeatureIndex > startPoint) & (lastFeatureIndex < (startPoint+rawData2.getFeaturesSize())))
        visited[lastFeatureIndex - startPoint] = 1;

    //Enviamos cada caso que corresponde
    free(probVector);
    claseStudy = lastFeatureIndex;
    proceso = (int)floor((float)claseStudy/sizePerProcessGlobal);
    if (proceso == mpi_size)
        proceso --;
    claseStudy -=startPoint;
    if (mpi_rank == proceso) {
        probVector = (double*) malloc(rawData2.getValuesRange(claseStudy)*sizeof(double));
        for (uint i = 0; i < dataSizes; i++)
            vectInteres[i] = rawData2.getFeature(claseStudy)[i];
        vectInteres[dataSizes] = rawData2.getValuesRange(claseStudy);
        for (uint i = 0; i < rawData2.getValuesRange(claseStudy); i++)
            probVector[i] = prob2.getProbability(claseStudy,i);
    }

    MPI_Bcast(vectInteres, dataSizes+1, MPI_CHAR, proceso, MPI_COMM_WORLD);
    if (mpi_rank != proceso)
        probVector = (double *) malloc(vectInteres[dataSizes] * sizeof(double));
    MPI_Bcast(probVector, vectInteres[dataSizes], MPI_DOUBLE, proceso, MPI_COMM_WORLD);


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
    int threadsMaximos = omp_get_max_threads();
    double* acumThread = (double*) malloc(threadsMaximos*sizeof(double));
    int* numFeature = (int*)malloc(threadsMaximos*sizeof(int));

    while (etapas < opts.selectedFeatures) {
        //acum = infinito inicializacion simple y logica para meter el primer valor directamente.
        acum = -std::numeric_limits<double>::infinity();
        for (int z = 0; z < threadsMaximos; ++z)
            acumThread[z] = -std::numeric_limits<double>::infinity();
        //Se divide en bloques de tamaño (en la medida de lo posible) homogeneos
        #pragma omp parallel for private(mrmr) schedule(dynamic,20)
        	for (j = startPoint; j < startPoint+rawData2.getFeaturesSize(); ++j) {
	            if (visited[j-startPoint] == 0 && j != opts.classIndex) {
        	        redundancesParcialesProcess[j-startPoint]+=mutualInfoProcess.get(vectInteres,j-startPoint,probVector);
                	mrmr = relevancesParcialesProcess[j-startPoint] - (redundancesParcialesProcess[j-startPoint] / (double) etapas);
                	if (mrmr > acumThread[omp_get_thread_num()]) {
                	    acumThread[omp_get_thread_num()] = mrmr;
                	    numFeature[omp_get_thread_num()] = j;
                	}
            	    }
        	}
        for( int z = 0; z < threadsMaximos; ++z){
            if (acumThread[z] > acum){
                acum = acumThread[z];
                newFeatureIndex = numFeature[z];
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
            lastFeatureIndex = newFeatureIndex;
        }
        int featurePass = (int) lastFeatureIndex;
        MPI_Bcast(&featurePass, 1, MPI_INT, 0, MPI_COMM_WORLD);
        lastFeatureIndex = featurePass;
        if ((lastFeatureIndex > startPoint) & (lastFeatureIndex < (startPoint+rawData2.getFeaturesSize())))
            visited[lastFeatureIndex - startPoint] = 1;

        //Se envia todo nuevamente
        free(probVector);
        claseStudy = lastFeatureIndex;
        proceso = (int)floor((float)claseStudy/sizePerProcessGlobal);
        if (proceso == mpi_size)
            proceso --;
        claseStudy -=startPoint;
        if (mpi_rank == proceso) {
            probVector = (double*) malloc(rawData2.getValuesRange(claseStudy)*sizeof(double));
            for (uint i = 0; i < dataSizes; i++)
                vectInteres[i] = rawData2.getFeature(claseStudy)[i];
            vectInteres[dataSizes] = rawData2.getValuesRange(claseStudy);
            for (uint i = 0; i < rawData2.getValuesRange(claseStudy); i++)
                probVector[i] = prob2.getProbability(claseStudy,i);
        }
        MPI_Bcast(vectInteres, dataSizes+1, MPI_CHAR, proceso, MPI_COMM_WORLD);
        if (mpi_rank != proceso)
            probVector = (double *) malloc(vectInteres[dataSizes] * sizeof(double));
        MPI_Bcast(probVector, vectInteres[dataSizes], MPI_DOUBLE, proceso, MPI_COMM_WORLD);

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
    free(acumThread);
    free(numFeature);
    free(probVector);
    free(vectInteres);

    free(MrmrFeaturesProcess);
    free(MrmrFeature);

    free(visited);
    free(relevancesParcialesProcess);free(redundancesParcialesProcess);

    rawData2.destroy();prob2.destroy();
    printf("\n");
    MPI_Finalize();
    return (0);
}
