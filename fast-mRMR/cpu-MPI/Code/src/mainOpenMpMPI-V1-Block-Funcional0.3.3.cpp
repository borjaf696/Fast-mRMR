
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
uint getMaxFeature(double MrmrFeature[],int mpi_size){
    uint newFeatureIndex = -5;
    double max = -std::numeric_limits<double>::infinity();
    for (int i = 0; i < 2*mpi_size; i+=2)
        if (MrmrFeature[i] > max) {
            max = MrmrFeature[i];
            newFeatureIndex = MrmrFeature[i + 1];
        }
    return newFeatureIndex;
}

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
void getMaxRelevance2(double* classRelevances, uint classIndex, uint rawDataSize, uint globalInteresValue, int mpi_rank, double *result, int suma) {
    uint i = 0;
    uint newFeature = -1;
    double relevance = -1242424234;
    for (i = 0; i < rawDataSize ; ++i) {
        if (classRelevances[i] > relevance && i != classIndex) {
            relevance = classRelevances[i];
            newFeature = globalInteresValue*mpi_rank+i+suma;
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

    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    options opts;
    uint j = 0;
    uint newFeatureIndex = 0;
    uint lastFeatureIndex = 0;
    uint sizePerProcess = 0, sizePerProcessGlobal = 0;
    double mrmr = 0;
    double acum = 0;

    double * relevancesParciales, *probVector;
    double * redundancesParciales;
    std::vector<int> selectedFeatures;
    int* visited;

    Timer tm;
    opts = parseOptions(argc, argv);
    tm.start();
    RawData rawData = RawData(opts.file, mpi_size,mpi_rank);
    ProbTable prob = ProbTable(rawData);
    MutualInfo mutualInfo = MutualInfo(rawData, prob);
    visited = (int*) malloc(rawData.getFeaturesSizeReal()*sizeof(int));
    sizePerProcessGlobal = (int) ceil(rawData.getFeaturesSize()/mpi_size);
    sizePerProcess = (int) ceil(rawData.getFeaturesSize()/mpi_size);
    uint suma = 0, diff =rawData.getFeaturesSize()-sizePerProcessGlobal * mpi_size ;
    if ( diff > mpi_rank) {
        sizePerProcess++;
        suma = (mpi_rank > 0) ? (mpi_rank):0; ;
    }else
        suma = rawData.getFeaturesSize()-sizePerProcessGlobal*mpi_size;
 
    int claseStudy = lastFeatureIndex;
    int proceso = (int)floor((float)claseStudy/sizePerProcessGlobal);
    int sumaInteres = (proceso > diff) ? diff:proceso;
    while (proceso*sizePerProcessGlobal+sumaInteres > claseStudy) {
        proceso--;
        sumaInteres = (proceso > diff) ? diff:proceso;
    }
    int dataSizes = rawData.getDataSize();
    claseStudy -=sizePerProcessGlobal*mpi_rank;
    t_feature vectInteres2 = (t_feature)malloc((dataSizes+2+sizeof(double)*256)*sizeof(t_data));
    if (proceso > mpi_size)
        proceso =mpi_size - 1;
    if (mpi_rank == proceso) {
        probVector = (double*) malloc(rawData.getValuesRange(claseStudy)*sizeof(double));
        for (uint i = 0; i < dataSizes; i++)
            vectInteres2[i] = rawData.getFeature(claseStudy)[i];
        vectInteres2[dataSizes] = rawData.getValuesRange(claseStudy) > 255 ? 255:rawData.getValuesRange(claseStudy);
        vectInteres2[dataSizes+1] = rawData.getValuesRange(claseStudy) > 255? 1:0;
        int pos = dataSizes+2;
        for (uint i = 0; i < rawData.getValuesRange(claseStudy); i++) {
            probVector[i] = prob.getProbability(claseStudy, i);
            double a =prob.getProbability(claseStudy,i);
            unsigned char *p = (unsigned char*)&a;
            for (int j = 0; j != sizeof(double); ++j)
                vectInteres2[pos+j] = p[j];
            pos+=8;
        }
    }
    MPI_Bcast(vectInteres2,dataSizes+2+sizeof(double)*256,MPI_CHAR, proceso, MPI_COMM_WORLD);
    int tamanho = vectInteres2[dataSizes]+vectInteres2[dataSizes+1];
    rawData.getFeatureStudio(vectInteres2);
    rawData.setValuesFeatureStudio(tamanho);
    mutualInfo.update(rawData);
    if (mpi_rank != proceso)
        probVector = (double *) malloc(tamanho * sizeof(double));
    int count = 0;
    for (int j = dataSizes+2; j < (dataSizes+2+sizeof(double)*tamanho); j+=8) {
        double out;
        memcpy(&out, &(vectInteres2[j]), sizeof(double));
        probVector[count] = out;
        count++;
    }

    for (uint z = 0;z < rawData.getFeaturesSizeReal();z++)
        visited[z] = 0;

    relevancesParciales = (double*) malloc(sizePerProcess*sizeof(double));
    redundancesParciales = (double*) malloc(sizePerProcess*sizeof(double));

 
    int iterator = 0;
    Timer timeRel;
    timeRel.start();

    int startPoint = mpi_rank*sizePerProcessGlobal+suma;
    #pragma omp parallel for schedule(dynamic)
        for (uint l = startPoint; l < startPoint + sizePerProcess; ++l) {
            relevancesParciales[l-startPoint] = mutualInfo.get(l-startPoint, probVector);
            redundancesParciales[l-startPoint] = 0;
        }
    double timeMeasureRel = timeRel.stop();
 
    double * parInteres = (double*)malloc(2*sizeof(double)), * vectorInteres;
    getMaxRelevance2(relevancesParciales, opts.classIndex, sizePerProcess, sizePerProcessGlobal, mpi_rank, parInteres,suma);

    vectorInteres = (double*)malloc(mpi_size*2* sizeof(double));

    MPI_Allgather(parInteres, 2, MPI_DOUBLE,
               vectorInteres, 2, MPI_DOUBLE, MPI_COMM_WORLD);
    double max = vectorInteres[1];
    newFeatureIndex = vectorInteres[0];
    for (int i = 0; i < mpi_size * 2; i += 2)
        if (vectorInteres[i + 1] > max) {
            max = vectorInteres[i + 1];
            newFeatureIndex = (int) vectorInteres[i];
        }

    if (!mpi_rank)
        selectedFeatures.push_back(newFeatureIndex);
    lastFeatureIndex = newFeatureIndex;
    if ((lastFeatureIndex > startPoint) & (lastFeatureIndex < (startPoint+sizePerProcess)))
        visited[lastFeatureIndex-startPoint] = 1;

 
    free(probVector);
    claseStudy = lastFeatureIndex;
    proceso = (int)floor((float)claseStudy/sizePerProcessGlobal);
    sumaInteres = (proceso > diff) ? diff:proceso;
    while (proceso*sizePerProcessGlobal+sumaInteres > claseStudy) {
        proceso--;
        sumaInteres = (proceso > diff) ? diff:proceso;
    }
    if (proceso > mpi_size)
        proceso = mpi_size - 1;
    claseStudy -=startPoint;
    if (mpi_rank == proceso) {
        probVector = (double*) malloc(rawData.getValuesRange(claseStudy)*sizeof(double));
        for (uint i = 0; i < dataSizes; i++)
            vectInteres2[i] = rawData.getFeature(claseStudy)[i];
        vectInteres2[dataSizes] = rawData.getValuesRange(claseStudy) > 255 ? 255:rawData.getValuesRange(claseStudy);
        vectInteres2[dataSizes+1] = rawData.getValuesRange(claseStudy) > 255? 1:0;
        int pos = dataSizes+2;
        for (uint i = 0; i < rawData.getValuesRange(claseStudy); i++) {
            probVector[i] = prob.getProbability(claseStudy, i);
            double a =prob.getProbability(claseStudy,i);
            unsigned char *p = (unsigned char*)&a;
            for (int j = 0; j != sizeof(double); ++j)
                vectInteres2[pos+j] = p[j];
            pos+=8;
        }
    }
    MPI_Bcast(vectInteres2,dataSizes+2+sizeof(double)*256,MPI_CHAR, proceso, MPI_COMM_WORLD);
    tamanho = vectInteres2[dataSizes]+vectInteres2[dataSizes+1];
    rawData.getFeatureStudio(vectInteres2);
    rawData.setValuesFeatureStudio(tamanho);
    mutualInfo.update(rawData);
    if (mpi_rank != proceso)
        probVector = (double *) malloc(tamanho * sizeof(double));
    count = 0;
    for (int j = dataSizes+2; j < (dataSizes+2+sizeof(double)*tamanho); j+=8) {
        double out;
        memcpy(&out, &(vectInteres2[j]), sizeof(double));
        probVector[count] = out;
        count++;
    }
 
    free(parInteres);
    if (!mpi_rank)
        free(vectorInteres);

 
    Timer timeMrmr;
    timeMrmr.start();
    double * MrmrFeature, * MrmrFeaturesProcess;
    uint  etapas = 1;
    MrmrFeaturesProcess = (double*) malloc(2*sizeof(double));
    MrmrFeature = (double *) malloc(2 * mpi_size * sizeof(double));
    int threadsMaximos = omp_get_max_threads();
    double* acumThread = (double*) malloc(threadsMaximos*sizeof(double));
    int* numFeature = (int*)malloc(threadsMaximos*sizeof(int));

    while (etapas < rawData.getFeaturesSize() - 1 /*-1 because class is discarded*/
           and etapas < opts.selectedFeatures) {
 
        acum = -std::numeric_limits<double>::infinity();
	mrmr = -std::numeric_limits<double>::infinity();
        for (int z = 0; z < threadsMaximos; ++z)
            acumThread[z] = -std::numeric_limits<double>::infinity();
 
        #pragma omp parallel for schedule(dynamic) private(mrmr)
            for (j = startPoint; j < startPoint+sizePerProcess; ++j) {
                if (visited[j-startPoint] == 0 && j != opts.classIndex) {
                    //redundancesParciales[j-startPoint] += mutualInfo.get(lastFeatureIndex, j);
                    //redundancesParciales[j-startPoint] += mutualInfo.get(vectInteres2,j-startPoint,probVector);
                    redundancesParciales[j-startPoint] += mutualInfo.get(j-startPoint,probVector);
                    mrmr = relevancesParciales[j-startPoint] - (redundancesParciales[j-startPoint] / (double) etapas);
                    if (acumThread[omp_get_thread_num()] < mrmr) {
                        acumThread[omp_get_thread_num()] = mrmr;
                        numFeature[omp_get_thread_num()] = j;
                    }
                }
            }
        for( int z = 0; z < threadsMaximos; ++z)
            if (acumThread[z] > acum) {
                acum = acumThread[z];
                newFeatureIndex = numFeature[z];
            }
        MrmrFeaturesProcess[0] = acum;
        MrmrFeaturesProcess[1] = (double)newFeatureIndex;
 
        MPI_Allgather(MrmrFeaturesProcess, 2, MPI_DOUBLE,
                   MrmrFeature, 2, MPI_DOUBLE, MPI_COMM_WORLD);
        double max = -std::numeric_limits<double>::infinity();
        for (int i = 0; i < 2*mpi_size; i+=2) {
            if (MrmrFeature[i] > max) {
                max = MrmrFeature[i];
                newFeatureIndex = MrmrFeature[i + 1];
            }
        }
        selectedFeatures.push_back(newFeatureIndex);
        lastFeatureIndex = newFeatureIndex;
        if ((lastFeatureIndex >= startPoint) & (lastFeatureIndex < (startPoint+sizePerProcess)))
            visited[lastFeatureIndex-startPoint] = 1;

 
        free(probVector);
        claseStudy = lastFeatureIndex;
        proceso = (int)floor((float)claseStudy/sizePerProcessGlobal);
        sumaInteres = (proceso > diff) ? diff:proceso;
        while (proceso*sizePerProcessGlobal+sumaInteres > claseStudy) {
            proceso--;
            sumaInteres = (proceso > diff) ? diff:proceso;
        }
        if (proceso > mpi_size)
            proceso = mpi_size - 1;
        claseStudy -=startPoint;
        if (mpi_rank == proceso) {
            probVector = (double*) malloc(rawData.getValuesRange(claseStudy)*sizeof(double));
            for (uint i = 0; i < dataSizes; i++)
                vectInteres2[i] = rawData.getFeature(claseStudy)[i];
            vectInteres2[dataSizes] = rawData.getValuesRange(claseStudy) > 255 ? 255:rawData.getValuesRange(claseStudy);
            vectInteres2[dataSizes+1] = rawData.getValuesRange(claseStudy) > 255? 1:0;
            int pos = dataSizes+2;
            for (uint i = 0; i < rawData.getValuesRange(claseStudy); i++) {
                probVector[i] = prob.getProbability(claseStudy, i);
                double a =prob.getProbability(claseStudy,i);
                unsigned char *p = (unsigned char*)&a;
                for (int j = 0; j != sizeof(double); ++j)
                    vectInteres2[pos+j] = p[j];
                pos+=8;
            }
        }
        MPI_Bcast(vectInteres2,dataSizes+2+sizeof(double)*256,MPI_CHAR, proceso, MPI_COMM_WORLD);
        tamanho = vectInteres2[dataSizes]+vectInteres2[dataSizes+1];
        rawData.getFeatureStudio(vectInteres2);
        rawData.setValuesFeatureStudio(tamanho);
        mutualInfo.update(rawData);
        if (mpi_rank != proceso)
            probVector = (double *) malloc(tamanho * sizeof(double));
        count = 0;
        for (int j = dataSizes+2; j < (dataSizes+2+sizeof(double)*tamanho); j+=8) {
            double out;
            memcpy(&out, &(vectInteres2[j]), sizeof(double));
            probVector[count] = out;
            count++;
        }
        etapas++;
    }
 
    double timeMeasureMrmr = timeMrmr.stop();
    double timeTotal = tm.stop();
    if (!mpi_rank){
        printf("Time measure: Relevance Time: %f, MrmrTime %f\n",timeMeasureRel, timeMeasureMrmr);
        printf("Final Time: %f\n",timeTotal);}

    if (!mpi_rank) {
        printf("Features selected: ");
        for (uint i = 0; i < selectedFeatures.size(); i++)
            printf("%d,", selectedFeatures[i]);
        printf("\n");
    }

    free(acumThread);free(numFeature);
    free(probVector);
    free(vectInteres2);

    free(MrmrFeaturesProcess);
    free(MrmrFeature);

    free(visited);
    free(relevancesParciales);free(redundancesParciales);

    rawData.destroy();prob.destroy();
    MPI_Finalize();
    return (0);
}
