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
#include <cmath>
#include <stdlib.h>

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
//	MPI_Init(&argc,&argv);
	options opts;
	uint i = 0;
	uint j = 0;
	uint newFeatureIndex = 0;
	uint lastFeatureIndex = 0;
    int threadsMaximos = omp_get_max_threads();
	double mrmr = 0;
	double acum = 0;
    double *acumThread = (double*)malloc(threadsMaximos*sizeof(double));
    int *numFeature = (int*)malloc(threadsMaximos*sizeof(int));
	/*vector<double> relevances;
	vector<double> redundances;*/
	double * relevances;
	double * redundances;
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


	//Get relevances between all features and class.
	Timer timeRel;
	timeRel.start();
	relevances = (double*) malloc (rawData.getFeaturesSize()*sizeof(double));
	redundances = (double*) malloc (rawData.getFeaturesSize()*sizeof(double));
	#pragma omp parallel for schedule(dynamic,10)
	for (i = 0; i < rawData.getFeaturesSize(); ++i) {
		/*relevances.push_back(mutualInfo.get(opts.classIndex, i));
		redundances.push_back(0);*/
		relevances[i] = mutualInfo.get(opts.classIndex,i);
		redundances[i] = 0;
	}

	/*for (uint i = 0; i < rawData.getFeaturesSize(); i++)
		printf("Relevances %d %f\n",i,relevances[i]);*/

	double timeMeasureRel = timeRel.stop();

	// Max relevance feature is added because no redundancy is possible.
	//newFeatureIndex = getMaxRelevance(relevances, opts.classIndex);
	newFeatureIndex = getMaxRelevance(relevances, opts.classIndex, rawData.getFeaturesSize());
	selectedFeatures.push_back(newFeatureIndex);
	visited[newFeatureIndex] = 1;
	lastFeatureIndex = newFeatureIndex;

	std::cout << newFeatureIndex << ",";

	Timer timeMrmr;
	timeMrmr.start();
	while (selectedFeatures.size() < rawData.getFeaturesSize() - 1 /*-1 because class is discarded*/ and selectedFeatures.size() < opts.selectedFeatures) {

		acum = -std::numeric_limits<double>::infinity();
        	for (int z = 0; z < threadsMaximos; ++z)
            		acumThread[z] = -std::numeric_limits<double>::infinity();
		#pragma omp parallel for schedule(dynamic,10) private(mrmr)
		for (j = 0; j < rawData.getFeaturesSize(); ++j){ 
			if (visited[j] == 0 && j != opts.classIndex) {
				redundances[j] += mutualInfo.get(lastFeatureIndex, j);
		                mrmr = relevances[j] - (redundances[j] / selectedFeatures.size());
		                if (mrmr > acumThread[omp_get_thread_num()]) {
                		    acumThread[omp_get_thread_num()] = mrmr;
		                    numFeature[omp_get_thread_num()] = j;
                		}
			}
		}
	        for (int z = 0; z < threadsMaximos;++z){
        	    if (acumThread[z] > acum){
                	acum = acumThread[z];
	                newFeatureIndex = numFeature[z];
        	    }
        	}
		//Last feature doesn't prints comma.
		if ( (selectedFeatures.size() == (opts.selectedFeatures - 1)) or (selectedFeatures.size() == (rawData.getFeaturesSize() -2)) ){
			std::cout << newFeatureIndex << "\n";
		}else{
			std::cout << newFeatureIndex << ",";
		}
		selectedFeatures.push_back(newFeatureIndex);
		visited[newFeatureIndex] = 1;
		lastFeatureIndex = newFeatureIndex;
	}
	double timeMeasureMrmr = timeMrmr.stop();
	printf("Measure Times: relevance Time: %f, mRMR-Time: %f\n",timeMeasureRel, timeMeasureMrmr);

    double timeTotal = tm.stop();
    printf("Final Total time: %f\n",timeTotal);

    free(acumThread);free(numFeature);
    free(visited);free(relevances);free(redundances);
	rawData.destroy();prob.destroy();
	printf("\n");
	return (0);
}

