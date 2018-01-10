#include "MutualInfo.h"
#include <math.h>
#include "omp.h"

MutualInfo::MutualInfo(RawData & rd, ProbTable & pt) :
		 rawData(rd), probTable(pt) {
}

MutualInfo::~MutualInfo() {

}

void MutualInfo::update(RawData rd) {
	rawData.setValuesFeatureStudio(rd.getValuesFeatureStudio());
	rawData.getFeatureStudio(rd.getFeatureStudio());
}

//Calculates the mutual information between the given features.
double MutualInfo::get(uint featureIndex1, uint featureIndex2) {
	uint range1 = rawData.getValuesRange(featureIndex1);
	uint range2 = rawData.getValuesRange(featureIndex2);
	uint i, j;
	t_prob mInfo = 0;
	t_prob jointProb = 0;
	t_prob marginalX = 0;
	t_prob marginalY = 0;
	t_prob division = 0;
	JointProb jointProbTable = JointProb(rawData, featureIndex1, featureIndex2);
	for (i = 0; i < range1; i++) {
		for (j = 0; j < range2; j++) {
			jointProb = jointProbTable.getProb(i, j);
			if (jointProb != 0) {
				marginalX = probTable.getProbability(featureIndex1, i);
				marginalY = probTable.getProbability(featureIndex2, j);
				division = jointProb / (marginalX * marginalY);
				mInfo += jointProb * log2(division);
			}
		}
	}
	return mInfo;
}

//Another way to work
double MutualInfo::get(t_feature featureInt, uint featureIndex2, double* prob1) {
	uint range1 = featureInt[rawData.getDataSize()]+featureInt[rawData.getDataSize()+1];
	uint range2 = rawData.getValuesRange(featureIndex2);
	uint i, j;
	t_prob mInfo = 0;
	t_prob jointProb = 0;
	t_prob marginalX = 0;
	t_prob marginalY = 0;
	t_prob division = 0;
	JointProb jointProbTable = JointProb(rawData, featureInt, featureIndex2);
	for (i = 0; i < range1; i++) {
		for (j = 0; j < range2; j++) {
			jointProb = jointProbTable.getProb(i, j);
			if (jointProb != 0) {
				marginalX = prob1[i];
				marginalY = probTable.getProbability(featureIndex2, j);
				division = jointProb / (marginalX * marginalY);
				mInfo += jointProb * log2(division);
			}
		}
	}
	jointProbTable.destroy();
	return mInfo;
}

double MutualInfo::get(uint featureIndex2, double* prob1) {
	uint range1 = rawData.getValuesFeatureStudio();
	uint range2 = rawData.getValuesRange(featureIndex2);
	uint i, j;
	t_prob mInfo = 0;
	t_prob jointProb = 0;
	t_prob marginalX = 0;
	t_prob marginalY = 0;
	t_prob division = 0;
	JointProb jointProbTable = JointProb(rawData, featureIndex2);
	for (i = 0; i < range1; i++) {
		for (j = 0; j < range2; j++) {
			jointProb = jointProbTable.getProb(i, j);
			if (jointProb != 0) {
				marginalX = prob1[i];
				marginalY = probTable.getProbability(featureIndex2, j);
				division = jointProb / (marginalX * marginalY);
				mInfo += jointProb * log2(division);
			}
		}
	}
	jointProbTable.destroy();
	return mInfo;
}
