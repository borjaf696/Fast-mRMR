#include "omp.h"
#include "JointProb.h"

JointProb::JointProb(RawData & rd, uint index1, uint index2) :
		rawData(rd) {
	this->index1 = index1;
	this->index2 = index2;
	this->valuesRange1 = rawData.getValuesRange(index1);
	this->valuesRange2 = rawData.getValuesRange(index2);
	this->datasize = rawData.getDataSize();
	this->data = (t_histogram) calloc(valuesRange1 * valuesRange2,
			sizeof(uint));

	calculate();
}

JointProb::JointProb(RawData & rd, t_feature featureOld, uint index2) :
		rawData(rd) {
	this->index2 = index2;
	this->valuesRange1 = featureOld[rawData.getDataSize()]+featureOld[rawData.getDataSize()+1];
	this->valuesRange2 = rawData.getValuesRange(index2);
	this->datasize = rawData.getDataSize();
	this->data = (t_histogram) calloc(valuesRange1 * valuesRange2,
									  sizeof(uint));
	this->feature = (t_feature)malloc((datasize)*sizeof(t_data));
	for (int i = 0; i < datasize;i++)
		feature[i] = featureOld[i];

	//calculate();
	calculate2();
}

JointProb::JointProb(RawData & rd, uint index2) :
		rawData(rd) {
	this->index2 = index2;
	this->valuesRange1 = rawData.getValuesFeatureStudio();
	this->valuesRange2 = rawData.getValuesRange(index2);
	this->datasize = rawData.getDataSize();
	this->data = (t_histogram) calloc(valuesRange1 * valuesRange2,
									  sizeof(uint));

	//calculate();
	calculate3();
}

JointProb::~JointProb() {
	free(data);
}
//Calculates the joint probability between the given features.
void JointProb::calculate() {
	t_feature h_vector1 = rawData.getFeature(index1);
	t_feature h_vector2 = rawData.getFeature(index2);

	//Calculate histogram in CPU
	for (uint i = 0; i < datasize; i++) {
		data[h_vector1[i] * valuesRange2 + h_vector2[i]]++;
	}
}
t_prob JointProb::getProb(t_data valueFeature1, t_data valueFeature2) {
	return (t_prob) data[valueFeature1 * valuesRange2 + valueFeature2]
			/ (t_prob) datasize;
}

void JointProb::calculate2() {
	t_feature h_vector1 = this->feature;
	t_feature h_vector2 = rawData.getFeature(index2);

	for (uint i = 0; i < datasize; i++) {
		data[h_vector1[i] * valuesRange2 + h_vector2[i]]++;
	}
}

void JointProb::calculate3() {
	t_feature h_vector1 = rawData.getFeatureStudio();
	t_feature h_vector2 = rawData.getFeature(index2);

	for (uint i = 0; i < datasize; i++) {
		data[h_vector1[i] * valuesRange2 + h_vector2[i]]++;
	}
}


void JointProb::destroy(){

}

