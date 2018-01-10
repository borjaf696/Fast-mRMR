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
	calculateDSandFS();
	loadData();
	calculateVR();
}

RawData::~RawData() {

}

size_t result;
/**
 *
 */
void RawData::destroy() {
	free(valuesRange);
	free(data);
}

/**Calculates
 *	DataSize: Number of patterns or samples
 *  FeaturesSize: Number of features
 */
void RawData::calculateDSandFS() {
	uint featuresSizeBuffer[1];
	uint datasizeBuffer[1];
	result = fread(datasizeBuffer, sizeof(uint), 1, dataFile);
	if (result != 1) 
		{fputs ("Reading error bufferDatafile\n",stderr); exit (3);}
	result = fread(featuresSizeBuffer, sizeof(uint), 1, dataFile);
	if (result != 1) 
		{fputs ("Reading error featuresDatafile\n",stderr); exit (3);}
	datasize = datasizeBuffer[0];
	featuresSize = featuresSizeBuffer[0];
}

void RawData::loadData() {
	uint i, j;
	t_data buffer[1];
	data = (t_data*) calloc(featuresSize, sizeof(t_data) * datasize);
	fseek(dataFile, 8, 0); 
	for (i = 0; i < datasize; i++) {
		for (j = 0; j < featuresSize; j++) {
			result = fread(buffer, sizeof(t_data), 1, dataFile);
			if (result != 1) 
				{fputs ("Reading error loadFile\n",stderr); exit (3);}
			long int index = (long)j*datasize+i;
			data[index] = buffer[0];
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
			long int index = (long) i*datasize+j;
			dataReaded = data[index];
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
 * Returns a vector containing a feature. Realmente devuelve un puntero
 */
t_feature RawData::getFeature(int index) {
	long int index2 = (long) index * datasize;
	return data + index2;
}

