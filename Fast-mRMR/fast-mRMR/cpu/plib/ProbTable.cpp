#include "ProbTable.h"
#include "omp.h"
#include <typeinfo>

ProbTable::ProbTable(RawData & rd): rawData(rd) {
	this->datasize = rawData.getDataSize();
	this->featuresSize = rawData.getFeaturesSize();
	this->valuesRange = rawData.getValuesRangeArray();
	this->table = (t_prob_table) malloc(featuresSize * sizeof(t_prob *));
	calculate();
}

ProbTable::~ProbTable() {

}

//Calculates the marginal probability table for each posible value in a feature.
//This table will be cached in memory to avoid repeating calculations.
void ProbTable::calculate() {
	Histogram histogram = Histogram(rawData);
	uint i = 0;
	uint j = 0;
	t_prob value = 0;
	t_histogram hist_data;
	for (i = 0; i < featuresSize; ++i) {
		table[i] = (t_prob*) malloc(valuesRange[i] * sizeof(t_prob));
		hist_data = histogram.getHistogram(i);
		for (j = 0; j < valuesRange[i]; ++j) {
			value = (t_prob) hist_data[j] / (t_prob) datasize;
			table[i][j] = value;
		}
		free(hist_data);
	}
}

/**
 * @param featureIndex The index of the feature you want to get the probability.
 * @return The selected feature prob value.
 */
double ProbTable::getProbability(uint index, t_data value) {
	return table[index][value];
}

void ProbTable::destroy() {
	for (uint i = 0; i < featuresSize; i++) {
		free(table[i]);
	}
	free(table);
}
