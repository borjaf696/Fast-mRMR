#include "Histogram.h"
#include "omp.h"

Histogram::Histogram(RawData & rd): rawData(rd) {
}

Histogram::~Histogram() {
}

//Calculates the histogram for the given feature index.
t_histogram Histogram::getHistogram(uint index) {
	uint vr = rawData.getValuesRange(index);
	t_feature data = rawData.getFeature(index);
	t_histogram h_acum = (t_histogram) malloc(vr * sizeof(uint));
	// Initialize to zero
	for (uint i = 0; i < vr; i++) {
		h_acum[i] = 0;
	}
	for (uint i = 0; i < rawData.getDataSize(); i++) {
		h_acum[data[i]]++;
	}
	return h_acum;
}

double* Histogram::getHistogram2(uint index) {
	uint vr = rawData.getValuesRange(index);
	t_feature data = rawData.getFeature(index);
	double *salida = (double*) malloc(vr * sizeof(double));
	double step = 1.0 / (double) rawData.getDataSize();
    for (uint i = 0; i < vr; i++) {
		salida[i] = 0;
	}
	for (uint i = 0; i < rawData.getDataSize(); i++) {
		salida[data[i]]+=step;
	}
	return salida;
}
