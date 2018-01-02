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
#include "omp.h"
#include "JointProb.h"
#include "mpi.h"
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

JointProb::JointProb(RawData & rd, t_feature featureOld, uint index2,uint index1) :
		rawData(rd) {
	this->index2 = index2;
	this->index1 = index1;
	this->valuesRange1 = featureOld[rawData.getDataSize()]+featureOld[rawData.getDataSize()+1];
	this->valuesRange2 = rawData.getValuesRange(index2);
	this->datasize = rawData.getDataSize();
	this->data = (t_histogram) calloc(valuesRange1 * valuesRange2, sizeof(uint));
	//this->feature = (t_feature)calloc(datasize,sizeof(t_data));
	//this->feature = featureOld;
	if ( this->data == NULL){
		printf("Memory limit exceeded\n");
		exit(3);
	}
	/*for (int i = 0; i < datasize;++i){
		if (i % 1000000 == 0)	
			printf("Value %d\n",this->feature[i]);
		int a = this->feature[i];
	}*/
	/*for (uint i = 0; i < datasize;i++)
		feature[i] = featureOld[i];*/

	//calculate();
	calculate2(featureOld);
}

JointProb::~JointProb() {
	//free(feature);
	free(data);
}
//Calculates the joint probability between the given features.
/*Cuenta ocurrencias de distintos valores asociados a cada valor de otra variable*/
void JointProb::calculate() {
	t_feature h_vector1 = rawData.getFeature(index1);
	t_feature h_vector2 = rawData.getFeature(index2);

	//Calculate histogram in CPU
	/*Hace una matriz de valuesRange1xvaluesRange2 -> con ocurrencias de h_vector1[n]xh_vector2[m]
	,n E Vr1 y m E Vr2. Matriz de coocurencias practicamente. Cuenta cuantas veces ocurren a la vez
	para poder determinar P(y|x) para cada valor de x e y.*/
	for (uint i = 0; i < datasize; i++) {
		data[h_vector1[i] * valuesRange2 + h_vector2[i]]++;
	}
}
//> La P(x,y) -> P(x)*P(y) en caso de ser independientes.
/*P(x,y) = P(y | x)P(x) = P(x | y)*P(y)*/
t_prob JointProb::getProb(t_data valueFeature1, t_data valueFeature2) {
	return (t_prob) data[valueFeature1 * valuesRange2 + valueFeature2]
			/ (t_prob) datasize;
}

void JointProb::calculate2(t_feature feature) {
	t_feature h_vector2 = rawData.getFeature(this->index2);
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
/*	if (rank == 3)
	for (int i = 0; i < datasize;++i){
                if (i % 1000000 == 0)
                        printf("Value %d\n",this->feature[i]);
                int a = this->feature[i];
        }*/
	if (data != NULL && h_vector2 != NULL && feature != NULL)
		for (uint i = 0; i < datasize; i++){
			data[feature[i] * valuesRange2 + h_vector2[i]]++;
		}
	else
		printf("Ahi hay algo que es nulo\n");
}


