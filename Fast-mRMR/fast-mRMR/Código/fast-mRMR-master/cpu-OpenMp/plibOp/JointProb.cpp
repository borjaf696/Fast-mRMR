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

JointProb::~JointProb() {
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

