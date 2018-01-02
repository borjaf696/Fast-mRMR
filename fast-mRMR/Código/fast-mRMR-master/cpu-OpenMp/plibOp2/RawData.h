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

#ifndef RAWDATA_H_
#define RAWDATA_H_
#include "utils.h"
#include <fstream>
#include <iostream>
#include <vector>

class RawData {
public:
	RawData(char * filename); //Constructor
	RawData();
	virtual ~RawData(); //Constructor abstracto

	//Metodos variados
	void calculateVR();
	void calculateDSandFS();
	void destroy();
	void setValues(Envio send);
	void setValues(t_data *data);
	uint getValuesRange(uint index);
	uint* getValuesRangeArray();
	t_data * getDataArray();
	void setValuesDataArray(t_data *ranges);
	void setValuesRangeArray(uint* ranges);
	uint getDataSize();
	uint getFeaturesSize();
	t_feature getFeature(int index);
	t_feature getVectGet();
	void setVectGet(int index);
	void setVectGet2(t_feature feature, uint range);

private:

	void loadData();
	t_data* data;
	uint featuresSize;
	uint datasize;
	uint* valuesRange;
	FILE * dataFile;
	t_feature vectGet;
};

#endif /* RAWDATA_H_ */
