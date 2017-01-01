#ifndef LOG_LIKLIHOODS_H
#define LOG_LIKLIHOODS_H

#include <float.h>
#include <algorithm>
#include "mcmc.h"


enum modelType
{
	noModel = 0,
	ModelA = 1,
	ModelB = 2,
	ModelC = 3,
	ModelD = 4,
	ModelE = 5,
	ModelBD = 6
};

std::string conv(const modelType m);

//  Selects a model based on a string, exits if the string is not valid
modelType modelFromString(const std::string & desc);

//	Allows the model to be output to a stream such as std::out as appropriate text
bool printVal(outputDataFile * f, modelType m);


//	Allows genomicPositions to be placed into StringEx's
namespace std
{
	std::string to_string(const modelType & m);
}
std::ostream& operator<< (std::ostream &out, const modelType & m);


struct bestResult
{
	bestResult() :minLL(DBL_MAX), minLL_dev(0), run(0), pos(0) {};

	//	For checking if paraeters have been estimated
	operator bool() const { return params.size(); };

	VEC_DATA_TYPE minLL, minLL_dev;
	size_t run, pos;
	dataVec params;
	dataVec param_dev[2];
};

void getBias(modelType m, dataVec & params, const dataVec & l,dataVec & bias);

#ifdef _DEBUG
#define VERIFY_SPEEDUP
#endif

double FLL_ModelA(const dataVec & param, const dataType & data);
double FLL_ModelB(const dataVec & param, const dataType & data);
double FLL_ModelC(const dataVec & param, const dataType & data);
double FLL_ModelD(const dataVec & param, const dataType & data);
double FLL_ModelE(const dataVec & param, const dataType & data);
double FLL_ModelBD(const dataVec & param, const dataType & data);

#endif

