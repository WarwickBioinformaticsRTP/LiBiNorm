#ifndef LOG_LIKLIHOODS_H
#define LOG_LIKLIHOODS_H

#include <float.h>
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

inline std::string conv(const modelType m)
{
	switch (m)
	{
	case noModel: return "No Model";
	case ModelA: return "Model A";
	case ModelB: return "Model B";
	case ModelC: return "Model C";
	case ModelD: return "Model D";
	case ModelE: return "Model E";
	case ModelBD: return "Model BD";
	}
	return "";
}

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

