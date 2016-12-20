#ifndef LOG_LIKLIHOODS_H
#define LOG_LIKLIHOODS_H

#include <float.h>
#include "mcmc.h"


struct bestResult
{
	bestResult() :minLL(DBL_MAX), minLL_dev(0), run(0), pos(0) {};
	VEC_DATA_TYPE minLL, minLL_dev;
	size_t run, pos;
	dataVec params;
	std::vector<dataVec> param_dev;
	dataVec norm;
	operator bool() const { return params.size(); };
};

void normaliseExpression(size_t m, bestResult & result, const dataVec & l);

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

