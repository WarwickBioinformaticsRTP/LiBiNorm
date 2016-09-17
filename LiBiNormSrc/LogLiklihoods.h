#ifndef LOG_LIKLIHOODS_H
#define LOG_LIKLIHOODS_H

#include "mcmc.h"
#include "dataVec.h"

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

