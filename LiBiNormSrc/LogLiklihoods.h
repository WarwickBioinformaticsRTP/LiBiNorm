#ifndef LOG_LIKLIHOODS_H
#define LOG_LIKLIHOODS_H

#include "mcmc.h"
#include "dataVec.h"


double FLL_ModelBD(const dataVec & param, const dataType & data);
double FLL_ModelB(const dataVec & param, const dataType & data);



#endif

