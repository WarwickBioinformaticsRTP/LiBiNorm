#ifndef MCMC_H
#define MCMC_H

#include <string>
#include "dataVec.h"
#include "params.h"

/*
% MODEL   model options structure
%    model.ssfun    -2*log(likelihood) function
%    model.priorfun -2*log(pior) prior function
%    model.sigma2   initial error variance
%    model.N        total number of observations
%    model.S20      prior for sigma2
%    model.N0       prior accuracy for S20
%    model.nbatch   number of datasets
*/
class modelType
{
public:
	double (*ssfun)(dataVec & param, dataType & data);

};



class optionsType
{
public:
	optionsType():nsimu(100){};
	size_t nsimu;
	double updatesigma;
	std::string method;
	std::vector<double> qcov;
};



class mcmc
{
public:
	mcmc(void);
	~mcmc(void);

	void mcmcrun(const modelType & model, const dataType & data,const paramSet & params,const optionsType & options);
};


#endif // !MCMC_H

