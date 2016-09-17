#ifndef MCMC_H
#define MCMC_H

#include <string>
#include "dataVec.h"
#include "params.h"


double rand(double a);

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
	double sigma2;
	double (*ssfun)(const dataVec & param, const dataType & data);
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
	std::vector<dataVec> _chain;
	dataVec _sschain;

	double rej,reju,ii,rejl;
	size_t nsimu;
public:

	const std::vector<dataVec> & chain() {return _chain;};
	const dataVec & sschain() {return _sschain;};
	double rejected() {return reju/nsimu;};
	mcmc(void);
	~mcmc(void);

	void mcmcrun(const modelType & model, const dataType & data,const paramSet & params,const optionsType & options);
	double priorfun(const dataVec & th, const dataVec & mu, const dataVec & sig);

};


#endif // !MCMC_H

