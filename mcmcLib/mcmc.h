#ifndef MCMC_H
#define MCMC_H

#include <string>
#include "dataVec.h"
#include "params.h"


#define MAX_DOUBLE std::numeric_limits<double>::max()

class dataType
{
public:
	dataVec fragData;
	std::vector<int> geneIndex;
	dataVec geneData[2];
};


double rand(double a);

/*
% options structure
%    options.ssfun    -2*log(likelihood) function
%    options.priorfun -2*log(pior) prior function
%    options.sigma2   initial error variance
%    options.Nsimu    number of mcmc cycles in a run
%    options.Nruns    total number of runs/observations
%    model.S20      prior for sigma2
%    model.N0       prior accuracy for S20
%    model.nbatch   number of datasets
*/

class optionsType
{
public:
	optionsType():nsimu(100){};

	double(*ssfun)(const dataVec & param, const dataType & data);
	double sigma2;
	size_t nsimu,Nruns;
	double jumpSize;
	dataVec qcov;
};


//	The class that manages the Monte Carlo Markov Chain determintaion of parameters
class mcmc
{
public:

	const std::vector<dataVec> & chain() {return _chain;};
	const dataVec & sschain() {return _sschain;};

	void mcmcrun(const dataType & data,const paramSet & params,const optionsType & options);
	double priorfun(const dataVec & th, const dataVec & mu, const dataVec & sig);

private:
	std::vector<dataVec> _chain;
	dataVec _sschain;

	double rej, reju, ii, rejl;
	size_t nsimu;

};


#endif // !MCMC_H

