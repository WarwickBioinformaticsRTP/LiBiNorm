#include <iostream>
#include "mcmc.h"

using namespace std;

//	The default prior function used in teh chain.  The original code allowed a custom function
//	to be defined but this was not being used for this application.
double mcmc::priorfun(const dataVec & th, const dataVec & mu, const dataVec & sig)
{
	  return sum(((th-mu)/sig)^2);
}

//	Performs a monte carl markov chain run
void mcmc::mcmcrun(const mcmcGeneData & data,const paramSet & params,const optionsType & options)
{
	const dataVec & thetamu = params.getMus();
	const dataVec & thetasig = params.getSigmas();

	dataVec oldpar = params.getvalues();

	dataVec R = options.qcov.diagchol();
	double sigma2 = options.sigma2;

	double ss = options.ssfun(oldpar,data);

	double ss1 = ss;
	double ss2 = ss;

	double oldprior = priorfun(oldpar,thetamu,thetasig);

	_chain.resize(options.nsimu);
	_sschain.resize(options.nsimu);
	_chain[0] = oldpar;
	_sschain[0] = ss;

	bool accept;
	double newprior,tst;
	int chainind = 0;
	for (size_t isimu = 1; isimu < options.nsimu; isimu++)
	{
		chainind++;
		dataVec u = randn(params.size());
		dataVec newpar = oldpar + u*R;

		//	If the parameters go out of range, reject them
		if (!params.isValid(newpar))
		{
			accept = false;
			newprior = 0;
			ss1 = MAX_DOUBLE;
			ss2 = ss;
		}
		else
		{
			newprior = priorfun(newpar, thetamu, thetasig);
			ss1 = options.ssfun(newpar, data);
			ss2 = ss;             //old ss
			tst = exp(-0.5*((ss1 - ss2) / sigma2) + newprior - oldprior); 
			if (tst <= 0)
				accept = false;
			else if (tst >= 1)
				accept = true;
			else if (tst > rand(1))
				accept = true;
			else
				accept = false;
		}

		if (accept)
		{
			_chain[chainind] = newpar;
			oldpar = newpar;
			oldprior = newprior;
			ss = ss1;
		}
		else
		{
			_chain[chainind] = oldpar;
		}

		_sschain[chainind] = ss;

	}
}
