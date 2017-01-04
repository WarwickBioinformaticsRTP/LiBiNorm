#include <iostream>
#include "mcmc.h"

using namespace std;


double mcmc::priorfun(const dataVec & th, const dataVec & mu, const dataVec & sig)
{
	  return sum(((th-mu)/sig)^2);
}
void mcmc::mcmcrun(const dataType & data,const paramSet & params,const optionsType & options)
{

	dataVec qcov = options.qcov;
	nsimu = options.nsimu;

	const dataVec & thetamu = params.getMus();
	const dataVec & thetasig = params.getSigmas();

	dataVec oldpar = params.getvalues();

	dataVec R = qcov.diagchol();


	double ss = options.ssfun(oldpar,data);

	double ss1 = ss;
	double ss2 = ss;

	double oldprior = priorfun(oldpar,thetamu,thetasig);

	double sigma2 = options.sigma2;

	_chain.resize(options.nsimu);
	_sschain.resize(options.nsimu);
	_chain[0] = oldpar;
	_sschain[0] = ss;

	rej=0; reju=0; ii=1; rejl = 0;
	bool accept,outbound;
	double newprior,tst;
	int chainind = 0;
	for (size_t isimu = 1; isimu < options.nsimu; isimu++)
	{
		chainind++;
		dataVec u = randn(params.size());
		dataVec newpar = oldpar + u*R;

		if (!params.isValid(newpar))
		{
			accept = false;
			newprior = 0;
			tst = 0;
			ss1 = MAX_DOUBLE;
			ss2 = ss;
			outbound = true;

		}
		else
		{
			outbound = false;
			newprior = priorfun(newpar, thetamu, thetasig);

			ss2 = ss;             //old ss
			ss1 = options.ssfun(newpar, data);

			tst = exp(-0.5*((ss1 - ss2) / sigma2) + newprior - oldprior); //???????????????
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
			rej++;
			reju++;
			if (outbound)
				rejl++;
		}

		_sschain[chainind] = ss;

	}
}
