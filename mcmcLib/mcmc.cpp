#include <iostream>
#include "mcmc.h"

using namespace std;



//	Performs a monte carl markov chain run
void mcmc::mcmcrun(const mcmcGeneData & data,const paramSet & params,const optionsType & options)
{
	dataVec oldpar = params.getvalues();

	dataVec R = options.qcov.diagchol();

	double ss = options.ssfun(oldpar,data);

	double ss1 = ss;
	double ss2 = ss;

	_chain.resize(options.nsimu);
	_sschain.resize(options.nsimu);
	_chain[0] = oldpar;
	_sschain[0] = ss;

	bool accept;
	double tst;
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
			ss1 = MAX_DOUBLE;
			ss2 = ss;
		}
		else
		{
			ss1 = options.ssfun(newpar, data);
			ss2 = ss;             //old ss
			tst = exp(-0.5*((ss1 - ss2) / options.sigma2));
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
			ss = ss1;
		}
		else
		{
			_chain[chainind] = oldpar;
		}

		_sschain[chainind] = ss;

	}
}
