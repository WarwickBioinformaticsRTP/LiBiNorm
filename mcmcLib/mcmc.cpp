#include <random>
#include "mcmc.h"

using namespace std;

double rand(double a)
{
	static std::random_device rd;
	static std::mt19937 gen(rd());
	static std::uniform_real_distribution<> dis(0, 1);
	return dis(gen) * a;
}


dataVec randn(size_t x)
{
	static	std::default_random_engine generator;
	static	std::normal_distribution<double> distribution;

	dataVec retVal(x);
	for (auto & i : retVal)
		i = distribution(generator);
	return retVal;
}


mcmc::mcmc(void)
{
}


mcmc::~mcmc(void)
{
}

double mcmc::priorfun(const dataVec & th, const dataVec & mu, const dataVec & sig)
{
	  return sum(((th-mu)/sig)^2);
}
void mcmc::mcmcrun(const modelType & model,const dataType & data,const paramSet & params,const optionsType & options)
{

	dataVec qcov = options.qcov;
	size_t npar = params.size();
	dataVec oldpar = params.getvalues();
	dataVec thetamu = params.getMus();
	dataVec thetasig = params.getSigmas();

/*	if (qcov.empty())
	{
		qcov = params.getSigmas()^2.0;

		for (size_t i = 0;i < qcov.size();i++)
		{
			if (qcov[i] ==0)
				qcov[i] = 1;
		}
	}

	double qcov_scale = 2.4 / sqrt(npar);*/
	dataVec R = qcov.chol();


	double ss = model.ssfun(oldpar,data);
//	ss = sseval(ssfun,ssstyle,oldpar,parind,value,local,data,modelfun);
	double ss1 = ss;
	double ss2 = ss;

	double oldprior = priorfun(oldpar,thetamu,thetasig);

	double sigma2 = model.sigma2;

	_chain.resize(options.nsimu);
	_sschain.resize(options.nsimu);
	_chain[0] = oldpar;
	_sschain[0] = ss;

	rej=0; reju=0; ii=1; rejl = 0;


/*
if isempty(qcov)
  qcov = thetasig.^2;
  ii = isinf(qcov)|isnan(qcov);
%  qcov(ii) = [abs(par0(ii))*0.05].^2; % default is 5% std
  qcov(ii) = [abs(value(ii))*0.05].^2; % default is 5% std
  qcov(qcov==0) = 1; % .. or one if we start from zero
  qcov = diag(qcov);
end

if isempty(adascale)||adascale<=0
  qcov_scale = 2.4 / sqrt(npar) ; % scale factor in R
else
  qcov_scale = adascale;
end

[cm,cn]=size(qcov);
if min([cm cn]) == 1 % qcov contains variances!
  s = sqrt(qcov(parind));
  R = diag(s); % *qcov_scale; % do NOT scale the initial qcov
  qcovorig = diag(qcov); % save qcov
  qcov = diag(qcov(parind));
else %  qcov has covariance matrix in it
  qcovorig = qcov; % save qcov
  qcov = qcov(parind,parind);
  R    = chol(qcov); % *qcov_scale;
end

*/



/*
Line 464 onwards
chainind = 1; % where we are in chain
for isimu=2:nsimu % simulation loop
  ii = ii+1; % local adaptation index (?)
  chainind = chainind+1;

  % waitbar
  if wbarupd;
	status = wbar('',isimu,nsimu);
	if status == -1 % waitbar killed, cancel the run and keep
					% the chain so far
	  message(verbosity,1,'Cancelling...\n');
	  chainind = chainind-1;
	  nsimu = isimu;
	  chain = chain(1:chainind,:);
	  sschain = sschain(1:chainind,:);
	  if updatesigma
		s2chain = s2chain(1:chainind,:);
	  end
	  if size(hchain,1)>1
		hchain = hchain(1:chainind,:);
	  end
	  break % break the nsimu loop
	end
  end
  message(verbosity,100,'i:%d/%d\n',isimu,nsimu);
*/

	bool accept,outbound;
	double newprior,tst;
	int chainind = 0;
	for (size_t isimu = 1; isimu < options.nsimu; isimu++)
	{
		chainind++;



/*
  % sample new candidate from Gaussian proposal
  u = randn(1,npar);
  newpar=oldpar+u*R;

  % reject points outside boundaries
  if any(newpar<low(parind)) || any(newpar>upp(parind))
	accept = 0;
	newprior = 0;
	tst      = 0;
	ss1      = Inf;
	ss2      = ss;
	outbound = 1;
  else
	outbound = 0;
	% prior SS for the new theta
	newprior = feval(priorfun,newpar,thetamu(parind),thetasig(parind));

	% calculate ss
	ss2 = ss;             % old ss
	ss1 = sseval(ssfun,ssstyle,newpar,parind,value,local,data,modelfun);

	tst = exp(-0.5*( sum((ss1-ss2)./sigma2) + newprior-oldprior) );

	if tst <= 0
	  accept = 0;
	elseif tst >= 1
	  accept = 1;
	elseif tst > rand(1,1)
	  accept = 1;
	else
	  accept = 0;
	end
	if shdebug && fix(isimu/shdebug) == isimu/shdebug
	  fprintf('%d: pri: %g, tst: %g, ss: %g\n',isimu, newprior,tst, ss1);
	end
  end

*/
	dataVec u = randn(params.size());
	dataVec newpar = oldpar + u*R;

	if (!params.isValid(newpar))
	{
		accept = false;
		newprior = 0;
		tst      = 0;
		ss1      = std::numeric_limits<double>::max();
		ss2      = ss;
		outbound = true;

	}
	else
	{
		outbound = false;
		newprior = priorfun(newpar,thetamu,thetasig);

		ss2 = ss;             //old ss
		ss1 = model.ssfun(newpar,data);

		tst = exp(-0.5*( (ss1-ss2)/sigma2) + newprior-oldprior); //???????????????
		if (tst <= 0)
			accept = false;
		else if (tst >= 1)
			accept = true;
		else if (tst > rand(1))
			accept = true;
		else
			accept = false;

	}


/*

  %%% DR -----------------------------------------------------
  if dodram == 1 && accept == 0 % & outbound == 0
	% we do a new try according to delayed rejection
	x.p   = oldpar;
	x.ss  = ss2;
	x.pri = oldprior;
	x.s2  = sigma2;

	y.p   = newpar;
	y.ss  = ss1;
	y.pri = newprior;
	y.s2  = sigma2;
	y.a   = tst;

	trypath = {x,y};
	itry    = 1;
	while accept == 0 && itry < Ntry
	  itry = itry+1;
	  z.p  = x.p + randn(1,npar)*RDR{itry};
	  z.s2 = sigma2;
	  if any(z.p<low(parind)) || any(z.p>upp(parind))
		z.a   = 0;
		z.pri = 0;
		z.ss  = Inf;
		trypath = {trypath{:},z};
		outbound = 1;
		continue
	  end

	  outbound = 0;
	  z.ss = sseval(ssfun,ssstyle,z.p,parind,value,local,data,modelfun);
	  z.pri = feval(priorfun,z.p,thetamu(parind),thetasig(parind));
	  trypath = {trypath{:},z};
	  alpha = alphafun(trypath{:});
	  trypath{end}.a = alpha;
	  if alpha >= 1 || rand(1,1) < alpha     %  accept
		accept   = 1;
		newpar   = z.p;
		ss1      = z.ss;
		newprior = z.pri;
		iacce(itry) = iacce(itry) + 1;
	  end
	  if shdebug && fix(isimu/shdebug) == isimu/shdebug
		fprintf('try %d: pri: %g, alpha: %g\n',itry, z.pri, alpha);
		fprintf(' p: %g\n',z.p);
	  end
	end
	if dostats2
	  evalchain(chainind) = itry;
	end
  end % DR --------------------------------------------------------
  */

/*
	  if accept
	%%% accept
	chain(chainind,:) = newpar;
	oldpar     = newpar;
	oldprior   = newprior;
	ss         = ss1;
	if dostats2
	  accechain(chainind) = 1;
	end
  else
	%%%% reject
	chain(chainind,:) = oldpar;
	rej        = rej + 1;
	reju       = reju + 1;
	if outbound
	  rejl     = rejl + 1;
	end
  end
  */

  if (accept)
  {
	  _chain[chainind] = newpar;
	  oldpar     = newpar;
	  oldprior   = newprior;
	  ss         = ss1;
  }
  else
  {
	  _chain[chainind] = oldpar;
	  rej++;
	  reju++;
	  if (outbound)
		  rejl++;
  }

  _sschain[chainind]=ss;

}
}