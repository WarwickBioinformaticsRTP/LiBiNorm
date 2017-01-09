#include <random>
#include <algorithm>
#include "params.h"

using namespace std;


//  Returns a random number, uniformly distributed between 0 and a.
double rand(double a)
{
	// Static members generate distribution in range 0 to 1.   This is scaled by input parameter for each call
	static std::random_device rd;
	static std::mt19937 gen(rd());
	static std::uniform_real_distribution<double> dis(0, 1);

	return dis(gen) * a;
}

//	Sets up a vector of values with a normal distribution
dataVec randn(size_t x)
{
	static std::random_device rd;
	static std::mt19937 gen(rd());
	static	std::normal_distribution<double> distribution;

	dataVec retVal(x);
	for (auto & i : retVal)
		i = distribution(gen);
	return retVal;
}

//	Sets up a parameter within the mcmc chain, using the min and max values to set an initial value
//	that is randomnly placed somewhere betweeen the min and the max values
paramType::paramType(std::string name,double min, double max):
	name(name),min(min),max(max), pri_mu(0), 
	pri_sig(MAX_DOUBLE),
	targetflag(true), 
	localflag(false)
{
	value = rand(max - min) + min;
};


//	Use a paramset to verify if a set of values lie within the prescribed range
bool paramSet::isValid(const dataVec & data) const
{
	if (data.size() != size())
		return false;

	for (size_t i = 0; i < size();i++)
		if ((data[i] < at(i).min) || (data[i] > at(i).max))
			return false;
	return true;
}

void paramSet::setValues(const dataVec & vals)
{
	for (size_t i = 0; i < min(size(), vals.size()); i++)
			at(i).value = vals[i];
}

dataVec paramSet::getvalues() const
{
	dataVec retVal(size());
	for (size_t i = 0;i < size();i++)
		retVal[i] = at(i).value;
	return retVal;
}

dataVec paramSet::getSigmas() const
{
	dataVec retVal(size());
	for (size_t i = 0;i < size();i++)
		retVal[i] = at(i).pri_sig;
	return retVal;
}

dataVec paramSet::getMus() const
{
	dataVec retVal(size());
	for (size_t i = 0;i < size();i++)
		retVal[i] = at(i).pri_mu;
	return retVal;
}
