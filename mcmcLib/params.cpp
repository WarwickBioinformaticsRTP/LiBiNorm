#include <algorithm>
#include "rand.h"
#include "params.h"

using namespace std;


//	Sets up a parameter within the mcmc chain, using the min and max values to set an initial value
//	that is randomnly placed somewhere betweeen the min and the max values
paramType::paramType(std::string name,double min, double max):
	name(name),min(min),max(max),initial(0),
	targetflag(true), 
	localflag(false)
{
	value = rand(max - min) + min;
};

paramType::paramType(std::string name, double min, double max, double initial) :
	name(name), min(min), max(max), initial(initial),
	targetflag(true),
	localflag(false)
{
	value = initial;
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

