#include "params.h"

using namespace std;

paramType::paramType(std::string name,double initial,double min, double max):
	name(name),value(initial),min(min),max(max), pri_mu(0), pri_sig(std::numeric_limits<double>::max()),targetflag(true), localflag(false){};


bool paramSet::isValid(const dataVec & data) const
{
	if (data.size() != size())
		return false;

	for (size_t i = 0; i < size();i++)
		if ((data[i] < at(i).min) || (data[i] > at(i).max))
			return false;
	return true;
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
