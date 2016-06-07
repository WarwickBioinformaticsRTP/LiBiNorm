#include "params.h"

using namespace std;

paramType::paramType(std::string name,double initial,double min, double max):
	name(name),value(initial),min(min),max(max), pri_mu(0), pri_sig(0),targetflag(true), localflag(false){};


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
