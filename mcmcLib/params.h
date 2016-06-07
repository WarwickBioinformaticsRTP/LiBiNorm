#ifndef PARAMS_H
#define PARAMS_H

#include <limits> 
#include <string>
#include "containerEx.h"
#include "dataVec.h"

class paramType
{
public:
	paramType(std::string name,double initial,double min = std::numeric_limits<double>::min(),double max = std::numeric_limits<double>::max());
	std::string name;
	double value, min, max, pri_mu, pri_sig;
	bool targetflag, localflag;

};

class paramSet: public vectorEx<paramType>
{
public:
	template<typename... P>paramSet(const P & ... params) : vectorEx<paramType>(params...){};

	dataVec getvalues() const;
	dataVec getSigmas() const;
};

#endif