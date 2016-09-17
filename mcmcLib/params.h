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
	//  template constructors do not require explicit declaration of the template type when they are used, as this is taken from the types
	//	in the constructor
	paramSet() {};
	template<typename... P>paramSet(const paramType & first,const P & ... rest) : vectorEx<paramType>(first,rest...){};
	bool isValid(const dataVec & data) const;

	dataVec getvalues() const;
	dataVec getSigmas() const;
	dataVec getMus() const;

};

#endif