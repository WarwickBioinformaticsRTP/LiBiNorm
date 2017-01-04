#ifndef PARAMS_H
#define PARAMS_H

#include <limits> 
#include <string>
#include "dataVec.h"


#define MAX_DOUBLE std::numeric_limits<double>::max()
#define MIN_DOUBLE std::numeric_limits<double>::min()
double rand(double a);
dataVec randn(size_t x);

class paramType
{
public:
	paramType(std::string name,double min = MIN_DOUBLE, double max = MAX_DOUBLE);
	std::string name;
	double value, min, max, pri_mu, pri_sig;
	bool targetflag, localflag;
};

class paramSet :public std::vector<paramType>
{
public:
	paramSet() {};
	//	This = operator allows the values to be assigned using an aggregate initiailiser for the base class
	paramSet & operator = (std::vector<paramType>(a)) { std::vector<paramType>::operator = (a); return *this; };

	bool isValid(const dataVec & data) const;

	dataVec getvalues() const;
	dataVec getSigmas() const;
	dataVec getMus() const;

};

#endif