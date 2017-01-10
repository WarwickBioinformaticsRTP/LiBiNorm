#ifndef PARAMS_H
#define PARAMS_H

#include <limits> 
#include "dataVec.h"


#define MAX_DOUBLE std::numeric_limits<double>::max()
#define MIN_DOUBLE std::numeric_limits<double>::min()
double rand(double a);
dataVec randn(size_t x);

//	Holds the subset of the read information associated with a gene that is used by the functions for 
//	calculating log liklyhoods with the mcmc
class mcmcGeneData
{
public:
	//	A list of all the read positions
	dataVec fragData;
	//	The list of the indexes to the genes associated with the fragments
	std::vector<int> geneIndex;
	//	The data associated with each gene, namely the length and the frequency with which it occurs
	dataVec geneLengths;
	dataVec geneFrequencies;
};


class paramType
{
public:
	paramType(std::string name,double min = MIN_DOUBLE, double max = MAX_DOUBLE);
	std::string name;
	double value, min, max;
	bool targetflag, localflag;
};

class paramSet :public std::vector<paramType>
{
public:
	paramSet() {};
	//	This = operator allows the values to be assigned using an aggregate initiailiser for the base class
	paramSet & operator = (std::vector<paramType>(a)) { std::vector<paramType>::operator = (a); return *this; };

	bool isValid(const dataVec & data) const;
	void setValues(const dataVec & vals);

	dataVec getvalues() const;

};

#endif