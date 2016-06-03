#ifndef DATAVEC_H
#define DATAVEC_H

#include <vector>
#include <string>

class transcriptDataMap;

class dataVec: public std::vector<double>
{
public:
	dataVec(size_t s = 0): std::vector<double>(s){};

	vector<double> & values() { return *this;};

	dataVec operator > (double a) 
	{
		dataVec retVal(size());
		for (size_t i = 0; i < size();i++)
			retVal[i] = at(i)>a?1:0;
		return retVal;
	}
	dataVec operator - (double a) 
	{
		dataVec retVal(size());
		for (size_t i = 0; i < size();i++)
			retVal[i] = at(i) - a;
		return retVal;
	}
	dataVec operator + (double a) 
	{
		dataVec retVal(size());
		for (size_t i = 0; i < size();i++)
			retVal[i] = at(i) + a;
		return retVal;
	}
	dataVec operator / (double a) 
	{
		dataVec retVal(size());
		for (size_t i = 0; i < size();i++)
			retVal[i] = at(i) / a;
		return retVal;
	}
	dataVec operator * (double a)
	{
		dataVec retVal(size());
		for (size_t i = 0; i < size();i++)
			retVal[i] = at(i) * a;
		return retVal;
	}
	dataVec operator ^ (double a)
	{
		dataVec retVal(size());
		for (size_t i = 0; i < size();i++)
			retVal[i] = pow(at(i),a);
		return retVal;
	}

	dataVec operator * (const dataVec & a)
	{
		dataVec retVal(size());
		for (size_t i = 0; i < size();i++)
			retVal[i] = at(i) * a[i];
		return retVal;
	}
	dataVec operator / (const dataVec & a)
	{
		dataVec retVal(size());
		for (size_t i = 0; i < size();i++)
			retVal[i] = at(i) / a[i];
		return retVal;
	}
	dataVec operator + (const dataVec & a)
	{
		dataVec retVal(size());
		for (size_t i = 0; i < size();i++)
			retVal[i] = at(i) + a[i];
		return retVal;
	}
	dataVec operator - (const dataVec & a)
	{
		dataVec retVal(size());
		for (size_t i = 0; i < size();i++)
			retVal[i] = at(i) - a[i];
		return retVal;
	}
	dataVec operator - ()
	{
		dataVec retVal(size());
		for (size_t i = 0; i < size();i++)
			retVal[i] = -at(i);
		return retVal;
	}

	bool does_not_contain_null()
	{
		for (auto & i : *this)
			if (i == 0)
				return false;
		return true;
	}

	dataVec & append(const dataVec a)
	{
		insert(end(),a.begin(),a.end());
		return *this;
	}
	dataVec & append(size_t N,double value)
	{
		insert(end(),N,value);
		return *this;
	}

	void selectAtMost(size_t s)
	{
		if (s > size())
			return;
		iterator i = begin();
		for (size_t j = 0;j < s;j++)
		{
			std::swap (*(i),*(begin() + rand() % size()));
		}
		resize(s);
	}

	dataVec & removeInvalidValues(double maxVal)
	{
		//	Sort in place for maximum efficiency, if we find an invalid value, replace with one from the end;
		iterator i = begin(), j = end();
		while (i != j)
		{
			if ((*i < 0) || (*i > maxVal))
				std::swap(*i,*--j);
			else
				i++;
		}
		resize(j - begin());
		return *this;
	}
};

inline double sum(const dataVec & a)
{
	double retVal = 0;
	for (auto & i : a)
		retVal += i;
	return retVal;
}

inline dataVec exp(const dataVec & a)
{
	dataVec retVal(a.size());
	for (size_t i = 0; i < a.size();i++)
		retVal[i] = exp(a[i]);
	return retVal;
}
inline dataVec log(const dataVec & a)
{
	dataVec retVal(a.size());
	for (size_t i = 0; i < a.size();i++)
		retVal[i] = log(a[i]);
	return retVal;
}

inline dataVec operator * (double a, const dataVec & b)
{
	dataVec retVal(b.size());
	for (size_t i = 0; i < b.size();i++)
		retVal[i] = a * b[i];
	return retVal;

}

inline dataVec operator < (double a, const dataVec & b)
{
	dataVec retVal(b.size());
	for (size_t i = 0; i < b.size();i++)
		retVal[i] = a < b[i]?1:0;
	return retVal;

}
inline dataVec operator - (double a, const dataVec & b)
{
	dataVec retVal(b.size());
	for (size_t i = 0; i < b.size();i++)
		retVal[i] = a - b[i];
	return retVal;

}
inline dataVec operator - (double a, dataVec && b)
{
	for (size_t i = 0; i < b.size();i++)
		b[i] = a - b[i];
	return b;
}





typedef std::vector<double> paramType;
class dataType 
{
	std::vector<dataVec> data;
public:
	dataVec & operator [] (int i) { return data[i];};

	dataType() : data(3){};
};

#endif