#ifndef DATAVEC_H
#define DATAVEC_H

#include <stdlib.h>
#include <vector>
#include <string>
#include <valarray>
#include "printEx.h"


inline bool printVal(outputDataFile * f,long double value)
{
	if (value != 0) 
		fprintf(f->fout,"%f",value);
	return true;
};


class dataVec: public std::vector<double>
{
public:
	dataVec(size_t s = 0): std::vector<double>(s){};
	dataVec(const std::vector<double> & a): std::vector<double>(a){};

#ifdef _DEBUG
	virtual ~dataVec()
	{
	}
#endif
	vector<double> & values() { return *this;};

	dataVec operator > (double a) const
	{
		size_t s = size();
		dataVec retVal(s);
		for (size_t i = 0; i < s;i++)
			retVal[i] = at(i)>a?1:0;
		return retVal;
	}
	dataVec operator < (const dataVec & a) const
	{
		size_t s = size();
		dataVec retVal(s);
		for (size_t i = 0; i < s;i++)
			retVal[i] = at(i)<a[i]?1:0;
		return retVal;
	}
	dataVec operator < (dataVec && a) const
	{
		size_t s = size();
		for (size_t i = 0; i < s;i++)
			a[i] = at(i)<a[i]?1:0;
		return a;
	}
	dataVec operator * (const dataVec & a) const
	{
		size_t s = size();
		dataVec retVal(s);
		for (size_t i = 0; i < s;i++)
			retVal[i] = at(i) * a[i];
		return retVal;
	}
	dataVec operator * (dataVec && a) const
	{
		size_t s = size();
		for (size_t i = 0; i < s;i++)
			a[i] *= at(i);
		return a;
	}
	dataVec operator / (const dataVec & a) const 
	{
		size_t s = size();
		dataVec retVal(s);
		for (size_t i = 0; i < s;i++)
			retVal[i] = at(i) / a[i];
		return retVal;
	}
	dataVec operator / (dataVec && a) const 
	{
		size_t s = size();
		for (size_t i = 0; i < s;i++)
			a[i] = at(i) / a[i];
		return a;
	}
	dataVec operator + (const dataVec & a) const
	{
		size_t s = size();
		dataVec retVal(s);
		for (size_t i = 0; i < s;i++)
			retVal[i] = at(i) + a[i];
		return retVal;
	}
	dataVec operator + (dataVec && a) const
	{
		size_t s = size();
		for (size_t i = 0; i < s;i++)
			a[i] += at(i);
		return a;
	}

	dataVec operator - (const dataVec & a) const 
	{
		size_t s = size();
		dataVec retVal(s);
		for (size_t i = 0; i < s;i++)
			retVal[i] = at(i) - a[i];
		return retVal;
	}
	dataVec operator - () const
	{
		size_t s = size();
		dataVec retVal(s);
		for (size_t i = 0; i < s;i++)
			retVal[i] = -at(i);
		return retVal;
	}

	bool does_not_contain_null()
	{
		for (double & i : *this)
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

#define RANDOM
	void selectAtMost(size_t s)
	{
		if (s > size())
			return;
		//	Swap the first s entries with the entry at some other position, then resize to just have the s entries
#ifdef RANDOM
		iterator i = begin();
		for (size_t j = 0;j < s;j++)
		{
			std::swap (*(i++),*(begin() + rand() % size()));
		}
#endif
		resize(s);
	}

	dataVec & removeInvalidValues(double maxVal)
	{
		//	Sort in place for maximum efficiency, if we find an invalid value, replace with one from the end;
		//	Note that if we swap with a value from the end we have to check it as well to see if it is invalid
		iterator i = begin(), j = end();
		while (i != j)
		{
			if ((*i < 0) || (*i >= maxVal))
				std::swap(*i,*--j);
			else
				i++;
		}
		resize(j - begin());
		return *this;
	}

	dataVec chol();

	dataVec operator () (const vector<int> & i) const
	{
		dataVec retVal(i.size());
		for (size_t j = 0;j < i.size();j++)
			retVal[j] = at(i[j]);
		return retVal;
	}

};

class dataArray : public std::vector<dataVec>
{
public:
	dataArray(size_t s): std::vector<dataVec>(s,dataVec(s)){};
};



inline double sum(const dataVec & a)
{
	double retVal = 0;
	for (auto & i : a)
		retVal += i;
	return retVal;
}

inline dataVec exp(dataVec && a)
{
	for (double & i : a)
		i = exp(i);
	return a;
}
inline dataVec log(dataVec && a)
{
	for (double & i : a)
		i = log(i);
	return a;
}

inline dataVec operator * (double a, const dataVec & b)
{
	dataVec retVal(b.size());
	for (size_t i = 0; i < b.size();i++)
		retVal[i] = b[i] * a;
	return retVal;
}

inline dataVec operator * (double a, dataVec && b)
{
	for (double & i : b)
		i *= a;
	return b;
}

inline dataVec operator < (double a, const dataVec & b)
{
	dataVec retVal(b.size());
	for (size_t i = 0; i < b.size();i++)
		retVal[i] = a < b[i]?1:0;
	return retVal;
}
inline dataVec operator < (double a, dataVec && b)
{
	for (double & i : b)
		i = a < i?1:0;
	return b;
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
	for (double & i : b)
		i = a - i;
	return b;

}
inline dataVec operator + (dataVec && a,double b)
{
	for (double & i : a)
		i += b;
	return a;
}
inline dataVec operator + (const dataVec & a ,double b)
{
	dataVec retVal(a.size());
	for (size_t i = 0; i < a.size();i++)
		retVal[i] = a[i] + b;
	return retVal;
}

inline dataVec operator - (dataVec && a,double b)
{
	for (double & i : a)
		i -= b;
	return a;
}
inline dataVec operator - (const dataVec & a ,double b)
{
	dataVec retVal(a.size());
	for (size_t i = 0; i < a.size();i++)
		retVal[i] = a[i] - b;
	return retVal;
}
inline dataVec operator * (dataVec && a,double b)
{
	for (double &  i : a)
		i *= b;
	return a;
}
inline dataVec operator * (const dataVec & a ,double b)
{
	dataVec retVal(a.size());
	for (size_t i = 0; i < a.size();i++)
		retVal[i] = a[i] * b;
	return retVal;
}

inline dataVec operator / (dataVec && a,double b)
{
	for (double &  i : a)
		i /= b;
	return a;
}
inline dataVec operator / (const dataVec & a ,double b)
{
	dataVec retVal(a.size());
	for (size_t i = 0; i < a.size();i++)
		retVal[i] = a[i] / b;
	return retVal;
}

inline dataVec operator ^ (const dataVec & a,int b)
{
	dataVec _ret(a);
	for (int p = 0;p < (b - 1);p++)
		for (double & i : _ret)
			i *= i;
	return _ret;
}

inline dataVec operator ^ (dataVec && a,int b)
{
	for (int p = 0;p < (b - 1);p++)
		for (double & i : a)
			i *= i;
	return a;
}

inline dataVec operator > (dataVec && a,double b)
{
	for (double &  i : a)
		i = i>b?1:0;
	return a;
}
/*
double sum(const dataVec & a);
dataVec exp(dataVec && a);
dataVec log(dataVec && a);
dataVec operator * (double a, const dataVec & b);
dataVec operator * (double a, dataVec && b);
dataVec operator < (double a, const dataVec & b);
dataVec operator < (double a, dataVec && b);
dataVec operator - (double a, const dataVec & b);
dataVec operator - (double a, dataVec && b);
dataVec operator + (dataVec && a,double b);
dataVec operator + (const dataVec & a ,double b);
dataVec operator - (dataVec && a,double b);
dataVec operator - (const dataVec & a ,double b);
dataVec operator * (dataVec && a,double b);
dataVec operator * (const dataVec & a ,double b);
dataVec operator / (dataVec && a,double b);
dataVec operator / (const dataVec & a ,double b);
dataVec operator ^ (dataVec && a,int b);
dataVec operator > (dataVec && a,double b);

*/

inline bool printVal(outputDataFile * f,const dataVec & value)
{
	printVal(f,(const std::vector<double> &) value);
	return true;
};

class dataType 
{
public:
	dataVec fragData;
	std::vector<int> geneIndex;
	dataVec geneData[2];
};

#endif
