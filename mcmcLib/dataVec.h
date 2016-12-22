#ifndef DATAVEC_H
#define DATAVEC_H

#include <stdlib.h>
#include <vector>
#include <string>
#include <valarray>
#include "libCommon.h"
#include "printEx.h"

#define VEC_DATA_TYPE double

class dataVec: public std::vector<VEC_DATA_TYPE>
{
public:
	void resize(size_t s = 0);
	dataVec(size_t s = 0,VEC_DATA_TYPE v = 0);
	dataVec(const std::vector<VEC_DATA_TYPE> & a);
	dataVec(const dataVec & a);
	dataVec(std::initializer_list<VEC_DATA_TYPE> a) : std::vector<VEC_DATA_TYPE>(a) {};

	static void clearCache();

	~dataVec();

	vector<VEC_DATA_TYPE> & values() { return *this;};

	dataVec operator > (VEC_DATA_TYPE a) const
	{
		size_t s = size();
		dataVec retVal(s);
		for (size_t i = 0; i < s;i++)
			retVal.at(i) = at(i)>a?1:0;
		return retVal;
	}
	dataVec operator < (const dataVec & a) const
	{
		size_t s = size();
		dataVec retVal(s);
		for (size_t i = 0; i < s;i++)
			retVal.at(i) = at(i)<a.at(i)?1:0;
		return retVal;
	}
	dataVec operator < (dataVec && a) const
	{
		size_t s = size();
		for (size_t i = 0; i < s;i++)
			a.at(i) = at(i)<a.at(i)?1:0;
		return a;
	}
	dataVec operator * (const dataVec & a) const
	{
		size_t s = size();
		dataVec retVal(s);
		for (size_t i = 0; i < s;i++)
			retVal.at(i) = at(i) * a.at(i);
		return retVal;
	}
	dataVec operator * (dataVec && a) const
	{
		size_t s = size();
		for (size_t i = 0; i < s;i++)
			a.at(i) *= at(i);
		return a;
	}
	dataVec operator / (const dataVec & a) const 
	{
		size_t s = size();
		dataVec retVal(s);
		for (size_t i = 0; i < s;i++)
			retVal.at(i) = at(i) / a.at(i);
		return retVal;
	}
	dataVec operator / (dataVec && a) const 
	{
		size_t s = size();
		for (size_t i = 0; i < s;i++)
			a.at(i) = at(i) / a.at(i);
		return a;
	}
	dataVec operator + (const dataVec & a) const
	{
		size_t s = size();
		dataVec retVal(s);
		for (size_t i = 0; i < s;i++)
			retVal.at(i) = at(i) + a.at(i);
		return retVal;
	}
	dataVec operator + (dataVec && a) const
	{
		size_t s = size();
		for (size_t i = 0; i < s;i++)
			a.at(i) += at(i);
		return a;
	}

	dataVec operator - (const dataVec & a) const 
	{
		size_t s = size();
		dataVec retVal(s);
		for (size_t i = 0; i < s;i++)
			retVal.at(i) = at(i) - a.at(i);
		return retVal;
	}
	dataVec operator - () const
	{
		size_t s = size();
		dataVec retVal(s);
		for (size_t i = 0; i < s;i++)
			retVal.at(i) = -at(i);
		return retVal;
	}

	bool does_not_contain_null()
	{
		for (VEC_DATA_TYPE & i : *this)
			if (i == 0)
				return false;
		return true;
	}

	dataVec & append(const dataVec a)
	{
		insert(end(),a.begin(),a.end());
		return *this;
	}
	dataVec & append(size_t N,VEC_DATA_TYPE value)
	{
		insert(end(),N,value);
		return *this;
	}

	dataVec diagchol();

	dataVec operator () (const vector<int> & i) const
	{
		size_t s(i.size());
		dataVec retVal(s);
		for (size_t j = 0;j < s;j++)
			retVal.at(j) = at(i.at(j));
		return retVal;
	}

};

class dataArray : public std::vector<dataVec>
{
public:
	dataArray(size_t s): std::vector<dataVec>(s,dataVec(s)){};
};



inline VEC_DATA_TYPE sum(const dataVec & a)
{
	VEC_DATA_TYPE retVal = 0;
	for (auto & i : a)
		retVal += i;
	return retVal;
}

/*

This is some code for a fast version of the exponent function that 
was found at.   It did not prove to be make much difference so is not currently being used
https://codingforspeed.com/using-faster-exponential-approximation/

inline
double exp1(double x) {
  x = 1.0 + x / 256.0;
  x *= x; x *= x; x *= x; x *= x;
  x *= x; x *= x; x *= x; x *= x;
  return x;
}
*/
inline dataVec exp(dataVec && a)
{
	for (VEC_DATA_TYPE & i : a)
	{
		i = exp(i);
	}
	return a;
}
inline dataVec log(dataVec && a)
{
	for (VEC_DATA_TYPE & i : a)
		i = log(i);
	return a;
}
inline dataVec sqrt(dataVec && a)
{
	for (VEC_DATA_TYPE & i : a)
		i = sqrt(i);
	return a;
}
inline dataVec operator * (VEC_DATA_TYPE a, const dataVec & b)
{
	size_t s(b.size());
	dataVec retVal(s);
	for (size_t i = 0; i < s;i++)
		retVal.at(i) = b.at(i) * a;
	return retVal;
}

inline dataVec operator * (VEC_DATA_TYPE a, dataVec && b)
{
	for (VEC_DATA_TYPE & i : b)
		i *= a;
	return b;
}

inline dataVec operator < (VEC_DATA_TYPE a, const dataVec & b)
{
	size_t s(b.size());
	dataVec retVal(s);
	for (size_t i = 0; i < s;i++)
		retVal.at(i) = a < b.at(i)?1:0;
	return retVal;
}
inline dataVec operator < (VEC_DATA_TYPE a, dataVec && b)
{
	for (VEC_DATA_TYPE & i : b)
		i = a < i?1:0;
	return b;
}

inline dataVec operator - (VEC_DATA_TYPE a, const dataVec & b)
{
	size_t s(b.size());
	dataVec retVal(s);
	for (size_t i = 0; i < s;i++)
		retVal.at(i) = a - b.at(i);
	return retVal;

}
inline dataVec operator - (VEC_DATA_TYPE a, dataVec && b)
{
	for (VEC_DATA_TYPE & i : b)
		i = a - i;
	return b;

}
inline dataVec operator + (dataVec && a,VEC_DATA_TYPE b)
{
	for (VEC_DATA_TYPE & i : a)
		i += b;
	return a;
}
inline dataVec operator + (const dataVec & a ,VEC_DATA_TYPE b)
{
	size_t s(a.size());
	dataVec retVal(s);
	for (size_t i = 0; i < s;i++)
		retVal.at(i) = a.at(i) + b;
	return retVal;
}

inline dataVec operator += (dataVec & a, const dataVec & b)
{
	_ASSERT_EXPR(a.size() == b.size(), "/= vector sizes do not match");
	for (size_t i = 0; i < a.size(); i++)
		a.at(i) += b.at(i);
	return a;
}

inline dataVec operator - (dataVec && a,VEC_DATA_TYPE b)
{
	for (VEC_DATA_TYPE & i : a)
		i -= b;
	return a;
}
inline dataVec operator - (const dataVec & a ,VEC_DATA_TYPE b)
{
	size_t s(a.size());
	dataVec retVal(s);
	for (size_t i = 0; i < s;i++)
		retVal.at(i) = a.at(i) - b;
	return retVal;
}
inline dataVec operator * (dataVec && a,VEC_DATA_TYPE b)
{
	for (VEC_DATA_TYPE &  i : a)
		i *= b;
	return a;
}
inline dataVec operator * (const dataVec & a ,VEC_DATA_TYPE b)
{
	size_t s(a.size());
	dataVec retVal(s);
	for (size_t i = 0; i < s;i++)
		retVal.at(i) = a.at(i) * b;
	return retVal;
}

inline dataVec operator / (dataVec && a,VEC_DATA_TYPE b)
{
	for (VEC_DATA_TYPE &  i : a)
		i /= b;
	return a;
}
inline dataVec operator /= (dataVec & a, VEC_DATA_TYPE b)
{
	for (VEC_DATA_TYPE & i : a)
		i /= b;
	return a;
}
inline dataVec operator /= (dataVec & a, const dataVec & b)
{
	_ASSERT_EXPR(a.size() == b.size(), "/= vector sizes do not match");
	for (size_t i = 0;i < a.size();i++)
		a.at(i) /= b.at(i);
	return a;
}
inline dataVec operator / (const dataVec & a ,VEC_DATA_TYPE b)
{
	size_t s(a.size());
	dataVec retVal(s);
	for (size_t i = 0; i < s;i++)
		retVal.at(i) = a.at(i) / b;
	return retVal;
}

inline dataVec operator ^ (const dataVec & a,int b)
{
	dataVec _ret(a);
	for (int p = 0;p < (b - 1);p++)
		for (VEC_DATA_TYPE & i : _ret)
			i *= i;
	return _ret;
}

inline dataVec operator ^ (dataVec && a,int b)
{
	for (int p = 0;p < (b - 1);p++)
		for (VEC_DATA_TYPE & i : a)
			i *= i;
	return a;
}

inline dataVec operator > (dataVec && a,VEC_DATA_TYPE b)
{
	for (VEC_DATA_TYPE &  i : a)
		i = i>b?1:0;
	return a;
}
/*
VEC_DATA_TYPE sum(const dataVec & a);
dataVec exp(dataVec && a);
dataVec log(dataVec && a);
dataVec operator * (VEC_DATA_TYPE a, const dataVec & b);
dataVec operator * (VEC_DATA_TYPE a, dataVec && b);
dataVec operator < (VEC_DATA_TYPE a, const dataVec & b);
dataVec operator < (VEC_DATA_TYPE a, dataVec && b);
dataVec operator - (VEC_DATA_TYPE a, const dataVec & b);
dataVec operator - (VEC_DATA_TYPE a, dataVec && b);
dataVec operator + (dataVec && a,VEC_DATA_TYPE b);
dataVec operator + (const dataVec & a ,VEC_DATA_TYPE b);
dataVec operator - (dataVec && a,VEC_DATA_TYPE b);
dataVec operator - (const dataVec & a ,VEC_DATA_TYPE b);
dataVec operator * (dataVec && a,VEC_DATA_TYPE b);
dataVec operator * (const dataVec & a ,VEC_DATA_TYPE b);
dataVec operator / (dataVec && a,VEC_DATA_TYPE b);
dataVec operator / (const dataVec & a ,VEC_DATA_TYPE b);
dataVec operator ^ (dataVec && a,int b);
dataVec operator > (dataVec && a,VEC_DATA_TYPE b);

*/

inline bool printVal(outputDataFile * f,const dataVec & value)
{
	printVal(f,(const std::vector<VEC_DATA_TYPE> &) value);
	return true;
};

#endif
