#ifndef DATAVEC_H
#define DATAVEC_H

#include <stdlib.h>
#include <vector>
#include <string>
#include <valarray>
#include "printEx.h"

#define VEC_DATA_TYPE double

class dataVec: public std::vector<VEC_DATA_TYPE>
{
public:
	dataVec(size_t s = 0);
	void resize(size_t s = 0);
	dataVec(size_t s,VEC_DATA_TYPE v);
	dataVec(const std::vector<VEC_DATA_TYPE> & a);
	dataVec(const dataVec & a);

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

//	Test code takes the first N samples rather than randomly picks samples, and uses the same 
//	algorithm as the MATLAB code for excluding invalid calues.   Used for comparing the two outputs
// #define TEST_CODE
	void selectAtMost(size_t s)
	{
		if (s > size())
			return;
		//	Swap the first s entries with the entry at some other position, then resize to just have the s entries
#ifndef TEST_CODE
		iterator i = begin();
		for (size_t j = 0;j < s;j++)
		{
			std::swap (*(i++),*(begin() + rand() % size()));
		}
#endif
		resize(s);
	}

	dataVec & removeInvalidValues(VEC_DATA_TYPE maxVal)
	{
		//	Sort in place for maximum efficiency, if we find an invalid value, replace with one from the end;
		//	Note that if we swap with a value from the end we have to check it as well to see if it is invalid
#ifdef TEST_CODE
		size_t offset = 0;
		for (size_t i = 0;i+offset < size();)
		{
			if ((at(i+offset) < 0) || (at(i+offset) >= maxVal))
				offset++;
			else
			{
				if (offset)
					at(i) = at(i+offset);
				i++;
			}
		}
		if (offset)
			resize(size()-offset);
#else
		iterator i = begin(), j = end();
		while (i != j)
		{
			if ((*i < 0) || (*i >= maxVal))
				std::swap(*i,*--j);
			else
				i++;
		}
		resize(j - begin());
#endif
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
//https://codingforspeed.com/using-faster-exponential-approximation/

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

class dataType 
{
public:
	dataVec fragData;
	std::vector<int> geneIndex;
	dataVec geneData[2];
};

#endif
