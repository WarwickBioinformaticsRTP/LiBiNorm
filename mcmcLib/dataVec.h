#ifndef DATAVEC_H
#define DATAVEC_H

#include <vector>
#include <string>

class dataVec: public std::vector<double>
{
public:
	dataVec(size_t s = 0): std::vector<double>(s){};
	dataVec(const std::vector<double> & a): std::vector<double>(a){};

	vector<double> & values() { return *this;};

	dataVec operator > (double a) const
	{
		dataVec retVal(size());
		for (size_t i = 0; i < size();i++)
			retVal[i] = at(i)>a?1:0;
		return retVal;
	}
	dataVec operator < (const dataVec & a) const
	{
		dataVec retVal(size());
		for (size_t i = 0; i < size();i++)
			retVal[i] = at(i)<a[i]?1:0;
		return retVal;
	}
	dataVec operator < (dataVec && a) const
	{
		for (size_t i = 0; i < size();i++)
			a[i] = at(i)<a[i]?1:0;
		return a;
	}
	dataVec operator * (const dataVec & a) const
	{
		dataVec retVal(size());
		for (size_t i = 0; i < size();i++)
			retVal[i] = at(i) * a[i];
		return retVal;
	}
	dataVec operator * (dataVec && a) const
	{
		for (size_t i = 0; i < size();i++)
			a[i] *= at(i);
		return a;
	}
	dataVec operator / (const dataVec & a) const 
	{
		dataVec retVal(size());
		for (size_t i = 0; i < size();i++)
			retVal[i] = at(i) / a[i];
		return retVal;
	}
	dataVec operator / (dataVec && a) const 
	{
		for (size_t i = 0; i < size();i++)
			a[i] = at(i) / a[i];
		return a;
	}
	dataVec operator + (const dataVec & a) const
	{
		dataVec retVal(size());
		for (size_t i = 0; i < size();i++)
			retVal[i] = at(i) + a[i];
		return retVal;
	}
	dataVec operator + (dataVec && a) const
	{
		for (size_t i = 0; i < size();i++)
			a[i] += at(i);
		return a;
	}

	dataVec operator - (const dataVec & a) const 
	{
		dataVec retVal(size());
		for (size_t i = 0; i < size();i++)
			retVal[i] = at(i) - a[i];
		return retVal;
	}
	dataVec operator - () const
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

	dataVec chol();

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
	for (size_t i = 0; i < a.size();i++)
		a[i] = exp(a[i]);
	return a;
}
inline dataVec log(dataVec && a)
{
	for (size_t i = 0; i < a.size();i++)
		a[i] = log(a[i]);
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
	for (size_t i = 0; i < b.size();i++)
		b[i] *= a;
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
	for (size_t i = 0; i < b.size();i++)
		b[i] = a < b[i]?1:0;
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
	for (size_t i = 0; i < b.size();i++)
		b[i] = a - b[i];
	return b;

}
inline dataVec operator + (dataVec && a,double b)
{
	for (size_t i = 0; i < a.size();i++)
		a[i] += b;
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
	for (size_t i = 0; i < a.size();i++)
		a[i] -= b;
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
	for (size_t i = 0; i < a.size();i++)
		a[i] *= b;
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
	for (size_t i = 0; i < a.size();i++)
		a[i] /= b;
	return a;
}
inline dataVec operator / (const dataVec & a ,double b)
{
	dataVec retVal(a.size());
	for (size_t i = 0; i < a.size();i++)
		retVal[i] = a[i] / b;
	return retVal;
}


inline dataVec operator ^ (dataVec && a,int b)
{
	for (int p = 0;p < b;p++)
		for (size_t i = 0; i < a.size();i++)
			a[i] *= a[i];
	return a;
}

inline dataVec operator > (dataVec && a,double b)
{
	for (size_t i = 0; i < a.size();i++)
		a[i] = a[i]>b?1:0;
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

class dataType 
{
	dataVec data[3];

public:
	dataVec & operator [] (int i) { return data[i];};
	const dataVec & operator [] (int i) const { return data[i];};

};

#endif