#include "dataVec.h"
#include <cmath>

using namespace std;

//	Cholesky_Decomposition returns the Cholesky Decomposition Matrix. 
dataVec dataVec::chol()
{

	dataArray L(size());

	//	Initialize and populate matrix L which will be the lower Cholesky
	for (size_t i = 0; i < size(); i++)
		for (size_t j = 0; j < size(); j++)
		{
			double temp = 0, temp2 = 0;
			if (i > j)
			{
//				if (j > 0)
//			{
//					for (k = 1; k < j + 1; k++)
//						temp2 += (L[i][k - 1] * L[j][k - 1]);
//				}
//				L[i][j] = (p[i][j] - temp2) / L[j][j];
			}
			else if (i == j)
			{
				for (size_t k = 0; k < i; k++)
					temp += pow(L[i][k], 2);
				L[i][j] = sqrt(at(i) - temp);
			}
//			else
//				L[i][j] = 0;
		}
	dataVec M(size());
	for (size_t i = 0; i < size(); i++)
		M[i] = L[i][i];

	return M;
};

/*
double sum(const dataVec & a)
{
	double retVal = 0;
	for (auto & i : a)
		retVal += i;
	return retVal;
}

dataVec exp(dataVec && a)
{
	for (size_t i = 0; i < a.size();i++)
		a[i] = exp(a[i]);
	return a;
}
dataVec log(dataVec && a)
{
	for (size_t i = 0; i < a.size();i++)
		a[i] = log(a[i]);
	return a;
}

dataVec operator * (double a, const dataVec & b)
{
	dataVec retVal(b.size());
	for (size_t i = 0; i < b.size();i++)
		retVal[i] = b[i] * a;
	return retVal;
}

dataVec operator * (double a, dataVec && b)
{
	for (size_t i = 0; i < b.size();i++)
		b[i] *= a;
	return b;
}

dataVec operator < (double a, const dataVec & b)
{
	dataVec retVal(b.size());
	for (size_t i = 0; i < b.size();i++)
		retVal[i] = a < b[i]?1:0;
	return retVal;
}
dataVec operator < (double a, dataVec && b)
{
	for (size_t i = 0; i < b.size();i++)
		b[i] = a < b[i]?1:0;
	return b;
}

dataVec operator - (double a, const dataVec & b)
{
	dataVec retVal(b.size());
	for (size_t i = 0; i < b.size();i++)
		retVal[i] = a - b[i];
	return retVal;

}
dataVec operator - (double a, dataVec && b)
{
	for (size_t i = 0; i < b.size();i++)
		b[i] = a - b[i];
	return b;

}
dataVec operator + (dataVec && a,double b)
{
	for (size_t i = 0; i < a.size();i++)
		a[i] += b;
	return a;
}
dataVec operator + (const dataVec & a ,double b)
{
	dataVec retVal(a.size());
	for (size_t i = 0; i < a.size();i++)
		retVal[i] = a[i] + b;
	return retVal;
}

dataVec operator - (dataVec && a,double b)
{
	for (size_t i = 0; i < a.size();i++)
		a[i] -= b;
	return a;
}
dataVec operator - (const dataVec & a ,double b)
{
	dataVec retVal(a.size());
	for (size_t i = 0; i < a.size();i++)
		retVal[i] = a[i] - b;
	return retVal;
}
dataVec operator * (dataVec && a,double b)
{
	for (size_t i = 0; i < a.size();i++)
		a[i] *= b;
	return a;
}
dataVec operator * (const dataVec & a ,double b)
{
	dataVec retVal(a.size());
	for (size_t i = 0; i < a.size();i++)
		retVal[i] = a[i] * b;
	return retVal;
}

dataVec operator / (dataVec && a,double b)
{
	for (size_t i = 0; i < a.size();i++)
		a[i] /= b;
	return a;
}
dataVec operator / (const dataVec & a ,double b)
{
	dataVec retVal(a.size());
	for (size_t i = 0; i < a.size();i++)
		retVal[i] = a[i] / b;
	return retVal;
}


dataVec operator ^ (dataVec && a,int b)
{
	for (int p = 0;p < b;p++)
		for (size_t i = 0; i < a.size();i++)
			a[i] *= a[i];
	return a;
}

dataVec operator > (dataVec && a,double b)
{
	for (size_t i = 0; i < a.size();i++)
		a[i] = a[i]>b?1:0;
	return a;
}

*/
