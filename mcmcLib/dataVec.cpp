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



