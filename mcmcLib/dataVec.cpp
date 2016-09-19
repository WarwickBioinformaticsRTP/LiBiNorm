#include "dataVec.h"
#include <cmath>
#include <thread>
#include <mutex>
#include <map>

using namespace std;
#define DATAVEC_CACHE



#ifdef DATAVEC_CACHE

#define CACHE_MINSIZE 500

map<thread::id,map<size_t,vector<vector<double> > > > cache;
mutex assign_vector_mtx; 

dataVec::dataVec(size_t s)
{
	if (s == 0)
		return;

	vector<vector<double> > & c = cache[this_thread::get_id()][s];
//	lock_guard<mutex> lock(assign_vector_mtx);
	if (c.size() && (s >= CACHE_MINSIZE))
	{
		swap(c.back());
		c.pop_back();
	}
	else
		resize(s);
};

dataVec::dataVec(const std::vector<double> & a)
{
//	lock_guard<mutex> lock(assign_vector_mtx);
	size_t s = a.size();
	vector<vector<double> > & c = cache[this_thread::get_id()][s];
	if (c.size() && (s >= CACHE_MINSIZE))
	{
		swap(c.back());
		c.pop_back();
		assign(a.begin(),a.end());
	}
	else
		vector<double>::operator=(a);
};
dataVec::dataVec(const dataVec & a)
{
//	lock_guard<mutex> lock(assign_vector_mtx);
	size_t s = a.size();
	vector<vector<double> > & c = cache[this_thread::get_id()][s];
	if (c.size() && (s >= CACHE_MINSIZE))
	{
		swap(c.back());
		c.pop_back();
		assign(a.begin(),a.end());
	}
	else
		vector<double>::operator=(a);
};


dataVec::~dataVec()
{
	size_t s = size();
	vector<vector<double> > & c = cache[this_thread::get_id()][s];
	if (s < CACHE_MINSIZE)
		return;
//	lock_guard<mutex> lock(assign_vector_mtx);
	c.emplace_back(move(*this));
}

#else
dataVec::dataVec(size_t s): std::vector<double>(s){};
dataVec::dataVec(const std::vector<double> & a): std::vector<double>(a){};
dataVec::~dataVec(){};
#endif


//	Cholesky_Decomposition returns the Cholesky Decomposition Matrix. 
dataVec dataVec::chol()
{

	dataArray L(size());

	//	Initialize and populate matrix L which will be the lower Cholesky
	for (size_t i = 0; i < size(); i++)
		for (size_t j = 0; j < size(); j++)
		{
			double temp = 0;// , temp2 = 0;
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

