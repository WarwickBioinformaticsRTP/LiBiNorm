#ifndef TRANSDATA_H
#define TRANSDATA_H

#include <map>
#include "mcmc.h"
#include "geneCountData.h"


using namespace std;
/*
class transcriptData
{
public:
	transcriptData(transcriptData && a) : positions(move(a.positions)),length(a.length),gene(a.gene)
	{
	};

	void remove_invalid_values()
	{
		for (auto & i : positions)
			i.removeInvalidValues(length);
	}

	vector<dataVec > positions;
	int length;
	std::string gene;
	int histoGram_ind;
	transcriptData() : positions(2), length(0) {	};

};



class transcriptDataMap : public vector<transcriptData>
{
public:
//	std::string loadData(const string filename, GeneCountData & countData, int Nlines = -1);
//	void remove_invalid_values();
//	void histc (const vector<int> E);
//	void transferTo(dataType & mcmcData,size_t maxLength,int maxTotReads);

//	vector<double> freq;
};
*/
#endif


