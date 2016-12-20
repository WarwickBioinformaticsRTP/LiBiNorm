#ifndef TRANSDATA_H
#define TRANSDATA_H

#include <map>
#include "mcmc.h"

using namespace std;

class transcriptData
{
public:
	vector<dataVec > positions;
	int length;
	std::string gene;
	int histoGram_ind;
	transcriptData(): positions(2),length(0) {	};

	transcriptData(transcriptData && a) : positions(move(a.positions)),length(a.length),gene(a.gene)
	{
	};

	void remove_invalid_values()
	{
		for (auto & i : positions)
			i.removeInvalidValues(length);
	}

};



class transcriptDataMap : public vector<transcriptData>
{
public:
	vector<double> freq;

	void remove_invalid_values();
	std::string loadData(const string filename,int Nlines = -1);
	void histc (const vector<int> E);
	void transferTo(dataType & mcmcData,size_t maxLength,int maxTotReads);


};

#endif


