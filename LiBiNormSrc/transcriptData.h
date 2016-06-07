#ifndef TRANSDATA_H
#define TRANSDATA_H

#include <map>
#include "dataVec.h"

using namespace std;

class transcriptData
{
public:
	int length;
	int histoGram_ind;
	vector<dataVec > counts;
	transcriptData(): counts(2),length(0) {	};

	transcriptData(transcriptData && a) : counts(move(a.counts)),length(a.length)
	{
	};

	void remove_invalid_values()
	{
		for (auto & i : counts)
			i.removeInvalidValues(length);
	}

};



class transcriptDataMap : public map<string,transcriptData>
{
public:
	vector<int> freq;

	void remove_invalid_values();
	int loadData(const string filename);
	void histc (const vector<int> E);
	void transferTo(dataType & mcmcData,size_t maxLength);


};

#endif


