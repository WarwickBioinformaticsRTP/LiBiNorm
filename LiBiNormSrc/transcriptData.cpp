
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include "parser.h"
#include "containerEx.h"
#include "transcriptData.h"

void transcriptDataMap::histc (const vector<int> E)
{
	freq.assign (E.size(),0);
	for (auto & i : *this)
	{
		double v = i.second.length;
		size_t l = 0;
		size_t h = E.size()-1;
		size_t k = l;
		while((h-l)>1)
		{
			k = (h+l)/2;
			if(v < E[k])
				h = k;
			else
				l = k;
		}
		if(v == E[h])
			k = h;
		else
			k = l;
		i.second.histoGram_ind = k;
		freq[k]++;
	}
}



void transcriptDataMap::remove_invalid_values()
{
	for (auto & i : *this)
		i.second.remove_invalid_values();
}

int transcriptDataMap::loadData(const string filename)
{
	std::ifstream f;
	f.open(filename);
	if (!f.is_open())
		exitFail("Unable to open file ", filename);

	string buffer,gene,direction;
	int a = 1;
	while (!f.eof())
	{
		//	Need to check for gene and length consistency and that the gene has not appeared before
		transcriptData tempData;

		getline(f,buffer);
		parseTsv(buffer,gene,tempData.length," ",direction,tempData.counts[0].values());
		getline(f,buffer);
		parseTsv(buffer,gene,tempData.length," ",direction,tempData.counts[1].values());
		emplace(gene,move(tempData));
		a++;
	}

	return EXIT_SUCCESS;
}

void transcriptDataMap::transferTo(dataType & mcmcData,size_t maxLength)
//void dataType::consolidateWith(transcriptDataMap & transData,size_t maxLength)
{
	vectorEx<int> bins(0,300);
	for (size_t i = 500;i <= 10000;i++)
		bins.push_back(i);
	bins.add(11000,12000,13000,15000,30000);


	histc(bins);

	for (auto & gene : *this)
	{
		for (auto & counts : gene.second.counts)
		{
			counts.selectAtMost(maxLength);

			mcmcData[0].append(counts);
			mcmcData[1].append(counts.size(),gene.second.length);
			mcmcData[1].append(counts.size(),freq[gene.second.histoGram_ind]);
		}

	}

}

