
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
		double v = i.length;
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
		i.histoGram_ind = k;
		freq[k]++;
	}
}



void transcriptDataMap::remove_invalid_values()
{
	for (auto & i : This)
		i.remove_invalid_values();
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
		if (buffer.size())
		{
			parseTsv(buffer,tempData.gene,tempData.length," ",direction,tempData.counts[0].values());
			getline(f,buffer);
			parseTsv(buffer,tempData.gene,tempData.length," ",direction,tempData.counts[1].values());
			push_back(move(tempData));

			a++;
		}
	}

	return EXIT_SUCCESS;
}

void transcriptDataMap::transferTo(dataType & mcmcData,size_t maxLength)
//void dataType::consolidateWith(transcriptDataMap & transData,size_t maxLength)
{
	vectorEx<int> bins(0,300);
	for (size_t i = 500;i <= 10000;i+=500)
		bins.push_back(i);
	bins.add(11000,12000,15000,30000);


	histc(bins);

	freq[0] = freq[0]*2;
	freq[1] = freq[1]*2;
	freq[21] = freq[21]/2;
	freq[22] = freq[22]/2;
	freq[23] = freq[23]/6;
	freq[24] = freq[24]/6;


	size_t geneIndex = 0;
	mcmcData.geneData[0].resize(size());
	mcmcData.geneData[1].resize(size());


	srand( (unsigned)time( NULL ) );

	for (auto & gene : *this)
	{
		for (auto & counts : gene.counts)
		{
			counts.selectAtMost(maxLength);

			mcmcData.fragData.append(counts);

			mcmcData.geneIndex.insert(mcmcData.geneIndex.end(),counts.size(),geneIndex);//gene.length);
		}
		mcmcData.geneData[0][geneIndex] = gene.length;
		mcmcData.geneData[1][geneIndex] = freq[gene.histoGram_ind];

		geneIndex++;
	}
}

