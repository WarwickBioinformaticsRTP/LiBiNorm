
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include "parser.h"
#include "containerEx.h"
#include "transcriptData.h"


#define MAX_LENGTH_OF_GENE_FOR_PARAM_ESTIMATION 20000

//	Find the number of read position values in 'this' (a vector) that sits within each bin of the
//	histogram defined by E.   The results go into the 'freq'vector
//	This is used as part of the liklyhood calculations
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

string transcriptDataMap::loadData(const string filename, int Nlines)
{
	std::ifstream f;
	f.open(filename);
	if (!f.is_open())
		exitFail("Unable to open file ", filename);

	string buffer,gene,direction;
	int a = 1;
	string lastGene;
	while (!f.eof() & (Nlines-- != 0))
	{
		//	Need to check for gene and length consistency and that the gene has not appeared before
		transcriptData tempData;

		getline(f,buffer);
		if (buffer.size())
		{
			parseTsv(buffer,tempData.gene,tempData.length," ",direction,tempData.positions[0].values());
			getline(f,buffer);
			parseTsv(buffer,tempData.gene,tempData.length," ",direction,tempData.positions[1].values());
			push_back(move(tempData));

			a++;
		}
		lastGene = tempData.gene;
	}

	progMessage("Last gene = ", lastGene);
	return lastGene;
}

//	Transfers information for up to maxLength reads from up to Ngenes genes or transcripts
//	into the form which can be used by the mcmc chain
void transcriptDataMap::transferTo(dataType & mcmcData,size_t maxLength, int maxTotReads)
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

	int Nreads = 0;

	for (auto & gene : This)
	{
		if (gene.length < MAX_LENGTH_OF_GENE_FOR_PARAM_ESTIMATION)
		{
			//	For the forward and the reverse counts
			for (auto & positions : gene.positions)
			{
				positions.selectAtMost(maxLength);

				//	fragData contains the count 
				mcmcData.fragData.append(positions);

				mcmcData.geneIndex.insert(mcmcData.geneIndex.end(), positions.size(), geneIndex);//gene.length);
				Nreads += positions.size();
			}
			mcmcData.geneData[0][geneIndex] = gene.length;
			mcmcData.geneData[1][geneIndex] = freq[gene.histoGram_ind];

			geneIndex++;
			if ((maxTotReads) && (Nreads > maxTotReads))
				break;
		}
	}
	optMessage(Nreads, " reads used for parameter determination");
}

