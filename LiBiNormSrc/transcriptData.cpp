
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include "libCommon.h"
#include "parser.h"
#include "containerEx.h"
#include "LiBiNorm.h"


//	Find the number of read position values in 'this' (a vector) that sits within each bin of the
//	histogram defined by E.   The results go into the 'freq'vector
//	This is used as part of the liklyhood calculations
void GeneCountData::histc (const vector<int> E)
{
	freq.assign (E.size(),0);
	histoGram_ind.assign(lengths.size(),-1);
	for (size_t i = 0;i < names.size();i++)
	{
		if (names[i].substr(0, 2) != "__")
		{
			double v = lengths[i];
			size_t l = 0;
			size_t h = E.size() - 1;
			size_t k = l;
			while ((h - l) > 1)
			{
				k = (h + l) / 2;
				if (v < E[k])
					h = k;
				else
					l = k;
			}
			if (v == E[h])
				k = h;
			else
				k = l;
			histoGram_ind[i] = k;
			freq[k]++;
		}
	}
}

void GeneCountData::remove_invalid_values()
{
	for (size_t i = 0;i < names.size();i++)
	{
		for (size_t j = 0;j < 2;j++)
			This[names[i]].positions[j].removeInvalidValues(lengths[i]);
	}
}

string GeneCountData::loadData(const string filename, int Ngenes)
{
	std::ifstream f;
	f.open(filename);
	if (!f.is_open())
		exitFail("Unable to open file ", filename);

	//The reference value
	addEntry("reference", DEFAULT_NORMALISATION_GENE_LENGTH);


	string buffer,gene,direction;
	int a = 1;
	string lastGene;
	while (!f.eof() & (Ngenes-- != 0))
	{
		//	Need to check for gene and length consistency and that the gene has not appeared before

		getline(f,buffer);
		string gene, n1, n2;
		if (buffer.size())
		{
			rnaPosVec posPositions,negPositions;

			parseTsv(buffer,gene,n1," ",direction,(std::vector<rna_pos_type> &)posPositions);
			getline(f,buffer);
			parseTsv(buffer,gene,n2," ",direction, (std::vector<rna_pos_type> &)negPositions);

			int length, count;

			if (strchr(n1.c_str(), ':') == NULL)
			{
				length = atoi(n1.c_str());
				count = posPositions.size() + negPositions.size();
			}
			else
				parser(n1, ":", length, count);

			addEntry(gene, length, count, posPositions, negPositions);
			a++;
		}
		lastGene = gene;
	}

	progMessage("Last gene = ", lastGene);
	return lastGene;
}

//	Transfers information for up to maxLength reads from up to Ngenes genes or transcripts
//	into the form which can be used by the mcmc chain
void GeneCountData::transferTo(dataType & mcmcData, size_t maxLength, int maxTotReads)
{
	vectorEx<int> bins{ { 0,300 } };
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

	//	Start at one because we dont include the reference gene (except there are no reads
	//	in the reference gene so this makes no difference
	for (size_t i = 1;i < names.size();i++)
	{
		if ((lengths[i] < MAX_LENGTH_OF_GENE_FOR_PARAM_ESTIMATION) && (names[i].substr(0,2) != "__"))
		{
			//	For the forward and the reverse counts
			for (size_t j = 0;j < 2;j++)
			{
				rnaPosVec & positions = at(names[i]).positions[j];
				positions.selectAtMost(maxLength);

				//	fragData contains the count 
				mcmcData.fragData.append(conv(positions));

				mcmcData.geneIndex.insert(mcmcData.geneIndex.end(), positions.size(), geneIndex);//gene.length);
				Nreads += positions.size();
			}
			mcmcData.geneData[0][geneIndex] = lengths[i];
			mcmcData.geneData[1][geneIndex] = freq[histoGram_ind[i]];

			geneIndex++;
			if ((maxTotReads) && (Nreads > maxTotReads))
				break;
		}
	}
	mcmcData.geneData[0].resize(geneIndex);
	mcmcData.geneData[1].resize(geneIndex);

	optMessage(Nreads, " reads used for parameter determination");
}

