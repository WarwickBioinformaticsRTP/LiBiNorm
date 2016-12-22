#include <fstream>
#include "LiBiNorm.h"
#include "GeneCountData.h"
#include "stringEx.h"
#include "containerEx.h"
#include "parser.h"

using namespace std;

VEC_DATA_TYPE & GeneCountData::operator()(const std::string gene)
{
	static VEC_DATA_TYPE dummy;
	auto i = find(gene);
	if (i == end())
	{
		i = errorCounts.find(gene);
		if (i == errorCounts.end())
		{
			_ASSERT_EXPR(false, "Looking for gene that was not in the reference genome");
			return dummy;
		}
		return errCounts.at((*i).second.index);
	}
	return rawCounts.at((*i).second.index);
}



void GeneCountData::useSelectedGenes(const std::string & filename)
{
	ifstream file;
	file.open(filename);

	if (!file.is_open())
	{
		progMessage("Unable to read gene list from ", filename);
		return;
	}

	string line, gene;
	while (!file.eof())
	{
		std::getline(file, line);
		parser(line, " \n\r", line, gene);
		addEntry(gene);
	};
}

void GeneCountData::outputGeneCount(TsvFile & output, const string name)
{
	output.print(name, rawCounts[find(name)->second.index]);
}


bool GeneCountData::outputGeneCounts(const string & filename, bool withDetails)
{
	TsvFile output;

	if (!output.open(filename))
		return false;

	if (withDetails && (bias.size()))
		output.print("Gene", "Normalised count","Raw count","RNA length","Raw count","Bias");

	if (bias.size())
	{
		//	Dont start at 0 as 0 is the reference for normalisation
		for (size_t i = 1; i < names.size(); i++)
		{
			output.printStart(names[i], (int)floor(rawCounts[i] / ((lengths[i] == 0) ? 1 : bias[i]) + 0.5));

			if (withDetails)
				output.printMiddle(rawCounts[i], lengths[i], (lengths[i] == 0) ? 1 : bias[i]);

			output.printEnd();
		}
	}
	else
	{
		for (auto i : This)
		{
			if (i.second.index != 0)
			{
				output.print(i.first, (int)rawCounts[i.second.index]);
			}
		}
	}
	for (size_t i = 0;i < errorNames.size();i++)
		output.print(errorNames[i], (int)errCounts[i]);
	return true;
}

bool GeneCountData::outputRNApositions(const string & filename)
{
	TsvFile output;

	if (!output.open(filename))
		exitFail("Unable to open ", filename, " for position data");

	for (size_t i = 1;i < names.size();i++)
	{
		long len = lengths[i];
		long count = rawCounts[i];
		output.print(names[i], _s(len,":",count," plus"), This[names[i]].positions[0]);
		output.print(names[i], _s(len,":", count," minus"), This[names[i]].positions[1]);
	}

	return true;
}

//	Find the number of read position values in 'this' (a vector) that sits within each bin of the
//	histogram defined by E.   The results go into the 'freq'vector
//	This is used as part of the liklyhood calculations
void GeneCountData::histc(const vector<int> E)
{
	freq.assign(E.size(), 0);
	histoGram_ind.assign(lengths.size(), -1);
	for (size_t i = 0; i < names.size(); i++)
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

void GeneCountData::remove_invalid_values()
{
	for (size_t i = 0; i < names.size(); i++)
	{
		for (size_t j = 0; j < 2; j++)
			This[names[i]].positions[j].removeInvalidValues(lengths[i]);
	}
}

void GeneCountData::addEntry(string name, VEC_DATA_TYPE length)
{
	auto geneData = find(name);
	if (geneData == end())
	{
		rawCounts.push_back(0);
		names.push_back(name);
		lengths.push_back(length);
		emplace(name, rawCounts.size() - 1);
	}
}
void GeneCountData::addEntry(string name, VEC_DATA_TYPE length, VEC_DATA_TYPE count, rnaPosVec & posPositions, rnaPosVec & negPositions)
{
	auto geneData = find(name);
	if (geneData == end())
	{
		rawCounts.push_back(count);
		names.push_back(name);
		lengths.push_back(length);
		emplace(name, geneAttribute(rawCounts.size() - 1, posPositions, negPositions));
	}
}

void GeneCountData::addErrorEntry(string name)
{
	auto geneData = errorCounts.find(name);
	if (geneData == errorCounts.end())
	{
		errCounts.push_back(0);
		errorNames.push_back(name);
		errorCounts.emplace(name, errCounts.size() - 1);
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


	string buffer, gene, direction;
	int a = 1;
	stringEx lastGene;
	while (!f.eof() & (Ngenes-- != 0))
	{
		//	Need to check for gene and length consistency and that the gene has not appeared before

		getline(f, buffer);
		string gene, n1, n2;
		if (buffer.size())
		{
			rnaPosVec posPositions, negPositions;

			parseTsv(buffer, gene, n1, " ", direction, (std::vector<rna_pos_type> &)posPositions);
			getline(f, buffer);
			parseTsv(buffer, gene, n2, " ", direction, (std::vector<rna_pos_type> &)negPositions);

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

	if(lastGene)
		progMessage("Last gene = ", lastGene);
	return lastGene;
}

//	Transfers information for up to maxLength reads from up to Ngenes genes or transcripts
//	into the form which can be used by the mcmc chain
void GeneCountData::transferTo(dataType & mcmcData, size_t maxLength, int maxTotReads)
{
	vectorEx<int> bins{ { 0,300 } };
	for (size_t i = 500; i <= 10000; i += 500)
		bins.push_back(i);
	bins.add(11000, 12000, 15000, 30000);

	histc(bins);

	freq[0] = freq[0] * 2;
	freq[1] = freq[1] * 2;
	freq[21] = freq[21] / 2;
	freq[22] = freq[22] / 2;
	freq[23] = freq[23] / 6;
	freq[24] = freq[24] / 6;

	size_t geneIndex = 0;
	mcmcData.geneData[0].resize(size());
	mcmcData.geneData[1].resize(size());

	srand((unsigned)time(NULL));

	int Nreads = 0;

	//	Start at one because we dont include the reference gene (except there are no reads
	//	in the reference gene so this makes no difference
	for (size_t i = 1; i < names.size(); i++)
	{
		if (lengths[i] < MAX_LENGTH_OF_GENE_FOR_PARAM_ESTIMATION)
		{
			//	For the forward and the reverse counts
			for (size_t j = 0; j < 2; j++)
			{
				rnaPosVec & positions = at(names[i]).positions[j];
				positions.selectAtMost(maxLength);

				//	fragData contains the count 
				mcmcData.fragData.append(positions);

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


void GeneCountData::calculateOtherExpressionMeasures()
{

}


