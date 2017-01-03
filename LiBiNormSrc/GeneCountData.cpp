#include <fstream>
#include "LiBiNorm.h"
#include "GeneCountData.h"
#include "stringEx.h"
#include "containerEx.h"
#include "parser.h"

using namespace std;


rnaPosVec & rnaPosVec::removeInvalidValues(rna_pos_type maxVal)
{
	//	Sort in place for maximum efficiency, if we find an invalid value, replace with one from the end;
	//	Note that if we swap with a value from the end we have to check it as well to see if it is invalid
#ifdef TEST_CODE
	size_t offset = 0;
	for (size_t i = 0; i + offset < size();)
	{
		if ((at(i + offset) < 0) || (at(i + offset) >= maxVal))
			offset++;
		else
		{
			if (offset)
				at(i) = at(i + offset);
			i++;
		}
	}
	if (offset)
		resize(size() - offset);
#else
	iterator i = begin(), j = end();
	while (i != j)
	{
		if ((*i < 0) || (*i >= maxVal))
			std::swap(*i, *--j);
		else
			i++;
	}
	resize(j - begin());
#endif
	return *this;
}

//	Test code takes the first N samples rather than randomly picks samples, and uses the same 
//	algorithm as the MATLAB code for excluding invalid calues.   Used for comparing the two outputs
// #define TEST_CODE
void rnaPosVec::selectAtMost(size_t s)
{
	if (s > size())
		return;
	//	Swap the first s entries with the entry at some other position, then resize to just have the s entries
#ifndef TEST_CODE
	iterator i = begin();
	for (size_t j = 0; j < s; j++)
	{
		std::swap(*(i++), *(begin() + rand() % size()));
	}
#endif
	resize(s);
}


//	Returns the address of the counter associated with a gene or error condition
//	such that the value can be incremened when an associated error if found

VEC_DATA_TYPE & GeneCountData::count(const std::string & gene)
{
	static VEC_DATA_TYPE dummy;
	auto i = readPositionData.find(gene);
	if (i == readPositionData.end())
	{
		i = errorCounts.find(gene);
		if (i == errorCounts.end())
		{
			_ASSERT_EXPR(false, "Looking for gene that was not in the reference genome");
			return dummy;
		}
		return errCounts.at((*i).second.index);
	}
	return counts.at((*i).second.index);
}

//	Returns the actual length of a gene (rather than the normalised length) given the name  
VEC_DATA_TYPE GeneCountData::length(const std::string & gene)
{
	auto i = readPositionData.find(gene);
	if (i == readPositionData.end())
	{
		return 0;
	}
	return lengths[0].at((*i).second.index);
}

//	Reads in a file contain a list of gene names.  Only these genes will then be used
//	in the subsequent analysis
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

//	Print the results to varying levels of detail
//	0 = htseq-count compatible
//	1 = basic normalised results
//	2 = As 1 but with header
//	3 = As 2 but more columns

bool GeneCountData::outputGeneCounts(const string & filename, int detailLevel, stringEx model)
{
	TsvFile output;

	model.replace(" ", "");

	if (!output.open(filename))
		return false;

	if ((bias.size()) && (detailLevel > 0))
	{
		lengths[1] = lengths[0] * bias;

		// Definitions taken from http://www.rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/
		for (size_t i = 0; i < ((bias.size()) ? 2 : 1); i++)
		{
			VEC_DATA_TYPE scalingFactor = sum(counts) / 1000000;

			RPM[i] = counts / scalingFactor;
			RPKM[i] = RPM[i] / lengths[i];

			RPK[i] = counts / lengths[i];
			scalingFactor = sum(RPK[i]) / 1000000;
			TPM[i] = RPK[i] / scalingFactor;
		}

		switch (detailLevel)
		{
		case 1:
			output.print("Gene", "count", "length", "Bias", _s("FPKM_", model), _s("TPM_", model));
			break;
		case 2:
			output.print("", "", "RNA", "Bias:", model + " +", model + " +");
			output.print("Gene", "count", "length", model, "FPKM", "TPM");
			output.print();
			break;
		case 3:
			output.print("", "", "RNA", "Bias:", model + " +", model + " +", model + " +", model + " +", "Raw", "Raw", "Raw", "Raw");
			output.print("Gene", "count", "length", model, "RPM", "RPKM", "RPK", "TPM", "RPM", "RPKM", "RPK", "TPM");
			output.print();
			break;
		}

		//	Dont start at 0 as 0 is the reference for normalisation
		for (size_t i = 1; i < names.size(); i++)
		{
			output.printStart(names[i], counts[i]);

			switch (detailLevel)
			{
			case 1:
			case 2:
				output.printMiddle(lengths[0][i], bias[i],RPKM[1][i], TPM[1][i]);
				break;
			case 3:
				output.printMiddle(lengths[0][i], bias[i],
					RPM[1][i], RPKM[1][i], RPK[1][i], TPM[1][i],
					RPM[0][i], RPKM[0][i], RPK[0][i], TPM[0][i]);
				break;
			}
			output.printEnd();
		}
	}
	else
	{
		//	We are just printing the counts in htseq-count mode, which means that they
		//	need to be in name order, and not in the order they are found 
		for (auto i : readPositionData)
		{
			if (i.second.index != 0)
				output.print(i.first, (int)counts[i.second.index]);
		}
	}

	//	Print the error counts at the end
	if (detailLevel != 1)
	{
		for (size_t i = 0; i < errorNames.size(); i++)
			output.print(errorNames[i], (int)errCounts[i]);
	}

	return true;
}

bool GeneCountData::outputRNApositions(const string & filename)
{
	TsvFile output;

	if (!output.open(filename))
		exitFail("Unable to open ", filename, " for position data");

	for (size_t i = 1;i < names.size();i++)
	{
		long len = lengths[0][i];
		long count = counts[i];
		string & name = names[i];
		output.print(name, _s(len,":",count," plus"), readPositionData[name].positions[0]);
		output.print(name, _s(len,":", count," minus"), readPositionData[name].positions[1]);
	}

	return true;
}

//	Find the number of read position values in 'this' (a vector) that sits within each bin of the
//	histogram defined by E.   The results go into the 'freq'vector
//	This is used as part of the liklyhood calculations
void GeneCountData::histc(const vector<int> E)
{
	freq.assign(E.size(), 0);
	histoGram_ind.assign(lengths[0].size(), -1);
	for (size_t i = 0; i < names.size(); i++)
	{
		double v = lengths[0][i];
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
			readPositionData[names[i]].positions[j].removeInvalidValues(lengths[0][i]);
	}
}

void GeneCountData::addEntry(string name, VEC_DATA_TYPE length)
{
	auto geneData = readPositionData.find(name);
	if (geneData == readPositionData.end())
	{
		counts.push_back(0);
		names.push_back(name);
		lengths[0].push_back(length);
		readPositionData.emplace(name, counts.size() - 1);
	}
}
void GeneCountData::addEntry(string name, VEC_DATA_TYPE length, VEC_DATA_TYPE count, rnaPosVec & posPositions, rnaPosVec & negPositions)
{
	auto geneData = readPositionData.find(name);
	if (geneData == readPositionData.end())
	{
		counts.push_back(count);
		names.push_back(name);
		lengths[0].push_back(length);
		readPositionData.emplace(name, geneAttribute(counts.size() - 1, move(posPositions), move(negPositions)));
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
	mcmcData.geneData[0].resize(readPositionData.size());
	mcmcData.geneData[1].resize(readPositionData.size());

	srand((unsigned)time(NULL));

	int Nreads = 0;

	//	Start at one because we dont include the reference gene (except there are no reads
	//	in the reference gene so this makes no difference
	for (size_t i = 1; i < names.size(); i++)
	{
		if (lengths[0][i] < MAX_LENGTH_OF_GENE_FOR_PARAM_ESTIMATION)
		{
			//	For the forward and the reverse counts
			for (size_t j = 0; j < 2; j++)
			{
				rnaPosVec & positions = readPositionData.at(names[i]).positions[j];
				positions.selectAtMost(maxLength);

				//	fragData contains the count 
				mcmcData.fragData.append(positions);

				mcmcData.geneIndex.insert(mcmcData.geneIndex.end(), positions.size(), geneIndex);//gene.length);
				Nreads += positions.size();
			}
			mcmcData.geneData[0][geneIndex] = lengths[0][i];
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


