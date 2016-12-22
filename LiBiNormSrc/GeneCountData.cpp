#include <fstream>
#include "GeneCountData.h"
#include "stringEx.h"
#include "parser.h"

using namespace std;

dataVec conv(const std::vector<rna_pos_type> a)
{
	dataVec retVal(a.size());
	for (size_t i = 0; i < a.size(); i++)
		retVal.at(i) = a.at(i);
	return retVal;
};


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

	if (norm.size())
	{
		//	Dont start at 0 as 0 is the reference for normalisation
		for (size_t i = 1; i < names.size(); i++)
		{
			output.printStart(names[i], (int)floor(rawCounts[i] / ((lengths[i] == 0) ? 1 : norm[i]) + 0.5));

			if (withDetails)
				output.printMiddle(rawCounts[i], lengths[i], (lengths[i] == 0) ? 1 : norm[i]);

			output.printEnd();
		}
	}
	else
	{
		for (auto i : This)
		{
			if ((i.first.substr(0, 2) != "__") && (i.second.index != 0))
			{
				output.print(i.first, rawCounts[i.second.index]);
			}
		}
		outputGeneCount(output, "__no_feature");
		outputGeneCount(output, "__ambiguous");
		outputGeneCount(output, "__too_low_aQual");
		outputGeneCount(output, "__not_aligned");
		outputGeneCount(output, "__alignment_not_unique");
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
		long len = lengths[i];
		long count = rawCounts[i];
		output.print(names[i], _s(len,":",count," plus"), This[names[i]].positions[0]);
		output.print(names[i], _s(len,":", count," minus"), This[names[i]].positions[1]);
	}

	return true;
}
