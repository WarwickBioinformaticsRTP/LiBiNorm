#include "GeneCountData.h"
#include "parser.h"
#include <fstream>

using namespace std;

dataVec conv(const std::vector<rna_pos_type> a)
{
	dataVec retVal(a.size());
	for (size_t i = 0; i < a.size(); i++)
		retVal.at(i) = a.at(i);
	return retVal;
};

GeneCountData::GeneCountData()
{
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
		if (names.empty() || (names.back() != gene))
		{
			names.emplace_back(gene);
		}
	};
}

bool GeneCountData::outputGeneCounts(const string & filename, bool withDetails)
{
	TsvFile output;

	if (!output.open(filename))
		return false;

	for (size_t i = 1; i < names.size(); i++)
	{
		if (withDetails)
			output.print(names[i], rawCounts[i] / norm[i], rawCounts[i], lengths[i], norm[i]);
		else
			output.print(names[i], rawCounts[i] / norm[i], rawCounts[i], lengths[i], norm[i]);

	}
}