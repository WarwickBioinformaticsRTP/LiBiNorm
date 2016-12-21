#ifndef GENE_COUNT_DATA_H
#define GENE_COUNT_DATA_H

#include <map>
#include <string>
#include "dataVec.h"

typedef long rna_pos_type;

dataVec conv(const std::vector<rna_pos_type> a);

static std::string blankString = "";

//	Result options
static std::string ambiguousString = "__ambiguous";
static std::string noFeatureString = "__no_feature";
static std::string lowQualString = "__too_low_aQual";
static std::string notAlignedString = "__not_aligned";
static std::string notUnique = "__alignment_not_unique";



class geneAttribute 
{
public:
	geneAttribute(size_t index = 0) :index(index) {};

	//	Number of associated reads
	size_t index;
	//	and their locations
	std::vector<rna_pos_type> posPositions, negPositions;

	void reset() {
		posPositions.clear();
		negPositions.clear();

	}

};


class GeneCountData : public std::map<const std::string, geneAttribute >
{
public:
	GeneCountData();

	VEC_DATA_TYPE & operator[](const std::string gene) 
	{
		static VEC_DATA_TYPE dummy;
		auto i = find(gene);
		if (i == end())
		{
			_ASSERT_EXPR(false, "Looking for gene that was not in the reference genome");
			return dummy;
		}
		else
			return rawCounts.at((*i).second.index); 
	}

	//	For printing out the list of counts.  entry.second is the count data one entry per attribute being investigated
	void print(TsvFile & output, const std::string & index, double norm, size_t length)
	{
		output.print(index, rawCounts[at(index).index] * norm, rawCounts[at(index).index],length,norm);
	};
	void print(TsvFile & output, const std::string & index, double norm = 1)
	{
		output.print(index, rawCounts[at(index).index] * norm, norm);
	};

	bool outputGeneCounts(const std::string & filename, bool withDetails);

	//	Needed if we decide the data is not name ordered and have to restart
	void reset() {
		for (auto & gene : This)
			gene.second.reset();
	};

	void addEntry(std::string name, VEC_DATA_TYPE length = 0)
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

	void useSelectedGenes(const std::string & filename);

	std::vector<std::string> names;
	dataVec lengths;
	dataVec rawCounts;
	dataVec norm;

};

#endif