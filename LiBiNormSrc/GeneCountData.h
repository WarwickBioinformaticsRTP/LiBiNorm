#ifndef GENE_COUNT_DATA_H
#define GENE_COUNT_DATA_H

#include <map>
#include "parser.h"
#include "mcmc.h"

typedef long rna_pos_type;

dataVec conv(const std::vector<rna_pos_type> a);

static std::string blankString = "";

//	Result options
static std::string ambiguousString = "__ambiguous";
static std::string noFeatureString = "__no_feature";
static std::string lowQualString = "__too_low_aQual";
static std::string notAlignedString = "__not_aligned";
static std::string notUnique = "__alignment_not_unique";


class rnaPosVec : public std::vector<rna_pos_type>
{
public:
	rnaPosVec & removeInvalidValues(rna_pos_type maxVal)
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
	void selectAtMost(size_t s)
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

};

inline bool printVal(outputDataFile * f, rnaPosVec value)
{
	return printVal(f, (std::vector<rna_pos_type>)value);
};


class geneAttribute 
{
public:
	geneAttribute(size_t index = 0) :index(index) {};
	geneAttribute(size_t index, rnaPosVec & posPos, rnaPosVec & negPos) :index(index)
	{
		swap(posPos, positions[0]);
		swap(negPos, positions[1]);
	};

	//	Number of associated reads
	size_t index;
	//	and their locations
	rnaPosVec positions[2];

	void reset() {
		positions[0].clear();
		positions[1].clear();

	}

};


class GeneCountData : public std::map<const std::string, geneAttribute >
{
public:
	std::string loadData(const std::string filename, int Nlines = -1);
	void remove_invalid_values();
	void histc (const std::vector<int> E);
	void transferTo(dataType & mcmcData,size_t maxLength,int maxTotReads);

	std::vector<double> freq;
	std::vector<int> histoGram_ind;


	VEC_DATA_TYPE & operator()(const std::string gene) 
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




	bool outputGeneCounts(const std::string & filename, bool withDetails = false);
	bool outputRNApositions(const std::string & filename);

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
	void addEntry(std::string name, VEC_DATA_TYPE length,VEC_DATA_TYPE count,rnaPosVec & posPositions,rnaPosVec & negPositions)
	{
		auto geneData = find(name);
		if (geneData == end())
		{
			rawCounts.push_back(count);
			names.push_back(name);
			lengths.push_back(length);
			emplace(name, geneAttribute(rawCounts.size() - 1,posPositions,negPositions));
		}
	}

	void useSelectedGenes(const std::string & filename);

	std::vector<std::string> names;
	dataVec lengths;
	dataVec rawCounts;
	dataVec norm;

};

#endif