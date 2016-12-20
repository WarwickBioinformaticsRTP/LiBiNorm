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

//	Conatins information relating to a specific gene/region type  combination
class geneTypeInfo
{
public:
	//	Number of associated reads
	size_t count;
	//	and their locations
	std::vector<rna_pos_type> posPositions, negPositions;

	geneTypeInfo() :count(0) {};
	void reset() {
		count = 0;
		posPositions.clear();
		negPositions.clear();

	}

};

//	Support function used by the print class for printing gtf tags
inline bool printVal(outputDataFile * f, const geneTypeInfo & gti)
{
	fprintf(f->fout, "%zu", gti.count);
	return true;
};


class geneAttribute : public geneTypeInfo
{
public:
	void print(TsvFile & output, const std::string & index, double norm)
	{
			output.print(index, count * norm);
	};
	void print(TsvFile & output, const std::string & index, double norm, size_t length)
	{
		output.print(index, count*norm, count, length, norm);
	};

	const size_t operator++(int) {
		//	Increments the count for the 'totals' counts for which there is no 'type' information 
		//	such as notUnique and loqQualString for which we use a dummy 'blank' entry 
		return count++;
	}

};


class GeneCountData : public std::map<const std::string, geneAttribute >
{
public:
	GeneCountData();

	//	For printing out the list of counts.  entry.second is the count data one entry per attribute being investigated
	void print(TsvFile & output, const std::string & index, double norm, size_t length)
	{
		at(index).print(output, index, norm, length);
	};
	void print(TsvFile & output, const std::string & index, double norm = 1)
	{
		at(index).print(output, index, norm);
	};

	//	Needed if we decide the data is not name ordered and have to restart
	void reset() {
		for (auto & gene : This)
			gene.second.reset();
	};



	void useSelectedGenes(const std::string & filename);

	std::vector<std::string> names;
	dataVec lengths;
	dataVec rawCounts;

};

#endif