#ifndef GENE_COUNT_DATA_H
#define GENE_COUNT_DATA_H

#include <map>
#include "parser.h"
#include "mcmc.h"

typedef long rna_pos_type;

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
	rnaPosVec & removeInvalidValues(rna_pos_type maxVal);

	//	Test code takes the first N samples rather than randomly picks samples, and uses the same 
	//	algorithm as the MATLAB code for excluding invalid calues.   Used for comparing the two outputs
	void selectAtMost(size_t s);
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


	void reset() {
		positions[0].clear();
		positions[1].clear();
	}

	//	Number of associated reads
	size_t index;
	//	and their locations
	rnaPosVec positions[2];


};

//	This class holds and processes all of the count information associated with the set of genes
//	or transcripts
class GeneCountData : public std::map<const std::string, geneAttribute >
{
public:

	VEC_DATA_TYPE & count(const std::string gene);
	const VEC_DATA_TYPE & length(const std::string gene);
	//	Needed if we decide the data is not name ordered and have to restart
	void reset() { for (auto & gene : This)	gene.second.reset(); };

	void addEntry(std::string name, VEC_DATA_TYPE length = 0);
	void addEntry(std::string name, VEC_DATA_TYPE length, VEC_DATA_TYPE count, rnaPosVec & posPositions, rnaPosVec & negPositions);
	void addErrorEntry(std::string name);


	std::string loadData(const std::string filename, int Nlines = -1);
	void remove_invalid_values();
	void histc (const std::vector<int> E);
	void transferTo(dataType & mcmcData,size_t maxLength,int maxTotReads);

	void calculateOtherExpressionMeasures();

	bool outputGeneCounts(const std::string & filename, const std::string & title = "", bool withDetails = false);
	void outputGeneCount(TsvFile & output, const std::string name);
	bool outputRNApositions(const std::string & filename);

	void useSelectedGenes(const std::string & filename);

	std::vector<std::string> names;
	dataVec lengths;
	dataVec bias;

	//	The first entry is for the raw data, and the second for the normalised data
	dataVec counts[2];
	dataVec RPM[2], RPKM[2], RPK[2], TPM[2];


	//	For data associated with reads that do not map
	std::map<const std::string, geneAttribute > errorCounts;
	std::vector<std::string> errorNames;
	dataVec errCounts;

	std::vector<double> freq;
	std::vector<int> histoGram_ind;
};

#endif