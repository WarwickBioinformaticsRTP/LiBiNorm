#ifndef GENE_COUNT_DATA_H
#define GENE_COUNT_DATA_H

#include <map>
#include <random>

#include "Options.h"
#include "parser.h"
#include "mcmc.h"



class intRandClass
{
public:
	intRandClass() 
	{
#ifdef SELECT_READS_SEED
		seed = SELECT_READS_SEED;
#else
		seed = rd();
#endif
		gen.seed(seed);
	};

	unsigned int value(unsigned int max) { return gen() % max; };
	void reseed(unsigned int value) {
		seed = value;  
		gen.seed(seed);
	};
	unsigned int theSeed() { return seed; };

private:
	unsigned int seed;
	std::random_device rd;
	std::mt19937 gen;
	std::uniform_int_distribution<> dis;
};

extern intRandClass intRand;

typedef long rna_pos_type;

static std::string blankString = "";

//	Result options
static std::string ambiguousString = "__ambiguous";
static std::string noFeatureString = "__no_feature";
static std::string lowQualString = "__too_low_aQual";
static std::string notAlignedString = "__not_aligned";
static std::string notUnique = "__alignment_not_unique";

//
//	The rnaPosVec class holds the set of rna positions associated with one strand direction of a gene
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

//	Holds the read position information associated with a specific gene.
class geneAttribute 
{
public:
	geneAttribute(size_t index = 0) :index(index) {};

	//	Used by the copy constructor in GeneCountData::addEntry.  Provides a more efficient way of copying the
	//	position vectors as they are just about to be discarded
	geneAttribute(size_t index, rnaPosVec && posPos, rnaPosVec && negPos) :index(index)
	{
		swap(posPos, positions[0]);
		swap(negPos, positions[1]);
	};

	void reset() {
		positions[0].clear();
		positions[1].clear();
	}

	//	Index to the position where the associated name length and count information is held for this gene
	size_t index;
	//	and their locations
	rnaPosVec positions[2];
};

//	This class holds and processes all of the count information associated with the set of genes
//	or transcripts
class GeneCountData
{
public:

	VEC_DATA_TYPE & count(const std::string & gene);
	VEC_DATA_TYPE length(const std::string & gene);
	//	Needed if we decide the data is not name ordered and have to restart
	void reset() { for (auto & gene : readPositionData)	gene.second.reset(); };

	void addEntry(std::string name, VEC_DATA_TYPE length = 0);
	void addEntry(std::string name, VEC_DATA_TYPE length, VEC_DATA_TYPE count, rnaPosVec & posPositions, rnaPosVec & negPositions);
	void addErrorEntry(std::string name);


	std::string loadData(const std::string filename, int Nlines = -1);
	void remove_invalid_values();
	void histc (const std::vector<int> E);
	void transferTo(mcmcGeneData & mcmcData,size_t maxLength,int maxTotReads);

	bool outputGeneCounts(const std::string & filename, int detailLevel = 0, stringEx model = "");
	bool outputLandscape(const std::string & filename);

	void useSelectedGenes(const std::string & filename);

	//	Names, bias, lengths and RPM data are held in a series of vectors sharing a common gene order
	std::vector<std::string> names;
	dataVec bias;
	//	The first entry is for the raw lengths, and the second for the normalised lengths
	dataVec lengths[2];
	dataVec counts;
	dataVec RPM[2], RPKM[2], RPK[2], TPM[2];

	//	Read position data for all of the genes
	std::map<const std::string, geneAttribute > readPositionData;

	//	For data associated with reads that do not map
	std::map<const std::string, geneAttribute > errorCounts;
	std::vector<std::string> errorNames;
	dataVec errCounts;

	std::vector<double> freq;
	std::vector<int> histoGram_ind;
};

#endif