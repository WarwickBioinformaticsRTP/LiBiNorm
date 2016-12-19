#ifndef LIBICOUNT_H
#define LIBICOUNT_H

#include "FeatureFileEx.h"
#include "dataVec.h"

using namespace BamTools;

enum mode {
	intersect_union,
	intersect_strict,
	intersect_nonempty,
	intersect_all
};


inline dataVec makeVec(const std::vector<rna_pos_type> & b)
{
	dataVec _ret(b.size());
	for (size_t i = 0; i < b.size(); i++)
		_ret.at(i) = b.at(i);
	return _ret;
};

class LiBiCount
{

	size_t bamCounter,maxCacheSize;
	int minqual;

	void incBamCounter(const BamAlignment * ba = 0,int size = -1);

public:
	bool useStrand,reverseStrand,htSeqCompatible,nameOrder;
	mode countMode;
	stringEx outputFilename,countsFilename,landscapeFilename, normaliseResultsFilename,tempDirectory;

	geneCountsClass geneCounts;

	featureFileEx genomeDef;

	BamReader reader;
	RefVector references;

	TsvFile outputFile;

	int main(int argc, char **argv);
	bool AReadIsMapped(const BamAlignment & ba);

	bool processNameOrderedBamData();
	bool processPositionOrderedBamData();
	void processCachedReads(size_t cacheFileCount);

	bool outputGeneCounts(const std::string & filename,bool withDetails = false);
	bool outputRNApositions(const stringEx & filename);
	
	void addRead(const regionLists & segments,const featureFileEx & gtfData);

	void fileCompare(int argc, char **argv);

};


#endif
