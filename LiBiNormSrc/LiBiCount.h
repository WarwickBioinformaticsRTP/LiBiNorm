#ifndef LIBICOUNT_H
#define LIBICOUNT_H

#include "GtfFileEx.h"

using namespace BamTools;

enum mode {
	intersect_union,
	intersect_strict,
	intersect_nonempty,
	intersect_all
};




class LiBiCount
{

	size_t bamCounter,maxCacheSize;
	int minqual;

	void incBamCounter(const BamAlignment * ba = 0,size_t size = -1);

public:
	bool useStrand,reverseStrand,verbose,htSeqCompatible;
	mode countMode;
	stringEx outputFilename,resultsFilename;

	geneCountsClass geneCounts;

	gtfFileEx genomeDef;

	BamReader reader;
	RefVector references;

	TsvFile outputFile;

	int main(int argc, char **argv);
	bool isValidAlignment(const BamAlignment & ba);

	bool processOrderedBamData();
	bool processUnorderedBamData();
	void processCachedReads(size_t cacheFileCount);

	bool outputGeneCounts(const std::string & filename);

	void addRead(const regionLists & segments,const gtfFileEx & gtfData);

	void fileCompare(int argc, char **argv);

};


#endif
