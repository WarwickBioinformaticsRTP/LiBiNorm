#ifndef LIBICOUNT_H
#define LIBICOUNT_H

#include "FeatureFileEx.h"
#include "LiBiNorm.h"
#include "dataVec.h"

using namespace BamTools;

//	The htseq-count modes
enum mode {
	intersect_union,
	intersect_strict,
	intersect_nonempty,
	intersect_all
};

class LiBiCount : private LiBiNorm
{
public:
	int main(int argc, char **argv);

private:
	//	For reading and processing bam data
	bool processNameOrderedBamData();
	bool processPositionOrderedBamData();
	void processCachedReads(size_t cacheFileCount);

	//Support functions for reading bam data
	bool AReadIsMapped(const BamAlignment & ba);
	void addRead(const regionLists & segments, const featureFileEx & gtfData);
	void incBamCounter(const BamAlignment * ba = 0, int size = -1);

	//	For comparing two files.  Not currently used
	void fileCompare(int argc, char **argv);
	
	//  For reading and persisting header data from a bam file
	BamReader reader;
	RefVector references;

	//	Config data for reading the ba, file
	bool useStrand, reverseStrand, nameOrder;
	mode countMode;
	int minqual;
	size_t bamCounter, maxCacheSize;

	//	Directory data for outputting results
	stringEx outputFilename, tempDirectory;
	TsvFile outputFile;

	//	For obtaining information from the gtf/gff file
	featureFileEx genomeDef;
};


#endif
