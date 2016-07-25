#ifndef LIBICOUNT_H
#define LIBICOUNT_H

#include "GtfFileEx.h"

using namespace BamTools;

enum mode {
	intersect_union,
	intersect_strict,
	intersect_nonempty
};



class LiBiCount
{

public:
	bool useStrand,reverseStrand,verbose;
	mode countMode;

	geneCountsClass geneCounts;

	gtfFileEx genomeDef;

	BamReader reader;
	RefVector references;

	TsvFile outputFile;

	int main(int argc, char **argv);
	void processBamData();
	void outputGeneCounts(const std::string & filename);
	void addRead(const regionLists & segments,const gtfFileEx & gtfData);



	void fileCompare(const std::string & maode);

};

#endif
