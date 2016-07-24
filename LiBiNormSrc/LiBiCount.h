#ifndef LIBICOUNT_H
#define LIBICOUNT_H

#include "GtfFileEx.h"

using namespace std;
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

	class geneCounts : public mapZeroDef<string,size_t>
	{
	public:
		void print(const string index,TsvFile & output)	{
			output.printEnd(index,This[index]);
		};
	} geneCounts;

	gtfFileEx genomeDef;

	BamReader reader;
	RefVector references;

	TsvFile outputFile;

	int main(int argc, char **argv);
	void processBamData();
	void outputGeneCounts(const string & filename);
	void addRead(const regionLists & segments,const gtfFileEx & gtfData);



	void fileCompare(const std::string & maode);

};

#endif
