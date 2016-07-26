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

	size_t bamCounter;

	void incBamCounter(const BamAlignment * ba = 0,size_t size = -1);

public:
	bool useStrand,reverseStrand,verbose;
	mode countMode;
	stringEx resultsFilename;

	geneCountsClass geneCounts;

	gtfFileEx genomeDef;

	BamReader reader;
	RefVector references;

	TsvFile outputFile;

	int main(int argc, char **argv);
	bool processOrderedBamData();
	bool processUnorderedBamData();
	void outputGeneCounts(const std::string & filename);
	void addRead(const regionLists & segments,const gtfFileEx & gtfData);


	void fileCompare(const std::string & maode);

};


class cacheRead : public regionLists
{
public:
	cacheRead() : file (0) {};
	~cacheRead() {delete(file);};

	std::ifstream * file;

	bool open(const std::string filename)
	{
		file = new std::ifstream();
		file ->open(filename);
		if (!file ->is_open()) return false;
		readNext();
		return true;
	}
	bool readNext()
	{
		if (file->eof())
			return false;
		data.clear();
		std::string line;
		getline(*file,line);
		parseTsv(line,name,data);
		return true;
	}
};


#endif
