#ifndef LIBICOUNT_H
#define LIBICOUNT_H

#include <mutex>
#include <atomic>
#include "GtfFileEx.h"

using namespace BamTools;

enum mode {
	intersect_union,
	intersect_strict,
	intersect_nonempty
};


class threadData
{
	geneCountsClass & sourcegeneCounts;

public:
	geneCountsClass geneCounts;
	std::vector<stringEx> samOutput;

	threadData(geneCountsClass & geneCounts) : sourcegeneCounts(geneCounts)
	{
	};
	~threadData() 
	{
		static std::mutex threadDataMutex;
		std::lock_guard<std::mutex> guard(threadDataMutex);
		for (auto i : geneCounts)
		{
			for (auto j : i.second)
				sourcegeneCounts[i.first][j.first] += j.second;
		}
			
	};

};
class LiBiCount
{

	std::atomic<size_t> bamCounter;
	int cacheFileCount;

	size_t cacheSize;


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
	bool processUnorderedBamData(int nThreads);
	bool processUnorderedBamDataThread();
	void processCachedReads();

	bool outputGeneCounts(const std::string & filename);
	void addRead(const regionLists & segments,const gtfFileEx & gtfData,threadData & threadData);


	void fileCompare(const std::string & maode);

};


class cacheRead : public regionLists
{
public:
	cacheRead() : file (0) {};
	~cacheRead(); 

	std::ifstream * file;

	bool open(const std::string filename);
	bool readNext();
	void close();
};


#endif
