
#ifndef GTFFILEEX_H
#define GTFFILEEX_H

#include "containerEx.h"
#include "gtfFile.h"
#include "Regions.h"

class gtfRegion;

class geneCountsClass : public std::map<std::string, mapZeroDef<std::string,size_t> >
{
public:
	void print(const std::string index,TsvFile & output)	{
		for (auto & entry: This[index])
		{
			output.printEnd(index,entry.first,entry.second);
		}
	};
	void reset()	{
		for (auto & gene: This)
		{
			for (auto & type: gene.second)
			{
				type.second = 0;
			}
		}
	};
};

typedef std::multimap<size_t, gtfRegion> chromosomeGtfData ;

struct gtfOverlap
{
	size_t start,finish;
	bool strict;
	const std::string & geneName;
	const std::string & type;
	gtfOverlap(size_t start,size_t finish,bool strict,const std::string & geneName,const std::string & type): 
		start(start),finish(finish),strict(strict),geneName(geneName),type(type){};
};

class gtfRegion
{
public:
	size_t start,finish;
	stringEx name;
	std::string type;
	char strand;
	chromosomeGtfData::iterator overlaps;
	gtfRegion(	size_t start, size_t finish,const std::string & name,char strand,const std::string & type ):start(start),finish(finish),name(name),strand(strand),type(type){};
	void checkOverlap(const region & segment,std::vector<gtfOverlap> & overlaps) const;
};


typedef std::map<std::string,chromosomeGtfData> genomeGtfRegions;
typedef std::multimap<size_t,chromosomeGtfData::iterator> chromosomeEndIndexMap;
typedef std::map<std::string,chromosomeEndIndexMap > genomeEndIndexMap;


class gtfFileEx : public gtfFile
{
public: 
	genomeGtfRegions genomeGtfData; 

	genomeEndIndexMap genomeEndIndex;

	void index(geneCountsClass & geneCounts);
	void outputChromData(const std::string & filename);

};

#endif
