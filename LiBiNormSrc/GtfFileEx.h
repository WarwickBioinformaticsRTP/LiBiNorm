
#ifndef GTFFILEEX_H
#define GTFFILEEX_H

#include "containerEx.h"
#include "gtfFile.h"
#include "Regions.h"

class gtfRegion;

typedef std::multimap<size_t, gtfRegion> chromosomeGtfData ;

struct gtfOverlap
{
	size_t start,finish;
	bool strict;
	const gtfRegion & gtfReg;
	gtfOverlap(size_t start,size_t finish,bool strict,const gtfRegion & gtfReg): start(start),finish(finish),strict(strict),gtfReg(gtfReg){};
};

class gtfRegion
{
public:
	size_t start,finish;
	stringEx name;
	setEx<std::string> type;
	char strand;
	chromosomeGtfData::iterator overlaps;
	gtfRegion(	size_t start, size_t finish,const std::string & name,char strand,setEx<std::string> && type ):start(start),finish(finish),name(name),strand(strand),type(type){};
	bool checkOverlap(const region & segment,std::vector<gtfOverlap> & overlaps) const;
};


typedef std::map<std::string,chromosomeGtfData> genomeGtfRegions;
typedef std::multimap<size_t,chromosomeGtfData::iterator> chromosomeEndIndexMap;
typedef std::map<std::string,chromosomeEndIndexMap > genomeEndIndexMap;


class gtfFileEx : public gtfFile
{
public: 
	genomeGtfRegions genomeGtfData; 

	genomeEndIndexMap genomeEndIndex;

	void index(mapZeroDef<std::string,size_t> & geneCounts);
	void outputChromData(const std::string & filename);

};

#endif
