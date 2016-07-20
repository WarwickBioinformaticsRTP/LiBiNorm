
#ifndef GTFFILEEX_H
#define GTFFILEEX_H

#include "containerEx.h"
#include "gtfFile.h"
#include "Regions.h"

class gtfRegion;

typedef std::multimap<size_t, gtfRegion> chromosomeData ;

class gtfRegion
{
public:
	size_t start,finish;
	std::string name;
	setEx<std::string> type;
	char strand;
	chromosomeData::iterator overlaps;
	gtfRegion(	size_t start, size_t finish,const std::string & name,char strand,setEx<std::string> && type ):start(start),finish(finish),name(name),strand(strand),type(type){};
	bool checkOverlap(const region & segment,bool & strict) const;
};




class gtfFileEx : public gtfFile
{
public: 
	std::map<std::string,chromosomeData> chromData; 

	std::map<std::string,std::multimap<size_t,chromosomeData::iterator> > chromEndIndex;

	void index(mapZeroDef<std::string,size_t> & geneCounts);
	void outputChromData(const std::string & filename);

};

#endif
