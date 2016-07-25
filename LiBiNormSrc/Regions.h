#ifndef REGIONS_H
#define REGIONS_H

#include<map>
#include<string>
#include "api/BamReader.h"
#include "libCommon.h"

#define _DEBUG 1

class region 
{
public:
	size_t start,end;
	char strand;
	region(size_t start,size_t end,bool revStrand) : start(start),end(end),strand(revStrand?'-':'+') {};
};


class regionList : public std::map<size_t,region>
{
public:
//	bool standardPair;
	void GetRegions(const BamTools::BamAlignment & ba);
	bool combineRegion(size_t start,size_t end,bool revStrand);
	void add(size_t start,size_t end,bool revStrand);
};

class regionLists : public std::map<int,regionList>
{
public:
	regionLists(){};
	std::string name;
	void GetRegions(const BamTools::BamAlignment & ba);
//	regionLists(regionLists & rl): std::map<int,regionList>(rl){};
	regionLists(regionLists && rl): std::map<int,regionList>(move(rl))
	{
	};
};

#endif

