#ifndef REGIONS_H
#define REGIONS_H

#include<map>
#include<string>
#include "printEx.h"
#include "api/BamReader.h"
#include "libCommon.h"

#define _DEBUG 1

class region 
{
public:
	int start,end;
	char strand;
	region(){};
	region(int start,int end,bool revStrand) : start(start),end(end),strand(revStrand?'-':'+') {};
};



class regionList 
{
public:
	std::map<int,region> data;
	void GetRegions(const BamTools::BamAlignment & ba);

	void combine(const regionList & rl);
	void combineRegion(const region & r);
	void add(int start,int end,bool revStrand);
};


class regionLists
{
public:
	std::map<int,regionList> data;
	std::string name;

	regionLists(){};
	regionLists(const std::string & line);
	regionLists(regionLists && rl): data(move(rl.data)) {};

	void GetRegions(const BamTools::BamAlignment & ba);

	void combine(const regionLists & rl);
	void print(TsvFile & file);

};

namespace parserInternal
{
	void parseval(const char *& start,regionList & rl,size_t & len);
	void parseval(const char *& start,region & r,size_t & len);
}

inline bool printVal(outputDataFile * f,const region & value)
{
	f->printStart(value.start,value.end,value.strand);
	return false;
};


inline bool printVal(outputDataFile * f,const regionList & value)
{
	f->printStart(value.data.size());
	for (auto i: value.data) 
		f->printMiddle(i.first,i.second);
	return false;
};

inline bool printVal(outputDataFile * f,const regionLists & value)
{
	f->printStart(value.data.size());
	for (auto i: value.data) 
		f->printMiddle(printZero(i.first),i.second);
	return true;
};


#endif

