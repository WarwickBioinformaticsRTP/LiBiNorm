#ifndef LIBICOUNT_H
#define LIBICOUNT_H

#include <vector>
#include <utility>
#include "containerEx.h"
#include "api/BamReader.h"
#include "gtfFile.h"

using namespace std;
using namespace BamTools;

enum mode {
	intersect_union,
	intersect_strict,
	intersect_nonempty
};


class region 
{
public:
	size_t start,end;
	char strand;
	region(size_t start,size_t end,bool revStrand) : start(start),end(end),strand(revStrand?'-':'+') {};
};


class regionList : public map<size_t,region>
{
public:
	void GetRegions(const BamAlignment & ba);
	bool combineRegion(size_t start,size_t end,bool revStrand);
	void add(size_t start,size_t end,bool revStrand);
};

class regionLists : public map<int,regionList>
{
public:
	_DBG( string name;)
	void GetRegions(const BamAlignment & ba);
};

class gtfRegion;

class chromosomeData : public multimap<size_t, gtfRegion> 
{

};


class gtfRegion
{
public:
	size_t start,finish;
	string name;
	setEx<string> type;
	char strand;
	multimap<size_t,chromosomeData::iterator> overlaps;
	gtfRegion(	size_t start, size_t finish,const string & name,char strand,setEx<string> && type ):start(start),finish(finish),name(name),strand(strand),type(type){};
	bool checkOverlap(const region & segment,bool & strict) const;
};




class gtfFileEx : public gtfFile
{
public: 
	map<string,chromosomeData> chromData; 

	map<string,multimap<size_t,chromosomeData::iterator> > chromEndIndex;

	void index(mapZeroDef<string,size_t> & geneCounts);
	void outputChromData(const string & filename);

};


class LiBiCount
{

public:
	bool useStrand,reverseStrand;
	mode countMode;

	mapZeroDef<string,size_t> geneCounts;
	RefVector references;


	int main(int argc, char **argv);
	void outputGeneCounts(const string & filename);
	void addRead(const regionLists & regions,const gtfFileEx & gtfData);

};

#endif
