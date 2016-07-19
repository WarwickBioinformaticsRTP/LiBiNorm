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



class gtfRegion
{
public:
	size_t start,finish;
	string name;
	setEx<string> type;
	char strand;
	vector<gtfRegion *> overlaps;
	gtfRegion(	size_t start, size_t finish,const string & name,char strand,setEx<string> && type ):start(start),finish(finish),name(name),strand(strand),type(type){};
	bool checkOverlap(const region & segment,bool & strict);
};


class chromosomeData : public multimap<size_t, gtfRegion> 
{

};


class gtfFileEx : public gtfFile
{
	map<string,chromosomeData> chromData; 
public: 
	bool useStrand,reverseStrand;

	mapZeroDef<string,size_t> geneCounts;

	void index();
	void outputChromData(const string & filename);
	void outputGeneCounts(const string & filename);
	void addRead(const string & chromosome,const regionList regions,mode countMode);

};


class LiBiCount
{

public:
	int main(int argc, char **argv);


};

#endif
