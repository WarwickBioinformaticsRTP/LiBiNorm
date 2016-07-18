#ifndef LIBICOUNT_H
#define LIBICOUNT_H

#include <vector>
#include <utility>
#include "api/BamReader.h"
#include "gtfFile.h"

using namespace std;
using namespace BamTools;

class region : public pair<size_t,size_t>
{
public:
	region(size_t first,size_t second) : pair<size_t,size_t>(first,second){};
};


class regionList : public vector<region>
{
public:
	bool combine(size_t start,size_t end);
	void add(size_t start,size_t end);
};

void GetRegions(const BamAlignment & ba,regionList & regions);


class gtfRegion
{
public:
	size_t start,finish;
	string name;
	char strand;
	vector<gtfRegion *> overlaps;
	gtfRegion(	size_t start, size_t finish,const string & name,char strand):start(start),finish(finish),name(name),strand(strand){};
	bool checkOverlap(const region & segment);
};


class chromosomeData : public multimap<size_t, gtfRegion> 
{

};


class gtfFileEx : public gtfFile
{
	map<string,chromosomeData> chromData; 
public: 
	void index();
	void outputChromData(const string & filename);
	void addRead(const string & chromosome,const regionList regions);

};


class LiBiCount
{

public:
	int main(int argc, char **argv);


};

#endif
