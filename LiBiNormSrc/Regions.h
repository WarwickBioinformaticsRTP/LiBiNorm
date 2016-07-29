#ifndef REGIONS_H
#define REGIONS_H

#include<map>
#include<string>
#include "printEx.h"
#include "api/BamReader.h"
#include "libCommon.h"

//#define _DEBUG 1

class Cigar : public std::vector<BamTools::CigarOp>
{
public:
	Cigar(){};
	Cigar(const std::vector<BamTools::CigarOp> & co): std::vector<BamTools::CigarOp>(co) {};
};

class cacheEntry
{
	public:
		int refId;
		int position;
		char strand;
		Cigar cigar;
		cacheEntry(void){}
		cacheEntry(const BamTools::BamAlignment & ba):
			refId(ba.RefID),
			position(ba.Position+1),
			strand((ba.IsReverseStrand() == ba.IsFirstMate())?'-':'+'),
			cigar(ba.CigarData)
		{};
};



class cacheRead : public cacheEntry
{
public:
	std::string name;
	cacheRead() : file (0) {};
	~cacheRead(); 

	std::ifstream * file;

	bool open(const std::string filename);
	bool readNext();
	void close();
};

class region 
{
public:
	int start,end;
	char strand;
	region(int start,int end,char strand) : start(start),end(end),strand(strand) {};
};



class regionList 
{
public:
	std::map<int,region> data;

	regionList(const cacheEntry & read);
	regionList(){};

	void combine(const regionList & rl);
	void combineRegion(const region & r);

};


class regionLists
{
public:
	std::map<int,regionList> data;
	const std::string & name;

	regionLists(const cacheEntry & read,std::string name) :name(name) {
		data[read.refId].combine(regionList(read));
	};

	void combine(const cacheEntry & read){
			data[read.refId].combine(regionList(read));
	};
};

namespace parserInternal
{
	void parseval(const char *& start, Cigar & cigar,size_t & len);
}

inline bool printVal(outputDataFile * f,const Cigar & cigar)
{
	for (auto & i: cigar)
		fprintf(f->fout,"%c%i",i.Type,i.Length);
	return true;
};

inline bool printVal(outputDataFile * f,const cacheEntry & read)
{
	f->printStart(printZero(read.refId),read.position,read.strand,read.cigar);
	return true;
};

#endif

