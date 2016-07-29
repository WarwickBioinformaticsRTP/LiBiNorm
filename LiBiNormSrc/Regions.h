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
//		std::string name;
		int refId;
		int position;
		char strand;
		Cigar cigar;
		cacheEntry(void){}
		cacheEntry(const BamTools::BamAlignment & ba):/*name(ba.Name),*/refId(ba.RefID),position(ba.Position+1),
			strand((ba.IsReverseStrand() == ba.IsFirstMate())?'-':'+'),cigar(ba.CigarData){};
//			strand(ba.IsReverseStrand()?'-':'+'),cigar(ba.CigarData){};
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
	region(){};
	region(int start,int end,char strand) : start(start),end(end),strand(strand) {};
};



class regionList 
{
public:
	std::map<int,region> data;


	regionList(size_t start,char strand,const std::vector<BamTools::CigarOp> & co);
	regionList(const cacheEntry & read) : regionList(read.position,read.strand,read.cigar){};
	regionList(){};

	void GetRegions(const BamTools::BamAlignment & ba);
	void combine(const regionList & rl);
	void combineRegion(const region & r);
	void add(int start,int end,char strand);
};


class regionLists
{
public:
	std::map<int,regionList> data;
	std::string name;

	regionLists(){};
//	regionLists(const std::string & line);
	regionLists(const cacheEntry & read){
		data[read.refId].combine(regionList(read));
	};

	regionLists(regionLists && rl): data(move(rl.data)) {};

	void GetRegions(const BamTools::BamAlignment & ba);

	void combine(const regionLists & rl);
	void combine(const cacheRead & read){
			data[read.refId].combine(regionList(read));
	};
};

namespace parserInternal
{
//	void parseval(const char *& start,regionList & rl,size_t & len);
//	void parseval(const char *& start,region & r,size_t & len);
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


/*
inline bool printVal(outputDataFile * f,const region & value)
{
	f->printStart(value.start,value.end,value.strand);
	return true;
};


inline bool printVal(outputDataFile * f,const regionList & value)
{
	f->printStart(value.data);
	return true;
};

inline bool printVal(outputDataFile * f,const regionLists & value)
{
	f->printStart(value.data);
	return true;
};

template<class _Kty,class _Ty>
inline bool printVal(outputDataFile * f,const std::map<_Kty,_Ty> & data)
{
	f->printStart(data.size());
	for (auto i: data) 
		f->printMiddle(printZero(i.first),i.second);
	return true;
};
*/


#endif

