
#ifndef GTFFILEEX_H
#define GTFFILEEX_H

#include "containerEx.h"
#include "libCommon.h"
#include "gtfFile.h"
#include "Regions.h"

class gtfRegion;

static std::string blankString = "";

//	Result options
static std::string ambiguousString = "__ambiguous";
static std::string noFeatureString = "__no_feature";
static std::string lowQualString = "__too_low_aQual";
static std::string notAlignedString = "__not_aligned";
static std::string notUnique = "__alignment_not_unique";


class gtfGeneAttribute : public mapZeroDef<std::string,size_t> 
{
public:
	void print(const std::string & index,TsvFile & output)	
	{
		for (auto & entry: This)
			output.print(index,entry.first,entry.second);
	};
	void reset()	
	{
		for (auto & entry: This)
			entry.second = 0;
	};

	const size_t operator++(int){
		size_t _R = at(blankString);
		at(blankString)++;
		return _R;
	}

};




//	A nested string map for holding counts for each identifier for each gene  
class geneCountsClass : public std::map<std::string,gtfGeneAttribute >
{
public:
	//	For printing out the list of counts.  entry.second is the counts
	void print(const std::string & index,TsvFile & output)	
	{
		at(index).print(index,output);
	};
	void reset()	{
		for (auto & gene: This)
			gene.second.reset();
	};
};

typedef std::multimap<size_t, gtfRegion> chromosomeGtfData ;

struct gtfOverlap
{
	size_t start,finish;
	bool strict;
	const std::string & geneName;
	const std::string & geneAttribute;
	gtfOverlap(size_t start,size_t finish,bool strict,const std::string & geneName,const std::string & geneAttribute): 
		start(start),finish(finish),strict(strict),geneName(geneName),geneAttribute(geneAttribute){};
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
