
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


class gtfGeneAttribute : public mapZeroDef<const std::string,size_t> 
{
public:
	void print(const std::string & index,TsvFile & output)	
	{
		//Print entries for all of the attributes being considered
		if (size() > 1)
			for (auto & entry: This)
				output.print(index,entry.first,_Z(entry.second));
		else
			for (auto & entry: This)
				output.print(index,_Z(entry.second));
	};
	void reset()	
	{
		for (auto & entry: This)
			entry.second = 0;
	};

	const size_t operator++(int){
	//	Increments the count for the 'totals' counts for which there is no 'type' information so 
	//	we use a dummy 'blank' entry 

		//	For consistency with definition of post operator, return velu before increment
		size_t _R = at(blankString);
		at(blankString)++;
		return _R;
	}

};




//	A nested string map for holding counts for each identifier for each gene  
class geneCountsClass : public std::map<const std::string,gtfGeneAttribute >
{
public:
	//	For printing out the list of counts.  entry.second is the count data one entry per attribute being investigated
	void print(const std::string & index,TsvFile & output)	
	{
		at(index).print(index,output);
	};

	//	Needed if we decide the data is not name ordered and have to restart
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
	//	Store references here as these are created and deleted lots and these removes the need to allocate and
	//	deallocate on the heap.
	const std::string & geneName;
	const std::string & featType;
	gtfOverlap(size_t start,size_t finish,bool strict,const std::string & geneName,const std::string & featType): 
		start(start),finish(finish),strict(strict),geneName(geneName),featType(featType){};
};

class gtfRegion
{
public:
	size_t start,finish;
	//	Store actual values here as these are only created once and then referenced lots, so this is more efficient
	const std::string name;
	const std::string type;
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
