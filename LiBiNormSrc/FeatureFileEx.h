
#ifndef FEATURE_FILE_EX_HEADER
#define FEATURE_FILE_EX_HEADER

#include "containerEx.h"
#include "libCommon.h"
#include "featureFile.h"
#include "Regions.h"
#include "dataVec.h"

class gtfRegion;

static std::string blankString = "";

//	Result options
static std::string ambiguousString = "__ambiguous";
static std::string noFeatureString = "__no_feature";
static std::string lowQualString = "__too_low_aQual";
static std::string notAlignedString = "__not_aligned";
static std::string notUnique = "__alignment_not_unique";

typedef long rna_pos_type;

dataVec conv(const std::vector<rna_pos_type> a);


//	Conatins information relating to a specific gene/region type  combination
class geneTypeInfo
{
public:
	//	Number of associated reads
	size_t count;
	//	and their locations
	std::vector<rna_pos_type> posPositions, negPositions;

	geneTypeInfo() :count(0) {};
	void reset() {
		count = 0;
		posPositions.clear();
		negPositions.clear();

	}
	//	For consistency with definition of post operator, return value before increment
	const size_t operator++(int) {
		size_t _R = count;
		count++;
		return _R;
	}

};

//	Support function used by the print class for printing gtf tags
inline bool printVal(outputDataFile * f, const geneTypeInfo & gti)
{
	fprintf(f->fout, "%zu", gti.count);
	return true;
};


class geneAttribute : public std::map<const std::string, geneTypeInfo>
{
public:
	void print(TsvFile & output,const std::string & index, double norm)
	{
		// Print entries for all of the attributes being considered.  Only output the type info if 
		// there are multiple region types being considered.  This is for consistency with
		// the htseq-count standard count output
		if (size() > 1)
			for (auto & entry: This)
				output.print(index,entry.first,entry.second);
		else
		{
			auto i = begin();
			output.print(index, (int)((*i).second.count * norm));
		}
	};
	void print(TsvFile & output, const std::string & index, double norm, size_t length)
	{
		auto i = begin();
		output.print(index, (int)((*i).second.count*norm), (*i).second.count, length,norm);
	};
	void reset()
	{
		for (auto & entry: This)
			entry.second.reset();
	};

	const size_t operator++(int){
	//	Increments the count for the 'totals' counts for which there is no 'type' information 
	//	such as notUnique and loqQualString for which we use a dummy 'blank' entry 
		return at(blankString)++;
	}

};




//	A nested string map for holding counts for each identifier (e.g. exon, gene) for each gene  
class geneCountsClass : public std::map<const std::string, geneAttribute >
{
public:
	//	For printing out the list of counts.  entry.second is the count data one entry per attribute being investigated
	void print(TsvFile & output,const std::string & index,double norm,size_t length)
	{
		at(index).print(output,index, norm,length);
	};
	void print(TsvFile & output, const std::string & index, double norm = 1)
	{
		at(index).print(output, index, norm);
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
	rna_pos_type RNAstartPos,RNAendPos;
	bool strict;
	//	Store references here as these are created and deleted lots and these removes the need to allocate and
	//	deallocate on the heap.
	const std::string & geneName;
	const std::string & featType;
	gtfOverlap(size_t start, size_t finish, rna_pos_type RNAstartPos, rna_pos_type RNAendPos,bool strict,const std::string & geneName,const std::string & featType):
		start(start),finish(finish), RNAstartPos(RNAstartPos), RNAendPos(RNAendPos),strict(strict),geneName(geneName),featType(featType)
	{
	};
};

class gtfRegion
{
public:
	size_t start,finish;
	rna_pos_type RNAstart;
	//	Store actual values here as these are only created once and then referenced lots, so this is more efficient
	const std::string name;
	const std::string type;
	char strand;
	chromosomeGtfData::iterator * overlaps;

	gtfRegion(gtfRegion && gtf) : start(gtf.start),finish(gtf.finish), RNAstart(gtf.RNAstart),name(std::move(gtf.name)),type(std::move(gtf.type)),strand(gtf.strand),
		overlaps(gtf.overlaps)
	{
		gtf.overlaps = 0;
	}

	gtfRegion(	size_t start, size_t finish,const std::string & name,char strand,const std::string & type );
	~gtfRegion();
	void checkOverlap(const region & segment,std::vector<gtfOverlap> & overlapList) const;
};


typedef std::map<std::string,chromosomeGtfData> genomeGtfRegions;
typedef std::multimap<size_t,chromosomeGtfData::iterator> chromosomeEndIndexMap;

typedef std::map<std::string,chromosomeEndIndexMap > genomeEndIndexMap;

typedef std::vector<gtfRegion *> gtfRegionList;

class geneData
{
public:
	gtfRegionList regions;
	bool overlaps;
	//	The total length of the regions associated with the gene
	rna_pos_type length;
	char strand;
	double normFactor;
	void addRegion(gtfRegion * newRegion,bool ol)
	{
		if (ol)
			overlaps = true;

		newRegion->RNAstart = length + 1;
		if (regions.size() == 0)
			strand = newRegion->strand;

		regions.push_back(newRegion);
		length += (newRegion->finish - newRegion->start + 1);

	};
	geneData() : overlaps(false),length(0), strand(' '), normFactor(1){};
};


class featureFileEx : public featureFile
{
public: 
	//	A container of all teh consolidated gtf regions
	genomeGtfRegions genomeGtfData; 
	//	A map of the ends of the gtf regions.   Used for finding overlaps
	genomeEndIndexMap genomeEndIndex;

	//	For doing comparison run with a specific set of genes
	setEx<std::string> geneSet;
	std::vector<std::string> geneList;


	//	A map of the regions associated with a gene
	std::map<std::string,geneData> genes;

	void index(geneCountsClass & geneCounts);
	void outputChromData(const std::string & filename);
	void useSelectedGenes(const std::string & filename);

};

#endif
