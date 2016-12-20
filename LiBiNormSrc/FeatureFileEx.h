
#ifndef FEATURE_FILE_EX_HEADER
#define FEATURE_FILE_EX_HEADER

#include "containerEx.h"
#include "libCommon.h"
#include "Regions.h"
#include "featureFile.h"
#include "dataVec.h"
#include "GeneCountData.h"

class featureRegion;


/*
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
//	This is used particularly when counts are being taken for more than one attribute

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
*/
typedef std::multimap<size_t, featureRegion> chromosomeFeatureData ;

struct featureOverlap
{
	featureOverlap(size_t start, size_t finish, rna_pos_type RNAstartPos, rna_pos_type RNAendPos, bool strict, const std::string & geneName, const std::string & featType) :
		start(start), finish(finish), RNAstartPos(RNAstartPos), RNAendPos(RNAendPos), strict(strict), geneName(geneName), featType(featType)
	{};

	size_t start,finish;
	rna_pos_type RNAstartPos,RNAendPos;
	bool strict;
	//	Store references here as these are created and deleted lots and these removes the need to allocate and
	//	deallocate on the heap.
	const std::string & geneName;
	const std::string & featType;
};

class featureRegion
{
public:
	featureRegion(featureRegion && gtf) : start(gtf.start), finish(gtf.finish), RNAstart(gtf.RNAstart), name(std::move(gtf.name)), type(std::move(gtf.type)), strand(gtf.strand),
		overlaps(gtf.overlaps)
	{
		gtf.overlaps = 0;
	}
	featureRegion(size_t start, size_t finish, const std::string & name, char strand, const std::string & type);
	~featureRegion();

	void checkOverlap(const region & segment, std::vector<featureOverlap> & overlapList) const;

	size_t start,finish;
	rna_pos_type RNAstart;
	//	Store actual values here as these are only created once and then referenced lots, so this is more efficient
	const std::string name;
	const std::string type;
	char strand;
	chromosomeFeatureData::iterator * overlaps;
};


typedef std::map<std::string,chromosomeFeatureData> genomeFeatureRegions;
typedef std::multimap<size_t,chromosomeFeatureData::iterator> chromosomeEndIndexMap;

typedef std::map<std::string,chromosomeEndIndexMap > genomeEndIndexMap;

typedef std::vector<featureRegion *> featureRegionList;

class geneData
{
public:
	geneData() : overlaps(false), length(0), strand(' '), normFactor(1) {};

	void addRegion(featureRegion * newRegion,bool ol)
	{
		if (ol)
			overlaps = true;

		newRegion->RNAstart = length + 1;
		if (regions.size() == 0)
			strand = newRegion->strand;

		regions.push_back(newRegion);
		length += (newRegion->finish - newRegion->start + 1);

	};

	featureRegionList regions;
	bool overlaps;
	//	The total length of the regions associated with the gene
	rna_pos_type length;
	char strand;
	double normFactor;
};


class featureFileEx : public featureFile
{
public: 
	void index(GeneCountData & geneCounts);
	void outputChromData(const std::string & filename);


	//	A container of all the consolidated feature regions
	genomeFeatureRegions genomeGtfData; 
	//	A map of the ends of the feature regions.   Used for finding overlaps
	genomeEndIndexMap genomeEndIndex;

	//	For doing comparison run with a specific set of genes
	setEx<std::string> geneSet;
	std::vector<std::string> geneList;


	//	A map of the regions associated with a gene
	std::map<std::string,geneData> genes;


};

#endif
