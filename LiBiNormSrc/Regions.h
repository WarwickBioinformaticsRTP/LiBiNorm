#ifndef REGIONS_H
#define REGIONS_H

#include<map>
#include<string>
#include "libCommon.h"
#include "printEx.h"
#include "api/BamReader.h"

//	Wrap the cigar data so that the parser and printVal methods will recognise it
class Cigar : public std::vector<BamTools::CigarOp>
{
public:
	Cigar(){};
	Cigar(const std::vector<BamTools::CigarOp> & co): std::vector<BamTools::CigarOp>(co) {};
};

//	This holds the data from a bam entry that we are actually interested in.  It is the basis
//	of 'in program' persisted data and also data that is persisted to cache files.  Does not include the
//	 read name as this is stored elesewhere (e.g as the index of a map of readData
class readData
{
	public:
		int refId;
		int position;
		char strand;
		Cigar cigar;
		int NH;
		int qual;
		readData(void){}

		//	This constructor creates the readData from the bam file entry.  This means that methos expecting 
		//	readData can be passed a bamAlignment.  Use an rValue constructor so that we can 'swallow up' the cigar data 
		//	rather than making a copy of it as once the readData has been created we will have no further use
		//	for the cigar data
		readData(BamTools::BamAlignment && ba):
			refId(ba.RefID),
			position(ba.Position+1),
			qual(ba.MapQuality),
			cigar(move(ba.CigarData))
		{
			//	This simulates line 155 in count py.  If there is no optional NH field in the first read
			//	then the pythin throws an error, which we simulate by setting NH to zero
			if (!ba.GetTag("NH",NH))
				NH = ba.IsFirstMate()?0:-1;
			if (ba.IsPaired())
			{
				if (ba.IsReverseStrand() == ba.IsFirstMate())
					strand = '-';
				else
					strand = '+';
			}
			else
				strand = ba.IsReverseStrand() ? '-' : '+';
		};
		readData(const BamTools::BamAlignment & ba):
			refId(ba.RefID),
			position(ba.Position+1),
			qual(ba.MapQuality),
			cigar(ba.CigarData)
		{
			//	This simulates line 155 in count py.  If there is no optional NH field in the first read
			//	then the pythin throws an error, which we simulate by setting NH to zero
			if (!ba.GetTag("NH",NH))
				NH = ba.IsFirstMate()?0:-1;
			if (ba.IsPaired())
				strand = (ba.IsReverseStrand() == ba.IsFirstMate())?'-':'+';
			else
				strand = ba.IsReverseStrand() ? '-' : '+';

		};
};

//	Used for reading back cached read information from cache files
class cacheEntry : public readData
{
	std::string fname;
	std::ifstream * file;
public:
	//	Holds the name of the read, which is not in the readData class
	std::string name;

	cacheEntry() : file (0) {};
	~cacheEntry(); 


	bool open(const std::string filename);
	bool readNext();
	void close();
};

//	A reagion within a chromosome
class region 
{
public:
	int start,end;
	char strand;
	region(int start,int end,char strand) : start(start),end(end),strand(strand) {};
};


//	A set of regions on one chromosome
class regionList 
{
public:
	regionList(){};

	//	Which is normally created from a read
	regionList(const readData & read);

	//	Map of the regions, indexed by the location on teh chromosome
	//	Separate entries for and -ve strands so needs to be a multimap to cater 
	//	for + and - entries starting at the same location (usually an artefact)
	std::multimap<int,region> data;


	//	For combining data from a second read
	void combine(const regionList & rl);

private:
		//	For combining each individual read within a regionList
	void combineRegion(const region & r);

};

//	The list of regions associated with a read pair
class regionLists
{
public:
	//	Contains the regions themselves in a map indexed by chromosome (indicated by refId, the chromosome identifier
	//	in the bam file.  Done this way to cater for a read pair where the reads are on difference chromosomes
	std::map<int,regionList> data;
	//	The name of the read
	const std::string & name;
	int NH;
	int qual;
	std::vector<char> strands;

	//	Creates a regionList from one of the reads, either from a bam entry or from cachedData.  Use emplace so that the
	//	regionList can be efficiently placed straight into the map.
	regionLists(const readData & read,const std::string & name) :name(name),NH(read.NH),qual(read.qual), strands(1,read.strand){
		data.emplace(read.refId,regionList(read));
//		strands.emplace_back(read.strand);
	};

	//	Adds the information associated with the second read, which will be placed in the existing chromosome
	//  or added to a new.
	void combine(const readData & read){
			data[read.refId].combine(regionList(read));
			if ((read.NH == 0) || (NH == 0))
				NH = 0;
			else if (read.NH > NH)
				NH = read.NH;
			if (read.qual < qual)
				qual = read.qual;
			strands.emplace_back(read.strand);
	};
};

//	Declare the availability of methods that are used by parser for parser cigar strings
//	and also the methods used for printing a cacheentry and the cigar data
//	These are both used for the temporary cache data that is placed on disk
namespace parserInternal
{
	void parseval(const char *& start, Cigar & cigar,size_t & len);
}

bool printVal(outputDataFile * f,const Cigar & cigar);
bool printVal(outputDataFile * f,const readData & read);

#endif

