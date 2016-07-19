#include <stdlib.h>
#include <crtdbg.h>
#include <ctime>
#include "libCommon.h"
#include "stringEx.h"
#include "containerEx.h"
#include "LiBiCount.h"


bool regionList::combineRegion(size_t start,size_t end)
{
	bool combined = false;
	size_t i;

	for (i = 0;i < size();i++)
	{
		if (start < at(i).second)
		{
			if (end > at(i).second)
				at(i).second = end;
			if (start < at(i).first)
				at(i).first = start;
			combined = true;
			break;
		}
	}
	if (combined)
	{
		if (i >= (size() -1))
			return false;
	}
	else
	{
		emplace_back(start,end);
		return false;
	}
	return true;
}

void regionList::add(size_t start,size_t end)
{
	if (size())
	{
		if (end < at(0).second)
			insert(begin(),region(start,end));
		else
			emplace_back(start,end);
	}
	else
		emplace_back(start,end);
}

bool gtfRegion::checkOverlap(const region & segment,bool & strict)
{
	if ((segment.first >= start) && (segment.second <= finish))
	{
		strict = true;
		return true;
	}
	else if (((segment.first < finish) && (segment.first >= start)) ||
				((segment.second < finish) && (segment.second >= start)))
	{
		strict = false;
		return true;
	}

	return false;
}


void regionList::GetRegions(const BamAlignment & ba) {

	//	If we already have some regions then we need to combine them
	bool combine = size();
	strand = ba.IsReverseStrand()?'-':'+';


	auto & CigarData = ba.CigarData;

	// initialize alignment end to starting position

	size_t start = ba.Position +1;
	size_t end = start;

	// iterate over cigar operations
	vector<CigarOp>::const_iterator cigarIter = CigarData.begin();
	vector<CigarOp>::const_iterator cigarEnd  = CigarData.end();
	for ( ; cigarIter != cigarEnd; ++cigarIter) {
		const CigarOp& op = (*cigarIter);

		switch ( op.Type ) {

			// increase end position on CIGAR chars [DMXN=]
			case Constants::BAM_CIGAR_DEL_CHAR      :
			case Constants::BAM_CIGAR_MATCH_CHAR    :
			case Constants::BAM_CIGAR_MISMATCH_CHAR :
			case Constants::BAM_CIGAR_SEQMATCH_CHAR :
				end += op.Length;
				break;

			case Constants::BAM_CIGAR_INS_CHAR :
				break;

			case Constants::BAM_CIGAR_REFSKIP_CHAR  :
				{
					if (combine)
						combine = combineRegion(start,end-1);
					else
						emplace_back(start,end-1);
					end = start = (end + op.Length);
					break;
				}

			default :
				break;
		}
	}
	if (combine)
		combine = combineRegion(start,end-1);
	else
		emplace_back(start,end-1);
}



void gtfFileEx::index()
{

	for (auto & chrom : entryMap)
	{
		//	In each chromosome
		chromosomeData & thisChromData = chromData[chrom.first];
		for (auto i = chrom.second.begin(); i != chrom.second.end();i++)
		{
			_DBG(size_t start = i->second.start;)
			size_t finish = i->second.finish;
			setEx<string> type(i->second.type);
			for (auto j = next(i,1);(j != chrom.second.end()) && (j->first < finish);)
			{
				auto k = j++;
				if (i->second.tags[0].val == k->second.tags[0].val)
				{
					if (k->second.finish > finish)
					{
						finish = k->second.finish;
						type.add(k->second.type);
						chrom.second.erase(k);
					}
					else if (k->second.finish <= finish)
					{
						type.add(k->second.type);
						chrom.second.erase(k);
					}
				}
			}


			string & attName = i->second.tags[0].val;
			thisChromData.emplace(i->first,gtfRegion(i->second.start,finish,i->second.tags[0].val,i->second.strand,move(type)));

		}
		//	And now produce overlap list:  The list of all regions that start before this region has ended.
		for (chromosomeData::iterator i = thisChromData.begin(); i != thisChromData.end();i++)
		{
			for (chromosomeData::iterator j = next(i,1);(j != thisChromData.end()) && (j->first < i->second.finish);j++)
				j->second.overlaps.push_back(&i->second);
		}
	}
}

void gtfFileEx::outputChromData(const string & filename)
{
	TsvFile output;
	output.open(filename);

	for(auto i : chromData)
	{
		for (auto j : i.second)
			output.printEnd(i.first,j.first,j.second.finish,j.second.strand,j.second.name,j.second.type);
	}
}

void gtfFileEx::outputGeneCounts(const string & filename)
{
	TsvFile output;
	output.open(filename);

	for(auto i : geneCounts)
	{
		output.printEnd(i.first,i.second);
	}

}


void gtfFileEx::addRead(const string & chromosome,const regionList regions,mode countMode)
{
	chromosomeData & thisChromData = chromData[chromosome];

	auto i = thisChromData.lower_bound(regions[0].first);

	if (i == thisChromData.end())
	{
		i--;
		if(i->second.start > regions.back().second)
		{
			geneCounts["__no_feature"]++;
			return;
		}
	}

	while ((i->second.finish > regions[0].first) && (i != thisChromData.begin()))
		i--;

	map<string,bool> genes;

	if (regions.back().second > i->second.start)
	{
		//It overlaps
		for (auto & segment : regions)
		{
			auto j = i;
			while ((j != thisChromData.end()) && (j->first < segment.second))
			{
				bool strict;
				if (!useStrand || ((j->second.strand == regions.strand) != reverseStrand))
				{
					if (j->second.checkOverlap(segment,strict))
					{
						auto g = genes.find(j->second.name);
						if (g != genes.end())
						{
							if (!strict)
								g->second = false;
						}
						else
							genes.emplace(j->second.name,strict);
						i = j;
						break;
					}
				}
				j++;
			}
		}
	}

	switch(countMode)
	{
	case intersect_union:
		if (genes.size() == 1)
			geneCounts[genes.begin()->first]++;
		else if (genes.size() > 1)
			geneCounts["__ambiguous"]++;
		else
			geneCounts["__no_feature"]++;
		break;
	case intersect_strict:
	case intersect_nonempty:
		{
			vector<string> strictNames,nonemptyNames;
			for (auto & gene : genes)
			{
				if (gene.second)
					strictNames.push_back(gene.first);
				else if (gene.second)
					nonemptyNames.push_back(gene.first);
			}
			if (strictNames.size() == 1)
				geneCounts[strictNames[0]]++;
			else if ((nonemptyNames.size() ==1) && (countMode == intersect_nonempty))
				geneCounts[nonemptyNames[0]]++;
			break;
		}
	}
}

int LiBiCount::main(int argc, char **argv)
{
	stringEx bamFileName,gtfFileName,
		id_attribute = "gene_name";

	setEx<string> feature_type("exon");

	gtfFileEx genomeDef;

	genomeDef.reverseStrand = false;
	genomeDef.useStrand = true;
	mode countMode = intersect_union;
//	setEx<string> feature_type("exon","five_prime_utr","three_prime_utr");

	if(argc < 1)
	{
		printf("Error: parameter wrong!\n");
		return EXIT_FAILURE;
	}
	else if(argc == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("     LiBiNorm:    RNA-seq library bias normalisation   \n");
		return EXIT_SUCCESS;
	}

	int ni = 1;
	while(ni < argc)
	{
		if(strcmp(argv[ni], "-b") == 0)
		{
			bamFileName = argv[++ni];
		}
		else if(strcmp(argv[ni], "-g") == 0)
		{
			gtfFileName = argv[++ni];
		}
		else
		{
			cout << "Invalid parameter";
		}
		ni++;
	}

	BamReader reader;
	if ( !reader.Open(bamFileName) ) 
		exitFail("Could not open input BAM files: ",bamFileName);

	// retrieve 'metadata' from BAM files.
	const SamHeader header = reader.GetHeader();
	const RefVector references = reader.GetReferenceData();



	clock_t begin = clock();



	genomeDef.open(gtfFileName,id_attribute,feature_type);
	cout << "Data read";

	genomeDef.index();
	genomeDef.outputChromData(gtfFileName.replaceSuffix(".txt"));


	BamAlignment ba1;
	BamAlignment ba2;
	size_t currStart;
	string currName;

	while (reader.GetNextAlignment(ba1))
	{
		regionList regions;
		bool OK = true;
		
		regions.GetRegions(ba1);

		if (ba1.IsFirstMate())
		{
			reader.GetNextAlignment(ba2);

			if (genomeDef.useStrand && (ba1.IsReverseStrand() != ba2.IsMateReverseStrand() != genomeDef.reverseStrand))
				OK = false;
			else
				regions.GetRegions(ba2);
		}

		if (OK)
			genomeDef.addRead(references[ba1.RefID].RefName,regions,countMode);
		else
			genomeDef.geneCounts["__no_feature"]++;

	}

	genomeDef.outputGeneCounts(bamFileName.replaceSuffix(".counts.txt"));


	clock_t end = clock();
	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

	cout << "Elapsed time " << elapsed_secs;

//	string test;
//	cin >> test;

	return EXIT_SUCCESS;
}
