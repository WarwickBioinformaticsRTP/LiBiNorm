#include <stdlib.h>
#include <crtdbg.h>
#include <ctime>
#include "libCommon.h"
#include "stringEx.h"
#include "containerEx.h"
#include "LiBiCount.h"


bool regionList::combine(size_t start,size_t end)
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

bool gtfRegion::checkOverlap(const region & segment)
{
	if ((segment.first >= start) && (segment.second <= finish))
		return true;
	else 
		return false;
}


void GetRegions(const BamAlignment & ba,regionList & regions) {

	//	If we already have some regions then we need to combine them
	bool combine = regions.size();


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
						combine = regions.combine(start,end-1);
					else
						regions.emplace_back(start,end-1);
					end = start = (end + op.Length);
					break;
				}

            default :
                break;
        }
    }
	if (combine)
		combine = regions.combine(start,end-1);
	else
		regions.emplace_back(start,end-1);
}



void gtfFileEx::index()
{

	for (auto & chrom : entryMap)
	{
		//	In each chromosome
		chromosomeData & thisChromData = chromData[chrom.first];
		for (auto i = chrom.second.begin(); i != chrom.second.end();i++)
		{
			size_t finish = i->second.finish;
			for (auto j = next(i,1);(j != chrom.second.end()) && (j->first < finish);)
			{
				auto k = j++;
				if (i->second.tags[0].val == k->second.tags[0].val)
				{
					if (k->second.finish > finish)
					{
						finish = k->second.finish;
						chrom.second.erase(k);
					}
					else if (k->second.finish <= finish)
					{
						chrom.second.erase(k);
					}
				}
			}


			string & attName = i->second.tags[0].val;
			thisChromData.emplace(i->first,gtfRegion(i->second.start,i->second.finish,i->second.tags[0].val,i->second.strand));

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
			output.printEnd(i.first,j.first,j.second.finish,j.second.strand,j.second.name);
	}

}

void gtfFileEx::addRead(const string & chromosome,const regionList regions)
{
		chromosomeData & thisChromData = chromData[chromosome];

		auto i = thisChromData.lower_bound(regions[0].first);

		if (i == thisChromData.end())
			return;

		while ((i->second.finish > regions[0].first) && (i != thisChromData.begin()))
			i--;

		set<string> genes;

		if (regions.back().second > i->second.start)
		{
			//It overlaps
			for (auto & segment : regions)
			{
				auto j = i;
				while ((j->first < segment.second) && (j != thisChromData.end()))
				{
					if (j->second.checkOverlap(segment))
					{
						genes.emplace(j->second.name);
						i = j;
						break;
					}
					j++;
				}
			}
		}

}



int LiBiCount::main(int argc, char **argv)
{
	stringEx bamFileName,gtfFileName,
		id_attribute = "gene_name";

	setEx<string> feature_type("exon");

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


	gtfFileEx genomeDef;

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
		
		GetRegions(ba1,regions);
		if (ba1.IsFirstMate())
		{
			reader.GetNextAlignment(ba2);
			GetRegions(ba2,regions);
		}
		genomeDef.addRead(references[ba1.RefID].RefName,regions);
	
	}


    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

	cout << "Elapsed time " << elapsed_secs;

	string test;
	cin >> test;

	return EXIT_SUCCESS;
}
