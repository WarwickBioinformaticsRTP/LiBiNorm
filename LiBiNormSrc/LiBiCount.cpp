#include <stdlib.h>
#include <crtdbg.h>
#include <ctime>
#include "libCommon.h"
#include "stringEx.h"
#include "containerEx.h"
#include "api/BamReader.h"
#include "LiBiCount.h"

using namespace std;
using namespace BamTools;

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
			thisChromData.emplace(i->first,region(i->second.finish,i->second.tags[0].val,i->second.strand));

		}
		//	And now produce overlap list
		for (chromosomeData::iterator i = thisChromData.begin(); i != thisChromData.end();i++)
		{
			for (chromosomeData::iterator j = i;(j != thisChromData.end()) && (j->first < i->second.finish);j++)
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



int LiBiCount::main(int argc, char **argv)
{
	stringEx bamFileName,gtfFileName,
		id_attribute = "gene_id";

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

/*	BamReader reader;
	if ( !reader.Open(bamFileName) ) 
		exitFail("Could not open input BAM files: ",bamFileName);

	// retrieve 'metadata' from BAM files.
	const SamHeader header = reader.GetHeader();
	const RefVector references = reader.GetReferenceData();
*/

	gtfFileEx genomeDef;

    clock_t begin = clock();



	genomeDef.open(gtfFileName,id_attribute,feature_type);
	cout << "Data read";

	genomeDef.index();

	genomeDef.outputChromData(gtfFileName.replaceSuffix(".txt"));

    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

	cout << "Elapsed time " << elapsed_secs;

	string test;
	cin >> test;

	return EXIT_SUCCESS;
}
