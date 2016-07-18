#include <stdlib.h>
#include <crtdbg.h>
#include <ctime>
#include "libCommon.h"
#include "containerEx.h"
#include "api/BamReader.h"
#include "LiBiCount.h"

using namespace std;
using namespace BamTools;

void gtfFileEx::index(const string & feature,const string & attribute)
{
	size_t attribute_id = gtfEntryTags::idMap[attribute];

	for (auto & chrom : entryMap)
	{
		chromosomeData & thisChromData = chromData[chrom.first];
		for (auto entry : chrom.second)
		{
			if (entry.type == feature)
			{
				string attName;
				for (auto & tag : entry.tags)
				{
					if (tag.type == attribute_id)
					{
						attName = tag.val;
						break;
					}
				}
				chromosomeData::iterator i = thisChromData.find(entry.start);
				bool duplicate = false;

				while (i != thisChromData.end() && (i->first == entry.start))
				{
					if ((i->second.finish == entry.finish) &&
						(i->second.name == attName))
					{
						duplicate = true;
						break;
					}
					i++;
				}
				if (!duplicate)
				{
					thisChromData.emplace(entry.start,region(entry.finish,attName));
				}
			}
		}
		//	And now produce overlap list
		for (chromosomeData::iterator i = thisChromData.begin(); i != thisChromData.end();i++)
		{
			for (chromosomeData::iterator j = i;(j != thisChromData.end()) && (j->first < i->second.finish);j++)
				j->second.overlaps.push_back(&i->second);
		}
	}
}




int LiBiCount::main(int argc, char **argv)
{
	string bamFileName,gtfFileName,
		id_attribute;// = "gene_id";

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

//	genomeDef.index(feature,attribute);

    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

	cout << "Elapsed time " << elapsed_secs;

	string test;
	cin >> test;

	return EXIT_SUCCESS;
}
