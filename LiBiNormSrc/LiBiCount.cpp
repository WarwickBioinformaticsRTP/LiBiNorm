#include "LiBiCount.h"
#include <stdlib.h>
#include <crtdbg.h>
#include "libCommon.h"
#include "containerEx.h"
#include "api/BamReader.h"

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
			if (entries[entry].type == feature)
			{
				string attName;
				for (auto & tag : entries[entry].tags)
				{
					if (tag.type == attribute_id)
					{
						attName = tag.val;
						break;
					}
				}
				chromosomeData::iterator i = thisChromData.find(entries[entry].start);
				bool duplicate = false;

				while (i != thisChromData.end() && (i->first == entries[entry].start))
				{
					if ((i->second.finish == entries[entry].finish) &&
						(i->second.name == attName))
					{
						duplicate = true;
						break;
					}
					i++;
				}
				if (!duplicate)
				{
					thisChromData.emplace(entries[entry].start,region(entries[entry].finish,attName));
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
		feature = "exon",
		attribute = "gene_id";

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

	genomeDef.open(gtfFileName);
	cout << "Data read";

	genomeDef.index(feature,attribute);

	string test;
	cout << "Finished";
	cin >> test;


	return EXIT_SUCCESS;
}
