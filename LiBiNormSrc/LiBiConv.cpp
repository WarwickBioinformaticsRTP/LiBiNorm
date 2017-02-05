#include "api/BamReader.h"
#include "libCommon.h"
#include "stringEx.h"
#include "fastaFile.h"
#include "LiBiConv.h"
#include "FeatureFileEx.h"

using namespace std;
using namespace BamTools;


int LiBiConv::main(int argc, char **argv)
{

	stringEx bamFileName, featureFileName;
	featureFileEx genomeDef;
	verbose = true;

	if (argc < 1)
	{
		printf("Error: parameter wrong!\n");
		return EXIT_FAILURE;
	}
	else if ((argc == 1) || ((argc == 2) && ((strcmp(argv[1], "-h") == 0) || (strcmp(argv[1], "--help") == 0))))
	{
		printf("Usage: LiBiNorm conv alignment_file gff_file\n");
		printf("This program takes an alignment file in BAM format and a feature file in\n");
		printf("GFF format and creates a feature file with the chromosome names converted to the names\n");
		printf("used for the chromosomes in the bam file\n");
		return EXIT_SUCCESS;
	}
	if (argc < 3)
		exitFail("Insufficient arguments");

	bamFileName = argv[argc - 2];
	featureFileName = argv[argc - 1];

	progMessage("Converting ", featureFileName);

	BamReader reader;
	RefVector references;

	if (!reader.Open(bamFileName))
		exitFail("Could not open input BAM files: ", bamFileName);


	// retrieve 'metadata' from BAM files.
	references = reader.GetReferenceData();
	for (auto & i : references)
	{
		if (strncasecmp(i.RefName.c_str(), "chr", 3) == 0)
			i.RefName = i.RefName.substr(3);
		genomeDef.addToChromosomeMap(i.RefLength, i.RefName);
	}

	if (!genomeDef.open(featureFileName, "","",true))
		exitFail("Could not open feature file: ", featureFileName);

	return EXIT_SUCCESS;

}

int LiBiConv::main2(int argc, char **argv)
{
	stringEx land_filename1, land_filename2;
	if (argc < 1)
	{
		printf("Error: parameter wrong!\n");
		return EXIT_FAILURE;
	}
	else if ((argc == 1) || ((argc == 2) && ((strcmp(argv[1], "-h") == 0) || (strcmp(argv[1], "--help") == 0))))
	{
		printf("Usage: LiBiNorm test landscapefile landscape file\n");
		printf("This program compares two landscape files\n");
		return EXIT_SUCCESS;
	}
	if (argc < 3)
		exitFail("Insufficient arguments");

	land_filename1 = argv[argc - 2];
	land_filename2 = argv[argc - 1];

	GeneCountData geneCounts1,geneCounts2;

	geneCounts1.loadData(land_filename1);
	geneCounts2.loadData(land_filename2);

	TsvFile resFile;
	resFile.open("commbined.txt");

	auto it1 = geneCounts1.readPositionData.begin();
	auto it2 = geneCounts2.readPositionData.begin();

	while ((it1 != geneCounts1.readPositionData.end()) && (it2 != geneCounts2.readPositionData.end()))
	{
		bool match = it1->first == it2->first;
		if (!match)
		{
			auto it1a = it1;
			for (int i = 0; (i < 10) && (!match) && (++it1a != geneCounts1.readPositionData.end()); i++)
			{
				if (it1a->first == it2->first)
				{
					for (int j = 0; j <= i; j++)
					{
						resFile.print(1, it1->first, it1->second.positions[0]);
						resFile.print(1, it1->first, it1->second.positions[1]);
						resFile.print();
						it1++;
					}
					match = true;
				}
			}
			if (!match)
			{
				auto it2a = it2;
				for (int i = 0; (i < 10) && (!match) && (++it2a != geneCounts2.readPositionData.end()); i++)
				{
					if (it2a->first == it1->first)
					{
						for (int j = 0; j <= i; j++)
						{
							resFile.print(2, it2->first, it2->second.positions[0]);
							resFile.print(2, it2->first, it2->second.positions[1]);
							resFile.print();
							it2++;
						}
						match = true;
					}
				}
			}
		}
		resFile.print(1, it1->first, it1->second.positions[0]);
		resFile.print(1, it1->first, it1->second.positions[1]);
		if (!match)
			resFile.print();
		resFile.print(2, it2->first, it2->second.positions[0]);
		resFile.print(2, it2->first, it2->second.positions[1]);
		resFile.print();

		it1++;
		it2++;

	}
	return EXIT_SUCCESS;

}

int LiBiConv::main3(int argc, char **argv)
{
	stringEx fastq_filename;
	if (argc < 1)
	{
		printf("Error: parameter wrong!\n");
		return EXIT_FAILURE;
	}
	else if ((argc == 1) || ((argc == 2) && ((strcmp(argv[1], "-h") == 0) || (strcmp(argv[1], "--help") == 0))))
	{
		printf("Usage: LiBiNorm genes fastqfile\n");
		printf("Extracts the gene names from a fasta file\n");
		return EXIT_SUCCESS;
	}
	if (argc < 2)
		exitFail("Insufficient arguments");

	fastq_filename = argv[argc - 1];

	fastaFileRead fasta;
	fasta.open(fastq_filename);

	TsvFile resFile;
	resFile.open(fastq_filename.replaceSuffix(".genes.txt"));

	bool OK = fasta.readEntry();
	while (OK)
	{
		resFile.print(fasta.NameStr());
		OK = fasta.readEntry();
	};
	return EXIT_SUCCESS;

}
