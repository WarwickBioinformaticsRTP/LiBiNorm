
#include "api/BamReader.h"

#include "LiBiTools.h"

#include "libCommon.h"
#include "stringEx.h"
#include "fastaFile.h"
#include "FeatureFileEx.h"
#include "GeneCountData.h"

using namespace BamTools;

#define LOOKAHEAD 1000

int LiBiTools::landMain(int argc, char **argv)
{
	stringEx land_filename1, land_filename2, gff_filename, bamFileName;

	if (argc < 1)
	{
		printf("Error: parameter wrong!\n");
		return EXIT_FAILURE;
	}
	else if ((argc == 1) || ((argc == 2) && ((strcmp(argv[1], "-h") == 0) || (strcmp(argv[1], "--help") == 0))))
	{
		printf("Usage: LiBiNorm test -g gff_file landscapefile1 landscapefile2 gff_file\n");
		printf("This program compares two landscape files\n");
		return EXIT_SUCCESS;
	}
	if (argc < 3)
		exitFail("Insufficient arguments");
	int ni = 1;
	while (ni < argc - 2)
	{
		if (strcmp(argv[ni], "-g") == 0)
		{
			gff_filename = argv[++ni];
		}
		else if (strcmp(argv[ni], "-b") == 0)
		{
			bamFileName = argv[++ni];
		}
		else
			exitFail("Unrecognised option", argv[ni]);
		ni++;
	}
	land_filename1 = argv[argc - 2];
	land_filename2 = argv[argc - 1];

	featureFileEx genomeDef;
	stringEx id_attribute = DEFAULT_GFF_ID_ATTRIBUTE,
		feature_type = DEFAULT_FEATURE_TYPE_EXON;

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

	if (!genomeDef.open(gff_filename, id_attribute, feature_type))
		exitFail("Could not open feature file: ", gff_filename);


	GeneCountData geneCounts1, geneCounts2;
	geneCounts1.loadData(land_filename1);
	geneCounts2.loadData(land_filename2);
	genomeDef.index(geneCounts1);

	TsvFile resFile;
	resFile.open("combined.txt");

	auto it1 = geneCounts1.readPositionData.begin();
	auto it2 = geneCounts2.readPositionData.begin();

	while ((it1 != geneCounts1.readPositionData.end()) && (it2 != geneCounts2.readPositionData.end()))
	{
		std::map<std::string, geneData>::iterator geneIterator;

		bool match = it1->first == it2->first;
		if (!match)
		{
			auto it1a = it1;
			for (int i = 0; (i < LOOKAHEAD) && (!match) && (++it1a != geneCounts1.readPositionData.end()); i++)
			{
				if (it1a->first == it2->first)
				{
					for (int j = 0; j <= i; j++)
					{
						if ((geneIterator = genomeDef.genes.find(it1->first)) != genomeDef.genes.end())
						{
							resFile.print(1, _s("chr", geneIterator->second.chromosome, ":", geneIterator->second.regions[0]->start), it1->first, it1->second.positions[0]);
							resFile.print(1, _s("chr", geneIterator->second.chromosome, ":", geneIterator->second.regions[0]->start), it1->first, it1->second.positions[1]);
							resFile.print();
						}
						it1++;
					}
					match = true;
				}
			}
			if (!match)
			{
				auto it2a = it2;
				for (int i = 0; (i < LOOKAHEAD) && (!match) && (++it2a != geneCounts2.readPositionData.end()); i++)
				{
					if (it2a->first == it1->first)
					{
						for (int j = 0; j <= i; j++)
						{
							if ((geneIterator = genomeDef.genes.find(it2->first)) != genomeDef.genes.end())
							{
								resFile.print(2, _s("chr", geneIterator->second.chromosome, ":", geneIterator->second.regions[0]->start), it1->first, it2->second.positions[0]);
								resFile.print(2, _s("chr", geneIterator->second.chromosome, ":", geneIterator->second.regions[0]->start), it2->first, it2->second.positions[1]);
								resFile.print();
							}
							it2++;
						}
						match = true;
					}
				}
			}
		}
		if ((geneIterator = genomeDef.genes.find(it1->first)) != genomeDef.genes.end())
		{
			resFile.print(1, _s("chr", geneIterator->second.chromosome, ":", geneIterator->second.regions[0]->start), it1->first, it1->second.positions[0]);
			resFile.print(1, _s("chr", geneIterator->second.chromosome, ":", geneIterator->second.regions[0]->start), it1->first, it1->second.positions[1]);
			if (!match)
				resFile.print();
		}
		if ((geneIterator = genomeDef.genes.find(it2->first)) != genomeDef.genes.end())
		{
			resFile.print(2, _s("chr",geneIterator->second.chromosome,":",geneIterator->second.regions[0]->start), it2->first, it2->second.positions[0]);
			resFile.print(2, _s("chr", geneIterator->second.chromosome, ":", geneIterator->second.regions[0]->start), it2->first, it2->second.positions[1]);
			resFile.print();
		}

		it1++;
		it2++;

	}
	return EXIT_SUCCESS;

}

int LiBiTools::geneMain(int argc, char **argv)
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
