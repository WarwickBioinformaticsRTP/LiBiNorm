
#include <fstream>
#include "api/BamReader.h"


#include "libCommon.h"
#include "stringEx.h"
#include "fastaFile.h"

#include "Options.h"
#include "FeatureFileEx.h"
#include "GeneCountData.h"
#include "LiBiTools.h"

using namespace std;
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
		printf("Usage: LiBiNorm land -g gff_file -b bamfile landscapefile1 landscapefile2\n");
		printf("This program compares two landscape files and produces a single file that combines the data from both\n");
		printf("The gff file is used to provide cooedinates for the genes so that they can be viewed easily on IGV\n");
		printf("The bam file allows the correct chromosome names to be used, the mapping being done based on chromosome length\n");
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
	genomeDef.index(geneCounts1,false);

	TsvFile resFile,missingGenesFile,extraGenesFile;
	resFile.open("combined.txt");
	missingGenesFile.open("missingGenes.txt");
	extraGenesFile.open("extraGenes.txt");

	readPositionDataClass::Iterator it1 = geneCounts1.readPositionData.begin();
	readPositionDataClass::Iterator it2 = geneCounts2.readPositionData.begin();

	while ((it1 != geneCounts1.readPositionData.end()) && (it2 != geneCounts2.readPositionData.end()))
	{
		std::map<std::string, geneData>::iterator geneIterator;

		bool match = it1->first == it2->first;
		if (!match)
		{
			readPositionDataClass::Iterator it1a = it1;
			for (int i = 0; (i < LOOKAHEAD) && (!match) && (++it1a != geneCounts1.readPositionData.end()); i++)
			{
				if (it1a.geneName() == it2.geneName())
				{
					for (int j = 0; j <= i; j++)
					{
						if ((geneIterator = genomeDef.genes.find(it1.geneName())) != genomeDef.genes.end())
						{
							resFile.print(1, _s("chr", geneIterator->second.chromosome, ":", geneIterator->second.regions[0]->start), it1.geneName(), it1.geneAttributes().positions[0]);
							resFile.print(1, _s("chr", geneIterator->second.chromosome, ":", geneIterator->second.regions[0]->start), it1.geneName(), it1.geneAttributes().positions[1]);
							resFile.print();
							missingGenesFile.print(it2.geneName(), _s(it2.geneAttributes().positions[0].size(), " plus"), it2.geneAttributes().positions[0]);
							missingGenesFile.print(it2.geneName(), _s(it2.geneAttributes().positions[1].size(), " minus"), it2.geneAttributes().positions[1]);
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
					if (it2a.geneName() == it1.geneName())
					{
						for (int j = 0; j <= i; j++)
						{
							if ((geneIterator = genomeDef.genes.find(it2.geneName())) != genomeDef.genes.end())
							{
								resFile.print(2, _s("chr", geneIterator->second.chromosome, ":", geneIterator->second.regions[0]->start), it1.geneName(), it2.geneAttributes().positions[0]);
								resFile.print(2, _s("chr", geneIterator->second.chromosome, ":", geneIterator->second.regions[0]->start), it2.geneName(), it2.geneAttributes().positions[1]);
								resFile.print();
								extraGenesFile.print(it1.geneName(), _s(it1.geneAttributes().positions[0].size(), " plus"), it1.geneAttributes().positions[0]);
								extraGenesFile.print(it1.geneName(), _s(it1.geneAttributes().positions[1].size(), " minus"), it1.geneAttributes().positions[1]);
							}
							it2++;
						}
						match = true;
					}
				}
			}
		}
		if ((geneIterator = genomeDef.genes.find(it1.geneName())) != genomeDef.genes.end())
		{
			resFile.print(1, _s("chr", geneIterator->second.chromosome, ":", geneIterator->second.regions[0]->start), it1.geneName(), it1.geneAttributes().positions[0]);
			resFile.print(1, _s("chr", geneIterator->second.chromosome, ":", geneIterator->second.regions[0]->start), it1.geneName(), it1.geneAttributes().positions[1]);
			if (!match)
				resFile.print();
		}
		if ((geneIterator = genomeDef.genes.find(it2.geneName())) != genomeDef.genes.end())
		{
			resFile.print(2, _s("chr",geneIterator->second.chromosome,":",geneIterator->second.regions[0]->start), it2.geneName(), it2.geneAttributes().positions[0]);
			resFile.print(2, _s("chr", geneIterator->second.chromosome, ":", geneIterator->second.regions[0]->start), it2.geneName(), it2.geneAttributes().positions[1]);
			resFile.print();
		}

		it1++;
		it2++;

	}
	return EXIT_SUCCESS;

}

int LiBiTools::landMain2(int argc, char **argv)
{
	stringEx landFilename, geneFilename;

	if (argc < 1)
	{
		printf("Error: parameter wrong!\n");
		return EXIT_FAILURE;
	}
	else if ((argc == 1) || ((argc == 2) && ((strcmp(argv[1], "-h") == 0) || (strcmp(argv[1], "--help") == 0))))
	{
		printf("Usage: LiBiNorm land2 landscapefile geneFile\n");
		printf("This program extracts data from a landscape file\n");
		return EXIT_SUCCESS;
	}
	if (argc < 3)
		exitFail("Insufficient arguments");

	landFilename = argv[argc - 2];
	geneFilename = argv[argc - 1];

	featureFileEx genomeDef;
	stringEx id_attribute = DEFAULT_GFF_ID_ATTRIBUTE,
		feature_type = DEFAULT_FEATURE_TYPE_EXON;


	GeneCountData geneCounts;
	geneCounts.loadData(landFilename);


	ifstream file;
	file.open(geneFilename);

	if (!file.is_open())
	{
		progMessage("Unable to read gene list from ", geneFilename);
		return EXIT_FAILURE;
	}

	setEx<string> genes;

	string line, gene;
	while (!file.eof())
	{
		std::getline(file, line);
		parser(line, " \n\r", gene);
		genes.emplace(gene);
	};

	TsvFile output;
	output.open(landFilename.replaceSuffix(".subset.txt"));
	TsvFile output2;
	output2.open(landFilename.replaceSuffix(".unused.txt"));
	TsvFile genelist;
	genelist.open(landFilename.replaceSuffix(".extraGenes.txt"));

	for (size_t i = 1; i < geneCounts.info.size(); i++)
	{
#ifdef COUNT_IN_LANDSCAPE
		long count = geneCounts.counts[i];
#endif
		string & name = geneCounts.info[i].name;
		if (genes.contains(name))
		{
			long len = geneCounts.lengths[0][i];
			string c;
#ifdef COUNT_IN_LANDSCAPE
			c = _s(":", count);
#endif
			output.print(name, _s(len, c, " plus"), geneCounts.readPositionData[name].positions[0]);
			output.print(name, _s(len, c, " minus"), geneCounts.readPositionData[name].positions[1]);
		}
		else
		{
			long len = geneCounts.lengths[0][i];
			string c;
#ifdef COUNT_IN_LANDSCAPE
			c = _s(":", count);
#endif
			output2.print(name, _s(len, c, " plus"), geneCounts.readPositionData[name].positions[0]);
			output2.print(name, _s(len, c, " minus"), geneCounts.readPositionData[name].positions[1]);
			genelist.print(name);
		}

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
