#include <random>
#include <chrono>
#include "MakeFastq.h"
#include "fastaFile.h"

#ifdef _DEBUG
#define BAMCACHESIZE 100
#else
#define BAMCACHESIZE 5000
#endif

bamRead::bamRead(const BamAlignment & ba) 
{
	operator = (ba);
}

bamRead & bamRead::operator = (const BamAlignment & ba)
{
	name = ba.Name;
	if (ba.IsMapped() && ba.IsReverseStrand())
	{
		readSeq = sequence(ba.QueryBases).complement();
		qualData = sequence(ba.Qualities).reversed();
	}
	else
	{
		readSeq = ba.QueryBases;
		qualData = ba.Qualities;
	}
	return This;
}

void bamRead::setName(const BamAlignment & ba)
{
	name = ba.Name;
}

void bamRead::setValues(const bamRead & br)
{
	readSeq = br.readSeq;
	qualData= br.qualData;
}


void bamRead::addSNP(double errorRate)
{
	double randVal = (double)rand()/((double)RAND_MAX+1);
	if (randVal > (errorRate*readSeq.size()))
		return;

	size_t p1 = randVal * readSeq.size();
	size_t p2 = (double)rand()/((double)RAND_MAX+1) * readSeq.size();

	readSeq[p1] = readSeq[p2];

}


void bamRead::output(FILE * f)
{
	fprintf(f,"@%s\n%s\n+\n%s\n",name.c_str(),readSeq.c_str(),qualData.c_str());
}

void bamReadCache::clear()
{
	for (auto & i: This)
		i.resetNames();
}


MakeFastq::MakeFastq()
{
}

bool MakeFastq::getNextAlignmentCore()
{
	bool OK = reader.GetNextAlignmentCore(ba);
	if (!OK)
	{
		cerr << "Rewinding" << endl;
		reader.Rewind();
		OK = reader.GetNextAlignmentCore(ba);
	}
	return OK;
}

bool MakeFastq::getNextAlignment()
{
	bool OK = reader.GetNextAlignment(ba);
	if (!OK)
	{
		cerr << "Rewinding" << endl;
		reader.Rewind();
		OK = reader.GetNextAlignment(ba);
	}
	return OK;
}
bool MakeFastq::getNextForeignAlignment()
{
	bool OK = foreignReader.GetNextAlignment(foreignBa);
	if (!OK)
	{
		cerr << "Rewinding foreign bam file" << endl;
		foreignReader.Rewind();
		OK = foreignReader.GetNextAlignment(foreignBa);
	}
	return OK;
}


int MakeFastq::main(int argc, char **argv)
{
	unsigned seed = chrono::system_clock::now().time_since_epoch().count();
	default_random_engine generator(seed);
    uniform_real_distribution<double> distribution(0.0,1.0);

	int size=1500000;
		
	double	mitoRate=1,unmappedRate=1,mappedRate=1,foreignBamDataRate = 0;;

	int	overamplified = -1;

	size_t start = 0;

	stringEx bamFileName,outputFileRoot,foreignBamData;

	stringEx skipChromosome;

	if(argc < 1)
	{
		printf("Error: parameter wrong!\n");
		return EXIT_FAILURE;
	}
	else if ((argc == 1) || ((argc == 2) && ((strcmp(argv[1],"-h") ==0) || (strcmp(argv[1],"--help")==0))))
	{
printf("Usage: LiBiNorm makefastq [options] alignment_file\n");
printf("Options:\n");
printf("  -h, --help				show this help message and exit\n");
printf("  -s N, --size=N			size of fastq file (1000)\n");
printf("  -v N, --overamplified=N	degree of overamplification\n");
printf("  -k ab, --skip=ab			skip chromosomes begining in ab\n");
printf("\n");
printf("Written by Nigel Dyer (nigel.dyer@warwick.ac.uk)\n");
		return EXIT_SUCCESS;
	}

	FILE * f1out,* f2out;

	int ni = 1;

	if (argc < 1)
		exitFail("Insufficient arguments");

	while(ni < argc-1)
	{
		bool opt2 = false;

		if((strcmp(argv[ni], "-l") == 0) || (opt2 = (strncmp(argv[ni], "--length=",9) == 0)))
		{
			size = atoi(opt2?argv[ni]+9:argv[++ni]);
		}
		else if((strcmp(argv[ni], "-o") == 0) || (opt2 = (strncmp(argv[ni], "--output=",9) == 0)))
		{
			outputFileRoot = opt2?argv[ni]+9:argv[++ni];
		}
		else if((strcmp(argv[ni], "-s") == 0) || (opt2 = (strncmp(argv[ni], "--start=",8) == 0)))
		{
			start = atoi(opt2?argv[ni]+8:argv[++ni]);
		}
		else if((strcmp(argv[ni], "-i") == 0) || (opt2 = (strncmp(argv[ni], "--mito=",7) == 0)))
		{
			mitoRate = atof(opt2?argv[ni]+9:argv[++ni]);
		}
		else if((strcmp(argv[ni], "-u") == 0) || (opt2 = (strncmp(argv[ni], "--unmapped=",11) == 0)))
		{
			unmappedRate = atof(opt2?argv[ni]+11:argv[++ni]);
		}
		else if((strcmp(argv[ni], "-m") == 0) || (opt2 = (strncmp(argv[ni], "--mapped=",9) == 0)))
		{
			mappedRate = atof(opt2?argv[ni]+9:argv[++ni]);
		}
		else if((strcmp(argv[ni], "-v") == 0) || (opt2 = (strncmp(argv[ni], "--overamplified=",9) == 0)))
		{
			overamplified = atoi(opt2?argv[ni]+9:argv[++ni]);
		}
		else if((strcmp(argv[ni], "-k") == 0) || (opt2 = (strncmp(argv[ni], "--skip=",6) == 0)))
		{
			skipChromosome = opt2?argv[ni]+7:argv[++ni];
		}
		else if((strcmp(argv[ni], "-f") == 0) || (opt2 = (strncmp(argv[ni], "--foreign=",10) == 0)))
		{
			foreignBamData = opt2?argv[ni]+10:argv[++ni];
			foreignBamDataRate = atof(opt2?argv[ni]+11+foreignBamData.size():argv[++ni]);
		}
		else
		{
			exitFail("Invalid parameter: ",string(argv[ni]));
		}
		ni++;
	}

	bamFileName = argv[argc-1];

	if (!outputFileRoot)
		outputFileRoot = bamFileName.removeSuffix();

	f1out = fopen((outputFileRoot + "_1.fastq").c_str(),"wb");
	f2out = fopen((outputFileRoot + "_2.fastq").c_str(),"wb");

	if ( !reader.Open(bamFileName) ) 
		exitFail("Could not open input BAM file: ",bamFileName);


	if (foreignBamData)
	{
		if ( !foreignReader.Open(foreignBamData) ) 
			exitFail("Could not open input BAM file: ",foreignBamData);
		for (size_t i = 0;i < start;i++)
			getNextForeignAlignment();
	}
	

	// retrieve 'metadata' from BAM files.
	RefVector references = reader.GetReferenceData();
	SamHeader header = reader.GetHeader();

	int mitoRef = -1;
	for (size_t i = 0;i < references.size();i++)
	{
		if ((references[i].RefName == "MT") || (references[i].RefName == "chrMT"))
		{
			mitoRef = i;
			break;
		}
	}
	if (mitoRef == -1)
		exitFail("No mitochondrial gene found");

	bool OK = getNextAlignment();
	
	size_t buffIndex = 0;

	int count = 0;

	srand(100);

	bamReadCache oaBuffer(BAMCACHESIZE);
	size_t overAmplifyRounds = 0;

	cerr << "Skipping reads" << endl;
	for (size_t i = 0;i < start;i++)
		OK = getNextAlignmentCore();

	ba.BuildCharData();

	cerr << "Creating fastq files" << endl;


	twoBamReads  readPair,readForeignPair;

	while ((OK) && (count < size))
	{

		bool useThis = true;

		int NH;
		if (ba.GetTag("NH",NH))
		{
			if ((NH > 1) && !ba.IsPrimaryAlignment())
				useThis = false;
		}

		if (skipChromosome && (ba.RefID != -1) && (references[ba.RefID].RefName.substr(0,skipChromosome.size()) == skipChromosome))
		{
			useThis = false;
		}

		if (ba.IsFirstMate())
		{
			double r2 = distribution(generator);

			if (!ba.IsMapped() && !ba.IsMateMapped())
			{
				if (r2 >= unmappedRate)
					useThis = false;
			}
			else if ((ba.IsMapped() && (ba.RefID == mitoRef)) || (ba.IsMateMapped() && (ba.MateRefID == mitoRef)))
			{
				if (r2 >= mitoRate)
					useThis = false;
			}
			else
			{
				if (r2 >= mappedRate)
					useThis = false;
			}
		}

		if (useThis)
		{
			if (overamplified > 1)
			{
				ba.BuildCharData();
				oaBuffer[buffIndex].addRead(ba);

				if (oaBuffer[buffIndex].namesMatch())
						buffIndex++;

				if (buffIndex >= BAMCACHESIZE)
				{
					size_t i = 0;

					size_t amp = 1;
					switch (overAmplifyRounds++ % 6)
					{
					case 0: amp = 2; break;
					default:amp = overamplified; break;
					}


					while ((i++ < (BAMCACHESIZE * amp)) && (count < size))
					{
						if ((++count % 10000) == 0)
							cerr << "Record " << count << endl;

						//Find a newname for the overamplified record						
						do {
							OK = getNextAlignment();
							NH = ba.GetTag("NH",NH);
						}
						while (!((NH <= 1) && ba.IsFirstMate()) && !((NH > 1) && ba.IsPrimaryAlignment() && ba.IsFirstMate()));

						//	Pick one of the records, rename it, add some errors and output it
						double SNPrate = 1.01/400;
						size_t pos = distribution(generator) * (BAMCACHESIZE-1);

						twoBamReads tbr(oaBuffer[pos]);
						tbr.setName(ba);
						tbr.addSNP(SNPrate);
						tbr.output(f1out,f2out);

					}
					oaBuffer.clear();
					buffIndex = 0;
				}

			}
			else
			{
				readPair.addRead(ba);

				if (readPair.namesMatch())
				{

					if ((++count % 10000) == 0)
						cerr << "Record " << count << endl;

					double r2 = distribution(generator);
					if ((foreignBamData) && (r2 < foreignBamDataRate))
					{
						while (!readForeignPair.namesMatch())
						{
							readForeignPair.addRead(foreignBa);
							getNextForeignAlignment();
						}
						readPair.setValues(readForeignPair);
						readForeignPair.resetNames();
					}
					readPair.output(f1out,f2out);
					readPair.resetNames();

				}
			}
		}

		OK = getNextAlignment();
	}

	fclose(f1out);
	fclose(f2out);

	return EXIT_SUCCESS;
}