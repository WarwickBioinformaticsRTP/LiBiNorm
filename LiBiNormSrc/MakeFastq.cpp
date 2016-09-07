
#include "MakeFastq.h"
#include "fastaFile.h"

#ifdef _DEBUG
#define BAMCACHESIZE 100
#else
#define BAMCACHESIZE 10000
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

void bamRead::output(FILE * f)
{
	fprintf(f,"@%s\n%s\n+\n%s\n",name.c_str(),readSeq.c_str(),qualData.c_str());
}

void bamReadCache::clear()
{
	for (auto & i: This)
		i.name.clear();
}


MakeFastq::MakeFastq(void)
{
}


int MakeFastq::main(int argc, char **argv)
{

	int size=1500000,
		
		mitoRate=1000,unmappedRate=1000,mappedRate=1000,

		overamplified = -1;

	stringEx bamFileName,outputFileRoot;

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

		if((strcmp(argv[ni], "-s") == 0) || (opt2 = (strncmp(argv[ni], "--size=",7) == 0)))
		{
			size = atoi(opt2?argv[ni]+7:argv[++ni]);
		}
		else if((strcmp(argv[ni], "-o") == 0) || (opt2 = (strncmp(argv[ni], "--output=",9) == 0)))
		{
			outputFileRoot = opt2?argv[ni]+9:argv[++ni];
		}
		else if((strcmp(argv[ni], "-i") == 0) || (opt2 = (strncmp(argv[ni], "--mito=",7) == 0)))
		{
			mitoRate = atoi(opt2?argv[ni]+9:argv[++ni]);
		}
		else if((strcmp(argv[ni], "-u") == 0) || (opt2 = (strncmp(argv[ni], "--unmapped=",11) == 0)))
		{
			unmappedRate = atoi(opt2?argv[ni]+1:argv[++ni]);
		}
		else if((strcmp(argv[ni], "-m") == 0) || (opt2 = (strncmp(argv[ni], "--mapped=",9) == 0)))
		{
			mappedRate = atoi(opt2?argv[ni]+1:argv[++ni]);
		}
		else if((strcmp(argv[ni], "-v") == 0) || (opt2 = (strncmp(argv[ni], "--overamplified=",9) == 0)))
		{
			overamplified = atoi(opt2?argv[ni]+1:argv[++ni]);
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

	BamReader reader;

	if ( !reader.Open(bamFileName) ) 
		exitFail("Could not open input BAM file: ",bamFileName);

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

	BamAlignment ba;

	bool OK = reader.GetNextAlignmentCore(ba);
	
	size_t buffIndex = 0;

	int count = 0;

	vector<bamReadCache> oaBuffer(2,BAMCACHESIZE);
	while ((OK) && (count < size))
	{

		bool useThis = true;

		int NH;
		if (ba.GetTag("NH",NH))
		{
			if ((NH > 1) && !ba.IsPrimaryAlignment())
				useThis = false;
		}


		if (useThis)
		{
			if (overamplified > 1)
			{
				ba.BuildCharData();
				if (ba.IsFirstMate())
					oaBuffer[0][buffIndex] = ba;
				else
					oaBuffer[1][buffIndex] = ba;

				if (oaBuffer[0][buffIndex].name && oaBuffer[1][buffIndex].name)
				{
					if (oaBuffer[0][buffIndex].name == oaBuffer[1][buffIndex].name)
						buffIndex++;
				}

				if (buffIndex >= BAMCACHESIZE)
				{
					size_t i = 0;
					while ((i++ < (BAMCACHESIZE * overamplified)) && (count < size))
					{
						if ((++count % 10000) == 0)
							cerr << "Record " << count << endl;
						size_t pos = ((long)rand() * (BAMCACHESIZE-1))/RAND_MAX;
						oaBuffer[0][pos].output(f1out);
						oaBuffer[1][pos].output(f2out);
					}
					oaBuffer[0].clear();
					oaBuffer[1].clear();
					buffIndex = 0;
				}

			}
			else
			{

				bamRead read(move(ba));

				auto r = rand();
				auto r2 = (r * 1000.0)/RAND_MAX;

				bool output = false;
				if (!ba.IsMapped())
				{
					if (r2 < unmappedRate)
						output = true;
				}
				else if (ba.RefID == mitoRef)
				{
					if (r2 < mitoRate)
						output = true;
				}
				else
				{
					if (r2 < mappedRate)
						output = true;
				} 

				if (output)
				{

					ba.BuildCharData();

					if (ba.IsFirstMate())
					{
						if ((++count % 10000) == 0)
							cerr << "Record " << count << endl;
						read.output(f1out);
					}
					else
						read.output(f2out);
				}
			}
		}

		OK = reader.GetNextAlignment(ba);
		if (!OK)
		{
			cerr << "Rewinding" << endl;
			reader.Rewind();
			OK = reader.GetNextAlignment(ba);
		}
	}

	fclose(f1out);
	fclose(f2out);

	return EXIT_SUCCESS;
}