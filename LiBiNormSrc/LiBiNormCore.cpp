#ifdef _WIN32
#define _CRTDBG_MAP_ALLOC
#include <crtdbg.h>
#endif

#include "LiBiNormCore.h"

using namespace std;

void LiBiNormCore::helpCommon()
{
	printf("  -n M, --normModel=M  Specifies that model M should be used rather than the default\n");
	printf("                        Model BD. M options: A or SMART,B or polyA,C,D,\n");
	printf("                        E or random,BD\n");
	printf("  -N FILENAME, --normFilename=FILENAME\n");
	printf("                        Normalise data trying all 6 models and output summary\n");
	printf("                        info to files with root FILENAME.  Best model selected\n");
	printf("                        unless overridden by -n\n");
	printf("  -p N, --threads=N     Number of threads for normalisation parameter\n");
	printf("                        determination (", DEF_THREADS, ")\n");
	printf("  -d N, --reads=N       Maximum number of reads using for normalisation\n");
	printf("                        parameter determination (", DEF_MAX_READS_FOR_PARAM_ESTIMATION, ")\n");
	printf("  -e N, --geneLength=N  Maximum length of transcripts used for normalisation\n");
	printf("                        parameter determination (", DEF_LENGTH_OF_GENE_FOR_PARAM_ESTIMATION, ")\n");
#ifdef USE_NELDER_MEAD_FOR_INITIAL_VALUES
	printf("  -o, --omit            Omit Nelder Mead parameter discovery. Use random initial values for MCMC\n");
#endif
#ifdef INITIAL_VALUES
	printf("  -i <filename>, --initial=FILENAME\n");
	printf("                        Set initial values for parameter discoverey from file\n");
#endif
	printf("  -q, --quiet           suppress progress report\n");
	printf("  -x, --debug           output debug messages\n");
#ifdef PAUSE_AT_END_OPTION
	printf("  -u	                pause at end rather than simply exiting\n");
#endif
	printf("  -c FILENAME, --counts=FILENAME\n");
	printf("                        Name of output file. default: writes to stdout\n");

}

bool LiBiNormCore::commandParseCommon(int & ni, int argc,char **argv)
{
		bool opt2 = false;
		if ((strcmp(argv[ni], "-N") == 0) || (opt2 = (strncmp(argv[ni], "--normFilename=", 15) == 0)))
		{
			normaliseResultsFilename = opt2 ? argv[++ni] + 15 : argv[++ni];
			return true;
		}
		if ((strcmp(argv[ni], "-n") == 0) || (opt2 = (strncmp(argv[ni], "--normModel=", 12) == 0)))
		{
			theModel = modelFromString(opt2 ? argv[++ni] + 12 : argv[++ni]);
			return true;
		}
		if ((strcmp(argv[ni], "-p") == 0) || (opt2 = (strncmp(argv[ni], "--threads=", 10) == 0)))
		{
			Nthreads = atoi(opt2 ? argv[ni] + 10 : argv[++ni]);
			if (Nthreads < 1)
				exitFail("At least 1 thread must be specified");
			return true;
		}
		if ((strcmp(argv[ni], "-d") == 0) || (opt2 = (strncmp(argv[ni], "--reads=", 8) == 0)))
		{
			maxReads = atoi(opt2 ? argv[ni] + 8 : argv[++ni]);
			if (maxReads < 1000)
				exitFail("At least 1000 reads must be specified");
			return true;
		}
		if ((strcmp(argv[ni], "-e") == 0) || (opt2 = (strncmp(argv[ni], "--geneLength=", 13) == 0)))
		{
			maxGeneLength = atoi(opt2 ? argv[ni] + 13 : argv[++ni]);
			return true;
		}
		if ((strcmp(argv[ni], "-c") == 0) || (opt2 = (strncmp(argv[ni], "--counts=", 9) == 0)))
		{
			countsFilename = opt2 ? argv[++ni] + 9 : argv[++ni];
			return true;
		}
		if ((strcmp(argv[ni], "-x") == 0) || (opt2 = (strncmp(argv[ni], "--debug", 7) == 0)))
		{
			debugPrint = true;
			return true;
		}
#ifdef PAUSE_AT_END_OPTION
		if (strcmp(argv[ni], "-u") == 0)
		{
			pauseAtEnd = true;
			return true;
		}
#endif

		if ((strcmp(argv[ni], "-q") == 0) || (opt2 = (strncmp(argv[ni], "--quiet", 7) == 0)))
		{
			verbose = false;
			return true;
		}
#ifdef USE_NELDER_MEAD_FOR_INITIAL_VALUES
		if ((strcmp(argv[ni], "-o") == 0) || (opt2 = (strncmp(argv[ni], "--omit", 6) == 0)))
		{
			nelderMead = false;
			return true;
		}
#endif
#ifdef INITIAL_VALUES
		if ((strcmp(argv[ni], "-i") == 0) || (opt2 = (strncmp(argv[ni], "--intial=", 9) == 0)))
		{
			parameterFilename = (opt2 ? argv[ni] + 8 : argv[++ni]);
			return true;
		}
#endif
		return false;
}

void LiBiNormCore::SetInitialParamsFromFile(const string & filename)
{
	parseTsvFile paramFile;
	if (!paramFile.open(filename))
		exitFail("Unable to read parameters from ", filename);
	paramFile.read(initialValues);

	//	Get rid of spurious values (possibly as a result of trailing tabs in the text file)
	for (modelType m : allModels())
		initialValues[m].resize(headers[m].size());
}

