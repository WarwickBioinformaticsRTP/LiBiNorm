#include "LiBiCount.h"
#include <stdlib.h>
#include <crtdbg.h>
#include "libCommon.h"
#include "containerEx.h"
#include "gtfFile.h"
#include "api/BamReader.h"

using namespace std;
using namespace BamTools;

int LiBiCount::main(int argc, char **argv)
{
	string bamFileName,gtfFileName;

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

	gtfFile genomeDef;

	genomeDef.open(gtfFileName);





	return EXIT_SUCCESS;
}
