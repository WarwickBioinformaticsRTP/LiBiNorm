#include "LiBiVariation.h"

using namespace std;

int LiBiVariation::main(int argc, char **argv)
{
	initClock();
	if (argc < 1)
	{
		printf("Error: parameter wrong!\n");
		return EXIT_FAILURE;
	}
	else if ((argc == 1) || ((argc == 2) && ((strcmp(argv[1], "-h") == 0) || (strcmp(argv[1], "--help") == 0))))
	{
		printf("/* ----------------------------- */\n");
		printf("     LiBiNorm variation:  Test of effect of parameter variation\n\n");
		printf("Options:\n");
		printf("  -h, --help            show this help message and exit\n");
		helpCommon();
		printf("  -i <filename>, --intial=<filename>  initial valuse\n");
		return EXIT_SUCCESS;
	}
	int ni = 1;
	while (ni < argc)
	{
		bool opt2 = false;
		if (commandParseCommon(ni, argc, argv))
		{
		}
		else if ((strcmp(argv[ni], "-i") == 0) || (opt2 = (strncmp(argv[ni], "--intial=", 9) == 0)))
		{
			parameterFilename = (opt2 ? argv[ni] + 8 : argv[++ni]);
		}
		else
		{
			exitFail("Invalid parameter: ", argv[ni]);
		}
		ni++;
	}

	if (!landscapeFilename)
		exitFail("Landscape file must be specified");
	geneCounts.loadData(landscapeFilename, -1);

	geneCounts.remove_invalid_values();
	geneCounts.transferTo(geneData, MAX_READS_GENE, maxReads);


	dataVecFile paramFile;
	vector<dataVec> paramSet = paramFile.read(parameterFilename);

	string filename(normaliseResultsFilename.replaceSuffix("_variation.txt"));
	TsvFile mcmcResult;
	//	Now output a table with the end points of each of the chains.
	if (!mcmcResult.open(filename))
		exitFail("Unable to open output File ", filename);

	//	First headers up to and including the maximum model that is run.   Always leave space
	//	for the intermediate models so the layout of the results is consistent
	mcmcResult.printStart("");
	for (modelType m : allModels())
	{
		mcmcResult.printMiddle(m);
		mcmcResult.printGaps(headers[m].size() + 1);
	}
	mcmcResult.printEnd();

	mcmcResult.printStart("");
	for (modelType m : allModels())
		mcmcResult.printMiddle(headers[m], "chain", "");
	mcmcResult.printEnd();


	mcmcResult.printStart("Log Opt");
	for (modelType m : allModels())
	{
		mcmcResult.printMiddle(paramSet[m], "" , "");
	}
	mcmcResult.printEnd();

	optionsType options;

	for (size_t p = 0; p < 2; p++)
	{
		for (double i = -1; i < 4; i += 0.2)
		{
			mcmcResult.printStart("");
			for (modelType m : allModels())
			{
				dataVec params = paramSet[m];
				params[p] += i;
				setSSfun(options, m);
				double ss1 = options.ssfun(params, geneData);
				mcmcResult.printMiddle(params, ss1, "");
			}
			mcmcResult.printEnd();
		}
		mcmcResult.print("");
	}
	return EXIT_SUCCESS;
}

