#include "parser.h"
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
		return EXIT_SUCCESS;
	}
	int ni = 1;
	while (ni < argc)
	{
		if (!commandParseCommon(ni, argc, argv))
		{
			exitFail("Invalid parameter: ", argv[ni]);
		}
		ni++;
	}

	//	Load up the initial values
	parseTsvFile paramFile;
	if (!paramFile.open(parameterFilename))
		exitFail("Unable to read parameters from ", parameterFilename);
	paramFile.read(initialValues);

	//	Load up the landscape file
	if (!landscapeFilename)
		exitFail("Landscape file must be specified");
	geneCounts.loadData(landscapeFilename, -1);
	geneCounts.remove_invalid_values();
	geneCounts.transferTo(geneData, MAX_READS_GENE, maxReads);

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
		mcmcResult.printMiddle(initialValues[m], "" , "");
	mcmcResult.printEnd();

	optionsType options;

	for (size_t p = 0; p < 2; p++)
	{
		for (double i = -1; i < 4; i += 0.2)
		{
			mcmcResult.printStart("");
			for (modelType m : allModels())
			{
				dataVec values = initialValues[m];
				values[p] += i;
				setSSfun(options, m);
				paramSet params = GetModelParams(m);
				double ss1 = options.ssfun(values, geneData,params);
				mcmcResult.printMiddle(values, ss1, "");
			}
			mcmcResult.printEnd();
		}
		mcmcResult.print("");
	}
	return EXIT_SUCCESS;
}

