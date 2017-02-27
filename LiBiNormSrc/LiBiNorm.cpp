#ifdef _WIN32
#define _CRTDBG_MAP_ALLOC
#include <crtdbg.h>
#endif

#include <mutex>
#include <thread>
#include "LiBiOptimiser.h"
#include "LiBiNorm.h"
#include "LiBiCount.h"
#include "LiBiDedup.h"
#include "LiBiConv.h"
#include "LiBiTools.h"
#include "LiBiVariation.h"
#include "MakeFastq.h"

using namespace std;

int main(int argc, char **argv)
{
#ifdef _WIN32
	_DBG( _CrtSetDbgFlag ( _CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF ));
#endif

	if ((argc == 1) || ((argc == 2) && ((strcmp(argv[1], "-h") == 0) || (strcmp(argv[1], "--help") == 0))))
	{
		printf("Usage: LiBiNorm <command> [options]\n");	
		printf("Commands:\n");	
		printf("     count            htseq-count replacement with optional bias correction\n");
		printf("     model            further bias correction analysis\n");
		printf("     conv	          renames chromosomes in a .gff3 file to match those in a bam file\n");
#ifdef INITIAL_VALUES
		printf("     variation        shows variation in Log Liklyhood with parameter\n");
#endif
#ifdef LIBITOOLS
		printf("     land             compares two different landscape files\n");
		printf("     land2            compares a landscape file and a gene list\n");
#endif
#ifdef DEDUP_MODE
		printf("     dedup            removes duplicates\n");
#endif
#ifdef MAKE_FASTQ_MODE
		printf("     makefastq        makes a fastq file from the bam file\n");
#endif
		printf("     -v, --version    version\n");
		printf("     -h, --help       this help\n");
	}
	else if (argc > 1)
	{
		string command(argv[1]);
		if (command == "count")
		{
			LiBiCount libiC;
			return libiC.main(argc-1,argv+1);
		}
		if (command == "model")
		{
			LiBiNorm norm;
			return norm.main(argc-1,argv+1);
		}
		if (command == "conv")
		{
			LiBiConv conv;
			return conv.main(argc - 1, argv + 1);
		}
#ifdef LIBITOOLS
		if (command == "variation")
		{
			LiBiVariation variation;
			return variation.main(argc - 1, argv + 1);
		}
		if (command == "land")
		{
			LiBiTools tools;
			return tools.landMain(argc - 1, argv + 1);
		}
		if (command == "land2")
		{
			LiBiTools tools;
			return tools.landMain2(argc - 1, argv + 1);
		}
		if (command == "genes")
		{
			LiBiTools tools;
			return tools.geneMain(argc - 1, argv + 1);
		}
#endif
#ifdef DEDUP_MODE
		if (command == "dedup")
		{
			LiBiDedup libiD;
			return libiD.main(argc - 1, argv + 1);
		}
#endif
#ifdef MAKE_FASTQ_MODE
		if (command == "makefastq")
		{
			MakeFastq makeFastq;
			return makeFastq.main(argc - 1, argv + 1);
		}
#endif
		if ((command == "--version") || (command == "-v"))
		{
			cout << "LiBiNorm version " << LIBINORM_VERSION << endl;
			return EXIT_SUCCESS;
		}
		exitFail("Invalid commmand:",command);
	}
	return EXIT_SUCCESS;
}


//	Multiple instances of this are called, each works through the requested mcmc runs as 
//	listed in threadLoopCounts which has a pair of integers associated with each 
void LiBiNorm::mcmcThread(optionsType options)
{
	static mutex mtx1;

	static vector<mutex> initValuesMutexes(allModels().size() + 1);
	if (nelderMead)
	{
		modelType m;
		bool modelsLeftToDo = true;
		while (modelsLeftToDo)
		{
			{
				lock_guard<mutex> lock(mtx1);

				while ((modelsLeftToDo = (nelderMeadCounter < allModels().size())))
				{
					m = allModels()[nelderMeadCounter++];
					if (threadLoopCounts[m].requested > 0)
						break;
				}
				if (!modelsLeftToDo)
					break;
				initValuesMutexes[m].lock();
				progMessage("Starting intial values ", m);
			}

			LiBiOptimiser  optimiser(geneData,m,bestResults[m].nelderMeadIterations);
			setSSfun(options, m);
			initialValues[m] = optimiser.getParams(m, options,initialValues[m]);

			initValuesMutexes[m].unlock();

			lock_guard<mutex> lock(mtx1);
			progMessage("Finishing intial values ", m);
		}
	}

	map<modelType, loop_counts >::iterator model_iterator = threadLoopCounts.begin();
	modelType currentModel;
	mcmcRunId mcmcRun;
	while (true)
	{
		//  Look for an iteration of a loop that is yet to be done.
		{
			lock_guard<mutex> lock(mtx1);
			//	Check to see if we have done all the runs associated with this model
			//	if so then move on to the next
			while (model_iterator->second.counter > model_iterator->second.requested)
			{
				if (++model_iterator == threadLoopCounts.end())
					//	All models have been done
					return;
			}
			mcmcRun = model_iterator->second.counter++;

			currentModel = model_iterator->first;
		}
		if (nelderMead)
		{
			//	Stop here if necessary to wait for the initial values to be completed.  
			//	Only likely to happen if there are more than 6 threads
			lock_guard<mutex> lock(initValuesMutexes[currentModel]);
		}
		{
			lock_guard<mutex> lock(mtx1);
			progMessage("Starting ", currentModel, ", iteration:", mcmcRun);
		}
		setSSfun(options, currentModel);

		paramSet params = GetModelParams(currentModel, &initialValues[currentModel], options.jumpSize);
		options.qcov = dataVec(params.size(), options.jumpSize);

		mcmc mcmcEngine;
		mcmcEngine.mcmcrun(geneData, params, options);

		//  Make sure only one thread at a time is outputting results
		lock_guard<mutex> lock(mtx1);

		progMessage("Finishing ", currentModel, ", iteration:", mcmcRun);

#ifdef STORE_ENDPOINTS
		Chain[currentModel].emplace(loop, mcmcEngine.chain().back());
		SSChain[currentModel].emplace(loop, mcmcEngine.sschain().back());
#endif

		//	Always store full set of results as these are needed to calculate the optimal parameters
		//	emplace/move them for efficiency
		fullResultChain[currentModel].emplace(mcmcRun, move(mcmcEngine._chain));
		fullResultSSChain[currentModel].emplace(mcmcRun, move(mcmcEngine._sschain));
	}
}



int LiBiNorm::main(int argc, char **argv)
{
	int Ngenes = -1;

	initClock();
	if(argc < 1)
	{
		printf("Error: parameter wrong!\n");
		return EXIT_FAILURE;
	}
	else if ((argc == 1) || ((argc == 2) && ((strcmp(argv[1], "-h") == 0) || (strcmp(argv[1], "--help") == 0))))
	{
		printf("/* ----------------------------- */\n");
		printf("     LiBiNorm model:    RNA-seq library bias normalisation using a landscape file   \n\n");
		printf("Usage: LiBiNorm model [options] landscapeFilename \n\n");
		printf("Options:\n");
		printf("  -h, --help            Show this help message and exit\n");
		helpCommon();
		printf("  -g N, --genes=N       Limit parameter discovery to using data from first N genes\n");
		printf("  -r N, --runs=N        Number of MCMC runs (", NUMBER_OF_MCMC_RUNS,")\n");
#ifdef USE_NELDER_MEAD_FOR_INITIAL_VALUES
		printf("  -s N, --mcmc=N        Length of each MCMC run (", NELDER_MCMC_ITERATIONS,")\n");
#else
		printf("  -s N, --mcmc=N        Length of each MCMC run (", MCMC_ITERATIONS, ")\n");
#endif
#ifdef SELECT_READ_SEED
		printf("  -y N, --seed=N        Set seed used for selecting subset of reads\n");
#endif
		printf("  -f, --full            Output complete set of MCMC run data\n");
		return EXIT_SUCCESS;
	}
	int ni = 1;
	while(ni < argc -1)
	{
		bool opt2 = false;
		if (commandParseCommon(ni, argc,argv))
		{
		}
		else if ((strcmp(argv[ni], "-g") == 0) || (opt2 = (strncmp(argv[ni], "--genes=", 8) == 0)))
		{
			Ngenes = atoi(opt2 ? argv[ni] + 8 : argv[++ni]);
			if (Ngenes < 10)
				exitFail("At least 10 genes must be specified");
		}
		else if ((strcmp(argv[ni], "-r") == 0) || (opt2 = (strncmp(argv[ni], "--runs=", 7) == 0)))
		{
			Nruns = atoi(opt2 ? argv[ni] + 7 : argv[++ni]);
			if ((Nruns < 1) || (Nruns > 200))
				exitFail("-r values must lie between 1 and 200");
		}
		else if ((strcmp(argv[ni], "-s") == 0) || (opt2 = (strncmp(argv[ni], "--mcmc=", 7) == 0)))
		{
			Nsimu = atoi(opt2 ? argv[ni] + 7 : argv[++ni]);
			if ((Nsimu < 1) || (Nsimu > 10000))
				exitFail("-s values must lie between 1 and 10000");
		}
#ifdef SELECT_READ_SEED
		else if ((strcmp(argv[ni], "-y") == 0) || (opt2 = (strncmp(argv[ni], "--seed=", 7) == 0)))
		{
			int seed = atoi(opt2 ? argv[ni] + 7 : argv[++ni]);
			intRand.reseed(seed);
		}
#endif
		else if ((strcmp(argv[ni], "-f") == 0) || (opt2 = (strncmp(argv[ni], "--full", 6) == 0)))
			outputFull = true;
		else
		{
			exitFail("Invalid parameter: ",argv[ni]);
		}
		ni++;
	}

	landscapeFilename = argv[argc - 1];

	if (!normaliseResultsFilename)
		normaliseResultsFilename = landscapeFilename;

	//	If we specifiy the model then run the other models just once 
	if (theModel == noModel)
		NrunsOtherModels = Nruns;
	else
		NrunsOtherModels = (Nruns ==1)?0:1;

	geneCounts.loadData(landscapeFilename, Ngenes);

	//	Load up the initial values
	if (parameterFilename)
	{
		SetInitialParamsFromFile(parameterFilename);
	}


	coreParameterEstimation();

	//	If we have explicitly specified the model then use it instead
	if (theModel == noModel)
	{
		bestModel = getBestModel();
		progMessage("Best model is ", bestModel);
		theModel = bestModel;
	}
	else
	{
		progMessage("Model selected by command line is ", theModel);
	}

	getBias(theModel, bestResults[theModel].params[logValue], geneCounts.lengths[0], geneCounts.bias);

	string filename = normaliseResultsFilename.replaceSuffix("_expression.txt");
	if (!geneCounts.outputGeneCounts(filename, outputFull ? 3 : 2, conv(theModel)))
		exitFail("Unable to output counts to :", filename);

	printResults();
	printBias();

	//********************************************************************************************
	//	This prints out all of the data for the full set of mcmc runs for each model

#ifdef PRINT_MCMC_RUN_DATA
	printAllMcmcRunData();
	printConsolidatedMcmcRunData();
#endif

	//	The basic count data in htseq-count format
	if (!geneCounts.outputGeneCounts(countsFilename, 1, conv(theModel)))
		exitFail("Unable to output counts to :", countsFilename);

	progMessage("Data modelled");
	elapsedTime();

	if (pauseAtEnd)
	{
		string test;
		cin >> test;
	}

	return EXIT_SUCCESS;
}

bool LiBiNorm::coreParameterEstimation()
{
	geneCounts.remove_invalid_values();
	geneCounts.transferTo(geneData, MAX_READS_GENE,maxReads,maxGeneLength);

	elapsedTime("Data loaded");

	optionsType options;

	options.jumpSize = MCMC_JUMP_SIZE;
	options.nsimu = Nsimu;
	options.Nruns = Nruns;

	options.sigma2 = MCMC_SIGMA;

	//	Set the number of iterations required of each of the models.
	for (modelType m : allModels())
	{
		threadLoopCounts[m].counter = 1;
		if (m == theModel)
			threadLoopCounts[m].requested = options.Nruns;
		else
			threadLoopCounts[m].requested = NrunsOtherModels;
	}

#ifdef FIXED_RESULTS
	theModel = M_FIXED_RESULTS;
	bestResults[theModel].params = dataVec{ FIXED_RESULTS };
#else

	nelderMeadCounter = 0;
	//	And then set the threads running
	if (Nthreads == 1)
	{
		//	Dont make additional threads if single threaded.  Makes debugging easier
		mcmcThread(options);
	}
	else
	{
		vector<thread> threads;
		for (size_t i = 0; i < Nthreads; i++)
			threads.emplace_back(&LiBiNorm::mcmcThread, this, options);

		for (auto & i : threads)
			i.join();
	}

#ifdef XXXX
	Unfinished code for
		vector<int> geneL{ 500,1000,2000,4000,8000 };
	dataVec TotalReads;

	for (size_t i = 0; i < geneL.size(); i++)
	{
		for (size_t j = 0; j < transData.size(); j++)
		{
			if (abs(transData[j].length - geneL[i]) < (0.1 * geneL[i]))
				TotalReads.append(transData[j].counts[0] / transData[j].length);
		}
	}
#endif

	//********************************************************************************************
	//	Find the optimal parameter values, which are associated with results found in the last
	//	iterations of all of the MCMC runs.
	for (modelType m : allModels())
	{
		multimap <double, dataVec *> & orderedResults = allOrderedResults[m];
		if (fullResultSSChain[m].size())
		{
			bestResult & br = bestResults[m];
			size_t Nparams = headers[m].size();
			for (size_t i = 0; i < 2; i++)
				br.params[i].resize(Nparams);
			for (size_t i = 0; i < 4; i++)
				br.param_dev[i].resize(Nparams);

			//	Use the second half of teh chain to find parameters and their variation
			for (mcmcRunId i = 1; i <= fullResultSSChain[m].size(); i++)
			{
				for (size_t j = fullResultSSChain[m][i].size() / 2; j < fullResultSSChain[m][i].size(); j++)
				{
					orderedResults.emplace(fullResultSSChain[m][i][j], &fullResultChain[m][i][j]);
				}
			}

#ifdef USE_PARAMS_FROM_LOWEST_LL
			//	Find the best LL.
			br.LLresult = orderedResults.begin()->first;
			//	Get the parameters from the most likly sample
			for (size_t p = 0; p < Nparams; p++)
			{
				VEC_DATA_TYPE v = (*orderedResults.begin()->second)[p];
				br.params[logValue][p] = v;
				if (p < 4)
					br.params[absValue][p] = pow(10, v);
				else
					br.params[absValue][p] = v;
			}
#else
			{
				//	Put the Log liklyhoods in order and find the median
				auto i = orderedResults.begin();
				size_t n = 0;
				for (; n < orderedResults.size() / 2; i++, n++) {};
				br.LLresult = i->first;
			}
			{
				//	Now put each of the params in order
				auto j = orderedResults.begin();
				vector<orderedVec> orderedParams(Nparams);
				for (; j != orderedResults.end(); j++)
				{
					for (size_t p = 0; p < Nparams; p++)
						orderedParams[p].add((*j->second)[p]);
				}
				//	and then find the median
				for (size_t p = 0; p < Nparams; p++)
				{
					VEC_DATA_TYPE v = orderedParams[p].median();
					br.params[logValue][p] = v;
					if (p < 4)
						br.params[absValue][p] = pow(10, v);
					else
						br.params[absValue][p] = v;
				}
			}
#endif
			//	Find the absolute distance from the selected 'result' LL in order
			orderedVec distanceFromOptimalLL;
			for (auto i = orderedResults.begin(); i != orderedResults.end(); i++)
				distanceFromOptimalLL.add(abs(i->first - br.LLresult));

			//	and then find the median = MAD
			br.LL_dev = distanceFromOptimalLL.median();

			//	Now go through the parameters creating ordered lists of the absolute distance from the selected param
			//	and also the positive and negative distances so that we can do single sided deviation measures
			typedef vector<orderedVec> diffList;
			vector<diffList> diffs(4, diffList(Nparams));

			size_t n = 0;
			auto j = orderedResults.begin();
			for (; n < orderedResults.size(); j++, n++)
			{
				for (size_t p = 0; p < Nparams; p++)
				{
					VEC_DATA_TYPE v = (*j->second)[p];
					if (v > br.params[logValue][p])
					{
						diffs[minLog][p].add(v - br.params[logValue][p]);
						diffs[logValue][p].add(v - br.params[logValue][p]);
					}
					else
					{
						diffs[maxLog][p].add(br.params[logValue][p] - v);
						diffs[logValue][p].add(br.params[logValue][p] - v);
					}
					if (p < 4)
						diffs[absValue][p].add(abs(br.params[absValue][p] - pow(10, v)));
					else
						diffs[absValue][p].add(abs(br.params[absValue][p] - v));
				}
			}
			//	And then find the medians
			for (size_t p = 0; p < Nparams; p++)
			{
				for (size_t t = 0; t < diffs.size(); t++)
					br.param_dev[t][p] = diffs[t][p].median();
			}
			progMessage(m);
		}
	}
/*	******************************************************************************************
	//	And then find the standard deviation
	for (modelType m :allModels())
	{
		bestResult & br = bestResults[m];
		VEC_DATA_TYPE LL_dev = 0;
		size_t size = br.params.size();
		vector<dataVec> param_dev(2, dataVec(size));
		vector<dataVec> p_N(2, dataVec(size));
		int N = 0;

		for (mcmcRunId i = 1; i <= fullResultSSChain[m].size(); i++)
		{
			for (size_t j = fullResultSSChain[m][i].size() - 1; j > fullResultSSChain[m][i].size() - 100 ; j--)
			{
				N++;
				VEC_DATA_TYPE diff = fullResultSSChain[m][i][j] - br.minLL;
				LL_dev += (diff * diff);
				for (size_t k = 0; k < size; k++)
				{
					VEC_DATA_TYPE diff = fullResultChain[m][i][j][k] - br.params[k];
					if (diff >= 0)
					{
						param_dev[0][k] += (diff * diff);
						p_N[0][k]++;
					}
					else
					{
						param_dev[1][k] += (diff * diff);
						p_N[1][k]++;
					}
				}
			}
		}
		br.minLL_dev = sqrt(LL_dev / N);
		br.param_dev[0] = sqrt(param_dev[0] / p_N[0]);
		br.param_dev[1] = sqrt(param_dev[1] / p_N[1]);
	}*/
#endif  // ifdef FIXED_RESULTS

	return true;
}

//	Find which model performed best based on the Log Liklihood
modelType LiBiNorm::getBestModel()
{
	modelType bestModel = noModel;
	double bestLL = MAX_DOUBLE;
	for (modelType m : allModels())
	{
		if (bestResults[m].LLresult < bestLL)
		{
			bestModel = m;
			bestLL = bestResults[m].LLresult;
		}
	}
	return bestModel;
}

//	The top level summary of the results, showing best LL and associated paremeters for each model
void LiBiNorm::printResults()
{
	string filename(normaliseResultsFilename.replaceSuffix("_results.txt"));
	TsvFile mcmcResult;
	//	Now output a table with the end points of each of the chains.
	if (!mcmcResult.open(filename))
		exitFail("Unable to open output File ", filename);

	//	First headers up to and including the maximum model that is run.   Always leave space
	//	for the intermediate models so the layout of the results is consistent
	mcmcResult.printStart("");
	for (modelType m : allModels())
	{
		string count;
		if (bestResults[m].nelderMeadIterations)
			count = _s("NM iterations:", bestResults[m].nelderMeadIterations);
		mcmcResult.printMiddle(m, count);
		mcmcResult.printGaps(headers[m].size());
	}
	mcmcResult.printEnd();

	mcmcResult.printStart("");
	for (modelType m : allModels())
		mcmcResult.printMiddle(headers[m], "", "");
	mcmcResult.printEnd();

	//	A row for the optimal parameters that were found for each model
	for (size_t i = 0; i < 2; i++)
	{
		mcmcResult.printStart((i==0)?"Log Opt":"Abs Opt");
		for (modelType m : allModels())
		{
			if (bestResults[m].params[i].size())
				mcmcResult.printMiddle(bestResults[m].params[i], bestResults[m].LLresult, "");
			else
				mcmcResult.printGaps(headers[m].size() + 2);
		}
		mcmcResult.printEnd();
	}

	//	And the initial Values
	if ((nelderMead) && (initialValues[ModelA].size()))
	{
		mcmcResult.printStart("Initial");
		for (modelType m : allModels())
			mcmcResult.printMiddle(initialValues[m], "");
	}
	mcmcResult.printEnd();

	//	A row for the deviations for each model
	for (size_t i = 0; i < 4; i++)
	{
		switch (i)
		{
		case 0:	mcmcResult.printStart("Spread Log"); break;
		case 1:	mcmcResult.printStart("Spread Abs"); break;
		case 2:	mcmcResult.printStart("pos Spread Log"); break;
		case 3:	mcmcResult.printStart("neg Spread Log"); break;
		}
		for (modelType m : allModels())
		{
			if (bestResults[m].param_dev[i].size())
				mcmcResult.printMiddle(bestResults[m].param_dev[i], bestResults[m].LL_dev, "");
			else
				mcmcResult.printGaps(headers[m].size() + 2);
		}
		mcmcResult.printEnd();
	}
	
	//	Print out the enpoints of all of the mcmcruns
	for (mcmcRunId i = 1; i <= Nruns; i++)
	{
		mcmcResult.printStart(_s("Chain end ", i));
		for (modelType m : allModels())
		{
			//	Was there a jth run of this model?  If so then print the end points.  *...rbegin() gets the last entry in the list.
			if (fullResultSSChain[m].size() && (i <= fullResultSSChain[m].rbegin()->first))
				mcmcResult.printMiddle(*fullResultChain[m][i].rbegin(), *fullResultSSChain[m][i].rbegin(), "");
			else
				mcmcResult.printGaps(headers[m].size() + 2);
		}
		mcmcResult.printEnd();
	}
	mcmcResult.close();
}



// #define BEST_RESULT_LOCATION  Use this to print out locations where best loglilihood is obtained in mcmc run
#ifdef BEST_RESULT_LOCATION
#define _BRL(A,B) A,B,
#else
#define _BRL(A,B)
#endif

void LiBiNorm::printBias()
{
	TsvFile mcmcResult;
	string filename(normaliseResultsFilename.replaceSuffix("_norm.txt"));
	if (!mcmcResult.open(filename))
		exitFail("Unable to open output File ", filename);

	//	Create a vector containing the list of frequencies that we are going to calculate the
	//	bias figures for
	dataVec lengths;
	lengths.push_back(DEFAULT_NORMALISATION_GENE_LENGTH);
	for (size_t i = 100; i <= MAX_GENE_LENGTH_FOR_NORM_PLOT; i += 100)
		lengths.push_back(i);

	map<size_t, dataVec> biases;

	for (modelType m : allModels())
	{
		if (bestResults[m])
			getBias(m, bestResults[m].params[logValue], lengths, biases[m]);
	}

	for (modelType m : allModels())
	{
		mcmcResult.print("","","Log likelihood",_BRL("mcmc run","mcmc iteration")headers[m]);
		mcmcResult.print(m,"Parameters",bestResults[m].LLresult,_BRL(bestResults[m].run, bestResults[m].pos) bestResults[m].params);
		if (biases[m].size())
		{
			mcmcResult.print("","Length", lengths);
			mcmcResult.print("","Bias",biases[m]);
		}
		else
		{
			mcmcResult.print("");
			mcmcResult.print("");
		}
		mcmcResult.print();
	}
	mcmcResult.close();
}

void LiBiNorm::printAllMcmcRunData()
{
	TsvFile mcmcResult;
	for (modelType modl : allModels())
	{
		string filename = normaliseResultsFilename.replaceSuffix("_", conv(modl,true), ".txt");
		if (!mcmcResult.open(filename))
			exitFail("Unable to open output file ", filename);

		//			mcmcResult.printMiddle(headers[modl], "chain", "");

		//	This ensures that at least one header is output, which ensures that there is something in the file
		//	even if no data were produced for this model
		for (size_t i = 0; i < fullResultChain[modl].size(); i++)
			mcmcResult.printMiddle(headers[modl], "LL", "");
		mcmcResult.printEnd();

		//	For each of the mcmc runs print out the results.  Each run is a column
		for (size_t i = 0; i < fullResultChain[modl][1].size(); i++)
		{
			for (mcmcRunId j = 1; j <= fullResultChain[modl].rbegin()->first; j++)
			{
				// fullResultChain[modl][j][i] is a vector of N values which are printed out as N tab separated values
				//	using the TsvFile support for printing vectors.
				mcmcResult.printMiddle(fullResultChain[modl][j][i], fullResultSSChain[modl][j][i], "");
			}
			mcmcResult.printEnd();
		}
		mcmcResult.close();
	}
}

void LiBiNorm::printConsolidatedMcmcRunData()
{
	TsvFile mcmcResult;
	string filename = normaliseResultsFilename.replaceSuffix("_model_cons.txt");
	if (!mcmcResult.open(filename))
		exitFail("Unable to open output File ", filename);

	map<size_t, multimap <double, dataVec *>::iterator > iterators;
	mcmcResult.printStart("");
	for (modelType m : allModels())
	{
		mcmcResult.printMiddle(headers[m], "chain", "");
		iterators[m] = allOrderedResults[m].begin();
	}
	mcmcResult.printEnd("");
	bool found = true;
	for (size_t i = 0; (i < 1000000) && found; i++)
	{
		found = false;
		mcmcResult.printStart(i);
		for (modelType m  : allModels())
		{
			if (iterators[m] != allOrderedResults[m].end())
			{
				found = true;
				mcmcResult.printMiddle(*(iterators[m]->second), iterators[m]->first, "");
				iterators[m]++;
			}
			else
				mcmcResult.printGaps(headers[m].size() + 2);
		}
		mcmcResult.printEnd();
	}
	mcmcResult.close();
}
