#ifdef _WIN32
#define _CRTDBG_MAP_ALLOC
#include <crtdbg.h>
#endif

#include <stdlib.h>
#include <vector>
#include <mutex>
#include <chrono>
#include <thread>
#include "libCommon.h"
#include "containerEx.h"
#include "mcmc.h"
#include "LiBiNorm.h"
#include "LiBiCount.h"
#include "LiBiDedup.h"
#include "LiBiConv.h"
#include "MakeFastq.h"

using namespace std;

//#define TEST_CODE
#ifdef TEST_CODE
extern map<thread::id,map<size_t,vector<vector<VEC_DATA_TYPE> > > > cache;
#endif

#ifdef _DEBUG
//	use this to test the calculations based on a specific result of the parameter derivation
// #define FIXED_RESULTS log10(3.9644),log10(147.39),log10(0.0079),log10(2.0128E-4),0
#define M_FIXED_RESULTS 5

#endif

int main(int argc, char **argv)
{


#ifdef _WIN32
	_DBG( _CrtSetDbgFlag ( _CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF ));
#endif

	if (argc == 1)
	{
		printf("Usage: LiBiNorm <command> [options]\n");	
		printf("Commands:\n");	
		printf("     count            htseq-count replacement\n");
		printf("     dedup            removes duplicates\n");
		printf("     conv	          renames chromosomes in a .gff3 file to match those in a bam file\n");
		printf("     makefastq        makes a fastq file from the bam file\n");
	}
	else if (argc > 1)
	{
		string command(argv[1]);
		if (command == "count")
		{
			LiBiCount libiC;
			return libiC.main(argc-1,argv+1);
		}
		else if (command == "dedup")
		{
			LiBiDedup libiD;
			return libiD.main(argc-1,argv+1);
		}
		else if (command == "makefastq")
		{
			MakeFastq makeFastq;
			return makeFastq.main(argc-1,argv+1);
		}
		else if (command == "model")
		{
			LiBiNorm norm;
			return norm.main(argc-1,argv+1);
		}
		else if (command == "conv")
		{
			LiBiConv conv;
			return conv.main(argc - 1, argv + 1);
		}
		if (command == "--version")
		{
			cout << "LiBiNorm version 1.0.3" << endl;
			return EXIT_SUCCESS;
		}
		else
			exitFail("Invalid commmand:",command);
	}
	return EXIT_SUCCESS;
}


void LiBiNorm::mcmcThread(paramSet params, optionsType options, modelType model)
{
	map<size_t,int>::iterator model_iterator = threadLoopCounts.begin();
	size_t loop;
	while (true)
	{
		//  Look for an iteration of a loop that is yet to be done.
		//	The iteration with the counter set to zero is not used for an mcmc run but is a dummy
		//	to ensure that the header information is set even if we are not running any iterations of this
		//	model
		{
			static mutex mtx; 
			lock_guard<mutex> lock(mtx);
			while(model_iterator->second == -1)
			{
				dataVec::clearCache();
				if (++model_iterator ==threadLoopCounts.end())
				{
					return;
				}
			}
			loop = model_iterator->second--;
			if (loop != 0)
				progMessage("Starting Model:",model_iterator->first," iteration:",loop);
			options.Model = model_iterator->first;
		}

		switch (options.Model)
		{
		case 1:
			model.ssfun = &FLL_ModelA;
			break;
		case 2:
			model.ssfun = &FLL_ModelB;
			break;
		case 3:
			model.ssfun = &FLL_ModelC;
			break;
		case 4:
			model.ssfun = &FLL_ModelD;
			break;
		case 5:
			model.ssfun = &FLL_ModelE;
			break;
		case 6:
			model.ssfun = &FLL_ModelBD;
			break;
		}



// #define _TEST
#ifdef _TEST
		vectorEx<double> p0{ {1.5, 1.6,-3.1, -3.2,0.6}};
#else
		vectorEx<double> p0{ {rand(3) - 1, rand(3), rand(4) - 5, rand(4) - 5, rand(1)} };
#endif

		switch (options.Model)
		{
		case 2: case 4: case 5:
			options.qcov = dataVec(4,options.jumpSize);

			params = { paramType("d", p0[0], -1 , 2)    // average length of fragments
			,paramType("h",  p0[1], 0 , 3)   // the minimum length of fragmenation
			,paramType("t1", p0[2], -5 , -1)   // theta1
			,paramType("t2", p0[3], -5, -1) // theta2
			//				,paramType("sig", p0[4], 0, 3) // sigma
			};

			break;
		case 3:
			options.qcov = dataVec(3,options.jumpSize);

/*			params = { paramType("d", p0[0], -1 , 2)    // average length of fragments
				,paramType("h",  p0[1], 0 , 3)   // the minimum length of fragmenation
				//				,paramType("t1", p0[2], -5 , -1)   // theta1
				,paramType("t2", p0[3], -5, -1) // theta2
				//				,paramType("sig", p0[4], 0, 3) // sigma
			};*/
			params = { {"d", p0[0], -1 , 2}    // average length of fragments
			, {"h",  p0[1], 0 , 3 }   // the minimum length of fragmenation
			//				,paramType("t1", p0[2], -5 , -1)   // theta1
			,{"t2", p0[3], -5, -1} // theta2
			//				,paramType("sig", p0[4], 0, 3) // sigma
			};
			break;
		case 1:
			options.qcov = dataVec(2,options.jumpSize);

			params = { paramType("d", p0[0], -1 , 2)    // average length of fragments
				,paramType("h",  p0[1], 0 , 3)   // the minimum length of fragmenation
				//				,paramType("t1", p0[2], -5 , -1)   // theta1
				//				,paramType("t2", p0[3], -5, -1) // theta2
			};
			break;
		case 6:
			options.qcov = dataVec(6,options.jumpSize);

			params = { paramType("d", p0[0], -1 , 2)    // average length of fragments
				,paramType("h",  p0[1], 0 , 3)   // the minimum length of fragmenation
				,paramType("t1", p0[2], -5 , -1)   // theta1
				,paramType("t2", p0[3], -5, -1) // theta2
				,paramType("a", p0[4], 0, 1) // alpha strength of model B
			};
			break;

		};

		if (loop == 0)
		{
			for (size_t i = 0; i < params.size(); i++)
				headers[options.Model].push_back(params[i].name);
		}
		else
		{
			mcmc mcmcEngine;
			mcmcEngine.mcmcrun(model, consData, params, options);

			static mutex mtx;
			lock_guard<mutex> lock(mtx);

			progMessage("Finishing Model:",model_iterator->first," iteration:",loop);

#ifdef STORE_ENDPOINTS
			Chain[options.Model].emplace(loop, mcmcEngine.chain().back());
			SSChain[options.Model].emplace(loop, mcmcEngine.sschain().back());
#endif

			//	Always store full set of results as these are needed to calculate the optimal parameters
			fullResultChain[options.Model].emplace(loop, mcmcEngine.chain());
			fullResultSSChain[options.Model].emplace(loop, mcmcEngine.sschain());

			RejectionRate[options.Model] += mcmcEngine.rejected();
		}
	}
}

void runThread(LiBiNorm * root,	paramSet params, optionsType options, modelType model)
{
	root->mcmcThread(params, options,model);
}


void LiBiNormCore::helpCommon()
{
	printf("  -l FILENAME, --landscape=FILENAME\n");
	printf("                        Name of file for landscape data)\n");
	printf("  -n, --normalise       Normalise rna-seq data using model 6 to correct for\n");
	printf("                        length related bias\n");
	printf("  -N FILENAME, --Normalise=FILENAME\n");
	printf("                        Normalise data trying all 6 models and output summary\n");
	printf("                        info to files with root FILENAME\n");
	printf("  -p N, --threads=N     Number of threads for normalisation parameter\n");
	printf(_s("                        determination (", DEF_THREADS, ")\n"));
	printf("  -d N, --reads=N       Maximum number of reads using for normalisation\n");
	printf(_s("                        parameter determination (", DEF_MAX_READS_FOR_PARAM_ESTIMATION, ")\n"));
	printf("  -q, --quiet           suppress progress report\n");
	printf("  -c FILENAME, --counts=FILENAME\n");
	printf("                        Name of output file. default: writes to stdout)\n");

}

bool LiBiNormCore::commandParseCommon(int & ni, char **argv)
{
		bool opt2 = false;
		if ((strcmp(argv[ni], "-l") == 0) || (opt2 = (strncmp(argv[ni], "--landscape=", 12) == 0)))
		{
			landscapeFilename = opt2 ? argv[++ni] + 12 : argv[++ni];
			return true;
		}
		if ((strcmp(argv[ni], "-N") == 0) || (opt2 = (strncmp(argv[ni], "--Normalise", 11) == 0)))
		{
			normalise = true;
			normaliseResultsFilename = opt2 ? argv[++ni] + 11 : argv[++ni];
			return true;
		}
		if ((strcmp(argv[ni], "-p") == 0) || (opt2 = (strncmp(argv[ni], "--threads=", 10) == 0)))
		{
			Nthreads = atoi(opt2 ? argv[ni] + 10 : argv[++ni]);
			return true;
		}
		if ((strcmp(argv[ni], "-d") == 0) || (opt2 = (strncmp(argv[ni], "--reads=", 8) == 0)))
		{
			maxReads = atoi(opt2 ? argv[ni] + 8 : argv[++ni]);
			return true;
		}
		if ((strcmp(argv[ni], "-c") == 0) || (opt2 = (strncmp(argv[ni], "--counts=", 9) == 0)))
		{
			countsFilename = opt2 ? argv[++ni] + 9 : argv[++ni];
			return true;
		}

		return false;
}

int LiBiNorm::main(int argc, char **argv)
{
	int Ngenes = -1;
	normalise = true;
	NrunsOtherModels = Nruns;

	initClock();
	if(argc < 1)
	{
		printf("Error: parameter wrong!\n");
		return EXIT_FAILURE;
	}
	else if ((argc == 1) || ((argc == 2) && ((strcmp(argv[1], "-h") == 0) || (strcmp(argv[1], "--help") == 0))))
	{
		printf("/* ----------------------------- */\n");
		printf("     LiBiNorm:    RNA-seq library bias normalisation   \n\n");
		printf("Options:\n");
		printf("  -h, --help            show this help message and exit\n");
		helpCommon();
		printf("  -g N, --genes=N       Number of genes to be included in calculations\n");
		printf(_s("  -r N, --runs=N        Number of mcmc runs (", NUMBER_OF_MCMC_RUNS,")\n"));
		printf(_s("  -s N, --mcmc=N        Length of each simulation (", MCMC_ITERATIONS,")\n"));
		printf("  -m N, --model=N       Just run for model N -r times.  All other models run once\n");
		printf("  -f, --full            Output complete set of montecarlo data\n");
		return EXIT_SUCCESS;
	}
	int ni = 1;
	while(ni < argc)
	{
		bool opt2 = false;
		if (commandParseCommon(ni, argv))
		{
		}
		else if ((strcmp(argv[ni], "-g") == 0) || (opt2 = (strncmp(argv[ni], "--genes=", 8) == 0)))
		{
			Ngenes = atoi(opt2 ? argv[ni] + 8 : argv[++ni]);
		}
		else if ((strcmp(argv[ni], "-r") == 0) || (opt2 = (strncmp(argv[ni], "--runs=", 7) == 0)))
		{
			Nruns = atoi(opt2 ? argv[ni] + 7 : argv[++ni]);
			NrunsOtherModels = (theModel) ? 1 : Nruns;
		}
		else if ((strcmp(argv[ni], "-s") == 0) || (opt2 = (strncmp(argv[ni], "--mcmc=", 7) == 0)))
		{
			Nsimu = atoi(opt2 ? argv[ni] + 7 : argv[++ni]);
		}
		else if ((strcmp(argv[ni], "-m") == 0) || (opt2 = (strncmp(argv[ni], "--model=", 8) == 0)))
		{
			theModel = atoi(opt2 ? argv[ni] + 8 : argv[++ni]);
			NrunsOtherModels = 1;
		}
		else if ((strcmp(argv[ni], "-f") == 0) || (opt2 = (strncmp(argv[ni], "--full", 6) == 0)))
			outputFull = true;
		else
		{
			exitFail("Invalid parameter: ",argv[ni]);
		}
		ni++;
	}

	if (!landscapeFilename)
		exitFail("Landscape file must be specified");

	if (!normaliseResultsFilename)
		normaliseResultsFilename = landscapeFilename;

	string lastGene = geneCounts.loadData(landscapeFilename, Ngenes);
	coreParameterEstimation();

	size_t bestModel = 0;
	double bestLL = 1E99;
	for (size_t m = 1; m <= N_MODELS; m++)
	{
		if (bestResults[m].minLL < bestLL)
		{
			bestModel = m;
			bestLL = bestResults[m].minLL;
		}
	}
	progMessage("Best model is model ", bestModel);

	normaliseExpression(bestModel, bestResults[bestModel].params, geneCounts.lengths, geneCounts.norm);
	//	And then the counts and the bias for the genes themselves
	string filename = normaliseResultsFilename.replaceSuffix("_expression.txt");
	if (!geneCounts.outputGeneCounts(filename, true))
		exitFail("Unable to output counts to :", filename);

	printResults(lastGene);
	printBias();

	//********************************************************************************************
	//	This prints out all of the data for the full set of mcmc runs for each model
	if (outputFull)
	{
		printAllMcmcRunData();
		printConsolidatedMcmcRunData();
	}

	//	The basic count data in htseq-count format
	if (!geneCounts.outputGeneCounts(countsFilename))
		exitFail("Unable to output counts to :", countsFilename);

	progMessage("Data modelled");
	elapsedTime();

#ifdef _WIN32
	string x;
	cin >> x;
#endif
	return EXIT_SUCCESS;
}

bool LiBiNorm::coreParameterEstimation()
{
	geneCounts.remove_invalid_values();
	geneCounts.transferTo(consData, MAX_READS_GENE,maxReads);

	elapsedTime("Data loaded");

	paramSet params;
	optionsType options;
	modelType model;

	options.jumpSize = 0.01;
	options.nsimu = Nsimu;
	options.Nruns = Nruns;

//	double drscale  = 0;
//	double adaptint = 0;

	options.updatesigma = 0;
	options.method = "mh";
	model.sigma2 = 1;

	fullResultSSChain.resize(N_MODELS + 1);
	fullResultChain.resize(N_MODELS + 1);

	RejectionRate.resize(N_MODELS + 1);

	//	Set the number of iterations required of each of the models.
	for (size_t m = 1; m <= N_MODELS; m++)
	{
		if (m == theModel)
			threadLoopCounts[m] = options.Nruns;
		else
			threadLoopCounts[m] = NrunsOtherModels;
	}


#ifdef FIXED_RESULTS
	maxModel = M_FIXED_RESULTS;
	bestResults[maxModel].params = dataVec{ FIXED_RESULTS };
#else
	//	And then set the threads running
	vector<thread> threads;
	for (size_t i = 0; i < Nthreads; i++)
		threads.emplace_back(thread(runThread, this, params, options, model));

	for (auto & i : threads)
		i.join();


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
	//	Find the optimal parameter values, which are associated with the lowest likelihood value found in the last
	//	1000 iterations of all of the runs.

	for (size_t m = 1; m <= N_MODELS; m++)
	{
		multimap <double, dataVec *> & orderedResults = allOrderedResults[m];

		bestResult & br = bestResults[m];
		for (size_t i = 1; i <= fullResultSSChain[m].size(); i++)
		{
			for (size_t j = fullResultSSChain[m][i].size() - 1; j > fullResultSSChain[m][i].size() / 2; j--)
			{
				if (outputFull)
					orderedResults.emplace(fullResultSSChain[m][i][j], &fullResultChain[m][i][j]);

				if (fullResultSSChain[m][i][j] < br.minLL)
				{
					br.minLL = fullResultSSChain[m][i][j];
					br.params = fullResultChain[m][i][j];
					br.run = i;
					br.pos = j;
				}
			}
		}
	}

	//******************************************************************************************
	//	And then find the standard deviation
	for (size_t m = 1; m <= N_MODELS; m++)
	{
		bestResult & br = bestResults[m];
		VEC_DATA_TYPE LL_dev = 0;
		size_t size = br.params.size();
		vector<dataVec> param_dev(2, dataVec(size));
		vector<dataVec> p_N(2, dataVec(size));
		int N = 0;

		for (size_t i = 1; i <= fullResultSSChain[m].size(); i++)
		{
			for (size_t j = fullResultSSChain[m][i].size() - 1; j > fullResultSSChain[m][i].size() / 2; j--)
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
	}
#endif

	return true;
}



void LiBiNorm::printResults(const string & lastGene)
{
	string filename(normaliseResultsFilename.replaceSuffix("_results.txt"));
	TsvFile mcmcResult;
	//	Now output a table with the end points of each of the chains.
	if (!mcmcResult.open(filename))
		exitFail("Unable to open output File ", filename);

	//	First headers up to and including the maximum model that is run.   Always leave space
	//	for the intermediate models so the layout of the results is consistent
	mcmcResult.printStart(lastGene);
	for (size_t m = 1; m <= N_MODELS; m++)
		mcmcResult.printMiddle(headers[m], "chain", "");
	mcmcResult.printEnd();

	//	A row for the optimal parameters that were found for each model
	mcmcResult.printStart("Best");
	for (size_t m = 1; m <= N_MODELS; m++)
	{
		if (bestResults[m].params.size())
			mcmcResult.printMiddle(bestResults[m].params, bestResults[m].minLL, "");
		else
			mcmcResult.printGaps(headers[m].size() + 2);
	}
	mcmcResult.printEnd();
	//	A row for the standard deviations for each model
	for (size_t i = 0; i < 2; i++)
	{
		mcmcResult.printStart((i == 0) ? "Pos Dev" : "Neg Dev");
		for (size_t m = 1; m <= N_MODELS; m++)
		{
			if (bestResults[m].param_dev[i].size())
				mcmcResult.printMiddle(bestResults[m].param_dev[i], bestResults[m].minLL_dev, "");
			else
				mcmcResult.printGaps(headers[m].size() + 2);
		}
		mcmcResult.printEnd();
	}
	//	The number of rows is set by the maximum number of runs, which may be for one or all models
	for (size_t i = 1; i <= Nruns; i++)
	{
		mcmcResult.printStart(_s("Chain end ", i));
		for (size_t m = 1; m <= N_MODELS; m++)
		{
			//	Was there a jth run of this model?  If so then print the end points
			if (fullResultSSChain[m].size() && (i <= fullResultSSChain[m].rbegin()->first))
				mcmcResult.printMiddle(*fullResultChain[m][i].rbegin(), *fullResultSSChain[m][i].rbegin(), "");
			else
				mcmcResult.printGaps(headers[m].size() + 2);
		}
		mcmcResult.printEnd();
	}
	mcmcResult.close();
}


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

	map<size_t, dataVec> norms;

	for (size_t m = 1; m <= N_MODELS; m++)
	{
		if (bestResults[m])
			normaliseExpression(m, bestResults[m].params, lengths,norms[m]);
	}

	for (size_t m = 1; m <= N_MODELS; m++)
	{
		mcmcResult.print(m, bestResults[m].minLL, bestResults[m].run, bestResults[m].pos, bestResults[m].params);
		if (norms[m].size())
		{
			mcmcResult.print(m, lengths);
			mcmcResult.print(m, norms[m]);
		}
		else
		{
			mcmcResult.print(m);
			mcmcResult.print(m);
		}
		mcmcResult.print();
	}
	mcmcResult.close();
}

void LiBiNorm::printAllMcmcRunData()
{
	TsvFile mcmcResult;
	for (size_t modl = 1; modl <= N_MODELS; modl++)
	{
		string filename = normaliseResultsFilename.replaceSuffix("_model_", modl, ".txt");
		if (!mcmcResult.open(filename))
			exitFail("Unable to open output file ", filename);

		//			mcmcResult.printMiddle(headers[modl], "chain", "");

		//	This ensures that at least one header is output, which ensures that there is something in the file
		//	even if no data were produced for this model
		for (size_t i = 0; i < fullResultChain[modl].size(); i++)
			mcmcResult.printMiddle(headers[modl], "chain", "");
		mcmcResult.printEnd();

		//	For each of the mcmc runs print out the results.  Each run is a column
		for (size_t i = 0; i < fullResultChain[modl][1].size(); i++)
		{
			for (size_t j = 1; j <= fullResultChain[modl].rbegin()->first; j++)
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
	for (size_t m = 1; m <= N_MODELS; m++)
	{
		mcmcResult.printMiddle(headers[m], "chain", "");
		iterators[m] = allOrderedResults[m].begin();
	}
	mcmcResult.printEnd("");
	bool found = true;
	for (size_t i = 0; (i < 1000) && found; i++)
	{
		found = false;
		mcmcResult.printStart(i);
		for (size_t m = 1; m <= N_MODELS; m++)
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