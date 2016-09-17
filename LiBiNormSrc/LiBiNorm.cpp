#ifdef _WIN32
#define _CRTDBG_MAP_ALLOC
#include <crtdbg.h>
#endif

#include <stdlib.h>
#include <vector>
#include <mutex>
#include <thread>
#include "libCommon.h"
#include "containerEx.h"
#include "mcmc.h"
#include "LogLiklihoods.h"
#include "LiBiNorm.h"
#include "LiBiCount.h"
#include "LiBiDedup.h"
#include "MakeFastq.h"

using namespace std;

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
		else
			exitFail("Invalid commmand:",command);
	}
	return EXIT_SUCCESS;
}


void LiBiNorm::mcmcThread(paramSet params, optionsType options, modelType model)
{
	while (true)
	{
		size_t loop;
		{
			static mutex mtx; 
			lock_guard<mutex> lock(mtx);
			loop = threadLoopCount++;
			cerr << "Starting " << loop << endl;

		}
#ifdef _DEBUG
//#define _TEST
#endif
#ifdef _TEST
		vectorEx<double> p0(1.5, 1.6,-3.1, -3.2,0.6);
#else
		vectorEx<double> p0(rand(3), rand(3), rand(4)-5, rand(4)-5, rand(1));
#endif
		switch (options.Model)
		{
		case 2: case 4: case 5:
			options.qcov = vector<double>(4,options.jumpSize);

			params = paramSet(paramType("d", p0[0], -1 , 2)    // average length of fragments
				,paramType("h",  p0[1], 0 , 3)   // the minimum length of fragmenation
				,paramType("t1", p0[2], -5 , -1)   // theta1
				,paramType("t2", p0[3], -5, -1) // theta2
				//				,paramType("sig", p0[4], 0, 3) // sigma
				);
			break;
		case 3:
			options.qcov = vector<double>(3,options.jumpSize);

			params = paramSet(paramType("d", p0[0], -1 , 2)    // average length of fragments
				,paramType("h",  p0[1], 0 , 3)   // the minimum length of fragmenation
				//				,paramType("t1", p0[2], -5 , -1)   // theta1
				,paramType("t2", p0[3], -5, -1) // theta2
				//				,paramType("sig", p0[4], 0, 3) // sigma
				);
			break;
		case 1:
			options.qcov = vector<double>(2,options.jumpSize);

			params = paramSet(paramType("d", p0[0], -1 , 2)    // average length of fragments
				,paramType("h",  p0[1], 0 , 3)   // the minimum length of fragmenation
				//				,paramType("t1", p0[2], -5 , -1)   // theta1
				//				,paramType("t2", p0[3], -5, -1) // theta2
				);
			break;
		case 6:
			options.qcov = vector<double>(6,options.jumpSize);

			params = paramSet(paramType("d", p0[0], -1 , 2)    // average length of fragments
				,paramType("h",  p0[1], 0 , 3)   // the minimum length of fragmenation
				,paramType("t1", p0[2], -5 , -1)   // theta1
				,paramType("t2", p0[3], -5, -1) // theta2
				,paramType("a", p0[4], 0, 1) // alpha strength of model B
				);
			break;

		};

		mcmc mcmcEngine;

		mcmcEngine.mcmcrun(model,consData,params,options);

		static mutex mtx; 

		lock_guard<mutex> lock(mtx);

		cerr << "Finishing " << loop << endl;

		Chain.push_back(mcmcEngine.chain().back());
		SSChain[options.Model].push_back(mcmcEngine.sschain().back());
		RejectionRate[options.Model] += mcmcEngine.rejected();
		if (threadLoopCount >= options.Nruns)
			break;
	}

}

void runThread(LiBiNorm * root,	paramSet params, optionsType options, modelType model)
{
	root->mcmcThread(params, options,model);
}

int LiBiNorm::main(int argc, char **argv)
{

	size_t Nthreads = 1;
	size_t minModel = 1;
	size_t maxModel = 6;
	size_t Nruns = 100;

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
		if(strcmp(argv[ni], "-c") == 0)
		{
			consFileName = argv[++ni];
		}
		else if(strcmp(argv[ni], "-p") == 0)
		{
			Nthreads = atoi(argv[++ni]);
		}
		else if(strcmp(argv[ni], "-n") == 0)
		{
			Nruns = atoi(argv[++ni]);
		}
		else if(strcmp(argv[ni], "-m") == 0)
		{
			minModel = atoi(argv[++ni]);
			maxModel = minModel;
		}
		else
		{
			cout << "Invalid parameter";
		}
		ni++;
	}



	transData.loadData(consFileName);
	transData.remove_invalid_values();
	transData.transferTo(consData,100);
	
	cerr << "Data loaded" << endl;


	string method = "mh";

	paramSet params;
	optionsType options;
	modelType model;

	options.jumpSize = 0.01;
	options.nsimu = 2000;
	options.Nruns = Nruns;

#ifdef _DEBUG
	options.Nruns = 6;
#else
//	size_t Nruns = 100;
#endif

//	double drscale  = 0;
//	double adaptint = 0;

	options.updatesigma = 0;

	options.method = method;

	model.sigma2 = 1;
	SSChain.resize(maxModel+1);
	RejectionRate.resize(maxModel+1);

	for (options.Model = minModel; options.Model < maxModel+1;options.Model++)
	{
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



		threadLoopCount = 0;


		vector<thread> threads;
		for (size_t i = 0;i < Nthreads;i++)
			threads.emplace_back(thread(runThread,this,params, options,model));

		for (auto & i : threads)
			i.join();


		TsvFile testResult;
		testResult.open(consFileName.replaceSuffix("_tempOut.txt"));

		for (size_t i = 0;i < Chain.size();i++)
		{
			for (size_t j = minModel;j < maxModel+1;j++) 
			{
				testResult.printMiddle(Chain[i],SSChain[j][i],"");
			}
			testResult.printEnd();
		}
		
	}

	cerr << "Data modelled" << endl;


	string x;
	cin >> x;

	return EXIT_SUCCESS;
}



