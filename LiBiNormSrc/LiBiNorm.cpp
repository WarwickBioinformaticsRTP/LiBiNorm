#ifdef _WIN32
#define _CRTDBG_MAP_ALLOC
#include <crtdbg.h>
#endif

#include <stdlib.h>
#include <vector>
#include <mutex>
#include <chrono>
#include <thread>
#include <float.h>
#include "libCommon.h"
#include "containerEx.h"
#include "mcmc.h"
#include "LogLiklihoods.h"
#include "LiBiNorm.h"
#include "LiBiCount.h"
#include "LiBiDedup.h"
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
		if (command == "--version")
		{
			cout << "LiBiNorm version 1.0.2" << endl;
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
			cerr << "Starting Model:" << model_iterator->first << " iteration:" << loop << endl;
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
//		vectorEx<double> p0(1.590355415,0.249897587,-3.760142979,-3.226093068);
//		vectorEx<double> p0(-0.458447029374,1.852235843328);

		vectorEx<double> p0(1.5, 1.6,-3.1, -3.2,0.6);

#else
		vectorEx<double> p0(rand(3)-1, rand(3), rand(4)-5, rand(4)-5, rand(1));
#endif
		
		switch (options.Model)
		{
		case 2: case 4: case 5:
			options.qcov = dataVec(4,options.jumpSize);

			params = paramSet(paramType("d", p0[0], -1 , 2)    // average length of fragments
				,paramType("h",  p0[1], 0 , 3)   // the minimum length of fragmenation
				,paramType("t1", p0[2], -5 , -1)   // theta1
				,paramType("t2", p0[3], -5, -1) // theta2
				//				,paramType("sig", p0[4], 0, 3) // sigma
				);
			break;
		case 3:
			options.qcov = dataVec(3,options.jumpSize);

			params = paramSet(paramType("d", p0[0], -1 , 2)    // average length of fragments
				,paramType("h",  p0[1], 0 , 3)   // the minimum length of fragmenation
				//				,paramType("t1", p0[2], -5 , -1)   // theta1
				,paramType("t2", p0[3], -5, -1) // theta2
				//				,paramType("sig", p0[4], 0, 3) // sigma
				);
			break;
		case 1:
			options.qcov = dataVec(2,options.jumpSize);

			params = paramSet(paramType("d", p0[0], -1 , 2)    // average length of fragments
				,paramType("h",  p0[1], 0 , 3)   // the minimum length of fragmenation
				//				,paramType("t1", p0[2], -5 , -1)   // theta1
				//				,paramType("t2", p0[3], -5, -1) // theta2
				);
			break;
		case 6:
			options.qcov = dataVec(6,options.jumpSize);

			params = paramSet(paramType("d", p0[0], -1 , 2)    // average length of fragments
				,paramType("h",  p0[1], 0 , 3)   // the minimum length of fragmenation
				,paramType("t1", p0[2], -5 , -1)   // theta1
				,paramType("t2", p0[3], -5, -1) // theta2
				,paramType("a", p0[4], 0, 1) // alpha strength of model B
				);
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

			cerr << "Finishing Model:" << model_iterator->first << " iteration:" << loop << endl;


#ifdef STORE_ENDPOINTS
			Chain[options.Model].emplace(loop, mcmcEngine.chain().back());
			SSChain[options.Model].emplace(loop, mcmcEngine.sschain().back());
#endif

			if (fullOutputMode || singleModel)
			{
				fullResultChain[options.Model].emplace(loop, mcmcEngine.chain());
				fullResultSSChain[options.Model].emplace(loop, mcmcEngine.sschain());
			}


			RejectionRate[options.Model] += mcmcEngine.rejected();
		}
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
	size_t Nsimu = 2000;
	fullOutputMode = false;
	singleModel = false;

	initClock();
	if(argc < 1)
	{
		printf("Error: parameter wrong!\n");
		return EXIT_FAILURE;
	}
	else if(argc == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("     LiBiNorm:    RNA-seq library bias normalisation   \n\n");
		printf("Options:\n");
		printf("  -h, --help            show this help message and exit\n");
		printf("  -c <filename>\n");
		printf("                        Location of consolidated location file data\n");
		printf("  -o <filename>\n");
		printf("                        Optional root for output file.  Default is that it\n");
		printf("                        is derived from the location file name\n");
		printf("  -p N\n");
		printf("                        Number of threads\n");
		printf("  -n N\n");
		printf("                        Number of mcmc iterations (100)\n");
		printf("  -s N\n");
		printf("                        Length of each simulation (100)\n");
		printf("  -m N\n");
		printf("                        Just run for model N -n times.  All other models run once\n");
		printf("  -M N\n");
		printf("                        Just run for model N -n times.  All other models not run\n");
		printf("  -f\n");
		printf("                        Output complete set of output filesN\n");
		return EXIT_SUCCESS;
	}

	int ni = 1;
	while(ni < argc)
	{
		if(strcmp(argv[ni], "-c") == 0)
		{
			consFileName = argv[++ni];
			if (!outputFileName)
				outputFileName = consFileName;
		}
		else if (strcmp(argv[ni], "-o") == 0)
		{
			outputFileName = argv[++ni];
		}
		else if(strcmp(argv[ni], "-p") == 0)
		{
			Nthreads = atoi(argv[++ni]);
		}
		else if(strcmp(argv[ni], "-n") == 0)
		{
			Nruns = atoi(argv[++ni]);
		}
		else if (strcmp(argv[ni], "-s") == 0)
		{
			Nsimu = atoi(argv[++ni]);
		}
		else if(strcmp(argv[ni], "-m") == 0)
		{
			minModel = atoi(argv[++ni]);
			maxModel = minModel;
		}
		else if (strcmp(argv[ni], "-M") == 0)
		{
			minModel = atoi(argv[++ni]);
			maxModel = minModel;
			singleModel = true;
		}
		else if (strcmp(argv[ni], "-f") == 0)
		{
			fullOutputMode = true;
		}
		else
		{
			exitFail("Invalid parameter: ",argv[ni]);
		}
		ni++;
	}



	transData.loadData(consFileName);
	transData.remove_invalid_values();
	transData.transferTo(consData,100);
	
	cerr << "Data loaded" << endl;
	elapsedTime();


	string method = "mh";

	paramSet params;
	optionsType options;
	modelType model;

	options.jumpSize = 0.01;
	options.nsimu = Nsimu;
	options.Nruns = Nruns;

#ifdef _DEBUG
//	options.Nruns = 6;
//	options.nsimu = 100;
#else
//	size_t Nruns = 100;
#endif

//	double drscale  = 0;
//	double adaptint = 0;

	options.updatesigma = 0;

	options.method = method;

	model.sigma2 = 1;
#ifdef STORE_ENDPOINTS
	SSChain.resize(maxModel+1);
	Chain.resize(maxModel+1);
#endif
	if (fullOutputMode || singleModel)
	{
		fullResultSSChain.resize(maxModel + 1);
		fullResultChain.resize(maxModel + 1);
	}
	RejectionRate.resize(maxModel+1);

	//	Set the number of iterations required of each of the models.
	for (size_t m = 1; m < maxModel+1;m++)
	{
		if (m >= minModel)
			threadLoopCounts[m] = options.Nruns;
		else
			threadLoopCounts[m] = singleModel?0:1;
	}

	struct bestResult
	{
		bestResult() :minLL(DBL_MAX), run(0), pos(0) {};
		VEC_DATA_TYPE minLL;
		size_t run, pos;
		dataVec params;
		dataVec norm;
		operator bool() const { return params.size(); };
	};

	map<size_t, bestResult> bestResults;

#ifdef FIXED_RESULTS
	maxModel = M_FIXED_RESULTS;
	bestResults[maxModel].params = dataVec{ FIXED_RESULTS };
#else
	//	And then set the threads running
	vector<thread> threads;
	for (size_t i = 0;i < Nthreads;i++)
		threads.emplace_back(thread(runThread,this,params, options,model));

	for (auto & i : threads)
		i.join();


	/*
		Unfinished code for 
	vector<int> geneL{ 500,1000,2000,4000,8000 };
	dataVec TotalReads;

	for (size_t i = 0; i < geneL.size(); i++)
	{
		for (size_t j = 0; j < transData.size(); j++)
		{
			if (abs(transData[j].length - geneL[i]) < (0.1 * geneL[i]))
				TotalReads.append(transData[j].counts[0]/ transData[j].length);
		}
	}
	*/


	//********************************************************************************************
	//	Find the optimal parameter values, which are associated with the lowest likelyhood value found in the last 
	//	1000 iterations of all of the runs.

	for (size_t m = 1; m <= maxModel; m++)
	{
		bestResult & br = bestResults[m];
		for (size_t i = 1; i <= fullResultSSChain[m].size(); i++)
		{
			for (size_t j = fullResultSSChain[m][i].size()-1; j > fullResultSSChain[m][i].size()/2; j--)
			{
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
#endif

	dataVec l;
	for (size_t i = 100; i <= 10000; i += 100)
		l.push_back(i);

	for (size_t m = 1; m <= maxModel; m++)
	{
		if (bestResults[m])
		{
			double d = pow(10, bestResults[m].params[0]);
			double h = pow(10, bestResults[m].params[1]);
			double t1, t2, a;
			switch (m)
			{
			case 2:
			case 4:
			case 5:
				t1 = pow(10, bestResults[m].params[2]);
				t2 = pow(10, bestResults[m].params[3]);
				break;
			case 3:
				t1 = 0;
				t2 = pow(10, bestResults[m].params[2]);
				break;
			case 6:
				t1 = pow(10, bestResults[m].params[2]);
				t2 = pow(10, bestResults[m].params[3]);
				a = bestResults[m].params[4];
				break;
			}

			dataVec & norm = bestResults[m].norm;
			norm.resize(l.size());

			switch (m)
			{
			case 1:
				norm = (2 * h<l)*(l - 2 * h) + l / d;
				break;
			case 2:
				norm = ((2 * h<l)*(t1*(exp(-2 * h*(t1 + t2)) - exp(-l*(t1 + t2))) + t2*(t1 + t2)*(l - 2 * h)*exp(-l*(t1 + t2))) / ((t1 + t2)*(t1 + t2)) +
					exp(-l*(t1 + t2))*(l*t2*t2 + t1*(exp(l*(t1 + t2)) + l*t2 - 1)) / ((t1 + t2)*(t1 + t2)) / d);
				break;
			case 3:
/*				for (size_t i = 0; i < l.size(); i++)
				{
					if (2 * h < l[i])
						norm[i] = (exp(-2 * h*t2 - l[i] * t1) - exp(-l[i] * (t1 + t2))) / t2 + (exp(-l[i] * t1) - exp(-l[i] * (t1 + t2))) / t2 / d;
					else
						norm[i] = (exp(-l[i] * t1) - exp(-l[i] * (t1 + t2))) / t2 / d;
				}
*/
			{
				dataVec exp_ml_t2 = exp(-l*(t2));
				norm = (2 * h < l)*(exp(-2 * h*t2) - exp_ml_t2) / t2 + (1 - exp_ml_t2) / t2 / d;
			}

				break;
			case 4:
				for (size_t i = 0; i < l.size(); i++)
				{
					if (2 * h < l[i])
						norm[i] = (exp(-2 * h*(t1 + t2)) - exp(-l[i] * (t1 + t2))) / (t1 + t2) + (1 - exp(-l[i] * (t1 + t2))) / (t1 + t2) / d;
					else
						norm[i] = (1 - exp(-l[i] * (t1 + t2))) / (t1 + t2) / d;
				}
				break;
			case 5:
/*  MATLAB
				if (2 * h<l(i))
					norm(i) = (exp(-l(i)*t1 - 2 * h*t2)*(t1 + t2) ^ 2 - exp(-l(i)*(t1 + t2))*t1 ^ 2 + t1*t2*exp(-2 * h*(t1 + t2))*(l(i)*t2 - 2 * h*t1 - 2 * h*t2 + l(i)*t1 - t2 / t1 - 2)) / (t1 + t2) ^ 2 / t1 ^ 2 / t2 + ...
					(l(i) - 1 / (t1 + t2) - 1 / t1 - t1 / t2 / (t1 + t2)*exp(-l(i)*(t1 + t2)) + (t1 + t2) / t1 / t2*exp(-l(i)*t1)) / (t1 + t2) / t1 / d;
				else
					norm(i) = (l(i) - 1 / (t1 + t2) - 1 / t1 - t1 / t2 / (t1 + t2)*exp(-l(i)*(t1 + t2)) + (t1 + t2) / t1 / t2*exp(-l(i)*t1)) / (t1 + t2) / t1 / d;
*/
/*				for (size_t i = 0; i < l.size(); i++)
				{
					if (2 * h < l[i])
						norm[i] = (exp(-l[i]*t1 - 2 * h*t2)*(t1 + t2)*(t1 + t2) - exp(-l[i]*(t1 + t2))*t1*t1 + t1*t2*exp(-2 * h*(t1 + t2))*(l[i]*t2 - 2 * h*t1 - 2 * h*t2 + l[i]*t1 - t2 / t1 - 2)) / ((t1 + t2)*(t1 + t2)) / (t1 *t1)/ t2 +
						(l[i] - 1 / (t1 + t2) - 1 / t1 - t1 / t2 / (t1 + t2)*exp(-l[i]*(t1 + t2)) + (t1 + t2) / t1 / t2*exp(-l[i]*t1)) / (t1 + t2) / t1 / d;
					else

						norm[i] = (l[i] - 1 / (t1 + t2) - 1 / t1 - t1 / t2 / (t1 + t2)*exp(-l[i]*(t1 + t2)) + (t1 + t2) / t1 / t2*exp(-l[i]*t1)) / (t1 + t2) / t1 / d;
				}
*/
				
/*	MATLAB
	norm =  (2*h<l).*(exp(-l*t1 - 2*h*t2)*(t1 + t2)^2 - exp(-l*(t1 + t2))*t1^2 + t1*t2*exp(-2*h*(t1 + t2))*(l*t2 -2*h*t1 -2*h*t2+l*t1 - t2/t1 - 2))/(t1 + t2)^2/t1^2/t2 + ...
	    (l-1/(t1 + t2) - 1/t1 - t1/t2/(t1+t2)*exp(-l*(t1 + t2))+(t1 + t2)/t1/t2*exp(-l*t1))/(t1 + t2)/t1/d;
*/

				norm = (2 * h<l)*(exp(-l*t1 - 2 * h*t2)*(t1 + t2)*(t1 + t2) - exp(-l*(t1 + t2))*t1*t1 + t1*t2*exp(-2 * h*(t1 + t2))*(l*t2 - 2 * h*t1 - 2 * h*t2 + l*t1 - t2 / t1 - 2)) / ((t1 + t2) * (t1 + t2) ) / (t1 * t1 ) / t2 +
					(l - 1 / (t1 + t2) - 1 / t1 - t1 / t2 / (t1 + t2)*exp(-l*(t1 + t2)) + (t1 + t2) / t1 / t2*exp(-l*t1)) / (t1 + t2) / t1 / d;

				norm /= t1;

				break;
			case 6:
			{
				norm = a*((2 * h < l)*(t1*(exp(-2 * h*(t1 + t2)) - exp(-l*(t1 + t2))) + t2*(t1 + t2)*(l - 2 * h)*exp(-l*(t1 + t2))) / ((t1 + t2) * (t1 + t2)) +
					(exp(-l*(t1 + t2))*(l*t2 *t2 + l*t2*t1 - t1) + t1) / ((t1 + t2) *(t1 + t2)) / d) +
					(1 - a)*((2 * h < l)*(exp(-2 * h*(t1 + t2)) - exp(-l*(t1 + t2))) / (t1 + t2) +
					(1 - exp(-l*(t1 + t2))) / (t1 + t2) / d);


				//					norm = a*((2 * h<l).*(t1.*(exp(-2 * h*(t1 + t2)) - exp(-l.*(t1 + t2))) + t2*(t1 + t2).*(l - 2 * h).*exp(-l.*(t1 + t2))) / (t1 + t2) ^ 2 + ...
				//						(exp(-l.*(t1 + t2)).*(l.*t2 ^ 2 + l.*t2*t1 - t1) + t1) / (t1 + t2) ^ 2 / d) + ...
				//						(1 - a)*((2 * h<l).*(exp(-2 * h*(t1 + t2)) - exp(-l.*(t1 + t2))) / (t1 + t2) + ...
				//						(1 - exp(-l.*(t1 + t2))) / (t1 + t2) / d);
				break;
			}
			}
			norm = norm * l[9] / norm[9];
			norm /= l;

		}
	}

	//********************************************************************************************
	TsvFile mcmcResult;

	stringEx filename(outputFileName.replaceSuffix("_norm.txt"));
	if (!mcmcResult.open(filename))
		exitFail("Unable to open output File ", filename);

	for (size_t m = 1; m <= maxModel; m++)
	{
		mcmcResult.print(m, bestResults[m].minLL, bestResults[m].run, bestResults[m].pos,bestResults[m].params);
		if (bestResults[m].norm.size())
		{
			mcmcResult.print(m, l);
			mcmcResult.print(m, bestResults[m].norm);
		}
		else
		{
			mcmcResult.print(m);
			mcmcResult.print(m);
		}
		mcmcResult.print();
	}
	mcmcResult.close();


#ifdef STORE_ENDPOINTS
	//	Now output a table with the end points of each of the chains.
	filename = outputFileName.replaceSuffix("_Chain.txt");
	if (!mcmcResult.open(filename))
		exitFail("Unable to open output File ", filename);

	//	First headers up to and including the maximum model that is run.   Always leave space
	//	for the intermediate models so the layout of the results is consistent
	for (size_t m = 1; m <= maxModel; m++)
		mcmcResult.printMiddle(headers[m], "chain", "");
	mcmcResult.printEnd();

	//	A row for the optimal parameters that were found for each model
	for (size_t m = 1; m <= maxModel; m++)
	{
		if (bestResults[m].params.size())
			mcmcResult.printMiddle(bestResults[m].params, bestResults[m].minLL, "");
		else
			mcmcResult.printGaps(headers[m].size() + 2);
	}
	mcmcResult.printEnd();
	//	The number of rows is set by the maximum number of runs, which may be for one or all models
	for (size_t i = 1; i <= options.Nruns; i++)
	{
		for (size_t m = 1; m <= maxModel; m++)
		{
			//	Was there a jth run of this model?  If so then print the end points
			if (SSChain[m].size() && (i <= SSChain[m].rbegin()->first))
				mcmcResult.printMiddle(Chain[m][i], SSChain[m][i], "");
			else
				mcmcResult.printGaps(headers[m].size() + 2);
		}
		mcmcResult.printEnd();
	}
	mcmcResult.close();
#endif
	//********************************************************************************************
	//	This prints out all of the data for the full set of mcmc runs for each model
	if (fullOutputMode)
	{

		for (size_t modl = 1; modl <= maxModel; modl++)
		{
			string filename = outputFileName.replaceSuffix("_model_", modl, ".txt");
			if(!mcmcResult.open(filename))
				exitFail("Unable to open output file ", filename);

//			mcmcResult.printMiddle(headers[modl], "chain", "");

			//	This ensures that at least one header is output, which ensures that there is something in the file
			//	even if no data were produced for this model
			for (size_t i = 0;i < fullResultChain[modl].size();i++)
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

	if (singleModel || fullOutputMode)
	{
		for (size_t modl = singleModel?minModel:1; modl <= maxModel; modl++)
		{
			multimap <double, dataVec *> orderedResults;

			for (size_t i = 0; i < fullResultChain[modl][1].size(); i++)
			{
				for (size_t j = 1; j <= fullResultChain[modl].rbegin()->first; j++)
				{
					orderedResults.emplace(fullResultSSChain[modl][j][i], &fullResultChain[modl][j][i]);
				}
			}

			string filename = outputFileName.replaceSuffix("_model_cons_", modl, ".txt");
			if (!mcmcResult.open(filename))
				exitFail("Unable to open output File ", filename);

			mcmcResult.print(headers[modl], "chain");
			for (multimap <double, dataVec *>::iterator j = orderedResults.begin(); j != orderedResults.end(); j++)
			{
				mcmcResult.print(*(j->second), j->first);
			}
			mcmcResult.close();

		}
	}

	cerr << "Data modelled" << endl;
	elapsedTime();

#ifdef _WIN32
	string x;
	cin >> x;
#endif
	return EXIT_SUCCESS;
}



