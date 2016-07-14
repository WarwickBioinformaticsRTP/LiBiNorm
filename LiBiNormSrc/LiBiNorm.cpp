#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>
#include <vector>
#include "libCommon.h"
#include "containerEx.h"
#include "mcmc.h"
#include "logLiklihoods.h"
#include "LiBiNorm.h"

using namespace std;

int main(int argc, char **argv)
{
	_CrtSetDbgFlag ( _CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF );

	//	New data management 100 runs 54 secs, no datavec memory management 50.4, ie leave it to the OS to sort out
	LiBiNorm libi;
	libi.main(argc,argv);

}

//#ifdef XXXX
int LiBiNorm::main(int argc, char **argv)
{


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
		else
		{
			cout << "Invalid parameter";
		}
		ni++;
	}

	transData.loadData(consFileName);
	transData.remove_invalid_values();
	transData.transferTo(consData,100);

	string method = "mh";

	vector<vector<dataVec> > Chain;
	paramSet params;
	optionsType options;
	modelType model;

	double JumpSize = 0.01;
	options.nsimu = 2000;
#ifdef _DEBUG
	size_t Nruns = 10;
#else
	size_t Nruns = 100;
#endif

	double drscale  = 0;
	double adaptint = 0;

	options.updatesigma = 0;

	options.method = method;

	model.sigma2 = 1;

	for (size_t Model = 2; Model < 3;Model++)
	{
		switch (Model)
		{
		case 2:
			model.ssfun = &FLL_ModelB;
			break;
		}

		for (size_t kk = 1; kk <= Nruns; kk++)
		{
			cout << kk << endl; 

#ifdef _DEBUG
#define _TEST
#endif
#ifdef _TEST
			vectorEx<double> p0(1.5, 1.6,-3.1, -3.2);
#else
			vectorEx<double> p0(rand(3), rand(3), rand(4)-5, rand(4)-5, rand(1));
#endif
			switch (Model)
			{
			case 2: case 4: case 5:
			options.qcov = vector<double>(4,JumpSize);

			params = paramSet(paramType("d", p0[0], -1 , 2)    // average length of fragments
				,paramType("h",  p0[1], 0 , 3)   // the minimum length of fragmenation
				,paramType("t1", p0[2], -5 , -1)   // theta1
				,paramType("t2", p0[3], -5, -1) // theta2
//				,paramType("sig", p0[4], 0, 3) // sigma
				);
			};

			mcmc mcmcEngine;
			
			mcmcEngine.mcmcrun(model,consData,params,options);

//			Chain.push_back(mcmcEngine.chain());

		}
		int a = 1;
	}



	cout << "Data loaded";
	string x;
	cin >> x;

	return EXIT_SUCCESS;
}



