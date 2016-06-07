#include <random>
#include "containerEx.h"
#include "libCommon.h"
#include "mcmc.h"
#include "logLiklihoods.h"
#include "LiBiNorm.h"

using namespace std;


double rand(double a)
{
	static std::random_device rd;
	static std::mt19937 gen(rd());
	static std::uniform_real_distribution<> dis(0, 1);
	return dis(gen) * a;
}


int main(int argc, char **argv)
{
	LiBiNorm lbn;
	return lbn.main(argc,argv);
}


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

	paramSet params;
	optionsType options;
	modelType model;

	double JumpSize = 0.01;
	options.nsimu = 2000;
	size_t Nruns = 100;

	double drscale  = 0;
	double adaptint = 0;

	options.updatesigma = 0;

	options.method = method;

	for (size_t Model = 2; Model < 3;Model++)
	{
		switch (Model)
		{
		case 2:
			model.ssfun = &FLL_ModelBD;
			break;
		}
		switch (Model)
		{
		case 2:
			break;
		}


		for (size_t kk = 1; kk <= Nruns; kk++)
		{
			//create input arguments for the dramrun function
			vectorEx<double> p0(rand(3), rand(3), rand(4)-5, rand(4)-5, rand(1));

			switch (Model)
			{
			case 2: case 4: case 5:
			options.qcov = vector<double>(4,JumpSize);

			params = paramSet(paramType("d", p0[1], -1 , 2),    // average length of fragments
				paramType("h",  p0[2], 0 , 3),   // the minimum length of fragmenation
				paramType("t1", p0[3], -5 , -1),   // theta1
				paramType("t2", p0[4], -5, -1), // theta2
				paramType("sig", p0[5], 0, 3)); // sigma
			};
		}


		mcmc mcmcEngine;


		mcmcEngine.mcmcrun(model,consData,params,options);
	}



	cout << "Data loaded";
	string x;
	cin >> x;

	return EXIT_SUCCESS;
}



