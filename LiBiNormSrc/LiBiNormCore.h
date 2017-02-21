#ifndef LIBINORMCORE_H
#define LIBINORMCORE_H

#include <map>
#include <mutex>
#include "stringEx.h"
#include "containerEx.h"
#include "GeneCountData.h"
#include "Options.h"
#include "ModelData.h"

//	For identifying the mcmc run for a specific model.  Numbered from 1.
typedef size_t mcmcRunId;


//	The class that contains common functionality for LiBiNorm count. model and variation

class LiBiNormCore
{
protected:
	LiBiNormCore() :normalise(true), pauseAtEnd(false),
#ifdef USE_NELDER_MEAD_FOR_INITIAL_VALUES
		nelderMead(true),
		Nsimu(NELDER_MCMC_ITERATIONS),
#else
		nelderMead(false),
		Nsimu(MCMC_ITERATIONS),
#endif
		theModel(noModel),
		maxReads(DEF_MAX_READS_FOR_PARAM_ESTIMATION),
		maxGeneLength(DEF_LENGTH_OF_GENE_FOR_PARAM_ESTIMATION),
		Nthreads(DEF_THREADS),
		Nruns(NUMBER_OF_MCMC_RUNS),
		NrunsOtherModels(0)
	{
		headers = getHeaders();
	};

	void helpCommon();
	bool commandParseCommon(int & ni, int argc, char **argv);
	void SetInitialParamsFromFile(const std::string & filename);

	bool normalise, pauseAtEnd,nelderMead;
	size_t Nsimu;
	modelType theModel;
	size_t maxReads,maxGeneLength,Nthreads;
	mcmcRunId Nruns,NrunsOtherModels;

	headerType headers;
	stringEx landscapeFilename, normaliseResultsFilename, countsFilename, parameterFilename;

	std::map<modelType, dataVec> initialValues;

	//	The specific data that will be used for the mcmc parameter determination
	mcmcGeneData geneData;
	GeneCountData geneCounts;
};

typedef std::map<modelType, std::vector<std::string> > headerType;


#endif


