#ifndef LIBIOPTIMISER_H
#define LIBIOPTIMISER_H
#include "nelderMeadOptimiser.h"
#include "ModelData.h"
#include "mcmc.h"

class LiBiOptimiser : public nelderMeadOptimiser
{
	optiDataType allOptiData;
	paramSet params;

public:

	optionsType * opts;
	const mcmcGeneData & geneData;
	const modelType currentModel;
	LiBiOptimiser(const mcmcGeneData & geneData, modelType currentModel) : geneData(geneData), currentModel(currentModel){};

	dataVec getParams(modelType m,optionsType & options);


	virtual ErrorPair ErrorFunc();
	virtual void SaveResults(bool toFile) {};

};

#endif

