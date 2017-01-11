#include "LiBiOptimiser.h"

ErrorPair LiBiOptimiser::ErrorFunc()
{
	ErrorPair ep;
	dataVec data;
	for (size_t i = 0; i < (*allOptiData[0]).size(); i++)
		data.push_back((*allOptiData[0])[i].value);

	ep = opts->ssfun(data, geneData);

	for (size_t i = 0; i < params.size(); i++)
	{
		if ((data[i] > params[i].max) || (data[i] < params[i].min))
			ep.WithWeightings += 10000;
	}

	return ep;
};


dataVec LiBiOptimiser::getParams(modelType m, optionsType & options)
{
	opts = & options;
	params  = GetModelParams(m);

	optiVector optiData;
	for (auto & p : params)
		optiData.push_back(optiItem((p.max+p.min)/2, true));

	allOptiData.push_back(&optiData);

	VEC_DATA_TYPE minError = optimise(allOptiData, 200,100, options.jumpSize);

	dataVec results;
	for (size_t i = 0; i < params.size(); i++)
		results.push_back(allOptiData[0]->at(i).value);
	return results;
}
