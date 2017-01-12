#include "Options.h"
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
		double diff = data[i] - params[i].max;
		if (diff > 0)
			ep.WithWeightings += (diff * diff ) * EDGE_PENALTY_MULTIPLIER;
		else
		{
			diff = params[i].min - data[i];
			if (diff > 0)
				ep.WithWeightings += (diff * diff) * EDGE_PENALTY_MULTIPLIER;
		}
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

	optimise(allOptiData, 50,100, options.jumpSize);

	dataVec results;
	for (size_t i = 0; i < params.size(); i++)
		results.push_back(allOptiData[0]->at(i).value);
	return results;
}
