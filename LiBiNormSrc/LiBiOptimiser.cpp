#include "Options.h"
#include "LiBiOptimiser.h"

ErrorPair LiBiOptimiser::ErrorFunc()
{
	ErrorPair ep;
	dataVec data;
	for (size_t i = 0; i < (*allOptiData[0]).size(); i++)
		data.push_back((*allOptiData[0])[i].value);

	ep = opts->ssfun(data, geneData);

	//	Stop parameter d drifting off
	for (size_t i = 0; i < params.size(); i++)
	{
		if ((i == 1) && (currentModel == ModelE))
			ep.WithWeightings += (data[i] * PARAMETER_WEIGHTING_SLOPE);

		VEC_DATA_TYPE diff = data[i] - params[i].max + 0.2;
		if (diff > 0)
			ep.WithWeightings += (diff * diff ) * EDGE_PENALTY_MULTIPLIER;
		else
		{
			diff = params[i].min - data[i] + 0.2;
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

	VEC_DATA_TYPE LL = optimise(allOptiData, NELDER_MEAD_ITERATIONS,20, options.jumpSize,conv(m));

	dataVec results;
	for (size_t i = 0; i < params.size(); i++)
		results.push_back(allOptiData[0]->at(i).value);
	results.push_back(LL);
	return results;
}
