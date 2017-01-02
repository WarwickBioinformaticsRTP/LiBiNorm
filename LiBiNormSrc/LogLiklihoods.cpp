
#include "LogLiklihoods.h"
#include "containerEx.h"
#include <map>
using namespace std;


/*

This code implements the six different models for the bias within an RNA transcript.

In each case the parameters (between 2 and 5 ) are passed in with param and the data, ie the information about the reads
and the genes are passed in in data.

data.fragData contains the information about the individual reads, one entry per read
data.geneIndex indeicates which gene a read is associated with, one entry per read
data.geneData gives inforation about each gene
data.genedata[0] is the lengths of the genes
data.geneData[1] are the frequencies

One modification arises from teh fact that the normalisation values are the same for all reads in the gene so only need
to be calculated per gene.  The (x) method uses the geneIndex to expand the one per gene vector to a one per read vector so that the normalisation 
can be performed for each read.

*/

//	Returns a list of all the models, which is used to iterate through the list
const std::vector<modelType> & allModels()
{
	static std::vector<modelType> list{ ModelA ,ModelB ,ModelC,ModelD,ModelE,ModelBD };
	return list;
};

std::string conv(const modelType m)
{
	switch (m)
	{
	case noModel: return "No Model";
	case ModelA: return "Model A";
	case ModelB: return "Model B";
	case ModelC: return "Model C";
	case ModelD: return "Model D";
	case ModelE: return "Model E";
	case ModelBD: return "Model BD";
	}
	return "";
}


namespace std
{
	string to_string(const modelType & m)
	{
		return conv(m);
	}
}
ostream& operator<< (ostream &out, const modelType & m)
{
	out << conv(m);
	return out;
}

//  Selects a model based on a string, exits if the string is not valid
modelType modelFromString(const string & desc)
{
	static map<string, modelType> mappings{
		{"A",ModelA },{"B",ModelB },{"C",ModelC },{"D",ModelD },{"E",ModelE },{"BD",ModelBD }, 
		{"a",ModelA },{"b",ModelB },{"c",ModelC },{"d",ModelD },{"e",ModelE },{"bd",ModelBD }, 
		{"SMART",ModelA},{ "smart",ModelA },{"polya",ModelB },{ "POLYA",ModelB },{ "polyA",ModelB },
		{"random",ModelE},{ "RANDOM",ModelE }
	};
	auto iter = mappings.find(desc);
	if (iter == mappings.end())
		exitFail("Unknown model description:", desc);
	return (*iter).second;
}

bool printVal(outputDataFile * f, modelType m)
{
	fputs(conv(m).c_str(), f->fout);
	return true;
};

paramSet GetModelParams(modelType model, optionsType & options)
{

	//	Use this to run the model with a specific set of parameters
	// #define _TEST
#ifdef _TEST
	vectorEx<double> p0{ { 1.5, 1.6,-3.1, -3.2,0.6 } };
#else
	vectorEx<double> p0{ { rand(3) - 1, rand(3), rand(4) - 5, rand(4) - 5, rand(1) } };
#endif

	paramSet params;

	switch (model)
	{
	case noModel:
		break;
	case ModelB: case ModelD: case ModelE:
		options.qcov = dataVec(4, options.jumpSize);
		params = { { "log d", p0[0], -1 , 2 }    // average length of fragments
			,{ "log h",  p0[1], 0 , 3 }   // the minimum length of fragmenation
			,{ "log t1", p0[2], -5 , -1 }   // theta1
			,{ "log t2", p0[3], -5, -1 } // theta2
		};

		break;
	case ModelC:
		options.qcov = dataVec(3, options.jumpSize);
		params = { { "log d", p0[0], -1 , 2 }    // average length of fragments
			,{ "log h",  p0[1], 0 , 3 }   // the minimum length of fragmenation
			,{ "log t2", p0[3], -5, -1 } // theta2
		};
		break;
	case ModelA:
		options.qcov = dataVec(2, options.jumpSize);
		params = { { "log d", p0[0], -1 , 2 }    // average length of fragments
			,{ "log h",  p0[1], 0 , 3 }   // the minimum length of fragmenation
		};
		break;
	case ModelBD:
		options.qcov = dataVec(5, options.jumpSize);
		params = { { "log d", p0[0], -1 , 2 }    // average length of fragments
			,{ "log h",  p0[1], 0 , 3 }   // the minimum length of fragmenation
			,{ "log t1", p0[2], -5 , -1 }   // theta1
			,{ "log t2", p0[3], -5, -1 } // theta2
			,{ "a", p0[4], 0, 1 } // alpha strength of model B
		};
		break;

	};

	return params;
}

headerType getHeaders()
{
	optionsType options;
	headerType _retVal;
	for (modelType m : allModels())
	{
		paramSet params = GetModelParams(m, options);
		for (auto i : params)
			_retVal[m].push_back(i.name);
	}
	return _retVal;
}




//	The log liklyhood calculations for each of the models
double FLL_ModelA(const dataVec & param, const dataType & data)
{
	double d = pow(10,param[0]);
	double h = pow(10,param[1]);

	const vector<int> & geneIndex = data.geneIndex;
	const dataVec L = data.geneData[0](geneIndex);
	const dataVec Freq_l = data.geneData[1](geneIndex);

	double last_l = 0;
	double norm=0;
	double LogL = 0;

	for (size_t i = 0;i < data.fragData.size(); i++)
	{
		const double & x = data.fragData[i];
		const double & l = L[i];
		const double & freq_l = Freq_l[i];

//	f_frag = (x> h).*(x < l-h) + 1/d;

		double f_frag = 1.0/d;
		if ((x> h) && (x < (l-h)))
			f_frag += 1.0;


		if (f_frag == 0)
				return 1E20;

		if (l != last_l)
		{


//	norm =  (2*h<l).*(l-2*h) + l/d;

			norm =  l/d;
			if (2*h<l)
				norm +=  (l-(2*h));
			last_l = l;
		}

		LogL += log(f_frag/norm)/freq_l;
	}


//LogL = -2*sum(log(f_frag./norm)./freq_l);

	LogL = -2*LogL;
	return LogL;
}


double FLL_ModelB(const dataVec & param, const dataType & data)
{

	double d = pow(10,param[0]);
	double h = pow(10,param[1]);
	double t1 = pow(10,param[2]);
	double t2 = pow(10,param[3]);

	const vector<int> & geneIndex = data.geneIndex;
	const dataVec L = data.geneData[0](geneIndex);
	const dataVec & Freq_l = data.geneData[1](geneIndex);

	double last_l = 0;
	double norm=0;
	double LogL = 0;

	//	f_frag = (x> h).*(x < l-h)./(t1+ t2) .* (t1.*exp(-2*l*(t1+t2)+(t1+t2)*(l-h+x)) + t2*exp(-l*(t1+t2))) + ...
	//         1./(t1+ t2) .* (t1.*exp(-2*l*(t1+t2)+(t1+t2)*(l+x)) + t2*exp(-l*(t1+t2)))/d;

	double t1_p_t2 = t1+t2;

	for (size_t i = 0;i < data.fragData.size(); i++)
	{
		const double & x = data.fragData[i];
		const double & l = L[i];
		const double & freq_l = Freq_l[i];

		double exp_l_t1_t2 = exp(-l * t1_p_t2);
		double t2_exp_l_t1_t2_full = t2*exp_l_t1_t2;
		double t1_exp_c_l_t1_t2 = t1*exp((x-l)*t1_p_t2);

		double f_frag = (t1_exp_c_l_t1_t2 + t2_exp_l_t1_t2_full)/d;
		if ((x> h) && (x < (l-h)))
			f_frag += (t1_exp_c_l_t1_t2/exp(h*t1_p_t2) + t2_exp_l_t1_t2_full);

		f_frag /= t1_p_t2;

		if (f_frag == 0)
			return 1E20;

		if (l != last_l)
		{

		//	    norm =  (2*h<l).*(t1.*(exp(-2*h*(t1+t2))-exp(-l*(t1+t2)))+t2*(t1+t2).*(l-2*h).*exp(-l*(t1+t2)))/(t1 + t2)^2 + ...
		//        (exp(-l*(t1+t2)).*(l.*t2^2+l.*t2*t1-t1)+t1)/(t1 + t2)^2/d;

			norm = ((exp_l_t1_t2*(l*t2*t1_p_t2-t1)+t1)/d);
			if ((2*h)<l)
				norm += (t1*(exp(-2*h*t1_p_t2)-exp_l_t1_t2)+t2*t1_p_t2*(l-(2*h))*exp_l_t1_t2);
			norm /= (t1_p_t2*t1_p_t2);

			last_l = l;
		}

		LogL += log(f_frag/norm)/freq_l;
	}

	return -2*LogL;
}

double FLL_ModelC(const dataVec & param, const dataType & data)
{
	//	function [LogL] = FLL_Deng(param, data)

	double d = pow(10,param[0]);
	double h = pow(10,param[1]);
	double t2 = pow(10,param[2]);

	const vector<int> & geneIndex = data.geneIndex;
	const dataVec L = data.geneData[0](geneIndex);
	const dataVec & Freq_l = data.geneData[1](geneIndex);

	double last_l = 0;
	double norm=0;
	double LogL = 0;

//    f_frag = (x> h).*(x < l-h).*exp(-t2*(x+h)) + ...
//    exp(-t2*(x))/d;

	for (size_t i = 0;i < data.fragData.size(); i++)
	{
		const double & x = data.fragData[i];
		const double & l = L[i];
		const double & freq_l = Freq_l[i];

		double f_frag = exp(-t2*(x))/d;
		if ((x> h) && (x < l-h))
			f_frag += exp(-t2*(x+h));

		if (f_frag == 0)
			return 1E20;


		if (l != last_l)
		{


//    norm =  (2*h<l).*(exp(-2*h*t2) - exp(-l*( t2)))/t2 + ...
//        (1-exp(-l*( t2)))/t2/d;

			double exp_ml_t2 = exp(-l*( t2));
//			dataVec normX =  (2*h<L)*(exp(-2*h*t2) - exp_mL_t2)/t2 + (1-exp_mL_t2)/t2/d;

			norm = (1-exp_ml_t2)/t2/d;
			if (2*h<l)
				norm += (exp(-2*h*t2) - exp_ml_t2)/t2;

//		LogL = -2*sum(log(f_frag/norm(geneIndex))/freq_l(geneIndex));

			last_l = l;
		}

		LogL += log(f_frag/norm)/freq_l;
	}

	return -2*LogL;
}


double FLL_ModelD(const dataVec & param, const dataType & data)
{

/*function [LogL] = FLL_Tang(param, data)

d = 10^(param(1));
h = 10^(param(2));
t1 = 10^(param(3));
t2 = 10^(param(4));
%sig = 10^(param(5));
x = data(1, :);
l = data(2, :);
freq_l = data(3, :);
*/
	double d = pow(10,param[0]);
	double h = pow(10,param[1]);
	double t1 = pow(10,param[2]);
	double t2 = pow(10,param[3]);
	double LogL = 1E20;

	const vector<int> & geneIndex = data.geneIndex;
	const dataVec L = data.geneData[0](geneIndex);
	const dataVec & Freq_l = data.geneData[1](geneIndex);

	double last_l = 0;
	double norm=0;
	double temp_l = 0;

//f_frag = (x> h).*(x < l-h)./(t1+ t2) .* (t1.*exp(-t1*(l-x) - 2*t2*h - t1*h) + t2.*exp(-t1*l-t2*(x+h))) + ...
//         1./(t1+ t2) .* (t1.*exp(-t1*(l-x)) + t2.*exp(-t1*l-t2*(x)))/d;

	double t1_p_t2 = t1+t2;
//	double t1_p_t2_sq = t1_p_t2*t1_p_t2;


	for (size_t i = 0;i < data.fragData.size(); i++)
	{
		const double & x = data.fragData[i];
		const double & l = L[i];
		const double & freq_l = Freq_l[i];

		double f_frag = 1/t1_p_t2 * (t1*exp(-t1*(l-x)) + t2*exp(-t1*l-t2*(x)))/d;

		if ((x > h) && (x < l-h))
			f_frag += (t1*exp(-t1*(l-x) - 2*t2*h - t1*h) + t2*exp(-t1*l-t2*(x+h)))/t1_p_t2;

		//if (isempty(find(f_frag == 0, 1)))
		if (f_frag == 0)
			return 1E20;


		if (l != last_l)
		{

			//    norm =  (2*h<l).*(exp(-2*h*(t1 + t2)) - exp(-l*(t1 + t2)))/(t1 + t2) + ...
			//    (1-exp(-l*(t1 + t2)))/(t1 + t2)/d;

			double exp_ml_t1_p_t2 = exp(-l*t1_p_t2);
			norm =  (2*h<l)*(exp(-2*h*t1_p_t2) - exp_ml_t1_p_t2)/t1_p_t2 + 
				(1-exp_ml_t1_p_t2)/t1_p_t2/d;

			last_l = l;


//			LogL = -2*sum(log(f_frag/norm(geneIndex))/Freq_l(geneIndex));
		}

		temp_l += log(f_frag/norm)/freq_l;
	}

	LogL = -2*temp_l;
	return LogL;
}

double FLL_ModelE(const dataVec & param, const dataType & data)
{
//	function [LogL] = FLL_Dan(param, data)

	double d = pow(10,param[0]);
	double h = pow(10,param[1]);
	double t1 = pow(10,param[2]);
	double t2 = pow(10,param[3]);

	const vector<int> & geneIndex = data.geneIndex;
	const dataVec L = data.geneData[0](geneIndex);
	const dataVec & Freq_l = data.geneData[1](geneIndex);

	double last_l = 0;
	double norm=0;
	double LogL = 0;

//	f_frag = (x> h).*(x < l-h)./t1/(t1+ t2).*(exp(-2*h*(t1+ t2))-exp(-(h +x)*(t1+t2)) - exp(-t1*h-2*h*t2-(l-x)*t1) + exp(-h*t2-l*t1-x*t2)) + ...
//     1./t1/(t1+ t2) .*(1 - exp(-x*(t1+t2)) - exp(-(l-x)*t1) + exp(-l*t1-x*t2))/d;

	double t1_p_t2 = t1+t2;
	double t1_p_t2_sq = t1_p_t2*t1_p_t2;
	double exp_m2_h_t1_p_t2 = exp(-2*h*t1_p_t2);

	for (size_t i = 0;i < data.fragData.size(); i++)
	{
		const double & x = data.fragData[i];
		const double & l = L[i];
		const double & freq_l = Freq_l[i];

		double	f_frag = 1/t1/t1_p_t2 *(1 - exp(-x*t1_p_t2) - exp(-(l-x)*t1) + exp(-l*t1-x*t2))/d;

		if ((x> h) && (x < l-h))
			f_frag += (exp_m2_h_t1_p_t2-exp(-(x + h)*t1_p_t2) - exp(-t1*h-2*h*t2-(l-x)*t1) + exp(-h*t2-l*t1-x*t2))/t1/t1_p_t2;

		if (f_frag == 0)
			return 1E20;

		if (l != last_l)
		{

//		   norm =  (2*h<l).*(exp(-l*t1 - 2*h*t2)*(t1 + t2)^2 - exp(-l*(t1 + t2))*t1^2 + t1*t2*exp(-2*h*(t1 + t2))*(l*t2 -2*h*t1 -2*h*t2+l*t1 - t2/t1 - 2))/(t1 + t2)^2/t1^2/t2 + ...
//        (l-1/(t1 + t2) - 1/t1 - t1/t2/(t1+t2)*exp(-l*(t1 + t2))+(t1 + t2)/t1/t2*exp(-l*t1))/(t1 + t2)/t1/d;

			norm = 	(l-1/t1_p_t2 - 1/t1 - t1/t2/t1_p_t2*exp(-l*t1_p_t2)+t1_p_t2/t1/t2*exp(-l*t1))/t1_p_t2/t1/d;

			if (2*h < l)
				norm += (exp(-l*t1 - 2*h*t2)*t1_p_t2_sq - exp(-l*t1_p_t2)*t1*t1 + t1*t2*exp(-2*h*t1_p_t2)*(l*t2 -2*h*t1 -2*h*t2+l*t1 - t2/t1 - 2))/(t1_p_t2*t1_p_t2*t1*t1*t2);

			last_l = l;
		}
		LogL += log(f_frag/norm)/freq_l;
	}
	return -2*LogL;
}


double FLL_ModelBD(const dataVec & param, const dataType & data)
{

	double d = pow(10,param[0]);
	double h = pow(10,param[1]);
	double t1 = pow(10,param[2]);
	double t2 = pow(10,param[3]);
	double a = param[4];

	const vector<int> & geneIndex = data.geneIndex;
	const dataVec L = data.geneData[0](geneIndex);
	const dataVec & Freq_l = data.geneData[1](geneIndex);

	double last_l = 0;
	double norm=0;
	double LogL = 0;

//	f_frag = a*( (x> h).*(x < l-h)./(t1+ t2) .* (t1.*exp(-2*l*(t1+t2)+(t1+t2)*(l-h+x)) + t2*exp(-l*(t1+t2))) + ...
//         1./(t1+ t2) .* (t1.*exp(-2*l*(t1+t2)+(t1+t2)*(l+x)) + t2*exp(-l*(t1+t2)))/d) + ...
//         (1-a)*( (x> h).*(x < l-h)./(t1+ t2) .* (t1.*exp(-t1*(l-x) - 2*t2*h - t1*h) + t2.*exp(-t1*l-t2*(x+h))) + ...
//         1./(t1+ t2) .* (t1.*exp(-t1*(l-x)) + t2.*exp(-t1*l-t2*(x)))/d);

	double t1_p_t2 = t1+t2;
	double t1_p_t2_sq = t1_p_t2*t1_p_t2;


	for (size_t i = 0;i < data.fragData.size(); i++)
	{
		const double & x = data.fragData[i];
		const double & l = L[i];
		const double & freq_l = Freq_l[i];


		double exp_ml_t1_p_t2 = exp(-l*t1_p_t2);
		double f_frag = 1/t1_p_t2 * (t1*exp(-2*l*t1_p_t2+t1_p_t2*(l+x)) + t2*exp_ml_t1_p_t2)/d;
		double f_frag_pt2 = (t1*exp(-t1*(l-x)) + t2*exp(-t1*l-t2*(x)))/t1_p_t2/d;

		if ((x > h) && (x < (l-h)))
		{
			f_frag += (t1 * exp(-2*l*t1_p_t2+t1_p_t2*(l-h+x)) + t2*exp_ml_t1_p_t2)/t1_p_t2;
			f_frag_pt2 += (t1*exp(-t1*(l-x) - 2*t2*h - t1*h) + t2*exp(-t1*l-t2*(x+h)))/t1_p_t2;
		}

		f_frag = a * f_frag + (1-a)*f_frag_pt2;

		if (f_frag == 0)
			return 1E20;

		if (l != last_l)
		{

//		    norm = a*( (2*h<l).*(t1.*(exp(-2*h*(t1+t2))-exp(-l*(t1+t2)))+t2*(t1+t2).*(l-2*h).*exp(-l*(t1+t2)))/(t1 + t2)^2 + ...
//        (exp(-l.*(t1+t2)).*(l.*t2^2+l.*t2*t1-t1)+t1)/(t1 + t2)^2/d) + ...
//        (1-a)*( (2*h<l).*(exp(-2*h*(t1 + t2)) - exp(-l*(t1 + t2)))/(t1 + t2) + ...
//        (1-exp(-l.*(t1 + t2)))/(t1 + t2)/d);

			norm = (exp_ml_t1_p_t2*(l*t2*t2+l*t2*t1-t1)+t1)/t1_p_t2_sq/d;
			double norm_pt2 = (1-exp_ml_t1_p_t2)/t1_p_t2/d;

			if (2*h<l)
			{
				double exp_m2_h_t1_p_t2 = exp(-2*h*t1_p_t2);
				norm += (t1*(exp_m2_h_t1_p_t2-exp_ml_t1_p_t2)+t2*t1_p_t2*(l-2*h)*exp_ml_t1_p_t2)/t1_p_t2_sq;
				norm_pt2 += (exp_m2_h_t1_p_t2 - exp_ml_t1_p_t2)/t1_p_t2;
			}

			norm = a*norm +(1-a) * norm_pt2;

		}

		LogL += log(f_frag/norm)/freq_l;
	}

	return -2*LogL;
}


void getBias(modelType m,dataVec & params,const dataVec & l, dataVec & bias)
{
	double d = pow(10, params[0]);
	double h = pow(10, params[1]);
	double t1, t2, a;
	switch (m)
	{
	case noModel:	//This should not happen
	case ModelA:
		break;
	case ModelB:
	case ModelD:
	case ModelE:
		t1 = pow(10, params[2]);
		t2 = pow(10, params[3]);
		break;
	case ModelC:
		t1 = 0;
		t2 = pow(10, params[2]);
		break;
	case ModelBD:
		t1 = pow(10, params[2]);
		t2 = pow(10, params[3]);
		a = params[4];
		break;
	}

	bias.resize(l.size());

	switch (m)
	{
	case noModel:
		break;
	case ModelA:
		bias = (2 * h < l)*(l - 2 * h) + l / d;
		break;
	case ModelB:
		bias = ((2 * h < l)*(t1*(exp(-2 * h*(t1 + t2)) - exp(-l*(t1 + t2))) + t2*(t1 + t2)*(l - 2 * h)*exp(-l*(t1 + t2))) / ((t1 + t2)*(t1 + t2)) +
			exp(-l*(t1 + t2))*(l*t2*t2 + t1*(exp(l*(t1 + t2)) + l*t2 - 1)) / ((t1 + t2)*(t1 + t2)) / d);
		break;
	case ModelC:
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
		bias = (2 * h < l)*(exp(-2 * h*t2) - exp_ml_t2) / t2 + (1 - exp_ml_t2) / t2 / d;
	}

	break;
	case ModelD:
		for (size_t i = 0; i < l.size(); i++)
		{
			if (2 * h < l[i])
				bias[i] = (exp(-2 * h*(t1 + t2)) - exp(-l[i] * (t1 + t2))) / (t1 + t2) + (1 - exp(-l[i] * (t1 + t2))) / (t1 + t2) / d;
			else
				bias[i] = (1 - exp(-l[i] * (t1 + t2))) / (t1 + t2) / d;
		}
		break;
	case ModelE:
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

		bias = (2 * h < l)*(exp(-l*t1 - 2 * h*t2)*(t1 + t2)*(t1 + t2) - exp(-l*(t1 + t2))*t1*t1 + t1*t2*exp(-2 * h*(t1 + t2))*(l*t2 - 2 * h*t1 - 2 * h*t2 + l*t1 - t2 / t1 - 2)) / ((t1 + t2) * (t1 + t2)) / (t1 * t1) / t2 +
			(l - 1 / (t1 + t2) - 1 / t1 - t1 / t2 / (t1 + t2)*exp(-l*(t1 + t2)) + (t1 + t2) / t1 / t2*exp(-l*t1)) / (t1 + t2) / t1 / d;

		bias /= t1;

		break;
	case ModelBD:
	{
		bias = a*((2 * h < l)*(t1*(exp(-2 * h*(t1 + t2)) - exp(-l*(t1 + t2))) + t2*(t1 + t2)*(l - 2 * h)*exp(-l*(t1 + t2))) / ((t1 + t2) * (t1 + t2)) +
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
	bias = bias * l[0] / bias[0];
	bias /= l;
}
