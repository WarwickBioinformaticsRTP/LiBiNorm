
#include <algorithm>
#include "LogLiklihoods.h"
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

This makes extensive use of the dataVec class which implements the required vector functions in a way that mimics matlab, 
allowing the original matlab equations to be used almost unaltered.

One modification arises from teh fact that the normalisation values are the same for all reads in the gene so only need
to be calculated per gene.  The (x) method uses the geneIndex to expand the one per gene vector to a one per read vector so that the normalisation 
can be performed for each read.

*/


double FLL_ModelA(const dataVec & param, const dataType & data)
{
	double d = pow(10,param[0]);
	double h = pow(10,param[1]);
	double LogL = 1E20;

	const vector<int> & geneIndex = data.geneIndex;
	const dataVec L = data.geneData[0](geneIndex);
	const dataVec & Freq_l = data.geneData[1](geneIndex);

	double last_l = 0;
	double norm=0;
	double temp_l = 0;

	for (size_t i = 0;i < data.fragData.size(); i++)
	{
		const double & x = data.fragData[i];
		const double & l = L[i];
		const double & freq_l = Freq_l[i];


//x = data(1, :);
//l = data(2, :);
//freq_l = data(3, :);

//f_frag = (x> h).*(x < l-h) + 1/d;
//	dataVec f_frag = ((x> h)*(x < (l-h))) + 1/d;
		double f_frag = 1.0/d;
		if ((x> h) && (x < (l-h)))
			f_frag += 1.0;


		if (f_frag == 0)
				return 1E20;

		if (l != last_l)
		{


//norm =  (2*h<l).*(l-2*h) + l/d;
//	dataVec norm =  (2*h<L)*(L-(2*h)) + L/d;

			norm =  l/d;
			if (2*h<l)
				norm +=  (l-(2*h));
			last_l = l;

//			LogL = -2*sum(log(f_frag/norm(geneIndex))/Freq_l(geneIndex));
		}

		temp_l += log(f_frag/norm)/freq_l;
	}

	LogL = -2*temp_l;


//LogL = -2*sum(log(f_frag./norm)./freq_l);
//	double LogL = -2*sum(log(f_frag/norm(geneIndex))/freq_l(geneIndex));

	return LogL;
}


double FLL_ModelB(const dataVec & param, const dataType & data)
{

	double d = pow(10,param[0]);
	double h = pow(10,param[1]);
	double t1 = pow(10,param[2]);
	double t2 = pow(10,param[3]);
//	double sig = pow(10,param[4]);
	double LogL = 1E20;

	const vector<int> & geneIndex = data.geneIndex;
	const dataVec L = data.geneData[0](geneIndex);
	const dataVec & Freq_l = data.geneData[1](geneIndex);

	double last_l = 0;
	double norm=0;
	double temp_l = 0;

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

	//	dataVec f_frag = ((x> h)*(x < (l-h)) * (t1*exp((x-l-h)*(t1+t2)) + t2*exp_l_t1_t2_full) + (t1*exp((x-l)*(t1+t2)) + t2*exp_l_t1_t2_full)/d)/(t1+t2);
//		dataVec f_frag = ((x> h)*(x < (l-h)) * (t1_exp_c_l_t1_t2/exp(h*t1_p_t2) + t2_exp_l_t1_t2_full) + (t1_exp_c_l_t1_t2 + t2_exp_l_t1_t2_full)/d)/t1_p_t2;

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

//			dataVec norm =  ((((2*h)<L)*(t1*(exp(-2*h*t1_p_t2)-exp_L_t1_t2)+t2*t1_p_t2*(L-(2*h))*exp_L_t1_t2)) + 
//					  ((exp_L_t1_t2*(L*t2*t1_p_t2-t1)+t1)/d))/(t1_p_t2*t1_p_t2);

			norm = ((exp_l_t1_t2*(l*t2*t1_p_t2-t1)+t1)/d);
			if ((2*h)<l)
				norm += (t1*(exp(-2*h*t1_p_t2)-exp_l_t1_t2)+t2*t1_p_t2*(l-(2*h))*exp_l_t1_t2);
			norm /= (t1_p_t2*t1_p_t2);

//		LogL = -2*sum(log(f_frag/norm(geneIndex))/freq_l(geneIndex));
			last_l = l;


//			LogL = -2*sum(log(f_frag/norm(geneIndex))/Freq_l(geneIndex));
		}

		temp_l += log(f_frag/norm)/freq_l;
	}

	LogL = -2*temp_l;
	return LogL;
}

double FLL_ModelC(const dataVec & param, const dataType & data)
{
	//	function [LogL] = FLL_Deng(param, data)
/*
d = 10^(param(1));
h = 10^(param(2));
%t1 = 10^(param(3));
t2 = 10^(param(3));
%sig = 10^(param(5));
x = data(1, :);
l = data(2, :);
freq_l = data(3, :);
*/
	double d = pow(10,param[0]);
	double h = pow(10,param[1]);
	double t2 = pow(10,param[2]);

	double LogL = 1E20;

	const vector<int> & geneIndex = data.geneIndex;
	const dataVec L = data.geneData[0](geneIndex);
	const dataVec & Freq_l = data.geneData[1](geneIndex);

	double last_l = 0;
	double norm=0;
	double temp_l = 0;

//    f_frag = (x> h).*(x < l-h).*exp(-t2*(x+h)) + ...
//    exp(-t2*(x))/d;

//	dataVec f_frag = (x> h)*(x < l-h)*exp(-t2*(x+h)) + exp(-t2*(x))/d;


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


//			LogL = -2*sum(log(f_frag/norm(geneIndex))/Freq_l(geneIndex));
		}

		temp_l += log(f_frag/norm)/freq_l;
	}

	LogL = -2*temp_l;
	return LogL;
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
	double t1_p_t2_sq = t1_p_t2*t1_p_t2;


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
	double LogL = 1E20;

	const vector<int> & geneIndex = data.geneIndex;
	const dataVec L = data.geneData[0](geneIndex);
	const dataVec & Freq_l = data.geneData[1](geneIndex);

	double last_l = 0;
	double norm=0;
	double temp_l = 0;


/*	d = 10^(param(1));
	h = 10^(param(2));
	t1 = 10^(param(3));
	t2 = 10^(param(4));
	%sig = 10^(param(5));
	x = data(1, :);
	l = data(2, :);
	freq_l = data(3, :);*/

/*	const dataVec & x = data.fragData;
	const vector<int> & geneIndex = data.geneIndex;
	const dataVec & L = data.geneData[0];
	const dataVec & freq_l = data.geneData[1];
	dataVec l = L(geneIndex);
*/

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

//	dataVec	f_frag = (x> h)*(x < l-h)/t1/t1_p_t2*(exp_m2_h_t1_p_t2-exp(-(x + h)*t1_p_t2) - exp(-t1*h-2*h*t2-(l-x)*t1) + exp(-h*t2-l*t1-x*t2)) + 
//		     1/t1/t1_p_t2 *(1 - exp(-x*t1_p_t2) - exp(-(l-x)*t1) + exp(-l*t1-x*t2))/d;

		double	f_frag = 1/t1/t1_p_t2 *(1 - exp(-x*t1_p_t2) - exp(-(l-x)*t1) + exp(-l*t1-x*t2))/d;

		if ((x> h) && (x < l-h))
			f_frag += (exp_m2_h_t1_p_t2-exp(-(x + h)*t1_p_t2) - exp(-t1*h-2*h*t2-(l-x)*t1) + exp(-h*t2-l*t1-x*t2))/t1/t1_p_t2;

		if (f_frag == 0)
			return 1E20;



		if (l != last_l)
		{

//		   norm =  (2*h<l).*(exp(-l*t1 - 2*h*t2)*(t1 + t2)^2 - exp(-l*(t1 + t2))*t1^2 + t1*t2*exp(-2*h*(t1 + t2))*(l*t2 -2*h*t1 -2*h*t2+l*t1 - t2/t1 - 2))/(t1 + t2)^2/t1^2/t2 + ...
//        (l-1/(t1 + t2) - 1/t1 - t1/t2/(t1+t2)*exp(-l*(t1 + t2))+(t1 + t2)/t1/t2*exp(-l*t1))/(t1 + t2)/t1/d;

//		dataVec norm =  (2*h<L)*(exp(-L*t1 - 2*h*t2)*t1_p_t2_sq - exp(-L*t1_p_t2)*t1*t1 + t1*t2*exp(-2*h*t1_p_t2)*(L*t2 -2*h*t1 -2*h*t2+L*t1 - t2/t1 - 2))/(t1_p_t2*t1_p_t2*t1*t1*t2) + 
//	       (L-1/t1_p_t2 - 1/t1 - t1/t2/t1_p_t2*exp(-L*t1_p_t2)+t1_p_t2/t1/t2*exp(-L*t1))/t1_p_t2/t1/d;


			norm = 	(l-1/t1_p_t2 - 1/t1 - t1/t2/t1_p_t2*exp(-l*t1_p_t2)+t1_p_t2/t1/t2*exp(-l*t1))/t1_p_t2/t1/d;

			if (2*h < l)
				norm += (exp(-l*t1 - 2*h*t2)*t1_p_t2_sq - exp(-l*t1_p_t2)*t1*t1 + t1*t2*exp(-2*h*t1_p_t2)*(l*t2 -2*h*t1 -2*h*t2+l*t1 - t2/t1 - 2))/(t1_p_t2*t1_p_t2*t1*t1*t2);

			last_l = l;
		}

//		LogL = -2*sum(log(f_frag/norm(geneIndex))/freq_l(geneIndex));

		temp_l += log(f_frag/norm)/freq_l;
	}

	LogL = -2*temp_l;
	return LogL;
}


double FLL_ModelBD(const dataVec & param, const dataType & data)
{

	double d = pow(10,param[0]);
	double h = pow(10,param[1]);
	double t1 = pow(10,param[2]);
	double t2 = pow(10,param[3]);
	double a = param[4];
	double LogL = 1E20;

	const vector<int> & geneIndex = data.geneIndex;
	const dataVec L = data.geneData[0](geneIndex);
	const dataVec & Freq_l = data.geneData[1](geneIndex);

	double last_l = 0;
	double norm=0;
	double temp_l = 0;

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

/*		dataVec f_fragX = a * ( (x > h)*(x < (l-h))/t1_p_t2 * (t1 * exp(-2*l*t1_p_t2+t1_p_t2*(l-h+x)) + t2*exp_mL_t1_p_t2(geneIndex)) +
			1/t1_p_t2 * (t1*exp(-2*l*t1_p_t2+t1_p_t2*(l+x)) + t2*exp_mL_t1_p_t2(geneIndex))/d) +
			(1-a)*( (x> h)*(x < l-h)/t1_p_t2 * (t1*exp(-t1*(l-x) - 2*t2*h - t1*h) + t2*exp(-t1*l-t2*(x+h))) +
			1/t1_p_t2 * (t1*exp(-t1*(l-x)) + t2*exp(-t1*l-t2*(x)))/d);
*/
		double f_frag = 1/t1_p_t2 * (t1*exp(-2*l*t1_p_t2+t1_p_t2*(l+x)) + t2*exp_ml_t1_p_t2)/d;

		if ((x > h) && (x < (l-h)))
			f_frag += (t1 * exp(-2*l*t1_p_t2+t1_p_t2*(l-h+x)) + t2*exp_ml_t1_p_t2)/t1_p_t2;

		double f_frag_pt2 = (t1*exp(-t1*(l-x)) + t2*exp(-t1*l-t2*(x)))/t1_p_t2/d;

		if ((x > h) && (x < (l-h)))
			f_frag_pt2 += (t1*exp(-t1*(l-x) - 2*t2*h - t1*h) + t2*exp(-t1*l-t2*(x+h)))/t1_p_t2;

		f_frag = a * f_frag + (1-a)*f_frag_pt2;

		if (f_frag == 0)
			return 1E20;

		if (l != last_l)
		{

//		    norm = a*( (2*h<l).*(t1.*(exp(-2*h*(t1+t2))-exp(-l*(t1+t2)))+t2*(t1+t2).*(l-2*h).*exp(-l*(t1+t2)))/(t1 + t2)^2 + ...
//        (exp(-l.*(t1+t2)).*(l.*t2^2+l.*t2*t1-t1)+t1)/(t1 + t2)^2/d) + ...
//        (1-a)*( (2*h<l).*(exp(-2*h*(t1 + t2)) - exp(-l*(t1 + t2)))/(t1 + t2) + ...
//        (1-exp(-l.*(t1 + t2)))/(t1 + t2)/d);

			double exp_m2_h_t1_p_t2 = exp(-2*h*t1_p_t2);
/*			dataVec normX = a*( (2*h<L)*(t1*(exp_m2_h_t1_p_t2-exp_mL_t1_p_t2)+t2*t1_p_t2*(L-2*h)*exp_mL_t1_p_t2)/t1_p_t2_sq + 
				(exp_mL_t1_p_t2*(L*t2*t2+L*t2*t1-t1)+t1)/t1_p_t2_sq/d) + 
				(1-a)*( (2*h<L)*(exp_m2_h_t1_p_t2 - exp_mL_t1_p_t2)/t1_p_t2 + 
				(1-exp_mL_t1_p_t2)/t1_p_t2/d);
*/
			norm = (exp_ml_t1_p_t2*(l*t2*t2+l*t2*t1-t1)+t1)/t1_p_t2_sq/d;
			if (2*h<l)
				norm += (t1*(exp_m2_h_t1_p_t2-exp_ml_t1_p_t2)+t2*t1_p_t2*(l-2*h)*exp_ml_t1_p_t2)/t1_p_t2_sq;

			double norm_pt2 = (1-exp_ml_t1_p_t2)/t1_p_t2/d;
			if (2*h<l)
				norm_pt2 += (exp_m2_h_t1_p_t2 - exp_ml_t1_p_t2)/t1_p_t2;

			norm = a*norm +(1-a) * norm_pt2;

		}

//		LogL = -2*sum(log(f_frag/norm(geneIndex))/freq_l(geneIndex));
		temp_l += log(f_frag/norm)/freq_l;
	}

	LogL = -2*temp_l;
	return LogL;
}

