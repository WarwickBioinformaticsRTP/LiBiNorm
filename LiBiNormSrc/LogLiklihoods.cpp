
#include <algorithm>
#include "LogLiklihoods.h"
using namespace std;


double FLL_ModelA(const dataVec & param, const dataType & data)
{
	double d = pow(10,param[0]);
	double h = pow(10,param[1]);

	const dataVec & x = data.fragData;
	const vector<int> & geneIndex = data.geneIndex;
	const dataVec & L = data.geneData[0];
	const dataVec & freq_l = data.geneData[1];
	dataVec l = L(geneIndex);

//x = data(1, :);
//l = data(2, :);
//freq_l = data(3, :);

//f_frag = (x> h).*(x < l-h) + 1/d;
	dataVec f_frag = ((x> h)*(x < (l-h))) + 1/d;

//norm =  (2*h<l).*(l-2*h) + l/d;
	dataVec norm =  (2*h<l)*(l-(2*h)) + l/d;

//LogL = -2*sum(log(f_frag./norm)./freq_l);
	double LogL = -2*sum(log(f_frag/norm)/freq_l(geneIndex));

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

	const dataVec & x = data.fragData;
	const vector<int> & geneIndex = data.geneIndex;
	const dataVec & L = data.geneData[0];
	const dataVec & freq_l = data.geneData[1];

	//	f_frag = (x> h).*(x < l-h)./(t1+ t2) .* (t1.*exp(-2*l*(t1+t2)+(t1+t2)*(l-h+x)) + t2*exp(-l*(t1+t2))) + ...
	//         1./(t1+ t2) .* (t1.*exp(-2*l*(t1+t2)+(t1+t2)*(l+x)) + t2*exp(-l*(t1+t2)))/d;

	dataVec exp_l_t1_t2 = exp(-L * (t1+t2));

	dataVec t2_exp_l_t1_t2_full = t2*exp_l_t1_t2(geneIndex);
	dataVec l = L(geneIndex);
	dataVec t1_exp_c_l_t1_t2 = t1*exp((x-l)*(t1+t2));

//	dataVec f_frag = ((x> h)*(x < (l-h)) * (t1*exp((x-l-h)*(t1+t2)) + t2*exp_l_t1_t2_full) + (t1*exp((x-l)*(t1+t2)) + t2*exp_l_t1_t2_full)/d)/(t1+t2);
	dataVec f_frag = ((x> h)*(x < (l-h)) * (t1_exp_c_l_t1_t2/exp(h*(t1+t2)) + t2_exp_l_t1_t2_full) + (t1_exp_c_l_t1_t2 + t2_exp_l_t1_t2_full)/d)/(t1+t2);

#ifdef VERIFY_SPEEDUP
	dataVec f_frag_full = (x> h)*(x < (l-h))/(t1+ t2) * (t1*exp(-2 * l *(t1+t2)+(l-h+x)*(t1+t2)) + t2*exp(-l * (t1+t2))) + 
			1/(t1+ t2) * (t1*exp(-2 * l *(t1+t2)+(l+x)*(t1+t2)) + t2*exp(-l * (t1+t2)))/d;
	dataVec verify = f_frag/f_frag_full;
#endif
	if (f_frag.does_not_contain_null())
	{

		//	    norm =  (2*h<l).*(t1.*(exp(-2*h*(t1+t2))-exp(-l*(t1+t2)))+t2*(t1+t2).*(l-2*h).*exp(-l*(t1+t2)))/(t1 + t2)^2 + ...
		//        (exp(-l*(t1+t2)).*(l.*t2^2+l.*t2*t1-t1)+t1)/(t1 + t2)^2/d;

		dataVec norm =  ((((2*h)<L)*(t1*(exp(-2*h*(t1+t2))-exp_l_t1_t2)+t2*(t1+t2)*(L-(2*h))*exp_l_t1_t2)) + 
				  ((exp_l_t1_t2*(L*t2*(t1+t2)-t1)+t1)/d))/((t1+t2)*(t1+t2));

#ifdef VERIFY_SPEEDUP
		dataVec norm_full =  (((2*h)<L)*(t1*(exp(-2*h*(t1+t2))-exp(-L*(t1+t2)))+t2*(t1+t2)*(L-(2*h))*exp(-L*(t1+t2)))/pow(t1 + t2,2)) + 
				  ((exp(-L*(t1+t2))*((L*t2*t2)+(L*t2*t1)-t1)+t1)/pow(t1 + t2,2)/d);
		verify = norm/norm_full;

#endif
		LogL = -2*sum(log(f_frag/norm(geneIndex))/freq_l(geneIndex));
	}
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

	const dataVec & x = data.fragData;
	const vector<int> & geneIndex = data.geneIndex;
	const dataVec & L = data.geneData[0];
	const dataVec & freq_l = data.geneData[1];
	dataVec l = L(geneIndex);

//    f_frag = (x> h).*(x < l-h).*exp(-t2*(x+h)) + ...
//    exp(-t2*(x))/d;

	dataVec f_frag = (x> h)*(x < l-h)*exp(-t2*(x+h)) + exp(-t2*(x))/d;
	
//	if (isempty(find(f_frag == 0, 1)))
	if (f_frag.does_not_contain_null())
	{

//    norm =  (2*h<l).*(exp(-2*h*t2) - exp(-l*( t2)))/t2 + ...
//        (1-exp(-l*( t2)))/t2/d;
		dataVec norm =  (2*h<l)*(exp(-2*h*t2) - exp(-l*( t2)))/t2 + (1-exp(-l*( t2)))/t2/d;

		LogL = -2*sum(log(f_frag/norm)/freq_l(geneIndex));
	}
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

	const dataVec & x = data.fragData;
	const vector<int> & geneIndex = data.geneIndex;
	const dataVec & L = data.geneData[0];
	const dataVec & freq_l = data.geneData[1];
	dataVec l = L(geneIndex);

//f_frag = (x> h).*(x < l-h)./(t1+ t2) .* (t1.*exp(-t1*(l-x) - 2*t2*h - t1*h) + t2.*exp(-t1*l-t2*(x+h))) + ...
//         1./(t1+ t2) .* (t1.*exp(-t1*(l-x)) + t2.*exp(-t1*l-t2*(x)))/d;
	dataVec f_frag = (x> h)*(x < l-h)/(t1+ t2) * (t1*exp(-t1*(l-x) - 2*t2*h - t1*h) + t2*exp(-t1*l-t2*(x+h))) + 
	         1/(t1+ t2) * (t1*exp(-t1*(l-x)) + t2*exp(-t1*l-t2*(x)))/d;

	//if (isempty(find(f_frag == 0, 1)))
	if (f_frag.does_not_contain_null())
	{

	//    norm =  (2*h<l).*(exp(-2*h*(t1 + t2)) - exp(-l*(t1 + t2)))/(t1 + t2) + ...
    //    (1-exp(-l*(t1 + t2)))/(t1 + t2)/d;
		dataVec norm =  (2*h<l)*(exp(-2*h*(t1 + t2)) - exp(-l*(t1 + t2)))/(t1 + t2) + 
			(1-exp(-l*(t1 + t2)))/(t1 + t2)/d;

		LogL = -2*sum(log(f_frag/norm)/freq_l(geneIndex));
	}
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

/*	d = 10^(param(1));
	h = 10^(param(2));
	t1 = 10^(param(3));
	t2 = 10^(param(4));
	%sig = 10^(param(5));
	x = data(1, :);
	l = data(2, :);
	freq_l = data(3, :);*/

	const dataVec & x = data.fragData;
	const vector<int> & geneIndex = data.geneIndex;
	const dataVec & L = data.geneData[0];
	const dataVec & freq_l = data.geneData[1];
	dataVec l = L(geneIndex);


//	f_frag = (x> h).*(x < l-h)./t1/(t1+ t2).*(exp(-2*h*(t1+ t2))-exp(-(h +x)*(t1+t2)) - exp(-t1*h-2*h*t2-(l-x)*t1) + exp(-h*t2-l*t1-x*t2)) + ...
//     1./t1/(t1+ t2) .*(1 - exp(-x*(t1+t2)) - exp(-(l-x)*t1) + exp(-l*t1-x*t2))/d;


	dataVec	f_frag = (x> h)*(x < l-h)/t1/(t1+ t2)*(exp(-2*h*(t1+ t2))-exp(-(x + h)*(t1+t2)) - exp(-t1*h-2*h*t2-(l-x)*t1) + exp(-h*t2-l*t1-x*t2)) + 
		     1/t1/(t1+ t2) *(1 - exp(-x*(t1+t2)) - exp(-(l-x)*t1) + exp(-l*t1-x*t2))/d;

//	if (isempty(find(f_frag == 0, 1)))
	if (f_frag.does_not_contain_null())
	{
//		   norm =  (2*h<l).*(exp(-l*t1 - 2*h*t2)*(t1 + t2)^2 - exp(-l*(t1 + t2))*t1^2 + t1*t2*exp(-2*h*(t1 + t2))*(l*t2 -2*h*t1 -2*h*t2+l*t1 - t2/t1 - 2))/(t1 + t2)^2/t1^2/t2 + ...
//        (l-1/(t1 + t2) - 1/t1 - t1/t2/(t1+t2)*exp(-l*(t1 + t2))+(t1 + t2)/t1/t2*exp(-l*t1))/(t1 + t2)/t1/d;
		dataVec norm =  (2*h<l)*(exp(-l*t1 - 2*h*t2)*(pow((t1 + t2),2)) - exp(-l*(t1 + t2))*pow(t1,2) + t1*t2*exp(-2*h*(t1 + t2))*(l*t2 -2*h*t1 -2*h*t2+l*t1 - t2/t1 - 2))/(pow((t1 + t2)*t1,2)*t2) + 
	       (l-1/(t1 + t2) - 1/t1 - t1/t2/(t1+t2)*exp(-l*(t1 + t2))+(t1 + t2)/t1/t2*exp(-l*t1))/(t1 + t2)/t1/d;
	
		LogL = -2*sum(log(f_frag/norm)/freq_l(geneIndex));
	}
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

	const dataVec & x = data.fragData;
	const vector<int> & geneIndex = data.geneIndex;
	const dataVec & L = data.geneData[0];
	const dataVec & freq_l = data.geneData[1];
	dataVec l = L(geneIndex);

//	f_frag = a*( (x> h).*(x < l-h)./(t1+ t2) .* (t1.*exp(-2*l*(t1+t2)+(t1+t2)*(l-h+x)) + t2*exp(-l*(t1+t2))) + ...
//         1./(t1+ t2) .* (t1.*exp(-2*l*(t1+t2)+(t1+t2)*(l+x)) + t2*exp(-l*(t1+t2)))/d) + ...
//         (1-a)*( (x> h).*(x < l-h)./(t1+ t2) .* (t1.*exp(-t1*(l-x) - 2*t2*h - t1*h) + t2.*exp(-t1*l-t2*(x+h))) + ...
//         1./(t1+ t2) .* (t1.*exp(-t1*(l-x)) + t2.*exp(-t1*l-t2*(x)))/d);


//	f_frag = a*( (x> h).*(x < l-h)./(t1+ t2) .* (t1.*exp(-2*l*(t1+t2)+(t1+t2)*(l-h+x)) + t2*exp(-l*(t1+t2))) + ...
//         1./(t1+ t2) .* (t1.*exp(-2*l*(t1+t2)+(t1+t2)*(l+x)) + t2*exp(-l*(t1+t2)))/d) + ...
//         (1-a)*( (x> h).*(x < l-h)./(t1+ t2) .* (t1.*exp(-t1*(l-x) - 2*t2*h - t1*h) + t2.*exp(-t1*l-t2*(x+h))) + ...
//         1./(t1+ t2) .* (t1.*exp(-t1*(l-x)) + t2.*exp(-t1*l-t2*(x)))/d);


	dataVec f_frag = a * ( (x > h)*(x < (l-h))/(t1+ t2) * (t1 * exp(-2*l*(t1+t2)+(t1+t2)*(l-h+x)) + t2*exp(-l*(t1+t2))) +
				1/(t1+ t2) * (t1*exp(-2*l*(t1+t2)+(t1+t2)*(l+x)) + t2*exp(-l*(t1+t2)))/d) +
				(1-a)*( (x> h)*(x < l-h)/(t1+ t2) * (t1*exp(-t1*(l-x) - 2*t2*h - t1*h) + t2*exp(-t1*l-t2*(x+h))) +
				1/(t1+ t2) * (t1*exp(-t1*(l-x)) + t2*exp(-t1*l-t2*(x)))/d);

	if (f_frag.does_not_contain_null())
	{

//		    norm = a*( (2*h<l).*(t1.*(exp(-2*h*(t1+t2))-exp(-l*(t1+t2)))+t2*(t1+t2).*(l-2*h).*exp(-l*(t1+t2)))/(t1 + t2)^2 + ...
//        (exp(-l.*(t1+t2)).*(l.*t2^2+l.*t2*t1-t1)+t1)/(t1 + t2)^2/d) + ...
//        (1-a)*( (2*h<l).*(exp(-2*h*(t1 + t2)) - exp(-l*(t1 + t2)))/(t1 + t2) + ...
//        (1-exp(-l.*(t1 + t2)))/(t1 + t2)/d);

//		    norm = a*( (2*h<l).*(t1.*(exp(-2*h*(t1+t2))-exp(-l*(t1+t2)))+t2*(t1+t2).*(l-2*h).*exp(-l*(t1+t2)))/(t1 + t2)^2 + ...
//        (exp(-l.*(t1+t2)).*(l.*t2^2+l.*t2*t1-t1)+t1)/(t1 + t2)^2/d) + ...
//        (1-a)*( (2*h<l).*(exp(-2*h*(t1 + t2)) - exp(-l*(t1 + t2)))/(t1 + t2) + ...
//        (1-exp(-l.*(t1 + t2)))/(t1 + t2)/d);

		dataVec norm = a*( (2*h<l)*(t1*(exp(-2*h*(t1+t2))-exp(-l*(t1+t2)))+t2*(t1+t2)*(l-2*h)*exp(-l*(t1+t2)))/pow((t1 + t2),2) + 
			(exp(-l*(t1+t2))*(l*pow(t2,2)+l*t2*t1-t1)+t1)/pow((t1 + t2),2)/d) + 
		(1-a)*( (2*h<l)*(exp(-2*h*(t1 + t2)) - exp(-l*(t1 + t2)))/(t1 + t2) + 
		(1-exp(-l*(t1 + t2)))/(t1 + t2)/d);

		LogL = -2*sum(log(f_frag/norm)/freq_l(geneIndex));
	}
	return LogL;
}

