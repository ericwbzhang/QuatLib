/*
 * simple_MCbasket_option.cpp
 *
 *  Created on: Nov 8, 2014
 *      Author: Eric
 */


#include "simple_MCbasket_option.h"
#include "LinearCongRand.h"
#include "BoxMuller.h"

double simple_MCbasket_option(simple_basket_option opt, long n)
{
	LinearCongRand u_generator= LinearCongRand(n);
	boost::numeric::ublas::vector<double> u= u_generator.rand();
	boost::numeric::ublas::vector<double> v;
	v.resize(u.size()/2);
	v.clear();

	long i=0;
	while (i<u.size())
	{
		BoxMuller bm= BoxMuller(u(i), u(i+1));
		u(i)= opt.S1*exp((opt.r-opt.q1-0.5*opt.sigma1*opt.sigma1)* opt.T+ opt.sigma1*sqrt(opt.T)* bm.z1);
		u(i+1)= opt.S2* exp( (opt.r- opt.q2-0.5*opt.sigma2*opt.sigma2)* opt.T+ opt.sigma2*sqrt(opt.T)*(opt.corr*bm.z1+sqrt(1-opt.corr*opt.corr)*bm.z2));
		v(i/2)= exp(-opt.r*opt.T)*std::max(u(i)+u(i+1)-opt.K,0.0);

		i=i+2;
	}

	if(v.size()!= u.size()/2) std::cout<< "Error in MC Basket Option!\n\n";

	i=0;
	double result=0;
	while (i<v.size())
	{
		result+= v(i);
		i++;
	}

	result=result/v.size();
	//cout<<v.size()<<std::endl;

	return result;
}

