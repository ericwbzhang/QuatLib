/*
 * MCEuro.cpp
 *
 *  Created on: Oct 31, 2014
 *      Author: Eric
 */

#include "MCEuro.h"
#include "LinearCongRand.h"
#include "inv_cdf_normal.h"
#include "boost/numeric/ublas/io.hpp"

MC_Euro::MC_Euro() {
	// TODO Auto-generated constructor stub

}

MC_Euro::~MC_Euro() {
	// TODO Auto-generated destructor stub
}

MC_Euro::MC_Euro(basic_option opt, long n)
{
	option=opt;
	N=n;

	this->load();
}

void MC_Euro::load()
{
	LinearCongRand unit_rand=LinearCongRand(N);
	boost::numeric::ublas::vector<double> normal=unit_rand.rand();
	for(long i=0; i<normal.size(); i++)
	{
		normal(i)= inv_cdf_normal(normal(i));
	}

	stock_path=normal;

	for (long i=0; i<stock_path.size(); i++)
	{
		stock_path(i)=option.S*exp((option.r-option.q-0.5*option.sigma*option.sigma)*option.T+option.sigma*sqrt(option.T)*normal(i));

	}


	price_path=stock_path;
	for(long i=0; i<price_path.size();i++)
	{
		double c=0;
		if(option.type==0)// call
		{
			c=std::max(stock_path(i)-option.K, 0.0);
		}else //put
		{
			c=std::max(option.K- stock_path(i), 0.0);
		}

		price_path(i)=exp(-option.r*option.T)*c;

	}


	delta_path=stock_path;
	for (long i=0; i<delta_path.size();i++)
	{
		if (option.type==0) //call
		{
			if(stock_path(i)>option.K)
				delta_path(i)=exp(-option.r*option.T)* stock_path(i)/option.S;
			else delta_path(i)=0;
		}else //put
		{
			if(stock_path(i)<option.K)
				delta_path(i)= -exp(-option.r*option.T) * stock_path(i)/option.S;
			else delta_path(i)=0;

		}

	}

	vega_path=stock_path;
	for (long i=0; i<vega_path.size();i++)
	{
		if(option.type==0)//call
		{
			if(stock_path(i)>option.K)
				vega_path(i)=stock_path(i)*exp(-option.r*option.T)* (-option.sigma*option.T+sqrt(option.T)*normal(i));
			else vega_path(i)=0;

		}else //put
		{
			if (stock_path(i)<option.K)
				vega_path(i)= -stock_path(i)*exp(-option.r*option.T)* (-option.sigma*option.T+sqrt(option.T)*normal(i));
			else vega_path(i)=0;
		}

	}
}


double MC_Euro::opt_price()
{
	double mean=0;
	for (long i=0 ; i<price_path.size(); i++)
		mean+=price_path(i);

	mean=mean/price_path.size();
	return mean;
}


double MC_Euro::opt_delta()
{
	double mean=0;
	for (long i=0 ; i<delta_path.size(); i++)
		mean+=delta_path(i);

	mean=mean/delta_path.size();
	return mean;

}

double MC_Euro::opt_vega()
{
	double mean=0;
	for (long i=0 ; i<vega_path.size(); i++)
		mean+=vega_path(i);

	mean=mean/vega_path.size();
	return mean;
}

boost::numeric::ublas::vector<double> MC_Euro::stockpath()
{
	return stock_path;
}




