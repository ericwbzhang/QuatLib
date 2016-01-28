/*
 * BS.cpp
 *
 *  Created on: Oct 11, 2014
 *      Author: ericmac
 */

#include "BS.h"
#include "normal_distribution.h"


BS::BS() {
	// TODO Auto-generated constructor stub

}

BS::~BS() {
	// TODO Auto-generated destructor stub
}

BS::BS(basic_option o)
{
	opt=o;
	opt_d1=(log(opt.S)-log(opt.K)+(opt.r-opt.q+0.5*opt.sigma*opt.sigma)*opt.T)/(opt.sigma*sqrt(opt.T));
	opt_d2=opt_d1- opt.sigma*sqrt(opt.T);
	if(opt.type==0)
	{
		opt_price= opt.S*exp(-opt.q*opt.T)*Normalcdf(opt_d1)-opt.K*exp(-opt.r*opt.T)*Normalcdf(opt_d2);

	}
	else
	{
		opt_price= -opt.S*exp(-opt.q*opt.T)*Normalcdf(-opt_d1)+opt.K*exp(-opt.r*opt.T)*Normalcdf(-opt_d2);

	}
}

double BS::price()
{
	return opt_price;
}

double BS::d1(){return opt_d1;}

double BS::d2() {return opt_d2;}

double BS::delta()
{
	double delta=0;
	if (opt.type==0)
	{
		delta= exp(-opt.q*opt.T)*Normalcdf(opt_d1);
	}
	else
	{
		delta= -exp(-opt.q*opt.T)*Normalcdf(-opt_d1);

	}
	return delta;
}

double BS::gamma()
{
	double  gamma= exp(-opt.q*opt.T)*Normalpdf(opt_d1)/(opt.S*opt.sigma*sqrt(opt.T));
	return gamma;
}

double BS::omega()
{
	double result=this->delta()*opt.S/opt_price;

	return result;

}

double BS::psi()
{
	double result=0;
	if (opt.type==0)
	{
		result=-Normalpdf(opt_d1)*sqrt(opt.T)/opt.sigma-opt.S*Normalcdf(opt_d1)*exp(-opt.q*opt.T)*opt.T+opt.K*exp(-opt.r*opt.T)*Normalpdf(opt_d2)*sqrt(opt.T)/opt.sigma;
	}
	else{
		result=-Normalpdf(opt_d1)*sqrt(opt.T)/opt.sigma-opt.S*Normalcdf(opt_d1)*exp(-opt.q*opt.T)*opt.T+opt.K*exp(-opt.r*opt.T)*Normalpdf(opt_d2)*sqrt(opt.T)/opt.sigma+opt.S*exp(-opt.q*opt.T)*opt.T;
	}
	return result;

}

double BS::rho()
{
	double result=0;
	if(opt.type==0)
	{
		result=opt.K*opt.T*exp(-opt.r*opt.T)*Normalcdf(opt_d2);

	}
	else{
		result=-opt.K*opt.T*exp(-opt.r*opt.T)*Normalcdf(-opt_d2);
	}
	return result;

}

double BS::theta()
{
	double result=0;
	if (opt.type==0)
	{
		result= -exp(-opt.q*opt.T)*opt.S*Normalpdf(opt_d1)*opt.sigma/(2*sqrt(opt.T))-opt.r*opt.K*exp(-opt.r*opt.T)*Normalcdf(opt_d2)+opt.q*opt.S*exp(-opt.q*opt.T)*Normalcdf(opt_d1);

	}
	else
	{
		result=-exp(-opt.q*opt.T)*opt.S*Normalpdf(opt_d1)*opt.sigma/(2*sqrt(opt.T))+opt.r*opt.K*exp(-opt.r*opt.T)*Normalcdf(-opt_d2)-opt.q*opt.S*exp(-opt.q*opt.T)*Normalcdf(-opt_d1);
	}
	return result;

}


double BS::vega()
{
	double result=opt.S*exp(-opt.q*opt.T)*Normalpdf(opt_d1)*sqrt(opt.T);
	return result;
}











