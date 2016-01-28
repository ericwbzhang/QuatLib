/*
 * TrinomialTree.cpp
 *
 *  Created on: Oct 25, 2014
 *      Author: ericmac
 */

#include "TrinomialTree.h"
#include "algorithm"

TrinomialTree::TrinomialTree() {
	// TODO Auto-generated constructor stub

}

TrinomialTree::~TrinomialTree() {
	// TODO Auto-generated destructor stub
}

TrinomialTree::TrinomialTree(basic_option opt, long n ,int Euro_Amer) :BinomialTree(opt, n, Euro_Amer)
{
	option_value.resize(2*N+1, N+1);
	option_value.clear();
}

void  TrinomialTree::load()
{
	t=option.T/N;
	u=exp(option.sigma*sqrt(3*t));
	d=1.0/u;
	m=1;
	p=1.0/6.0+(option.r-option.q-0.5*option.sigma*option.sigma)*sqrt(t/(12*option.sigma*option.sigma));
	q=1.0/6.0-(option.r-option.q-0.5*option.sigma*option.sigma)*sqrt(t/(12*option.sigma*option.sigma));
	v=2.0/3.0;
	disc=exp(-option.r*t);


	if(type==0)// Euro
	{
		if(option.type==0)//Euro call
		{
			for(long i=0; i<=2*N;i++)
				option_value(i,N)=std::max(option.S*pow(u,N-i)-option.K,0.0);
			for (long j=N-1;j>=0;j--)
				for (long i=0;i<=2*j;i++)
					option_value(i,j)=disc*(p*option_value (i,j+1)+v*option_value(i+1,j+1)+ q*option_value(i+2,j+1));


		}else //Euro put
		{
			for (long i=0;i<=2*N;i++)
				option_value(i,N)=std::max(option.K-option.S* pow(u,N-i),0.0);
			for (long j=N-1;j>=0; j--)
				for (long i=0;i<=2*j;i++)
					option_value(i,j)=disc*(p*option_value (i,j+1)+v*option_value(i+1,j+1)+ q*option_value(i+2,j+1));
		}

	}else //American
	{
		if(option.type==0)//American call
		{
			for (long i=0;i<=2*N;i++)
				option_value(i,N)=std::max(option.S* pow(u,N-i)-option.K,0.0);
			for (long j=N-1; j>=0;j--)
				for (long i=0;i<=2*j;i++)
					option_value(i,j)=std::max(disc* (p*option_value(i,j+1)+v*option_value(i+1,j+1)+ q*option_value(i+2,j+1)), option.S*pow(u, j-i)-option.K);


		}else //American put
		{

			for (long i=0;i<=2*N;i++)
				option_value(i,N)=std::max(option.K-option.S*pow( u, N-i),0.0);

				for (long j=N-1; j>=0;j--)
				for(long i=0;i<=2*j;i++)
					option_value(i,j)= std::max(disc* (p*option_value(i,j+1)+v*option_value(i+1,j+1)+ q*option_value(i+2,j+1)), option.K-option.S*pow(u, j-i));

		}
	}
}

double TrinomialTree::price()
{return option_value(0,0);}

double TrinomialTree::delta()
{
	return (option_value(0,1)-option_value(2,1))/(option.S*(u-d));

}

double TrinomialTree::gamma()
{
	double a=(option_value(0,2)-option_value(2,2))/(option.S * (u*u-1));
	double b=(option_value(2,2)-option_value(4,2))/(option.S * (1-d*d));
	return  (a-b)/(option.S*(u-d));
}


double TrinomialTree::theta()
{
	return (option_value(1,1)-option_value(0,0))/t;

}

