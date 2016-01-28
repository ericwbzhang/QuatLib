/*
 * BinomialTree.cpp
 *
 *  Created on: Oct 11, 2014
 *      Author: ericmac
 */

#include "BinomialTree.h"
#include <algorithm>

BinomialTree::BinomialTree() {
	// TODO Auto-generated constructor stub

}

BinomialTree::~BinomialTree() {
	// TODO Auto-generated destructor stub
}

BinomialTree::BinomialTree(basic_option opt, long n, int Euro_Amer)
{
	option= opt;
	N=n;
	type=Euro_Amer;
	option_value.resize(N+1, N+1);
	option_value.clear();
}



void BinomialTree::load()
{
	t=option.T/N;
	 u=exp(option.sigma*sqrt(t));
	 d=1/u;
	 p=(exp((option.r-option.q)*t)-d)/(u-d);
	 q = 1-p;
	 disc=exp(-option.r*t);

	if (type==0)//Euro
	{
		if(option.type==0) // euro call
		{
			for (long i=0; i<=N; i++)
				option_value(i,N)=std::max(option.S*pow(u,N-i)*pow(d,i)-option.K,0.0);
			for (long j=N-1;j>=0;j--)
				for (long i=0;i<=j;i++)
					option_value(i,j)=disc*(p*option_value(i,j+1)+ q*option_value(i+1,j+1));

		}else// euro put
		{
			for (long i=0; i<=N; i++)
				option_value(i,N)=std::max(option.K-option.S*pow(u,N-i)*pow(d,i),0.0);
			for (long j=N-1;j>=0;j--)
				for (long i=0;i<=j;i++)
					option_value(i,j)=disc*(p*option_value(i,j+1)+ q*option_value(i+1,j+1));
		}


	}else //American
	{
		if (option.type==0)// american call
		{
			for (long i=0; i<=N; i++)
				option_value(i,N)=std::max(option.S*pow(u,N-i)*pow(d,i)-option.K,0.0);
			for (long j=N-1;j>=0;j--)
				for (long i=0;i<=j;i++)
					option_value(i,j)=std::max(disc*(p*option_value(i,j+1)+ q*option_value(i+1,j+1)),option.S*pow(u, j-i)*pow(d,i)-option.K);

		}else //american put
		{
			for (long i=0; i<=N; i++)
				option_value(i,N)=std::max(option.K- option.S*pow(u,N-i)*pow(d,i),0.0);
			for (long j=N-1;j>=0;j--)
				for (long i=0;i<=j;i++)
					option_value(i,j)=std::max(disc*(p*option_value(i,j+1)+ q*option_value(i+1,j+1)),option .K- option.S*pow(u, j-i)*pow(d,i));


		}

	}
}

double BinomialTree::price()
{return option_value(0,0);}

double BinomialTree::average_price()
{
	BinomialTree BT1(option, N-1, type);
	BinomialTree BT2(option, N, type);

	return 0.5*(BT1.price()+BT2.price());
}

double BinomialTree::average_delta()
{
	BinomialTree BT1(option, N-1, type);
	BinomialTree BT2(option, N, type);


	return 0.5*(BT1.delta()+BT2.delta());

}
double BinomialTree::delta()
{
	double result=0;
	result= (option_value(0,1)-option_value(1,1))/ (option.S*(u-d));
	return result;

}

double BinomialTree::gamma()
{
	double result;
	double a,b,c;
	a=(option_value(0,2)-option_value(1,2))/ (option.S*(u*u-u*d));
	b=(option_value(1,2)-option_value(2,2))/ (option.S*(u*d-d*d));
	c=0.5*(option.S*(u*u-d*d));
	result=(a-b)/c;

	return result;
}

double BinomialTree::theta()
{
	double result;
	result= (option_value(1,2)-option_value(0,0))/(2*t);
	return result;

}

