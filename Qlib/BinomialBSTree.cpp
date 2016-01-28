/*
 * BinomialBSTree.cpp
 *
 *  Created on: Oct 12, 2014
 *      Author: ericmac
 */

#include "BinomialBSTree.h"
#include "BS.h"

BinomialBSTree::BinomialBSTree() {
	// TODO Auto-generated constructor stub

}

BinomialBSTree::~BinomialBSTree() {
	// TODO Auto-generated destructor stub
}

void BinomialBSTree::load()
{

	t=option.T/N;
	 u=exp(option.sigma*sqrt(t));
	 d=1/u;
	 p=(exp((option.r-option.q)*t)-d)/(u-d);
	 q = 1-p;
	 disc=exp(-option.r*t);


	if (type==0)//Euro
	{

			for (long i=0; i<=N-1; i++)
			{
				basic_option opt=option;
				opt.S*=pow(u,N-1-i)*pow(d,i);
				opt.T=t;
				BS bs=BS(opt);
				option_value(i,N-1)=bs.price();
			}
			for (long j=N-2;j>=0;j--)
				for (long i=0;i<=j;i++)
					option_value(i,j)=disc*(p*option_value(i,j+1)+ q*option_value(i+1,j+1));

	}else //American
	{
		if (option.type==0)// american call
		{
			for (long i=0; i<=N-1; i++)
			{
				basic_option opt=option;
				opt.S*=pow(u,N-1-i)*pow(d,i);
				opt.T=t;
				BS bs =BS(opt);
				option_value(i,N-1)=std::max(option.S*pow(u,N-1-i)*pow(d,i)-option.K, bs.price());
			}
			for (long j=N-2;j>=0;j--)
				for (long i=0;i<=j;i++)
					option_value(i,j)=std::max(disc*(p*option_value(i,j+1)+ q*option_value(i+1,j+1)),option.S*pow(u, j-i)*pow(d,i)-option.K);

		}else //american put
		{
			for (long i=0; i<=N-1; i++)
			{
				basic_option opt=option;
				opt.S*=pow(u,N-1-i)*pow(d,i);
				opt.T=t;
				BS bs =BS(opt);
				option_value(i,N-1)=std::max(option.K- option.S*pow(u,N-1-i)*pow(d,i),bs.price());
			}
			for (long j=N-2;j>=0;j--)
				for (long i=0;i<=j;i++)
					option_value(i,j)=std::max(disc*(p*option_value(i,j+1)+ q*option_value(i+1,j+1)),option .K- option.S*pow(u, j-i)*pow(d,i));


		}

	}

}

double BinomialBSTree::average_price()
{
	BinomialBSTree BBST1(option, N-1,type);
	BinomialBSTree BBST2(option, N, type);

	return 0.5*(BBST1.price()+BBST2.price());
}

double BinomialBSTree::average_delta()
{
	BinomialBSTree BBST1(option, N-1,type);
	BinomialBSTree BBST2(option, N, type);
	return 0.5*(BBST1.delta()+BBST2.delta());

}
