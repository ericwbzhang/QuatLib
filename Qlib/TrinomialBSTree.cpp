/*
 * TrinomialBSTree.cpp
 *
 *  Created on: Oct 25, 2014
 *      Author: ericmac
 */

#include "TrinomialBSTree.h"

TrinomialBSTree::TrinomialBSTree() {
	// TODO Auto-generated constructor stub

}

TrinomialBSTree::~TrinomialBSTree() {
	// TODO Auto-generated destructor stub
}


void TrinomialBSTree::load()
{
	t=option.T/N;
	u=exp(option.sigma*sqrt(3*t));
	d=1.0/u;
	m=1;
	p=1.0/6.0+(option.r-option.q-0.5*option.sigma*option.sigma)*sqrt(t/(12*option.sigma*option.sigma));
	q=1.0/6.0-(option.r-option.q-0.5*option.sigma*option.sigma)*sqrt(t/(12*option.sigma*option.sigma));
	v=2.0/3.0;
	disc=exp(-option.r*t);

	if (type==0)//Euro
	{
		for(long i=0 ; i<=2*N-2;i++)
		{
			basic_option opt=option;
			opt.S*=pow(u,N-1-i);
			opt.T=t;
			BS bs=BS(opt);
			option_value(i,N-1)=bs.price();
		}
		for (long j=N-2; j>=0;j--)
			for (long i=0;i<=2*j;i++)
				option_value(i,j)=disc*(p*option_value (i,j+1)+v*option_value(i+1,j+1)+ q*option_value(i+2,j+1));
	}else //American
	{
		if(option.type==0)// American Call
		{
			for (long i=0;i<=2*N-2;i++)
			{
				basic_option opt=option;
				opt.S=option.S*pow(u,N-1-i);
				opt.T=t;
				BS bs=BS(opt);

				option_value(i,N-1)=std::max(bs.price(), opt.S-opt.K);
			}
			for (long j=N-2; j>=0;j--)
				for (long i=0;i<=2*j;i++)
					option_value(i,j)=std::max(disc*(p*option_value (i,j+1)+v*option_value(i+1,j+1)+ q*option_value(i+2,j+1)), option.S*pow(u,j-i)-option.K);

		}else // american put
		{
			for (long i=0;i<=2*N-2;i++)
			{
				basic_option opt=option;
				opt.S=option.S*pow(u,N-1-i);
				opt.T=t;
				BS bs=BS(opt);

				option_value(i,N-1)=std::max(bs.price(), opt.K-opt.S);
			}
			for (long j=N-2; j>=0;j--)
				for (long i=0;i<=2*j;i++)
					option_value(i,j)=std::max(disc*(p*option_value (i,j+1)+v*option_value(i+1,j+1)+ q*option_value(i+2,j+1)), option.K- option.S*pow(u,j-i));

		}



	}
}
