/*
 * DiscreteDivdBT.cpp
 *
 *  Created on: Nov 27, 2014
 *      Author: Eric
 */

#include "DiscreteDivdBT.h"

DiscreteDivdBT_Amer::DiscreteDivdBT_Amer() {
	// TODO Auto-generated constructor stub

}

DiscreteDivdBT_Amer::~DiscreteDivdBT_Amer() {
	// TODO Auto-generated destructor stub
}

void DiscreteDivdBT_Amer::load(){

	t=option.T/N;
	u=exp(option.sigma*sqrt(t));
	d=1/u;
	p=(exp((option.r-option.q)*t)-d)/(u-d);
	q = 1-p;
	disc=exp(-option.r*t);
	for(std::vector<dividend>::iterator it=divd_sche.begin(); it!= divd_sche.end();it++)
	{
		if((*it).divd_type==0)// const divd
			option.S-=(*it).divd*exp(-(*it).t*option.r);
		else if ((*it).divd_type==1)	// proportional divd
			option.S*=(1-(*it).divd);

		(*it).t=round((*it).t/t);

	}
	std::vector<dividend>::reverse_iterator re_it= divd_sche.rbegin();
	for(long i=0; i<=N; i++)
	{
		option_value(i,N)=option.S* pow(u,N-i)* pow(d,i);
	}
	for (long j=N-1;j>=0; j--)
		for(long i=0;i<=j;i++)
		{
			if(j+1==(*re_it).t)
			{
				if((*re_it).divd_type==0)//const divd
				{
					option_value(i,j)=(option_value(i,j+1)+(*re_it).divd)/u;
				}else //prop divd
				{
					option_value(i,j)=(option_value(i,j+1))/(1-(*re_it).divd)/u;
				}
				re_it++;
			}else{
				option_value(i,j)= option_value(i,j+1)/u;
			}

		}



	if (option.type==0) //American call
	{
		for (long i=0 ; i<=N; i++)
		{
			option_value (i,N)=std::max(option_value(i,N)-option.K, 0.0);
		}
		for(long j=N-1; j>=0;j--)
			for(long i=0;i<=j; i++)
			{
				option_value(i,j)= std::max(disc*(p*option_value(i,j+1)+ q*option_value(i+1,j+1)), option_value(i,j)-option.K);
			}
	}else // American put
	{
		for(long i=0;i<=N;i++)
		{
			option_value(i,N)=std::max(option.K- option_value(i,N),0.0);
		}
		for(long j=N-1; j>=0; j--)
			for (long i=0;i<=j; i++)
			{
				option_value(i,j)=std::max(disc*(p*option_value(i,j+1)+ q*option_value(i+1,j+1)), option.K- option_value(i,j));
			}

	}

}


