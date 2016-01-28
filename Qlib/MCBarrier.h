/*
 * MCBarrier.h
 *
 *  Created on: Dec 5, 2014
 *      Author: Eric
 */

#ifndef MCBARRIER_H_
#define MCBARRIER_H_
#include "basic_option.h"
#include "boost/numeric/ublas/matrix.hpp"
#include "LinearCongRand.h"
#include "inv_cdf_normal.h"


struct barrier_option:public basic_option{
	double barrier=0;
	int barrier_type=0;
	// 0 if up_in; 1 if up_out; 2 if down_in; 3 if down_out;

	barrier_option(double stock, double strike, double rf, double div, double maturity, double vol, double bar, double corp, int bar_type): basic_option( stock,  strike,  rf,  div,  maturity,  vol,  corp){
		barrier=bar;
		barrier_type=bar_type;
	};
	barrier_option(basic_option opt, double bar, int bar_type){
		S=opt.S;
		K=opt.K;
		r=opt.r;
		q=opt.q;
		T=opt.T;
		sigma=opt.sigma;
		type=opt.type;
		barrier=bar;
		barrier_type= bar_type;
	};
};

class MCBarrier {

private:
	boost::numeric::ublas::matrix<double> path;
	std::vector<double> sample_value;


public:
	MCBarrier(){};
	virtual ~MCBarrier(){};
	MCBarrier(long N, long m, barrier_option option)
	{
		double t=option.T/(double)m;
		path.resize(N, m+1);
		path.clear();
		// clean the data container

		LinearCongRand rng(N*m);
		boost::numeric::ublas::vector<double> rand= rng.rand();
		long k=0;
		for(long i=0;i<N;i++)
			for(long j=1;j<m+1;j++)
			{
				path(i,j)= inv_cdf_normal(rand(k));
				k++;
			}
		// now matrix path is filled by normal rn

		for(long i=0; i<N;i++)
		{
			path(i,0)=option.S;
			for(long j=1;j<m+1;j++)
			{
				path(i,j)=path(i,j-1)*exp((option.r-option.q-0.5*option.sigma*option.sigma)*t+option.sigma*std::sqrt(t)*path(i,j));
			}
		}

		// now matrix path is the stock price movement path

		for(long i=0; i<N; i++)
		{
			bool cond=false;
			if (option.barrier_type==0) // up in
			{
				cond=false;
				for(long j=0;j <m+1; j++)
				{
					if (path(i,j)>option.barrier)
						cond=true;
				}
			}
			if (option.barrier_type==1) //up _out
			{
				cond=true;
				for(long j=0; j<m+1; j++)
				{
					if(path(i,j)>option.barrier)
						cond=false;
				}
			}
			if (option.barrier_type==2) // down in
			{
				cond=false;
				for(long j=0;j <m+1; j++)
				{
					if(path(i,j)<option.barrier)
						cond=true;
				}
			}
			if(option.barrier_type==3)// down out
			{
				cond=true;
				for(long j=0; j<m+1; j++)
				{
					if(path(i,j)< option.barrier)
						cond=false;
				}
			}

			// now for ith path we determin whether this path still needs consideration or not.
			if (cond)
			{
				double a;
				if (option.type==0) //call
				{
					a=std::max(path(i,m)-option.K,0.0);
				}else if(option.type==1) //put
				{
					a=std::max(option.K- path(i,m),0.0);
				}

				sample_value.push_back(a*exp(-option.r*option.T));
			}else{
				sample_value.push_back(0);
			}
		}

	};
	double option_price(){
		double result=0;
		for(std::vector<double>::iterator it=sample_value.begin();it!=sample_value.end();it++)
		{
			result+=*it;
		}
		result/=(double) sample_value.size();
		return result;
	}
};

#endif /* MCBARRIER_H_ */
