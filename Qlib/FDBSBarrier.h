/*
 * FDBSBarrier.h
 *
 *  Created on: Dec 15, 2014
 *      Author: Eric
 */

#ifndef FDBSBARRIER_H_
#define FDBSBARRIER_H_
#include "FDBSEuro.h"


class FD_BS_Barrier:public FD_BS_Euro  {
protected:
	bool left_bounded;
	bool right_bounded;
	double lbar,rbar ;
	single_var_fcn* gl;
	single_var_fcn* gr;
	single_var_fcn* gd;



public:
	FD_BS_Barrier(){};
	virtual ~FD_BS_Barrier(){delete gl; delete gr; delete gd;};

	FD_BS_Barrier(basic_option opt, double left_bar, double right_bar, double rate, long M, int method, bool iterative, double w=1.2, double tol=pow(10.0,-6.0), boost::numeric::ublas::vector<double> ini_guess=boost::numeric::ublas::vector<double>(0,0),bool resi_rule=false )
	{
		ini_param.load(opt, rate, M,method, iterative, w, tol,ini_guess, resi_rule);
		lbar=left_bar;
		rbar=right_bar;
		left_bounded=false;
		right_bounded=false;

		option=opt;
		x_left=log(option.S/option.K)+(option.r-option.q-0.5*option.sigma*option.sigma)*option.T-3*option.sigma*sqrt(option.T);
		x_right=x_left+6*option.sigma*sqrt(option.T);

		if(x_left<log(lbar/option.K))
		{
			x_left=log(lbar/option.K);
			left_bounded=true;
		}
		if(x_right>log(rbar/option.K))
		{
			x_right=log(rbar/option.K);
			right_bounded=true;
		}

		double tuo=option.T*option.sigma*option.sigma/2.0;
		double t=tuo/(double ) M;
		long N=std::floor((x_right-x_left)/sqrt(t/rate));
		double x=(x_right-x_left)/N;
		alpha=t/(x*x);// note alpha will be sightly less than rate
		a=(option.r-option.q)/(option.sigma*option.sigma)-.5;
		b=pow((option.r-option.q)/(option.sigma*option.sigma)+0.5,2.0)+2*option.q/(option.sigma*option.sigma);

		for (long i=0; i<N+1; i++)
		{
			S_grid.push_back(exp(x_left+i*x)*option.K);
		}
		for (long i=0; i<M+1; i++)
		{
			time_grid.push_back(option.T-2/(option.sigma*option.sigma)* (i*t));
		}


		FD_setting set;
		if(option.type==0)//call
		{
			if (left_bounded)// left barrier
			{
				gl=new g_left_barrier(x_left,a, b);
			}else {
				gl=new g_left_call;
			}


			if (right_bounded)// right barrier
			{
				gr=new g_right_barrier(x_right, a,b);
			}else{
				gr=new g_right_call(x_right,option.K, option.r, option.q, option.sigma, a, b);
			}

			gd=new g_down_call(option.K,a);

		}else //put
		{
			if (left_bounded)// left barrier
			{
				gl=new g_left_barrier(x_left, a ,b);
			}else {
				gl=new g_left_put(x_left, option.K, option.r, option.q, option.sigma, a, b);
			}


			if (right_bounded)// right barrier
			{
				gr=new g_right_barrier(x_right, a,b);
			}else{
				gr=new g_right_put;
			}

			gd=new g_down_put(option.K,a);

		}

		set.load(x_left,x_right,tuo,0,N,M, *gl, *gr, *gd, method, iterative, w, tol, ini_guess, resi_rule);

		FD=new FD_Heat_PDE (set);

		option_grid=FD->grid_matrix();

		for( long i=0; i<option_grid.size1();i++)
			for (long j=0; j<option_grid.size2();j++)
			{
				option_grid(i,j)*= exp(-a*(x_left+j*x)-b*i*t);
			}

		this->search_S0(); //initialize p;

	};

	double error_pointwise1(double exact)
	{
		return std::abs(this->price_approx1()-exact);
	};

	double error_pointwise2(double exact){
		return std::abs(this->price_approx2()-exact);
	};


private:
	struct g_left_call:public single_var_fcn{
		double operator ()(double t){
			return 0;
		};
	};

	struct g_left_barrier:public single_var_fcn{
		//double operator() (double t){return 0;}
		double x_l, c ,d;
		g_left_barrier(double X_l, double A, double B): x_l(X_l), c(A), d(B){std::cout<<"left barrier";};
		double operator() (double t){
			return 2*exp(c*x_l+d*t);
		};

	};


	struct g_left_put: public single_var_fcn{
		double x_l,K , r, q, sigma,c ,d;
		g_left_put(double X_l, double k, double R, double Q, double Sigma, double A, double B): x_l(X_l), K(k), r(R), q(Q), sigma(Sigma), c(A), d(B){};
		double operator() (double t){
			return K*exp(c*x_l+d*t)*(exp(-2*r*t/(sigma*sigma))- exp(x_l- 2*q*t/(sigma*sigma)));
		};
	};

	struct g_right_call:public single_var_fcn{
		double x_r, K, r, q,sigma, c, d;
		g_right_call(double X_r, double k, double R, double Q, double Sigma,double A, double B):x_r(X_r),K(k), r(R),q(Q), sigma(Sigma) ,c(A), d(B){};
		double operator() (double t){
			return K*exp(c*x_r+d*t)*(exp(x_r-2*q*t/(sigma*sigma))- exp(-2*r*t/(sigma*sigma)));
		};
	};

	struct g_right_put: public single_var_fcn{
		double operator() (double t){
			return 0;
		};
	};

	struct g_right_barrier:public single_var_fcn{
		//double operator( ) (double t) {return 0;}
		double x_r, c, d;
		g_right_barrier(double X_r, double A, double B): x_r(X_r), c(A), d(B){std::cout<<"right barrier";};
		double operator() (double t) {
			return 2*exp(c*x_r+d*t);
		};
	};

	struct g_down_call: public single_var_fcn{
		double K, c;
		g_down_call(double k, double A):K(k), c(A){};

		double operator ()(double x){
			return K*exp(c*x)*std::max(exp(x)-1,0.0);
		};
	};

	struct g_down_put: public single_var_fcn{
		double K, c;
		g_down_put(double k, double A): K(k), c(A){};
		double operator() (double x){
			return K*exp(c*x)*std::max(1-exp(x),0.0);
		};
	};
};

#endif /* FDBSBARRIER_H_ */
