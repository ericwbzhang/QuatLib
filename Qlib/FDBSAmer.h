/*
 * FDBSAmer.h
 *
 *  Created on: Dec 14, 2014
 *      Author: Eric
 */

#ifndef FDBSAMER_H_
#define FDBSAMER_H_
#include "FDBSEuro.h"
#include "FDHeatPDEConstr.h"



class FD_BS_Amer:public FD_BS_Euro {


public:

	std::vector<double> opt_exercise_price(){
		std::vector<double > result;
		boost::numeric::ublas::matrix<double> exe_mtx=this->exercise_mtx();
		for(long i=0; i<exe_mtx.size1(); i++)
		{
			long j=1;
			for(j=1; j<exe_mtx.size2();j++)
			{
				if(exe_mtx(i,j)!=exe_mtx(i,0))
					break;
			}
			double a=0.5*(S_grid[j]+S_grid[j-1]);
			result.push_back(a);

		}
		return result;
	};
	boost::numeric::ublas::matrix<double> exercise_mtx(){
		return FD->characteristic_mtx();
	};

	double error_pointwise1(double exact)
	{
		return std::abs(this->price_approx1()-exact);
	};

	double error_pointwise2(double exact){
		return std::abs(this->price_approx2()-exact);
	};


	FD_BS_Euro corresponding_FD_Euro(){
		FD_BS_Euro fd(ini_param.opt, ini_param.rate, ini_param.M, ini_param.method, ini_param.iterative, ini_param.w, ini_param.tol, ini_param.ini_guess, ini_param.resi_rule);
		return fd;
	};

	double calibrated_price(){
		FD_BS_Euro fd=this->corresponding_FD_Euro();
		BS bs(option);
		return this->price_approx1()+bs.price()-fd.price_approx1();
	};

	double calibrated_error_pointwise(double exact){
		return std::abs(exact-this->calibrated_price());
	};


	FD_BS_Amer(){};
	virtual ~FD_BS_Amer(){};
	FD_BS_Amer(basic_option opt, double rate, long M, int method, bool iterative, double w=1.2, double tol=pow(10.0, -6.0), boost::numeric::ublas::vector<double> ini_guess=boost::numeric::ublas::vector<double>(0,0),bool resi_rule=false)
	{
		ini_param.load(opt,rate, M, method, iterative,w, tol, ini_guess, resi_rule);

		option=opt;
		x_left=log(option.S/option.K)+(option.r-option.q-0.5*option.sigma*option.sigma)*option.T-3*option.sigma*sqrt(option.T);
		x_right=x_left+6*option.sigma*sqrt(option.T);
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
		double_var_fcn* low_con;
		if(option.type==0) //call
		{
			g_left_call gl;
			g_right_call gr(x_right, option.K, a, b);
			g_down_call gd(option.K, a);
			low_con=new low_bound_call(option.K, a, b);

			set.load(x_left, x_right, tuo,0,N,M, gl, gr, gd, method,iterative,w, tol, ini_guess, resi_rule);
		}else // put
		{
			g_left_put gl(x_left, option.K, a, b);
			g_right_put gr;
			g_down_put gd(option.K, a);
			low_con=new low_bound_put(option.K, a, b);

			set.load(x_left, x_right, tuo,0,N,M, gl, gr, gd, method,iterative,w, tol, ini_guess, resi_rule);
		}

		FD=new FD_Heat_PDE_Constr(set, low_con);
		option_grid=FD->grid_matrix();

		for( long i=0; i<option_grid.size1();i++)
			for (long j=0; j<option_grid.size2();j++)
			{
				option_grid(i,j)*= exp(-a*(x_left+j*x)-b*i*t);
			}
		this->search_S0(); //initialize p;



		delete low_con;
	};



private:
	struct g_down_call:public single_var_fcn{
		double K, a;
		g_down_call(double k, double A): K(k) , a(A){};
		double operator() (double x){
			return exp(a*x)*K*std::max(exp(x)-1, 0.0);
		};
	};

	struct g_left_call:public single_var_fcn{
		double operator() (double t){return 0;};
	};

	struct g_right_call:public single_var_fcn{
		double x_r, K, a ,b;
		g_right_call(double xr, double k, double A, double B): x_r(xr), K(k), a(A), b(B){};
		double operator() (double t){
			return K*(exp(x_r)-1)*exp(a*x_r+b*t);
		};
	};

	struct g_down_put:public single_var_fcn{
		double K, a;
		g_down_put(double k, double A): K(k), a(A){};
		double operator() (double x){
			return exp(a*x)*K*std::max(1-exp(x),0.0);
		};
	};

	struct g_left_put:public single_var_fcn{
		double x_l, K, a, b;
		g_left_put(double xl, double k, double A, double B): x_l(xl), K(k), a(A), b(B){};
		double operator() (double t){
			return K*(1-exp(x_l))*exp(a*x_l+b*t);
		}
	};

	struct g_right_put:public single_var_fcn{
		double operator() (double t){ return 0;};
	};


	struct low_bound_call:public double_var_fcn{
		double K, a, b;
		low_bound_call (double k, double A, double B):K(k), a(A), b(B){};
		double operator ( ) (double t, double x){
			return K*exp(a*x+b*t)*std::max(exp(x)-1,0.0);
		};
	};
	struct low_bound_put: public double_var_fcn{
		double K, a, b;
		low_bound_put(double k, double A, double B): K(k), a(A), b(B){};
		double operator() (double t, double x){
			return K * exp(a*x+b*t)* std::max(1-exp(x),0.0);
		};
	};




};

#endif /* FDBSAMER_H_ */
