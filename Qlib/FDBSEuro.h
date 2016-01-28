/*
 * FDBS.h
 *
 *  Created on: Dec 12, 2014
 *      Author: Eric
 */

#ifndef FDBSEURO_H_
#define FDBSEURO_H_
#include "basic_option.h"
#include "FDHeatPDE.h"


class FD_BS_Euro {
private:
	struct g_down_call: public single_var_fcn{
		double K, c;
		g_down_call(double k, double A):K(k), c(A){};

		double operator ()(double x)
		{
			return K*exp(c*x)*std::max(exp(x)-1,0.0);
		};
	};

	struct g_left_call:public single_var_fcn{
		double operator ()(double t){
			return 0;
		};
	};

	struct g_right_call:public single_var_fcn{
		double x_r, K, r, q,sigma, c, d;
		g_right_call(double X_r, double k, double R, double Q, double Sigma,double A, double B):x_r(X_r),K(k), r(R),q(Q), sigma(Sigma) ,c(A), d(B){};
		double operator() (double t){
			return K*exp(c*x_r+d*t)*(exp(x_r-2*q*t/(sigma*sigma))- exp(-2*r*t/(sigma*sigma)));
		};
	};


	struct g_down_put: public single_var_fcn{
		double K, c;
		g_down_put(double k, double A): K(k), c(A){};
		double operator() (double x){
			return K*exp(c*x)*std::max(1-exp(x),0.0);
		};
	};

	struct g_left_put: public single_var_fcn{
		double x_l,K , r, q, sigma,c ,d;
		g_left_put(double X_l, double k, double R, double Q, double Sigma, double A, double B): x_l(X_l), K(k), r(R), q(Q), sigma(Sigma), c(A), d(B){};
		double operator() (double t){
			return K*exp(c*x_l+d*t)*(exp(-2*r*t/(sigma*sigma))- exp(x_l- 2*q*t/(sigma*sigma)));
		};
	};

	struct g_right_put: public single_var_fcn{
		double operator() (double t){
			return 0;
		};
	};

protected:
	virtual void search_S0(){
		long i=S_grid.size()/2;
		long j;
		if (S_grid[i]>option.S){
			do{
				i--;
			}while(S_grid[i]>option.S);
			j=i+1;
		}else if(S_grid[i+1]<option.S){
			j=i+1;
			do{
				j++;
			}while(S_grid[j]<option.S);
			i=j-1;
		}else{
			j=i+1;
		}

		if (S_grid[j]==option.S)
		{
			p=j;
		}else
		{
			p=i;
		}
	};

protected:
	double alpha;
	double a, b;
	double x_left;
	double x_right;
	basic_option option;
	std::vector<double> S_grid; // grid of stock price
	std::vector<double> time_grid;// grid of time
	boost::numeric::ublas::matrix<double> option_grid;// grid of option price
	FD_Heat_PDE* FD;
	long p;// S_grid[p] is the largest stock price grid node which is no greater than S0;
	// ie S_grid[p] <S0 or = S0;


	struct param{
		basic_option opt;
		double rate;
		long M;
		int method;
		bool iterative;
		double w;
		double tol;
		boost::numeric::ublas::vector<double> ini_guess;
		bool resi_rule;

		void load(basic_option opt, double rate, long M, int method, bool iterative, double w, double tol, boost::numeric::ublas::vector<double> ini_guess,bool resi_rule){
			this->opt=opt;
			this->rate=rate;
			this->M=M;
			this->method=method;
			this->iterative=iterative;
			this->w=w;
			this->tol=tol;
			this->ini_guess=ini_guess;
			this->resi_rule=resi_rule;
		} ;
	};

	param ini_param;
public:
	FD_BS_Euro(){};
	virtual ~FD_BS_Euro(){delete FD;};

	virtual double get_alpha(){return alpha;};
	virtual double get_xleft(){return x_left;};
	virtual double get_xright(){return x_right;};
	virtual long get_N(){return S_grid.size()-1;};
	virtual long get_M(){return time_grid.size()-1;};
	virtual double get_deltax(){return (x_right-x_left)/this->get_N();};
	virtual double get_deltatuo(){return option.T*option.sigma*option.sigma/2.0/(double )this->get_M();};
	virtual double get_xcomp(){return log(option.S/option.K);};
	virtual double get_taufinal(){return option.T*option.sigma*option.sigma/2.0;};
	virtual long get_p(){return p;};

	FD_BS_Euro(basic_option opt, double rate, long M, int method, bool iterative, double w=1.2, double tol=pow(10.0,-6.0), boost::numeric::ublas::vector<double> ini_guess=boost::numeric::ublas::vector<double>(0,0),bool resi_rule=false)
	// method =0 if forward Euler; 1 backward Euler; 2 Crank Nicolson
	// we compute the knock-out barrier option
	// this constructor does not require that the initial price is on the grid
	{
		ini_param.load(opt, rate, M,method, iterative, w, tol,ini_guess, resi_rule);
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
		if (option.type==0) //call
		{
			g_left_call gl;
			g_right_call gr(x_right, option.K, option.r, option.q, option.sigma, a, b);
			g_down_call gd(option.K, a);
			set.load(x_left, x_right, tuo, 0, N, M, gl,gr, gd, method, iterative, w,tol,ini_guess,  resi_rule);
		}else // put
		{
			g_left_put gl(x_left, option.K, option.r, option.q, option.sigma, a, b);
			g_right_put gr;
			g_down_put gd(option.K, a);
			set.load(x_left, x_right, tuo, 0, N, M, gl,gr, gd, method, iterative, w,tol,ini_guess,  resi_rule);
		}


		FD= new FD_Heat_PDE (set);

		option_grid=FD->grid_matrix();
		for( long i=0; i<option_grid.size1();i++)
			for (long j=0; j<option_grid.size2();j++)
			{
				option_grid(i,j)*= exp(-a*(x_left+j*x)-b*i*t);
			}

		this->search_S0(); //initialize p;

	};

	virtual double error_pointwise1(){
		long i=p;
		long j=i+1;

		long m=option_grid.size1()-1;
		double appro= option_grid(m,i)*(S_grid[j]-option.S)+ option_grid(m,j)*(option.S-S_grid[i]);
		appro/=(S_grid[j]-S_grid[i]);

		BS bs(option);
		double exact=bs.price();

		return std::abs(appro-exact);

	};

	virtual double price_approx1(){
		long i=p;
		long j=i+1;

		long m=option_grid.size1()-1;
		double appro= option_grid(m,i)*(S_grid[j]-option.S)+ option_grid(m,j)*(option.S-S_grid[i]);
		appro/=(S_grid[j]-S_grid[i]);

		return appro;
	};

	virtual double error_pointwise2(){
		long i=p;
		long j=i+1;

		long m=option_grid.size1()-1;
		double appro= FD->grid_matrix()(m,i)*log(S_grid[j]/option.S)+FD->grid_matrix()(m,j)*log(option.S/ S_grid[i]);
		appro/=log(S_grid[j]/S_grid[i]);
		appro=appro*exp(-a*log(option.S/option.K)-b*option.T*option.sigma*option.sigma/2.0);

		BS bs(option);
		double exact= bs.price();
		return std::abs(appro-exact);
	};

	virtual double price_approx2(){
		long i=p;
		long j=i+1;

		long m=option_grid.size1()-1;
		double appro= FD->grid_matrix()(m,i)*log(S_grid[j]/option.S)+FD->grid_matrix()(m,j)*log(option.S/ S_grid[i]);
		appro/=log(S_grid[j]/S_grid[i]);
		appro=appro*exp(-a*log(option.S/option.K)-b*option.T*option.sigma*option.sigma/2.0);

		return appro;
	};

	virtual double error_RMS()
	{
		basic_option opt=option;
		long count=0;
		double sum=0;
		double result=0;

		for(long i=0; i<S_grid.size();i++){
			opt.S=S_grid[i];
			BS bs(opt);
			double temp= bs.price();
			if(temp/option.S< pow(10.0, -5.0))
				continue;

			count++;
			sum+= pow( temp- option_grid(option_grid.size1()-1, i),2.0)/ pow(temp, 2.0);
		}

		result=sqrt(sum/count);
		return result;

	};

	virtual boost::numeric::ublas::matrix<double> price_grid(){return option_grid;};

	FD_Heat_PDE* core_PDE(){
		return FD;
	};


	virtual double delta(){
		long m=option_grid.size1()-1;
		return (option_grid(m,p+1)-option_grid(m, p))/(S_grid[p+1]- S_grid[p]);
	};

	virtual double gamma(){
		long  m= option_grid.size1()-1;
		double u= (option_grid(m,p+2)-option_grid(m,p+1))/(S_grid[p+2]-S_grid[p+1]);
		double v=(option_grid(m, p)-option_grid(m,p-1))/ (S_grid[p]-S_grid[p-1]);
		double w= (S_grid[p+2]+S_grid[p+1])/2.0;
		double z= (S_grid[p]+S_grid[p-1])/2.0;

		return (u-v)/(w-z);
	};

	virtual double theta(){
		long m=option_grid.size1()-2;
		double a=option_grid(m,p)*(S_grid[p+1]-option.S)+option_grid(m,p+1)*(option.S-S_grid[p]);
		a/=(S_grid[p+1]-S_grid[p]);

		long k=time_grid.size()-1;
		return (this->price_approx1()-a)/(time_grid[k]-time_grid[k-1]);

	};

};

#endif /* FDBSEURO_H_ */
