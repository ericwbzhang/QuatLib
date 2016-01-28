/*
 * FDImpliedVol.h
 *
 *  Created on: Dec 15, 2014
 *      Author: Eric
 */

#ifndef FDIMPLIEDVOL_H_
#define FDIMPLIEDVOL_H_
#include "Secant.h"
#include "FDBSAmer.h"


struct price_sigma:public single_var_fcn{

	basic_option opt;
	double rate;
	long M;
	int method;
	bool iterative;
	double w;
	double tol;
	boost::numeric::ublas::vector<double> ini_guess;
	bool resi_rule;
	int type;
	double price;


	price_sigma(double price,int type,  basic_option opt, double rate, long M, int method, bool iterative, double w=1.2, double tol=pow(10.0,-6.0), boost::numeric::ublas::vector<double> ini_guess=boost::numeric::ublas::vector<double>(0,0),bool resi_rule=false){
		this->opt=opt;
		this->rate=rate;
		this->M=M;
		this->method=method;
		this->iterative=iterative;
		this->w=w;
		this->tol=tol;
		this->ini_guess=ini_guess;
		this->resi_rule=resi_rule;
		this->type=type;
		this->price=price;
	};

	double operator() (double sigma)
	{
		opt.sigma=sigma;
		double result;
		if(type==0) //Euro
		{
			FD_BS_Euro FD(opt,rate, M, method, iterative, w, tol, ini_guess, resi_rule);
			result=FD.price_approx1()-price;

		}else // Amer
		{
			FD_BS_Amer FD(opt,rate, M, method, iterative, w, tol, ini_guess, resi_rule);
			result=FD.price_approx1()-price;
		}

		return result;
	};
};


inline std::vector<double> FD_Implied_Vol(double price, double tol_secant, double guess1, double guess2, int type, basic_option opt, double rate, long M, int method, bool iterative, double w=1.2, double tol=pow(10.0, -6.0), boost::numeric::ublas::vector<double> ini_guess=boost::numeric::ublas::vector<double>(0,0),bool resi_rule=false)
//type=0 for Euro ; 1 for American
{
	price_sigma fcn(price,type, opt, rate, M, method, iterative, w, tol, ini_guess, resi_rule);
	return Secant(fcn, tol_secant, guess1, guess2);
};



#endif /* FDIMPLIEDVOL_H_ */
