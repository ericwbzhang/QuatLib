/*
 * Secant.h
 *
 *  Created on: Dec 15, 2014
 *      Author: Eric
 */

#ifndef SECANT_H_
#define SECANT_H_
#include "FDHeatPDE.h"

inline std::vector<double> Secant(single_var_fcn &fcn, double tol,double guess1, double guess2){
	double diff=0;
	std::vector<double> result;
	double x=guess2;
	double x_old=guess1;
	result.push_back(guess1);
	result.push_back(guess2);
	do{
		double temp=x;
		x= temp-(temp-x_old)*fcn(temp)/(fcn(temp)-fcn(x_old));
		x_old=temp;
		diff=std::abs(x-x_old);
		result.push_back(x);

	}while (diff>tol);

	return result;

};



#endif /* SECANT_H_ */
