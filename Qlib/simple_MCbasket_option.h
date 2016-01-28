/*
 * simple_MCbasket_option.h
 *
 *  Created on: Nov 8, 2014
 *      Author: Eric
 */

#ifndef SIMPLE_MCBASKET_OPTION_H_
#define SIMPLE_MCBASKET_OPTION_H_

struct simple_basket_option
{
	double S1, S2, K, T, r, sigma1, sigma2, corr;
	double q1=0; double q2=0;

	simple_basket_option(double s1, double s2, double t, double rf, double k, double vol1, double vol2, double correlation, double d1=0, double d2=0):S1(s1), S2(s2), T(t), r(rf), K(k), sigma1(vol1), sigma2(vol2), corr(correlation), q1(d1), q2(d2){};
	simple_basket_option(){};

};


double simple_MCbasket_option(simple_basket_option opt, long n);



#endif /* SIMPLE_MCBASKET_OPTION_H_ */
