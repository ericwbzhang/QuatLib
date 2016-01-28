/*
 * basic_option.h
 *
 *  Created on: Oct 11, 2014
 *      Author: ericmac
 */

#ifndef BASIC_OPTION_H_
#define BASIC_OPTION_H_

struct basic_option
{
	double S,K;
	double r=0;
	double q=0;
	double T;
	double sigma;
	double type=0;// 0 means call and 1 means put;

	basic_option(){};
	basic_option(double stock, double strike, double rf, double div, double maturity, double vol, double corp):S(stock), K(strike), r(rf), q(div), T(maturity), sigma(vol), type(corp) {};

};



#endif /* BASIC_OPTION_H_ */
