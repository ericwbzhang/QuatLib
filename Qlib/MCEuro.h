/*
 * MCEuro.h
 *
 *  Created on: Oct 31, 2014
 *      Author: Eric
 */

#ifndef MCEURO_H_
#define MCEURO_H_

#include "basic_option.h"
#include "boost/numeric/ublas/vector.hpp"


class MC_Euro {
protected:
	basic_option option;
	long N;
	boost::numeric::ublas::vector<double> price_path;
	boost::numeric::ublas::vector<double> delta_path;
	boost::numeric::ublas::vector<double> vega_path;
	boost::numeric::ublas::vector<double> stock_path;

	void load();
public:
	MC_Euro();
	virtual ~MC_Euro();
	MC_Euro(basic_option opt, long n);

	double opt_price();
	double opt_delta();
	double opt_vega();
	boost::numeric::ublas::vector<double> stockpath();

};

#endif /* MCEURO_H_ */
