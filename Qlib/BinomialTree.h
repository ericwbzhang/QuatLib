/*
 * BinomialTree.h
 *
 *  Created on: Oct 11, 2014
 *      Author: ericmac
 */

#ifndef BINOMIALTREE_H_
#define BINOMIALTREE_H_
#include "basic_option.h"
#include "boost/numeric/ublas/matrix.hpp"


class BinomialTree {
protected:
	basic_option option;
	long N;
	boost::numeric::ublas::matrix<double> option_value;
	int type;
	double t, u ,d, p ,q, disc;


public:
	BinomialTree(basic_option opt, long n, int Euro_Amer);
	BinomialTree();
	virtual ~BinomialTree();

	virtual void load();

	double price();
	double delta();
	double gamma();
	double theta();
	virtual double average_price();
	virtual double average_delta();

};

#endif /* BINOMIALTREE_H_ */
