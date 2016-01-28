/*
 * TrinomialTree.h
 *
 *  Created on: Oct 25, 2014
 *      Author: ericmac
 */

#ifndef TRINOMIALTREE_H_
#define TRINOMIALTREE_H_
#include "BinomialTree.h"


class TrinomialTree: public BinomialTree {
protected:
	double m;
	double v;

public:
	TrinomialTree(basic_option opt, long n, int Euro_Amer );
	TrinomialTree();
	virtual ~TrinomialTree();

	virtual void load();
	virtual double delta();
	virtual double gamma();
	virtual double theta();
	virtual double  price();


};

#endif /* TRINOMIALTREE_H_ */
