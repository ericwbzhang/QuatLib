/*
 * BBSRTree.h
 *
 *  Created on: Oct 25, 2014
 *      Author: ericmac
 */

#ifndef BBSRTREE_H_
#define BBSRTREE_H_
#include "BinomialBSTree.h"

class BBSRTree {

	BinomialBSTree * ptr1;
	BinomialBSTree * ptr2;

public:
	BBSRTree(basic_option opt, long n, int type);
	BBSRTree();
	virtual ~BBSRTree();

	double price();
	double delta();
	double gamma();
	double theta();

};

#endif /* BBSRTREE_H_ */
