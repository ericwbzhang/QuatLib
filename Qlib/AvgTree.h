/*
 * AvgTree.h
 *
 *  Created on: Oct 25, 2014
 *      Author: ericmac
 */

#ifndef AVGTREE_H_
#define AVGTREE_H_
#include "BinomialBSTree.h"

class AvgTree {
	// compute the average of n-1 and n

	BinomialTree * ptr1=NULL;
	BinomialTree* ptr2=NULL;

public:
	AvgTree(basic_option opt, long n, int method, int type);
	// return the average of n-1 and n.
		// method= 0 if BinomialTree; method =1 if BBS Tree;
		//type = 0 if Euro; type =1 if American

	AvgTree();
	virtual ~AvgTree();

	double price();
	double delta();
	double gamma();
	double theta();


};

#endif /* AVGTREE_H_ */
