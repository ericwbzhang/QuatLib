/*
 * TBSRTree.h
 *
 *  Created on: Oct 25, 2014
 *      Author: ericmac
 */

#ifndef TBSRTREE_H_
#define TBSRTREE_H_
#include "TrinomialBSTree.h"

class TBSRTree {
protected:
	TrinomialBSTree * ptr1;
	TrinomialBSTree * ptr2;

public:
	TBSRTree(basic_option opt, long n, int type);
	TBSRTree();
	virtual ~TBSRTree();

	double price();
	double delta();
	double gamma();
	double theta();


};

#endif /* TBSRTREE_H_ */
