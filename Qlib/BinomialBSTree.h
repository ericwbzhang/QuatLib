/*
 * BinomialBSTree.h
 *
 *  Created on: Oct 12, 2014
 *      Author: ericmac
 */

#ifndef BINOMIALBSTREE_H_
#define BINOMIALBSTREE_H_
#include "BinomialTree.h"

class BinomialBSTree:public BinomialTree {
public:
	BinomialBSTree(basic_option opt, long n, int Euro_Amer) :BinomialTree(opt, n, Euro_Amer){};
	BinomialBSTree();
	virtual ~BinomialBSTree();

	virtual void load();
	virtual double average_price();
	virtual double average_delta();

};

#endif /* BINOMIALBSTREE_H_ */
