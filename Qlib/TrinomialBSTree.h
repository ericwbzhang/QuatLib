/*
 * TrinomialBSTree.h
 *
 *  Created on: Oct 25, 2014
 *      Author: ericmac
 */

#ifndef TRINOMIALBSTREE_H_
#define TRINOMIALBSTREE_H_
#include "TrinomialTree.h"
#include "BS.h"

class TrinomialBSTree:public TrinomialTree {
public:
	TrinomialBSTree(basic_option opt, long n, int t): TrinomialTree(opt, n, t){};
	TrinomialBSTree();
	virtual ~TrinomialBSTree();

	virtual void load();


};



#endif /* TRINOMIALBSTREE_H_ */
