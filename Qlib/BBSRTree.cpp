/*
 * BBSRTree.cpp
 *
 *  Created on: Oct 25, 2014
 *      Author: ericmac
 */

#include "BBSRTree.h"

BBSRTree::BBSRTree() {
	// TODO Auto-generated constructor stub

}

BBSRTree::~BBSRTree() {
	// TODO Auto-generated destructor stub
	delete ptr1;
	delete ptr2;

}

BBSRTree::BBSRTree(basic_option opt, long n ,int type)
{
	ptr1=new BinomialBSTree ( opt, n/2, type);
	ptr2=new BinomialBSTree( opt, n, type);

	ptr1->load();
	ptr2->load();
}

double BBSRTree::price()
{
	return 2*ptr2->price()-ptr1->price();
}

double BBSRTree::delta()
{
	return 2*ptr2->delta()-ptr1->delta();
}

double BBSRTree::gamma()
{
	return 2*ptr2->gamma()-ptr1->gamma();
}

double BBSRTree::theta()
{
	return 2*ptr2->theta()-ptr1->theta();
}

