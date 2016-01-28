/*
 * TBSRTree.cpp
 *
 *  Created on: Oct 25, 2014
 *      Author: ericmac
 */

#include "TBSRTree.h"

TBSRTree::TBSRTree() {
	// TODO Auto-generated constructor stub

}

TBSRTree::~TBSRTree() {
	// TODO Auto-generated destructor stub
	delete ptr1;
	delete ptr2;

}

TBSRTree::TBSRTree(basic_option opt, long n, int type)
{
	ptr1= new TrinomialBSTree ( opt, n/2 , type);
	ptr2 = new TrinomialBSTree(opt, n, type);

	ptr1->load();
	ptr2->load();
}

double TBSRTree::price()
{return 2* ptr2->price()-ptr1->price();}

double TBSRTree::delta()
{return 2*ptr2->delta()-ptr1->delta();}

double TBSRTree::gamma()
{
	return 2*ptr2->gamma()-ptr1->gamma();
}

double TBSRTree::theta()
{return 2*ptr2->theta()- ptr1->theta();}


