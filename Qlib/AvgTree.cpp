/*
 * AvgTree.cpp
 *
 *  Created on: Oct 25, 2014
 *      Author: ericmac
 */

#include "AvgTree.h"

AvgTree::AvgTree() {
	// TODO Auto-generated constructor stub

}

AvgTree::~AvgTree() {
	// TODO Auto-generated destructor stub
	delete ptr1;
	delete ptr2;

}

AvgTree::AvgTree(basic_option opt, long n, int method, int type)
{
	// return the average of n-1 and n.
	// method= 0 if BinomialTree; method =1 if BBS Tree;
	//type = 0 if Euro; type =1 if American
	if (method==0)// Binomial Tree
	{
		ptr1=new BinomialTree(opt, n-1 , type);
		ptr2= new BinomialTree (opt, n, type);
	}
	else if (method ==1)// Binomial BS Tree
	{
		ptr1= new BinomialBSTree( opt, n-1, type);
		ptr2= new BinomialBSTree (opt, n, type);
	}


	ptr1->load();
	ptr2->load();
}

double AvgTree::price()
{return 0.5*(ptr1->price()+ptr2->price());}

double AvgTree::delta()
{
	return 0.5*(ptr1->delta()+ptr2->delta());
}

double AvgTree::gamma()
{return 0.5*(ptr1->gamma()+ptr2->gamma());}

double AvgTree::theta()
{
	return 0.5*(ptr1->theta()+ptr2->theta());
}



