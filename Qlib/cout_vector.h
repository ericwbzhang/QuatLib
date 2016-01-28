/*
 * cout_vector.h
 *
 *  Created on: Sep 1, 2014
 *      Author: ericmac
 */

#ifndef COUT_VECTOR_H_
#define COUT_VECTOR_H_

#include <iostream>
#include <vector>
#include <stdio.h>


template <class T>
ostream & operator << (ostream & os, std::vector<T>   vec)
{
	for(long i=0; i<vec.size(); i++)
	//for(vector< T>::iterator it=vec.begin(); it!= vec.end(); it++)
	{
		os<< vec[i]<<"\t";
	}
	return os;
}

#endif /* COUT_VECTOR_H_ */
