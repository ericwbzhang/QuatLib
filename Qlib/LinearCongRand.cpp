/*
 * LinearCongRand.cpp
 *
 *  Created on: Oct 30, 2014
 *      Author: Eric
 */

#include "LinearCongRand.h"

LinearCongRand::LinearCongRand() {
	// TODO Auto-generated constructor stub

}

LinearCongRand::~LinearCongRand() {
	// TODO Auto-generated destructor stub
}

LinearCongRand::LinearCongRand( long n)
{
	result.resize(n);
	long x=1;
	long k=pow(2.0, 31.0)-1;
	long a=39373;
	long c=0;
	for (long i=0; i<result.size(); i++)
	{
		x=(a*x+c)%k;
		result(i)=x*1.0/k;
	}
}

boost::numeric::ublas::vector<double> LinearCongRand::rand()
{
	return result;
}
