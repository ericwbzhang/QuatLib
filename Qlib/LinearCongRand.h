/*
 * LinearCongRand.h
 *
 *  Created on: Oct 30, 2014
 *      Author: Eric
 */

#ifndef LINEARCONGRAND_H_
#define LINEARCONGRAND_H_
#include <boost/numeric/ublas/vector.hpp>

class LinearCongRand {
protected:
	boost::numeric::ublas::vector<double> result;
public:
	LinearCongRand();
	LinearCongRand(long n);
	virtual ~LinearCongRand();

	boost::numeric::ublas::vector<double> rand();
};

#endif /* LINEARCONGRAND_H_ */
