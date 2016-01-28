/*
 * normal_distribution.h
 *
 *  Created on: Oct 11, 2014
 *      Author: ericmac
 */

#ifndef NORMAL_DISTRIBUTION_H_
#define NORMAL_DISTRIBUTION_H_

#include <boost/math/distributions/normal.hpp>

double Normalpdf(double a)
{
	boost::math::normal_distribution<double> normal(0,1.0);
	return pdf(normal,a);
}
double Normalcdf(double a)
{
	boost::math::normal_distribution<double> normal(0,1.0);
	return cdf(normal, a);
}



#endif /* NORMAL_DISTRIBUTION_H_ */
