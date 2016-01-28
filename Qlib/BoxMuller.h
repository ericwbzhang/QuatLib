/*
 * BoxMuller.h
 *
 *  Created on: Nov 8, 2014
 *      Author: Eric
 */

#ifndef BOXMULLER_H_
#define BOXMULLER_H_
#define _USE_MATH_DEFINES

#include <math.h>


using namespace std;

struct BoxMuller
{
	double u1=0;
	double u2=0;
	double z1=0;
	double z2=0;

	BoxMuller(double a, double b):u1(a), u2(b){
		z1=sqrt(-2*log(u1))*cos(2*M_PI*u2);
		z2=sqrt(-2*log(u1))*sin(2*M_PI*u2);
	};
};



#endif /* BOXMULLER_H_ */
