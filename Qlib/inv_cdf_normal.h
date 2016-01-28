/*
 * inv_cdf_normal.h
 *
 *  Created on: Oct 31, 2014
 *      Author: Eric
 */

#ifndef INV_CDF_NORMAL_H_
#define INV_CDF_NORMAL_H_


/* Generate normal random variable from uniforms
 * Beasley-Springer Moro algorithm for computing the inverse cumulative normal function
 * Copyright(C) Eric CrazyQuant 01/04/2012
 * key words: normal random variable, inverse normal cdf
 */

/* Simulation of inversed cumulative density function of normal distribution
 */


inline double inv_cdf_normal(double  x)
{
	if(x>0 && x<1)
	{
 // initialization
        const double a[4] = {2.50662823884, -18.61500062529, 41.39119773534, -25.44106049637};
        const double b[4] = {-8.47351093090, 23.08336743743, -21.06224101826, 3.13082909833};
        const double c[9] = {0.3374754822726147, 0.9761690190917186, 0.1607979714918209, 0.0276438810333863, 0.0038405729373609, 0.0003951896411919, 0.0000321767881768, 0.0000002888167364, 0.0000003960315187};
        double y = x - 0.5;
        double r,s,t,result;

        if (fabs(y) < 0.42){    //Beasley-Springer
            r = y * y;
            result = y * (a[0] + r * (a[1] + r * (a[2] + r * a[3]))) / (1.0 + r * (b[0] + r * (b[1] + r * (b[2] + r * b[3]))));
        }
        else{                   //Moro
            if(y <= 0) r = x;
            else r = 1 - x;

            s = log(-log(r));
            t = c[0] + s * (c[1] + s * (c[2] + s * (c[3] + s * (c[4] + s * (c[5] + s * (c[6] + s * (c[7] + s * c[8])))))));

     if(x > 0.5) result =  t;
     else result = -t;
        }
        return result;
	}else{
		std::cout<<"The input of inv_cdf_normal() is illegal.\n";
		exit(0);
	}

}



#endif /* INV_CDF_NORMAL_H_ */
