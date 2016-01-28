/*
 * cholesky.h
 *
 *  Created on: Sep 6, 2014
 *      Author: ericmac
 */

#ifndef CHOLESKY_H_
#define CHOLESKY_H_

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>

using namespace boost::numeric::ublas;

struct chole
{matrix<double> U; bool spd; };

chole cholesky(matrix<double> A);
chole cholesky_banded(matrix<double> A, long m);
chole cholesky_tridiagonal (matrix<double> A);

boost::numeric::ublas::vector<double> cholesky_solver(matrix<double> A, std::vector<double> b);
boost::numeric::ublas::vector<double> cholesky_banded_solver(matrix<double> A, long m, std::vector<double> b);
boost::numeric::ublas::vector<double> cholesky_tri_solver(matrix<double> A, std::vector<double> b);

boost::numeric::ublas::vector<double> cholesky_solver(matrix<double> A, boost::numeric::ublas::vector<double> b);

boost::numeric::ublas::vector<double> cholesky_banded_solver(matrix<double> A, long m, boost::numeric::ublas::vector<double> b);
boost::numeric::ublas::vector<double> cholesky_tri_solver(matrix<double> A, boost::numeric::ublas::vector<double> b);


#endif /* CHOLESKY_H_ */
