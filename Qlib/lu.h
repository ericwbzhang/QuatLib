/*
 * function.h
 *
 *  Created on: Aug 30, 2014
 *      Author: ericmac
 */

#ifndef FUNCTION_H_
#define FUNCTION_H_
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include<boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/triangular.hpp>

using namespace boost::numeric::ublas;

boost::numeric::ublas::vector<double> forward_subst (matrix<double> L, boost::numeric::ublas::vector<double> b);
boost::numeric::ublas::vector<double> backward_subst(matrix<double> U, boost::numeric::ublas::vector<double> b);

boost::numeric::ublas::vector<double> forward_subst (matrix<double> L, std::vector<double> b);
boost::numeric::ublas::vector<double> backward_subst(matrix<double> U, std::vector<double> b);


template <class T>
struct LU
{
	matrix <T> L; matrix <T> U;matrix <T> P;
};

LU<double> lu_no_pivoting (matrix<double> A);

matrix<double> row_pivot(matrix<double> &A,long a, long b);

LU<double> lu_row_pivoting(matrix<double>  A);

LU<double> lu_no_pivoting_banded(matrix<double> A, long m);

boost::numeric::ublas::vector<double> lu_solver(matrix<double> A, boost::numeric::ublas::vector<double> b);








#endif /* FUNCTION_H_ */
