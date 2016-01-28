/*
 * iterative_solver.h
 *
 *  Created on: Sep 19, 2014
 *      Author: ericmac
 */

#ifndef ITERATIVE_SOLVER_H_
#define ITERATIVE_SOLVER_H_

#include <vector>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>

using namespace boost::numeric::ublas;


struct Iter_sol
{boost::numeric::ublas::vector<double> x; long iter=0;boost::numeric::ublas::vector<double> resi;};


struct mtx_decomp
{
	matrix<double> orig;
	matrix<double> diag;
	matrix<double> low;
	matrix<double> up;
};

struct Iterative_Set{
	bool resi_rule=false;
	double tol=pow(10.0, -6.0);
	boost::numeric::ublas::vector<double> guess;
	double weight;

	Iterative_Set(){};
	Iterative_Set(double tolerance, double w, boost::numeric::ublas::vector<double> ini_guess, bool rule): tol(tolerance), weight(w),guess(ini_guess), resi_rule(rule){};
};

mtx_decomp mtx_decompose(matrix<double> A);


Iter_sol Jacobi_solver(matrix<double> A, boost::numeric::ublas::vector<double> b, double tol, boost::numeric::ublas::vector<double> guess, bool resi_rule=true);

inline Iter_sol Jacobi_solver(matrix<double>A, boost::numeric::ublas::vector<double> b, Iterative_Set set){
	return Jacobi_solver(A, b, set.tol, set.guess, set.resi_rule);
};

Iter_sol GS_solver(matrix<double> A, boost::numeric::ublas::vector <double > b, double tol, boost::numeric::ublas::vector<double> guess, bool resi_rule=true);

inline Iter_sol GS_solver(matrix<double> A, boost::numeric::ublas::vector<double> b, Iterative_Set set){
	return GS_solver(A,b,set.tol, set.guess, set.resi_rule);
};

Iter_sol SOR_solver(matrix<double> A, boost::numeric::ublas::vector<double > b,double weight ,double tol, boost::numeric::ublas::vector<double> guess, bool resi_rule=true);

inline Iter_sol SOR_solver(matrix<double> A, boost::numeric::ublas::vector<double> b, Iterative_Set set){
	return SOR_solver(A,b, set.weight, set.tol, set.guess, set.resi_rule);
};










#endif /* ITERATIVE_SOLVER_H_ */
