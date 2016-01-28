/*
 * interative_solver.cpp
 *
 *  Created on: Sep 19, 2014
 *      Author: ericmac
 */
#include "iterative_solver.h"
#include "lu.h"
#include <iomanip>
#include <iostream>
using namespace boost::numeric::ublas;


Iter_sol Jacobi_solver(matrix<double> A, boost::numeric::ublas::vector<double> b, double tol, boost::numeric::ublas::vector<double> guess, bool resi_rule)
// if resi_rule =true we use residual criterion
{

	mtx_decomp mtx_dec= mtx_decompose(A);

	long count=0;

	boost::numeric::ublas::vector<double> x=guess;
	boost::numeric::ublas::vector<double> resi= b- boost::numeric::ublas::prod(A,x);
	double stop_iter_resi= tol*norm_2(resi);

	boost::numeric::ublas::vector<double> b_new= b;
	for (long i=0; i<b_new.size(); i++)
		b_new(i)/=mtx_dec.diag(i,i);

	//std::cout<<"stop:"<<stop_iter_resi<<std::endl;

	if (resi_rule)
	while (norm_2(resi)> stop_iter_resi)
	{
		boost::numeric::ublas::vector<double> y= prod(mtx_dec.low+ mtx_dec.up, x);
		for (long i=0 ; i<y.size(); i++)
			y(i)/= mtx_dec.diag(i,i);

		x=b_new- y;
		resi=b-prod(A, x);
		if(count<=3)
			std::cout<<std::setprecision(15)<<x<<std::endl;
		count++;
		//std::cout<<x_mt<<"\n";
		//std::cout << count<<"\t\t"<<std::setprecision(9)<<vec_norm(mtx2vec(resi))<<std::endl;
	}
	else
	{
		double diff=1000;
		boost::numeric::ublas::vector<double> x_old=x;
		while(diff>tol)
		{
			boost::numeric::ublas::vector<double> y=prod(mtx_dec.low+ mtx_dec.up, x);
			for (long i=0 ; i<y.size(); i++)
				y(i)/= mtx_dec.diag(i,i);

			x=b_new- y;
			diff=norm_2(x-x_old);
			x_old=x;
			count++;

		}
	}
	resi=b-prod(A,x);
	Iter_sol result;
	result.x=x;
	result.iter= count;
	result.resi=resi;
	return result;

}

Iter_sol GS_solver(matrix<double> A, boost::numeric::ublas::vector<double> b, double tol, boost::numeric::ublas::vector<double> guess, bool resi_rule)
{
	mtx_decomp mtx_dec= mtx_decompose(A);

	matrix<double> M=mtx_dec.diag+mtx_dec.low;
	matrix<double> N=mtx_dec.up;

	long count=0;

	boost::numeric::ublas::vector<double> x=guess;
	boost::numeric::ublas::vector<double> resi= b-prod(A,x);
	double stop_iter_resi= tol*norm_2(resi);

	boost::numeric::ublas::vector<double> b_new=forward_subst(M,b);
	if (resi_rule)
	while (norm_2(resi)> stop_iter_resi)
	{
		boost::numeric::ublas::vector<double> y=prod(N,x);
		x=-forward_subst(M,y)+b_new;
		resi=b-prod(A,x);
		if(count<=3)
			std::cout<<std::setprecision(15)<<x<<std::endl;
		count++;
	}
	else
	{
		double diff=1000;
		boost::numeric::ublas::vector<double> x_old=x;
		while(diff>tol)
		{
			//std::cout<<count<<std::endl;
			boost::numeric::ublas::vector<double> y=prod(N,x);
			x=-forward_subst(M,y)+b_new;
			diff=norm_2(x-x_old);
			x_old=x;
			count++;
		}
	}
	resi=b-prod(A,x);
	Iter_sol result;
	result.x=x;
	result.iter= count;
	result.resi=resi;
	return result;


}

Iter_sol SOR_solver(matrix<double> A, boost::numeric::ublas::vector<double> b, double weight, double tol, boost::numeric::ublas::vector<double> guess, bool resi_rule)
{
	mtx_decomp mtx_dec= mtx_decompose(A);

	matrix<double> M=mtx_dec.diag+mtx_dec.low*weight;
	matrix<double> N=-weight* mtx_dec.up+(1-weight)*mtx_dec.diag;

	long count=0;

	boost::numeric::ublas::vector<double> x=guess;
	boost::numeric::ublas::vector<double> resi= b-prod(A,x);
	double stop_iter_resi= tol*norm_2(resi);
	boost::numeric::ublas::vector<double> b_new= weight* forward_subst(M,b);

	if (resi_rule)
	while(norm_2(resi)> stop_iter_resi)
	{
		boost::numeric::ublas::vector<double> y=prod(N,x);
		x=forward_subst(M, y)+b_new;
		resi=b-prod(A,x);
		if(count<=3)
			std::cout<<std::setprecision(15)<<x<<std::endl;
		count++;
	}
	else
	{
		double diff=1000;
		boost::numeric::ublas::vector<double> x_old=x;
		while(diff>tol)
		{
			boost::numeric::ublas::vector<double> y=prod(N,x);
			x=forward_subst(M,y)+b_new;
			diff=norm_2(x-x_old);
			x_old=x;
			count++;
		}
	}


	resi=b-prod(A,x);
	Iter_sol result;
	result.x=x;
	result.iter= count;
	result.resi=resi;
	return result;
}




mtx_decomp mtx_decompose(matrix<double> A)
{

	matrix<double> L(A);
	matrix<double> U(A);
	matrix<double> D(A);
	D.clear(); L.clear(); U.clear();
	for (long i=0; i<A.size1();i++)
		for (long j=0; j<A.size2(); j++)
		{
			if (i<j) // upper
				U(i,j)=A(i,j);
			else if (i==j) //diagnoal
				D(i,j)=A(i,j);
			else // lower
				L(i,j)=A(i,j);
		}
	mtx_decomp result;
	result.orig=A;
	result.diag=D;
	result.up=U;
	result.low=L;

	return result;

}





