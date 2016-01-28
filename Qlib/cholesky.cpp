/*
 * cholesky.cpp
 *
 *  Created on: Oct 29, 2014
 *      Author: Eric
 */

/*
 * cholesky.cpp
 *
 *  Created on: Sep 6, 2014
 *      Author: ericmac
 */
#include "cholesky.h"
#include "lu.h"
#include <boost/numeric/ublas/triangular.hpp>

using namespace boost::numeric::ublas;


chole cholesky(matrix<double> A)
{
	long n=A.size1();
	chole cho; cho.spd=true;
	triangular_matrix<double, upper> U(n,n);
	U.clear();
	for (long i=0; i<n-1; i++)
	{
		U(i,i)=std::sqrt(A(i,i));
		if(!(U(i,i)>0)) cho.spd=false;
		for(long k=i+1; k<n; k++)
		{
			U(i,k)= A(i,k)/U(i,i);
		}
		for(long j=i+1; j<n; j++)
			for(long k=j; k<n; k++)
				A(j,k)-=U(i,j)*U(i,k);
	}
	U(n-1,n-1)= std::sqrt(A(n-1,n-1));
	if(!(U(n-1,n-1)>0)) cho.spd=false;
	cho.U=U;
	return cho;
}

chole cholesky_banded (matrix<double> A, long m)
{
	long n=A.size1();
	chole cho; cho.spd=true;
	triangular_matrix<double, upper> U(n,n);
	U.clear();
		for (long i=0; i<n-1; i++)
		{
			U(i,i)=std::sqrt(A(i,i));
			if(!(U(i,i)>0)) cho.spd=false;
			for(long k=i+1; k<n; k++)
			{
				U(i,k)= A(i,k)/U(i,i);
			}
			for(long j=i+1; j<std::min(n,i+m+1); j++)
				for(long k=j; k<std::min(n,i+m+1); k++)
					A(j,k)-=U(i,j)*U(i,k);
		}
		U(n-1,n-1)= std::sqrt(A(n-1,n-1));
		if(!(U(n-1,n-1)>0)) cho.spd=false;
		cho.U=U;
		return cho;
}

chole cholesky_tridiagonal(matrix<double> A)
{
	long n=A.size1();
		chole cho; cho.spd=true;
		triangular_matrix<double, upper> U(n,n);
		U.clear();
			for (long i=0; i<n-1; i++)
			{
				U(i,i)=std::sqrt(A(i,i));
				if(!(U(i,i)>0)) cho.spd=false;
				U(i,i+1) = A(i,i+1)/U(i,i);
				A(i+1, i+1)-=std::pow(U(i,i+1),2);
			}
			U(n-1,n-1)= std::sqrt(A(n-1,n-1));
			if(!(U(n-1,n-1)>0)) cho.spd=false;
			cho.U=U;
			return cho;
}

vector<double > cholesky_solver(matrix<double> A,std::vector<double> b)
{
	chole cho=cholesky(A);
	vector<double> y= forward_subst(trans(cho.U), b);
	vector<double> x= backward_subst(cho.U,y);
	vector<double> result(x.size());

	return x;

}

vector<double> cholesky_banded_solver(matrix<double> A, long m, std::vector<double > b)
{
	chole cho=cholesky_banded(A,m);
		vector<double> y= forward_subst(trans(cho.U), b);
		vector<double> x= backward_subst(cho.U,y);

		return x;
}

vector<double> cholesky_tri_solver(matrix<double> A, std::vector<double> b){
	chole cho=cholesky_tridiagonal ( A);
	vector<double> y= forward_subst(trans(cho.U), b);
	vector<double> x= backward_subst(cho.U,y);

			return x;
}


vector<double> cholesky_solver(matrix<double> A, vector<double> b)
{
	std::vector<double> c;
	for(long i=0;i<b.size();i++)
		c.push_back(b(i));
	vector<double> result= cholesky_solver(A,c);
	return result;

}

vector<double> cholesky_banded_solver(matrix<double> A, long m, vector<double> b){

	std::vector<double> c;
	for(long i=0;i<b.size();i++)
			c.push_back(b(i));
	vector<double> result= cholesky_banded_solver(A,m, c);
	return result;


}
vector<double> cholesky_tri_solver(matrix<double> A, vector<double> b){
	std::vector<double> c;
	for(long i=0;i<b.size();i++)
			c.push_back(b(i));
	vector<double> result= cholesky_tri_solver(A, c);
	return result;

}



