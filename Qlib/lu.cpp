/*
 * function.cpp
 *
 *  Created on: Aug 30, 2014
 *      Author: ericmac
 */

#include "lu.h"
using namespace boost::numeric::ublas;

vector<double> forward_subst(matrix<double> L, std::vector<double> b)
{
	std::vector<double> x=b;
	x[0]= b[0]/L(0,0);
	for (long j=1; j<x.size(); j++)
	{
		double sum=0;
		for (long k=0; k<j;k++)
		{
			sum+=L(j,k)*x[k];
		}
		x[j]=(b[j]-sum)/L(j,j);
	}
	boost::numeric::ublas::vector<double>result(x.size());
	for (long i=0; i<result.size(); i++)
		result(i)=x[i];
	return result;
}

vector<double> backward_subst(matrix<double> U, std::vector<double> b)
{
	std::vector <double > x=b;
	long n=x.size()-1;
	x[n]=b[n]/U(n,n);
	for(long j=n-1;j>=0;j--)
	{
		double sum=0;
		for(long k=n; k>j;k--)
		{
			sum+=U(j,k)*x[k];
		}
		x[j]=(b[j]-sum)/U(j,j);
	}
	boost::numeric::ublas::vector<double>result(x.size());
	for (long i=0; i<result.size(); i++)
		result(i)=x[i];
	return result;
}


vector<double> forward_subst (matrix<double> L, vector<double> b)
{
	std::vector<double> c;
	for(long i=0; i<b.size();i++)
		c.push_back(b(i));
	vector<double> result=forward_subst(L, c);
	return result;
}

vector<double> backward_subst(matrix<double> U, vector<double> b){
	std::vector<double> c;
	for(long i=0; i<b.size();i++)
		c.push_back(b(i));
	vector<double> result=backward_subst(U, c);
	return result;
}


LU<double> lu_no_pivoting   (matrix<double> A)
{
	long n=A.size1();
	triangular_matrix<double,lower> L(n,n);
	triangular_matrix<double, upper> U(n,n);
	identity_matrix<double> P(n);
	L.clear();
	U.clear();
	for(long i=0;i<n-1;i++)
	{
		for(long k=i;k<n;k++)
		{
			U(i,k)=A(i,k);
			L(k,i)=A(k,i)/U(i,i);
		}
		for(long j=i+1;j<n;j++)
		{
			for(long k=i+1; k<n;k++)
			{
				A(j,k)-=L(j,i)*U(i,k);
			}
		}
	}
	L(n-1,n-1)=1; U(n-1,n-1)=A(n-1,n-1);
	LU<double> lu;
	lu.U=U; lu.L=L; lu.P=P;
	return lu;
}

LU<double> lu_no_pivoting_banded(matrix<double> A, long m)
{
	long n=A.size1();
	triangular_matrix<double,lower> L(n,n);L.clear();
	triangular_matrix<double, upper> U(n,n);U.clear();
	identity_matrix<double> P(n);

	for(long i=0;i<n-1;i++)
	{
		for (long k=i; k<n;k++)
		{
			U(i,k)=A(i,k);L(k,i)=A(k,i)/U(i,i);
		}
		for(long j=i+1; j<std::min(n,i+m+1);j++)
		{
			for(long k=i+1;k<std::min(n,i+m+1);k++)
			{
				A(j,k)-=L(j,i)*U(i,k);
			}
		}
	}
	L(n-1,n-1)=1; U(n-1,n-1)=A(n-1, n-1);

	LU<double> lu;
	lu.U=U; lu.L=L; lu.P=P;
	return lu;

}
matrix<double>  row_pivot(matrix<double> &A,long a, long b)
{
	long n=A.size2();
	for(long j=0;j<n;j++)
	{
		std::swap(A(a,j),A(b,j));
	}
	return A;
}



LU<double> lu_row_pivoting(matrix <double> A)
{
	long n=A.size1();
	matrix <double> L(n,n);L.clear();
	matrix< double> U(n,n);U.clear();
	identity_matrix <double> I(n,n);
	matrix<double> P(I);P.clear();
	//L=P;
	for(long i=0; i<n-1;i++)
	{
		std::vector<double> vec;
		vec.clear();
		for(long k=i;k<n;k++)
		{
			vec.push_back(std::abs(A(k,i)));
		}
		long a=std::max_element(vec.begin(),vec.end())-vec.begin()+i;
		if(a!=i)
		{
			row_pivot(A,a ,i);
			row_pivot(P,a,i);
			if(i>0) row_pivot(L, a,i);
		}
		for(long j=i;j<n;j++)
		{
			L(j,i)=A(j,i)/A(i,i);
			U(i,j)=A(i,j);
		}
		for(long j=i+1; j<n;j++)
		{
			for(long k=i+1; k<n;k++)
			{
				A(j,k)-=L(j,i)*U(i,k);
			}
		}
	}
	L(n-1,n-1)=1; U(n-1,n-1)=A(n-1,n-1);

	LU <double> lu;
	lu.P=P; lu.U=U;lu.L=L;
	return lu;
}

boost::numeric::ublas::vector<double>	lu_solver(matrix<double> A, boost::numeric::ublas::vector<double> b)
{
	LU<double> lu=lu_row_pivoting(A);
	// We solve PAx=Pb
	boost::numeric::ublas::vector<double> Pb=boost::numeric::ublas::prod(lu.P,b);
	boost::numeric::ublas::vector<double> y=forward_subst(lu.L,Pb);
	boost::numeric::ublas::vector<double> x=backward_subst(lu.U, y);
	return x;
}


