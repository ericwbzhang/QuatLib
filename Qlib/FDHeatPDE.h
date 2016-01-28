/*
 * FDHeatPDE.h
 *
 *  Created on: Dec 8, 2014
 *      Author: Eric
 */

#ifndef FDHEATPDE_H_
#define FDHEATPDE_H_
#include "boost/numeric/ublas/matrix.hpp"
#include "lu.h"
#include "iterative_solver.h"

struct single_var_fcn{
	virtual double operator() (double x)=0;
	virtual ~single_var_fcn(){};
};

struct double_var_fcn{
	virtual double operator() (double x, double y)=0;
	virtual ~double_var_fcn(){};
};


struct FD_setting{
	double left, right, up;
	double down=0;
	long N, M;
	single_var_fcn* g_left;
	single_var_fcn* g_right;
	single_var_fcn* f;

	int type; //type=0 forward; 1 backward; 2 Crank Nicolson
	bool iterative=false;
	double weight=0;
	double tol=pow(10.0, -6.0);
	boost::numeric::ublas::vector<double> guess;
	bool resi_rule=false;

	FD_setting(){};

	FD_setting(double l, double r, double u, double d, long n, long m, 	single_var_fcn &g_l , single_var_fcn& g_r, single_var_fcn& g_d , int _type, bool iter=false, double w=0, double _tol=pow(10.0,-6.0), boost::numeric::ublas::vector<double> ini_guess=boost::numeric::ublas::vector<double> (0,0),bool rule =false){
		left=l;
		right=r;
		up=u;
		down=d;
		N=n; M=m;
		g_left=&g_l;
		g_right=&g_r;
		f=&g_d;
		type=_type;
		iterative=iter;
		weight=w;
		tol=_tol;
		guess=ini_guess;
		resi_rule=rule;
	};

	void load(double l, double r, double u, double d, long n, long m, 	single_var_fcn &g_l , single_var_fcn& g_r, single_var_fcn& g_d , int _type, bool iter=false, double w=0, double _tol=pow(10.0,-6.0), boost::numeric::ublas::vector<double> ini_guess=boost::numeric::ublas::vector<double> (0,0),bool rule =false)
	{
		left=l;
		right=r;
		up=u;
		down=d;
		N=n; M=m;
		g_left=&g_l;
		g_right=&g_r;
		f=&g_d;
		type=_type;
		iterative=iter;
		weight=w;
		tol=_tol;
		guess=ini_guess;
		resi_rule=rule;

	};

};



class FD_Heat_PDE {
protected:
	FD_setting set;
	boost::numeric::ublas::matrix<double> grid;

public:
	FD_Heat_PDE(){};
	virtual ~FD_Heat_PDE(){};
	FD_Heat_PDE(FD_setting FD)
	// FD.type=0 if forward Euler;1 if backward Euler; 2 if Crank Nicolson
	{
		set=FD;
		double t=(FD.up-FD.down)/FD.M;
		double x=(FD.right-FD.left)/FD.N;
		double alpha= t/(x*x);
		grid.resize(FD.M+1, FD.N+1);


		for(long i=0; i<grid.size1(); i++)
		{
			grid(i,0)=(*FD.g_left) (t*i+FD.down);
			grid(i,grid.size2()-1)=(*FD.g_right)(t*i+FD.down);

		}// compute the left and right edge value and assign to the matrix grid.

		for(long j=1;j<grid.size2()-1;j++)
		{
			grid(0,j)=(*FD.f)(FD.left+x*j);
		}// compute the down edge value and assign to the matrix grid.




		if(FD.type==0) //forward Euler
		{
			boost::numeric::ublas::matrix<double> A(FD.N-1, FD.N-1);
			A.clear();
			boost::numeric::ublas::vector<double> U(FD.N-1);
			U.clear();
			for(long i=0; i<U.size();i++)
			{
				U(i)=(*FD.f)(FD.left+x*(i+1));
				A(i,i)=1-2*alpha;
			}
			for(long i=1;i<A.size1();i++)
			{
				A(i,i-1)=alpha;
				A(i-1,i)=alpha;
			}
			boost::numeric::ublas::vector<double> b(U.size());
			b.clear();
			// matrix A and U are ready

			for(long i=0; i<FD.M;i++)
			{

				b(0)=(*FD.g_left)(FD.down+i*t)*alpha;
				b(b.size()-1)=(*FD.g_right)(FD.down+i*t)*alpha;
				U=boost::numeric::ublas::prod(A,U)+b;
				for(long j=1;j<grid.size2()-1;j++)
				{

					grid(i+1,j)=U(j-1);
				}
				b.clear();
			}
		}else if (FD.type==1) //backward Euler
		{
			boost::numeric::ublas::matrix<double> A(FD.N-1, FD.N-1);
			A.clear();
			boost::numeric::ublas::vector<double> U(FD.N-1);
			U.clear();
			for(long i=0; i<U.size();i++)
			{
				U(i)=(*FD.f)(FD.left+x*(i+1));
				A(i,i)=1+2*alpha;
			}
			for(long i=1;i<A.size1();i++)
			{
				A(i,i-1)=-alpha;
				A(i-1,i)=-alpha;
			}
			boost::numeric::ublas::vector<double> b(U.size());
			b.clear();
			// matrix A and U are ready

			for(long i=0; i<FD.M;i++)
			{
				b(0)=(*FD.g_left)(FD.down+(i+1)*t)*alpha;
				b(b.size()-1)=(*FD.g_right)(FD.down+(i+1)*t)*alpha;

				if(FD.iterative) //SOR
				{
					if (FD.guess.size()!=U.size())
					{
						FD.guess.resize(U.size());
						FD.guess.clear();
					}
					Iter_sol result=SOR_solver(A, U+b, FD.weight, FD.tol, FD.guess, FD.resi_rule);
					U=result.x;
				}else// lu
				{
					LU<double> lu=lu_no_pivoting(A);
					boost::numeric::ublas::vector<double> y=forward_subst(lu.L, U+b);
					U=backward_subst(lu.U, y);
				}

				for(long j=1;j<grid.size2()-1;j++)
				{
					grid(i+1,j)=U(j-1);
				}
				b.clear();
			}

		}else if (FD.type==2) // Crank Nicolson
		{
			boost::numeric::ublas::matrix<double> A(FD.N-1, FD.N-1);
			A.clear();
			boost::numeric::ublas::matrix<double> B(FD.N-1, FD.N-1);
			B.clear();
			boost::numeric::ublas::vector<double> U(FD.N-1);
			U.clear();
			for(long i=0; i<U.size();i++)
			{
				U(i)=(*FD.f)(FD.left+x*(i+1));
				A(i,i)=1+alpha;
				B(i,i)=1-alpha;
			}
			for(long i=1;i<A.size1();i++)
			{
				A(i,i-1)=-0.5*alpha;
				A(i-1,i)=-0.5*alpha;
				B(i,i-1)=0.5*alpha;
				B(i-1,i)=0.5*alpha;
			}
			boost::numeric::ublas::vector<double> b(U.size());
			b.clear();
			// matrix A , B and U are ready

			for(long i=0; i<FD.M;i++)
			{
				b(0)=(*FD.g_left)(FD.down+i*t)*alpha*0.5+ (*FD.g_left)(FD.down+t*(i+1))*alpha* 0.5;
				b(b.size()-1)=(*FD.g_right)(FD.down+i*t)*alpha*0.5+(*FD.g_right)(FD.down+t*(i+1))*alpha*0.5;

				if(FD.iterative) //SOR
				{
					if (FD.guess.size()!=U.size())
					{
						FD.guess.resize(U.size());
						FD.guess.clear();
					}
					Iter_sol result=SOR_solver(A, boost::numeric::ublas::prod(B,U)+b, FD.weight, FD.tol, FD.guess, FD.resi_rule);
					U=result.x;
				}else// lu
				{
					LU<double> lu=lu_no_pivoting(A);
					boost::numeric::ublas::vector<double> y=forward_subst(lu.L, boost::numeric::ublas::prod(B,U)+b);
					U=backward_subst(lu.U, y);
				}

				for(long j=1;j<grid.size2()-1;j++)
				{
					grid(i+1,j)=U(j-1);
				}
				b.clear();
			}



		}

	};
	virtual boost::numeric::ublas::matrix<double> grid_matrix(){return grid;};

	virtual double max_pointwise_error(double_var_fcn &exa){
		double_var_fcn * exact=&exa;
		double max=0;
			for(long j=0; j<grid.size2();j++){
				double error= grid(grid.size1()-1,j)- (*exact)(set.up, set.left+(set.right-set.left)/(double) set.N*j);
				error=std::abs(error);
				if (max<error)
					max=error;

			}
		return max;
	};

	virtual double RMS_error(double_var_fcn & exa){
		double result=0;
		double_var_fcn* exact=&exa;
		for (long j=0 ; j<grid.size2();j++)
		{
			double exact_value=(*exact)(set.up, set.left+(set.right-set.left)/(double) set.N*j);
			double diff= exact_value-grid(grid.size1()-1, j);
			result+= pow(diff/ exact_value, 2.0);
		}
		result/=(grid.size2());
		result=std::sqrt(result);

		return result;

	};

	virtual boost::numeric::ublas::matrix<double> characteristic_mtx(){
		std::cout<< "FDHeatPDE has no characteristic matrix!\n\n";
        boost::numeric::ublas::matrix<double> result;

        return result;
	};



};

#endif /* FDHEATPDE_H_ */
