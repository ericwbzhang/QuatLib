/*
 * FDHeatPDEConstr.h
 *
 *  Created on: Dec 14, 2014
 *      Author: Eric
 */

#ifndef FDHEATPDECONSTR_H_
#define FDHEATPDECONSTR_H_
#include "FDHeatPDE.h"
#include "stdio.h"


class FD_Heat_PDE_Constr : public FD_Heat_PDE{
protected:
	boost::numeric::ublas::matrix<double> chara_mtx;

public:

	FD_Heat_PDE_Constr(){};
	virtual ~FD_Heat_PDE_Constr(){};
	boost::numeric::ublas::matrix<double> characteristic_mtx()	{
		return chara_mtx;
	}
	FD_Heat_PDE_Constr(FD_setting FD, double_var_fcn *low_con)
	// for Constrained FD Heat PDE solver we only implement forward and Crank Nicolson - SOR method
	// So type only take 0 or 2 with iterative=true

	{
		set=FD;

		double t=(FD.up-FD.down)/FD.M;
		double x=(FD.right-FD.left)/FD.N;
		double alpha= t/(x*x);
		grid.resize(FD.M+1, FD.N+1);
		chara_mtx.resize(FD.M+1, FD.N+1);
		grid.clear();
		chara_mtx.clear();


		for(long i=0; i<grid.size1(); i++)
		{
			grid(i,0)=(*FD.g_left) (t*i+FD.down);
			grid(i,grid.size2()-1)=(*FD.g_right)(t*i+FD.down);

		}// compute the left and right edge value and assign to the matrix grid.

		for(long j=1;j<grid.size2()-1;j++)
		{
			grid(0,j)=(*FD.f)(FD.left+x*j);
		}// compute the down edge value and assign to the matrix grid.




		if(FD.type==0)// forward Euler
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

					U(j-1)=std::max(U(j-1),(*low_con)(FD.down+(i+1)*t, FD.left+j*x));
					grid(i+1,j)=U(j-1);

				}
				b.clear();
			}

		}else if(FD.type==1){
			std::cout<< "Backward Euler method is not applicable!\n"<<std::endl;
		}
		else if(FD.type==2)// Crank Nicolson
			// we use SOR method
		{

			boost::numeric::ublas::matrix<double> B(FD.N-1, FD.N-1);
			B.clear();
			boost::numeric::ublas::vector<double> U(FD.N-1);
			U.clear();
			for(long i=0; i<U.size();i++)
			{
				U(i)=(*FD.f)(FD.left+x*(i+1));
				B(i,i)=1-alpha;
			}
			for(long i=1;i<B.size1();i++)
			{
				B(i,i-1)=0.5*alpha;
				B(i-1,i)=0.5*alpha;
			}
			boost::numeric::ublas::vector<double> b(U.size());
			b.clear();
			// matrix A , B and U are ready


			for( long i=0; i<FD.M;i++)
			{
				b(0)=(*FD.g_left)(FD.down+i*t)*alpha*0.5+ (*FD.g_left)(FD.down+t*(i+1))*alpha* 0.5;
				b(b.size()-1)=(*FD.g_right)(FD.down+i*t)*alpha*0.5+(*FD.g_right)(FD.down+t*(i+1))*alpha*0.5;
				b+=boost::numeric::ublas::prod(B,U);

				// AX=b the b is ready
				// Remember to update U in each iteration
				// and clear b at the end

				boost::numeric::ublas::vector<double> guess(b.size());
				guess.clear();
				for(long k=0; k<guess.size();k++)
				{
					guess(k)=(*low_con)(FD.down+(i+1)*t,FD.left+(k+1)*x);
				}
				// ******SOR******

				boost::numeric::ublas::vector<double> z=guess;
				boost::numeric::ublas::vector<double> z_old=guess;
				double diff=10000;

				while(diff>FD.tol)
				{
					z(0)=(1-FD.weight)*z_old(0)+ FD.weight*alpha/(2*(1+alpha))*(z_old(1)) +FD.weight/(1+alpha)*b(0);
					if(z(0)<guess(0))  z(0)=guess(0);
					for(long k=1; k<z.size()-1;k++)
					{
						z(k)=(1-FD.weight)*z_old(k)+ FD.weight*alpha/(2*(1+alpha))*(z(k-1)+z_old(k+1)) +FD.weight/(1+alpha)*b(k);
					if(z(k)<guess(k))  z(k)=guess(k);
					}
					z(FD.N-2)=(1-FD.weight)*z_old(FD.N-2)+ FD.weight*alpha/(2*(1+alpha))*(z(FD.N-3))+FD.weight/(1+alpha)*b(FD.N-2);
					if(z(FD.N-2)<guess(FD.N-2))  z(FD.N-2)= guess(FD.N-2);

					diff=boost::numeric::ublas::norm_2(z-z_old);
					z_old=z;

				}
				for(long j=1;j<chara_mtx.size2()-1;j++)
				{
					if(z(j-1)==guess(j-1))
						chara_mtx(i+1,j)=1;
				}

				U=z;
				b.clear();

				for(long j=1; j<grid.size2()-1;j++)
				{
					grid(i+1,j)=z(j-1);
				}

			}
			for(long j=1;j<grid.size2()-1;j++)
			{
				if (grid(0,j)>0)
					chara_mtx(0,j)=1;
			}
			for(long i=0;i<grid.size1();i++)
			{
				if(grid(i,0)>0)
					chara_mtx(i,0)=1;
				if(grid(i,FD.N)>0)
					chara_mtx(i,FD.N)=1;
			}

		}


	};

};

#endif /* FDHEATPDECONSTR_H_ */
