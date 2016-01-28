//============================================================================
// Name        : NMF.cpp
// Author      : Eric Z
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================


#include "csv2vector.h"
#include"vec2csv.h"
#include "cout_vector.h"
#include "iomanip"
#include "AvgTree.h"
#include "BBSRTree.h"
#include "BS.h"
#include "TBSRTree.h"
#include "MCEuro.h"
#include "boost/numeric/ublas/io.hpp"
#include "iostream"
#include "simple_MCbasket_option.h"
#include"LinearCongRand.h"
#include "inv_cdf_normal.h"
#include "iomanip"
#include "DiscreteDivdBT.h"
#include "FDBSEuro.h"
#include "MCBarrier.h"
#include "FDHeatPDE.h"
#include "FDBSAmer.h"
#include "FDBSBarrier.h"
#include "FDImpliedVol.h"
#include "Secant.h"





using namespace boost::numeric::ublas;

struct g_left:public single_var_fcn{
	double operator () (double x){
		return exp(x-2);
	};
};
struct g_right:public single_var_fcn{
	double operator() (double x){
		return exp(x+2);
	};
};
struct f:public single_var_fcn{
	double operator() (double x){
		return exp(x);
	};
};

struct exact: public double_var_fcn{
	double operator() (double t, double x){
		return  exp(x+t);
	};
};

int main()

{
	std::cout<<std::setprecision(10)<<std::endl;
	//basic_option option(41, 40, 0.04, 0.02,0.75,0.35,0);
	basic_option opt1(52,50,0.02, 0.015,11.0/12.0, 0.3,0);
	basic_option opt2(50,48,0.02,0.005,0.25,0.3,0);

	
	double L=40;
	double U=60;
	double xl=log(L/opt2.K);
	double xr=log(U/opt2.K);
	double xcomp=log(opt2.S/opt2.K);
	double tau=0.5*opt2.T*opt2.sigma*opt2.sigma;
	for(long M=4;M<300; M*=4)
	{
		double a=4;
		double t=tau/M;
		long N=std::floor((xr-xl)/sqrt(t/a));
		double x=(xr-xl)/N;
		double alpha=t/(x*x);

		std::cout<<alpha<<"\t"<<N<<"\t"<<xl<<"\t"<<xr<<"\t"<<xcomp<<"\t"<<tau<<"\t"<<t<<"\t"<<x<<std::endl;


	}





	for(int M=4; M<300;M=M*4)
	{
	FD_BS_Barrier FD_amer(opt2,40, 60, 0.5, M,0,false,1.2 ,pow(10.0, -8.0));
	std::cout<<FD_amer.get_alpha()<<"\t"<<FD_amer.get_N()<<"\t"<<FD_amer.get_xleft()<<"\t"<<FD_amer.get_xright()<<"\t"<<FD_amer.get_xcomp()<<"\t"<<FD_amer.get_taufinal()<< "\t"<<FD_amer.get_deltatuo()<<"\t"<<FD_amer.get_deltax()<<std::endl;

	std::cout<<"************\n\n\n";
	boost::numeric::ublas::matrix<double> u_mtx=FD_amer.core_PDE()->grid_matrix();
	//std::cout<<u_mtx<<"\n\n";

	std::cout<<u_mtx(u_mtx.size1()-1, FD_amer.get_p())<<"\t"<<u_mtx(u_mtx.size1()-1, FD_amer.get_p()+1)<<"\t";
	std::cout<<FD_amer.price_approx1()<<"\t";
	std::cout<< FD_amer.delta()<<"\t"<<FD_amer.gamma()<<"\t"<<FD_amer.theta()<<std::endl;

	}

	FD_BS_Barrier FD_amer(opt2,40, 60, 0.5, 4,2,false,1.2 ,pow(10.0, -8.0));
	std::cout<<FD_amer.core_PDE()->grid_matrix();

	std::cout<<FD_amer.price_grid();

	/*
	for(long M=4;M<270; M*=4)
	{
	FD_BS_Amer FD_amer(opt1, 0.5, M,0,true,1.2 ,pow(10.0, -8.0));
	std::cout<<FD_amer.calibrated_price()<<std::endl;
	}





	/*
	TBSRTree bt(opt1, 1280,  1);
	//bt.load();
	std::cout<<bt.price()<<std::endl;

	/*

	FD_BS_Barrier FD(option, 35,10000,5,4,2,true);

	std::cout<< std::setprecision(9)<<FD.core_PDE()->grid_matrix();

	std::cout<<FD.delta()<<std::endl<<FD.gamma()<<std::endl<<FD.theta()<<std::endl<<FD.price_approx1()<<std::endl;

	std::cout<<"*********\n\n";
	FD_BS_Amer FD_amer(option, 0.45, 256, 2, true);

	//FD_Heat_PDE *fd=FD_amer.core_PDE();
	//std::cout<<std::setprecision(10)<< fd->grid_matrix()<<std::endl;

	std::cout<< FD_amer.error_pointwise1(4.083817051176386)<< std::endl<< FD_amer.error_pointwise2(4.083817051176386)<<"\n\n\n";
	std::cout<<FD_amer.delta()<<std::endl<<FD_amer.gamma()<< std::endl<< FD_amer.theta()<<std::endl<<FD_amer.calibrated_price()<<std::endl<< FD_amer.calibrated_error_pointwise(4.083817051176386);


	std::cout<<FD_amer.exercise_mtx()<<std::endl<<FD_amer.opt_exercise_price();




	g_right gr;

	std::cout<<"\n\n**********\n";
	basic_option opt(38,40, 0.04, 0.01,1.0/12.0,0.3,1);
	std::cout<<FD_Implied_Vol(2.45,0.0001,0.1, 0.4,1,opt,0.4,16, 0, true)<<std::endl;
	//std::cout<<Secant( gr, 0.000001, 1,1.5);



/*
	simple_basket_option option(25,30, 0.5, 0.05, 50, 0.3, 0.2, 0.35);


	for (int k=0; k<=8; k++)
	{
		long N=pow(10.0,4.0)*pow(2.0,float(k));
		std::cout<<N<<"\n";
		std::cout<< simple_MCbasket_option( option, N)<<std::endl<<std::endl;

	}

*/
	/*
	matrix<double> A4(8,8);
	for (long i=0;i<A4.size1(); i++) A4(i,i)=9;
	for (long i=0; i<A4.size1()-2;i++) A4(i,i+2)=1;
	for(long i=2; i<A4.size1();i++)	 A4(i,i-2)=-2;
	for (long i=0;i<A4.size1()-3;i++) A4(i,i+3)=2;
	for (long i=3;i<A4.size1(); i++) A4(i,i-3)=-1;

	csv2vector c2v;
	vec2csv v2c;
	c2v.load("A4.csv");
	matrix<double> A4=c2v.data_mtx_double();
	std::cout<< A4<<std::endl;
	boost::numeric::ublas::vector<double> b4(8);
	for (long i=0;i<b4.size();i++)
		b4(i)= (3*i-4.0)/(i*i+1.0);
	//exact solution using lu


	boost::numeric::ublas::vector<double> guess(8);
	for (long i=0;i<guess.size(); i++) guess(i)=0;
	double tol= pow(10.0,-6.0);

	cout<<"*****************\nJacobi\n";
	// Jacobbi
	Iter_sol solution;
	solution= Jacobi_solver(A4,b4,tol,guess);
	std::cout<<std::setprecision(15)<<solution.x<<"\n"<<solution.iter<<std::endl;


	//GS
	cout<<"*****************\nGS\n";
	solution=GS_solver(A4,b4,tol, guess);
	cout<<solution.x<<"\n"<<solution.iter<<std::endl;


	///SOR
	cout<<"*****************\nSOR\n";
	double w1=0.9;
	double w2=1.15;
	solution=SOR_solver(A4, b4, w1, tol,guess);
	cout<<solution.x<<"\n"<<solution.iter<<std::endl;
	solution=SOR_solver(A4, b4, w2, tol,guess);
	cout<<solution.x<<"\n"<<solution.iter<<std::endl;










/*
	// Prepare
	vec2csv v2c;
	csv2vector c2v;
	c2v.load("9.csv");
	matrix<double> A=c2v.data_mtx_double();
	c2v.load("9_b.csv");
	std::vector<double> temp=c2v.data_double();


	boost::numeric::ublas::vector<double> guess (A.size2());
	guess.clear();
	boost::numeric::ublas::vector<double> b (temp.size());
	for (long i=0; i<b.size(); i++)
		b(i)=temp[i];
	double tol=pow(10,-6);
	cout<<A<<endl<<guess<<endl;
	cout<<b<<endl<<tol<<endl;

/*
	 // USE CHOLESKY

	cout<<"\nCholesky\n";
	boost::numeric::ublas::vector<double> exact_sol=cholesky_solver(A,b);
	std::cout<<exact_sol<<std::endl<<std::endl;
	std::cout<<prod(A,exact_sol)<<std::endl;
	v2c.load("Cholesky.csv", exact_sol,exact_sol.size(),1);
	v2c.export2csv();
	cout<<"**********************\n";
*/

/*
	Iter_sol solution;
	// Jacobi
	cout<<"\nJacobi_resi\n";
	solution= Jacobi_solver(A,b,tol,guess);
	cout<<solution.x<<endl<<solution.iter<<endl;
	cout<<prod(A,solution.x)<<std::endl;
	std::vector<double> result;
	result.clear();
	for(long i=0;i<solution.x.size();i++)
		result.push_back(solution.x(i));
	result.push_back(solution.iter);
	result.push_back(norm_2(solution.resi));
	v2c.load("Jacobi_resi.csv",result, result.size(),1);
	v2c.export2csv();
	result.clear();
	cout<<"*******************\n";

	cout<<"\nJacobi_cons\n";
	solution=Jacobi_solver(A,b,tol,guess,false);
	cout<<solution.x<<endl<<solution.iter<<endl;
	cout<<prod(A,solution.x)<<std::endl;

	result.clear();
	for(long i=0;i<solution.x.size();i++)
			result.push_back(solution.x(i));
	result.push_back(solution.iter);
	result.push_back(norm_2(solution.resi));
	v2c.load("Jacobi_cons.csv",result, result.size(),1);
	v2c.export2csv();
	cout<<"*******************\n";

	// GS
	cout<<"\nGS_resi\n";
	solution= GS_solver(A,b,tol,guess);
	cout<<solution.x<<endl<<solution.iter<<endl;
	cout<<prod(A,solution.x)<<std::endl;
	result.clear();
	for(long i=0;i<solution.x.size();i++)
			result.push_back(solution.x(i));
	result.push_back(solution.iter);
	result.push_back(norm_2(solution.resi));

	v2c.load("GS_resi.csv",result,result.size(),1);
	v2c.export2csv();
	cout<<"*******************\n";
	cout<<"\nGS_cons\n";
	solution= GS_solver(A,b,tol,guess,false);
	cout<<solution.x<<endl<<solution.iter<<endl;
	cout<<prod(A,solution.x)<<std::endl;

	result.clear();
	for(long i=0;i<solution.x.size();i++)
			result.push_back(solution.x(i));
	result.push_back(solution.iter);
	result.push_back(norm_2(solution.resi));

	v2c.load("GS_cons.csv",result,result.size(),1);
	v2c.export2csv();
	cout<<"*******************\n";

	// SOR
	cout<<"\nSOR_resi\n";
	double w=1.05;
	solution= SOR_solver(A,b,w,tol,guess);
	cout<<solution.x<<endl<<solution.iter<<endl;
	cout<<prod(A,solution.x)<<std::endl;

	result.clear();
	for(long i=0;i<solution.x.size();i++)
			result.push_back(solution.x(i));
	result.push_back(solution.iter);
	result.push_back(norm_2(solution.resi));


	v2c.load("SOR_resi.csv",result,result.size(),1);
	v2c.export2csv();
	cout<<"*******************\n";
	cout<<"\nSOR_cons\n";
	solution= SOR_solver(A,b,w,tol,guess,false);
	cout<<solution.x<<endl<<solution.iter<<endl;
	cout<<prod(A,solution.x)<<std::endl;

	result.clear();
	for(long i=0;i<solution.x.size();i++)
			result.push_back(solution.x(i));
	result.push_back(solution.iter);
	result.push_back(norm_2(solution.resi));

	v2c.load("SOR_cons.csv",result,result.size(),1);
	v2c.export2csv();
	cout<<"*******************\n";


/*
	//Repeat SOR
	double weight=1.02;

	result.clear();
	long r=0;
	while(weight<2)
	{

		solution=SOR_solver(A,b,weight,tol,guess);
		for(long i=0;i<solution.x.size();i++)
			result.push_back(solution.x(i));

		result.push_back(solution.iter);

		weight+=0.02;
		r++;
	}
	v2c.load("Repeat_SOR.csv", result,r,solution.x.size()+1);
	v2c.export2csv();


*/
	return 0;

}




