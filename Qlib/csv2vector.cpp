/*
 * csv2vector.cpp
 *
 *  Created on: Aug 31, 2014
 *      Author: ericmac
 */


#include "csv2vector.h"
#include <stdlib.h>
#include <fstream>
#include <boost/numeric/ublas/matrix.hpp>

using namespace std;



csv2vector::~csv2vector() {
	// TODO Auto-generated destructor stub
}


csv2vector::csv2vector()
{}

void csv2vector::load(string s)
{
	vector<long> a;
	vector<long> b;
	a.clear(); b.clear();
	this->load(s,a,b);
}

void csv2vector::load(string s, vector<long> skip_r, vector<long> skip_c)
// string s is the file name like "example.csv". Actually this class works for various file types as long as
//the data is separated by comma.
// skip_r and skip_c is vectors that store the indexes of row and column that you want to skip
// skip_r and skip_c must be in ascending order and no multiplicity, no exceeding the range of data.
{
	skip_row.clear(); skip_col.clear();
	vec_data.clear(); vec_str.clear();
	r=0;c=0;

	skip_row=skip_r;
	skip_col=skip_c;
	ifstream myfile(s.c_str());
	string line;
	if (myfile.is_open())
	{
		while (getline(myfile,line,'\n'))
		{
			r++;
			stringstream sstream(line);
			string s;
			while (getline(sstream,s, ',' ))
			{
				vec_data.push_back(s);
			}
		}

		c=vec_data.size()/r;
		this-> parse();
		myfile.close();
	}
	else {cout<<"The file is not open!"<<endl;}
}


void csv2vector::parse()
{
	vec_str=vec_data;
	vector<string > ::iterator it_str=vec_str.begin();
	vector<long> row_index = skip_row;
	row_index.insert(row_index.begin(), 0);
	for(vector<long> ::iterator it_r=row_index.begin()+1; it_r!= row_index.end(); it_r++)
	{
		it_str=vec_str.erase(it_str+ c*(*it_r-*(it_r-1)-1),it_str+ c*(*it_r-*(it_r-1)-1)+c);
	}

	it_str=vec_str.begin();
	vector<long> col_index=skip_col;
	col_index.insert(col_index.begin(),0);
	while (it_str!= vec_str.end())
	{
		for( vector<long> ::iterator it_c=col_index.begin()+1; it_c!= col_index.end(); it_c++)
		{
			it_str=vec_str.erase(it_str+*it_c-*(it_c-1)-1);
		}
		it_str+=c-*(col_index.end()-1);
	}
}

vector<string> csv2vector::data_str()
// return the result in vector<string> type

{
	return vec_str;
}
boost::numeric::ublas::matrix<string> csv2vector::data_mtx_str()
{
	boost::numeric::ublas::matrix<string> A(r-skip_row.size(),c-skip_col.size());
	std::vector<string>::iterator it=vec_str.begin();
	for(long i=0;i<A.size1();i++)
		for(long j=0;j<A.size2();j++)
			A(i,j)=*(it++);
	return A;
}
boost::numeric::ublas::matrix<double> csv2vector::data_mtx_double()
{
	boost::numeric::ublas::matrix<double> A(r-skip_row.size(),c-skip_col.size());
		std::vector<string>::iterator it=vec_str.begin();
		for(long i=0;i<A.size1();i++)
			for(long j=0;j<A.size2();j++)
				{stringstream ss(*(it++)); double d; ss>>d;
				A(i,j)=d;
				}
		return A;
}


long csv2vector::row()
//return the row of the vector

{
	return r-skip_row.size();
}
long csv2vector::col()
//return the col of the vector

{
	return c-skip_col.size();
}

vector<double> csv2vector::data_double()
//return the vector in vector<double> type

{
	vector<double > vec_d;
	double d=0;
	for(vector<string>::iterator it=vec_str.begin(); it!=vec_str.end(); it++)
	{
		stringstream ss(*it);
		ss>>d;
		vec_d.push_back(d);
	}

	return vec_d;
}


vector<string> csv2vector::title(vector<long> t_index, bool rORc)
//t_index is a vector containing the row/col indexes that you want read as string.
//rORc =true if it is row index, false if the it is in column index.
// the return is vector<string>

{
	vector<string> title;
	if(rORc==true)
	{
		for(vector<long>::iterator it=t_index.begin(); it!= t_index.end(); it++)
		{
			title.insert(title.end(),vec_data.begin()+ (*it-1)*c, vec_data.begin()+ (*it-1)*c+c);
		}
	}else
	{
		vector<long> t=t_index;
		t.insert(t.begin(),1);
		vector<string > ::iterator it_data= vec_data.begin();
		while (it_data!= vec_data.end())
		{
			for (vector<long > ::iterator it_t= t.begin()+1; it_t!= t.end(); it_t++	)
			{
				it_data+= (*it_t-*(it_t-1));
				title.push_back(*it_data);
			}
			it_data+= c-*(t.end()-1)+1;
		}
	}
	return title;
}

boost::numeric::ublas::matrix<string> csv2vector::title_mtx(vector<long> t_index, bool rORc)
{
	std::vector<string> vec= this-> title(t_index, rORc);
	boost::numeric::ublas::matrix<string> A;
	if(rORc==true)
	{
		A(t_index.size(), c);
		std::vector<string>::iterator it=vec.begin();
			for(long i=0;i<A.size1();i++)
			for(long j=0;j<A.size2();j++)
				A(i,j)=*(it++);
	}
	else{
		A(r, t_index.size());
		std::vector<string>::iterator it=vec.begin();
			for(long i=0;i<A.size1();i++)
			for(long j=0;j<A.size2();j++)
				A(i,j)=*(it++);
	}
	return A;

}




