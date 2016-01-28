/*
 * vec2csv.cpp
 *
 *  Created on: Sep 1, 2014
 *      Author: ericmac
 */

#include "vec2csv.h"
#include <string>
#include <sstream>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <fstream>


vec2csv::vec2csv()
{}

void vec2csv::load(string s, boost::numeric::ublas::matrix<string> m) {
	content.clear();
	filename =s;
	r=m.size1();c=m.size2();
	for(long i= 0; i<r; i++)
		for(long j=0; j<c; j++)
		{
			content.push_back(m(i,j));
		}
}

void vec2csv::load(string s, boost::numeric::ublas::matrix<double> m) {
	content.clear();
	filename =s;
	r=m.size1();c=m.size2();
	stringstream sstream;
	for(long i= 0; i<r; i++)
		for(long j=0; j<c; j++)
		{
			stringstream sstream;
			sstream<<m(i,j);
			content.push_back(sstream.str());
		}
}

void vec2csv::load(string s, boost::numeric::ublas::matrix<int> m) {
	content.clear();
	filename =s;
	r=m.size1();c=m.size2();
	stringstream sstream;
	for(long i= 0; i<r; i++)
		for(long j=0; j<c; j++)
		{
			stringstream sstream;
			sstream<<m(i,j);
			content.push_back(sstream.str());
		}
}

void vec2csv::load(string s, vector<int> vec, long row, long col)
{
	r=row; c=col;
	content.clear();
	filename=s;
	for(vector<int>::iterator it=vec.begin(); it!= vec.end(); it++)
	{
		stringstream sstream;
		sstream<<*it;
		content.push_back(sstream.str());
	}
}
void vec2csv::load(string s ,vector<double> vec, long row, long col){
	r=row; c=col;
	content.clear();
	filename=s;
	for(vector<double>::iterator it=vec.begin(); it!= vec.end(); it++)
	{
			stringstream sstream;
			sstream<<*it;
			content.push_back(sstream.str());
	}

}
void vec2csv::load(string s ,vector<string> vec, long row, long col){
	r=row; c=col;
	content.clear();
	filename=s;
	content=vec;
}

vec2csv::~vec2csv() {
	// TODO Auto-generated destructor stub
}

void vec2csv::add_title(vector<string> s, long t_row, long t_col, bool rORc)
{
	if(rORc==true)
	{
		title_h=s;
		tit_h_r=t_row;
		tit_h_c=t_col;
	}else
	{
		title_v=s;
		tit_v_r=t_row;
		tit_v_c=t_col;
	}
}

void vec2csv::export2csv()
{
	ofstream myfile(filename.c_str());
	if (myfile.is_open())
	{
		long i=1;
		for(vector<string> ::iterator it=title_h.begin(); it!= title_h.end(); it++)
		{
			myfile<<*it;
			if(div(i,tit_h_c).rem==0) myfile<<"\n";
			else myfile<<",";
			i++;
		}
		i=1;
		vector<string>::iterator it_titv= title_v.begin();
		for(vector<string> ::iterator it=content.begin(); it!= content.end(); it++)
		{
			if(div(i-1, c+ tit_v_c).rem==0)
			{
				for(long k=0; k< tit_v_c; k++)
				{
					it_titv+=k;
					myfile<<*(it_titv)<<",";
				}
				i+=tit_v_c;
			}
			myfile<<*it;
			if(div(i,c+tit_v_c).rem==0) myfile<<"\n";
			else myfile<<",";
			i++;
		}
		myfile.close();
	}else cout<< "Error! The file is closed.\n";
}

void vec2csv::load(string s, boost::numeric::ublas::vector<double> vec,long row, long col)
{
	std::vector<double> v;
	for(long i=0; i<vec.size(); i++)
		v.push_back(vec(i));
	this->load(s,v,row,col);
}


