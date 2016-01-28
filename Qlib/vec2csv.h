/*
 * vec2csv.h
 *
 *  Created on: Sep 1, 2014
 *      Author: ericmac
 */

#ifndef VEC2CSV_H_
#define VEC2CSV_H_
#include "vector"
#include <boost/numeric/ublas/matrix.hpp>
#include<boost/numeric/ublas/vector.hpp>

using namespace std;

class vec2csv {
private:
	vector<string> title_h, title_v;
	long tit_h_r=0; long tit_h_c=0;
	long tit_v_r=0; long tit_v_c=0;
	vector<string> content;
	long r=0; long c=0;
	string filename;
public:
	//The following constructors set that the content, which is input as matrix or vector,
	//is exported to csv by row x col. string s is the filename.
	vec2csv();
	void load(string s, boost::numeric::ublas::vector<double> vec,long row, long col);
	void load(string s, boost::numeric::ublas::matrix<string> m);
	void load(string s, boost::numeric::ublas::matrix<int> m);
	void load(string s, boost::numeric::ublas::matrix<double> m);
	void load(string s, vector<int> vec, long row, long col);
	void load(string s ,vector<double> vec, long row, long col);
	void load(string s, vector<string> vec, long row, long col);
	void add_title(vector<string> s, long t_row, long t_col, bool rORc);
	// add_title gives the opportunity to add horizontal or vertical title.
	//rORc is true if the title is horizontal, false if it is vertical.
	//t_row and t_col gives the number of rows and cols of the vector<string> s.
	//The title must have a size which can match content.

	void export2csv();
	//row dominated export.i.e. we first export horizontal title, by row, then export the vertical title and content, again, row by row.
	virtual ~vec2csv();
};

#endif /* VEC2CSV_H_ */

