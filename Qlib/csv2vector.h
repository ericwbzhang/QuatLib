/*
 * csv2vector.h
 *
 *  Created on: Aug 31, 2014
 *      Author: ericmac
 */

#ifndef CSV2VECTOR_H_
#define CSV2VECTOR_H_
#include <iostream>
#include <sstream>
#include <vector>
#include <boost/numeric/ublas/matrix.hpp>

using namespace std;

class csv2vector {
private:
	long   r=0;
	long   c=0;
	vector<long> skip_row;
	vector<long> skip_col;
	vector<string> vec_data;
	vector<string> vec_str;
	void parse();
public:
	csv2vector();
	virtual ~csv2vector();
	long int row();
	long int col();
	vector<string> data_str();
	vector<double> data_double();
	boost::numeric::ublas::matrix <string> data_mtx_str();
	boost::numeric::ublas::matrix<double> data_mtx_double();
	void load(string s, vector<long > skip_r, vector<long> skip_c);
	void load(string s);
	//string s is filename. skip_r is a vector which contains the index of rows that should be skipped. skip_c is the a vector which contains the index of cols that should be skipped.
	// skip_r and skip_c must be in ascending order and no multiplicity, no exceeding the range of data.

	vector<string> title(vector<long> t_index, bool rORc=true);
	//it returns the title, which is specified t_index that contains the index of title. rORc =true if the title is in row, false if the title is in column. In any case, the title vector is ordered by row. NEVER by COLUMN
	boost::numeric::ublas::matrix<string> title_mtx(vector<long> t_index, bool rORc=true);
};

#endif /* CSV2VECTOR_H_ */
