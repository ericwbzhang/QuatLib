/*
 * DiscreteDivdBT.h
 *
 *  Created on: Nov 27, 2014
 *      Author: Eric
 */

#ifndef DISCRETEDIVDBT_H_
#define DISCRETEDIVDBT_H_
#include "BinomialTree.h"
#include "vector"



struct dividend{
	double t, divd;
	int divd_type; // 0 const; 1 proportional

	dividend(double time, double div, int type): t(time), divd(div), divd_type(type){};
};


class DiscreteDivdBT_Amer :public BinomialTree{
protected:
	std::vector<dividend> divd_sche;
public:
	DiscreteDivdBT_Amer();
	DiscreteDivdBT_Amer(basic_option opt, long n, std::vector<dividend> divd_schedule):BinomialTree(opt, n, 1){
		divd_sche=divd_schedule;
	};
	virtual ~DiscreteDivdBT_Amer();

	virtual void load();
};


class DiscreteDivdBT_Euro{
protected:
	BinomialTree *bt;
public:
	DiscreteDivdBT_Euro(){};
	DiscreteDivdBT_Euro(basic_option opt, long n, std::vector<dividend> divd_sche){

		for(std::vector<dividend>::iterator it=divd_sche.begin(); it!= divd_sche.end();it++)
			{
				if((*it).divd_type==0)// const divd
					opt.S-=(*it).divd*exp(-(*it).t*opt.r);
				else if ((*it).divd_type==1)// proportional divd
				{
					opt.S*=(1-(*it).divd);
				}

			}
		bt=new BinomialTree(opt, n, 0);
		bt->load();
	};
	virtual ~DiscreteDivdBT_Euro(){delete bt;};

	inline double price(){return bt->price();};
	inline double delta(){return bt->delta();};
	inline double gamma(){return bt->gamma();};
	inline double theta(){return bt->theta();};

};
#endif /* DISCRETEDIVDBT_H_ */
