

#ifndef BS_H_
#define BS_H_

#include "basic_option.h"

class BS {
private:
	basic_option opt;
	double opt_d1;
	double opt_d2;
	double opt_price;

public:
	BS();
	BS(basic_option o);
	virtual ~BS();

	double price();
	double delta();
	double gamma();
	double vega();
	double theta();
	double rho();
	double psi();
	double omega();
	double d1();
	double d2();

	double option_S() {return opt.S;}
	double option_K() {return opt.K;}
	double option_T() {return opt.T;}
	double option_sigma() {return opt.sigma;}
	double option_r() {return opt.r;}
	double option_q() {return opt.q;}

	double option_type() {return opt.type;}


};

#endif /* BS_H_ */
