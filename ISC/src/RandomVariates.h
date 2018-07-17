#ifndef RANDOMVARIATES_H_
#define RANDOMVARIATES_H_
#include "RngStream.h"
#include <set>

double Normal(double mean, double variance);
double Normal(double mean, double variance, RngStream& indrng);
double Exponential(double arrivalRate);
double TruncatedNormal(double mean, double variance, double lowerbound, double upperbound);

// Gumbel  CDF:  exp(-exp(z/(mu + gamma))), where mu is parameter, gamma is Euler's constant,
// gamma = 0.5772. Mean 0, variance mu^2 * pi /6
double Gumbel(double mu);
double Gumbel(double mu, RngStream& indrng);

double ExtremeValue();
double ExtremeValue(RngStream& indrng);

int Poisson(double lammbda);
int Poisson(double lammbda, RngStream& indrng);

/* Generate unique random numbers from 0 to lim. To be used with a generate_n algorithm
 * Example: 
 * 
 * vector<int> x(SZ);
 * URandIntGen urg(MAX);
 * generate_n(x.begin(), SZ, urg);
 */

class URandIntGen 
{
	private :
		std::set<int> used;
  		int limit;
	public :
		URandIntGen(int lim) : limit(lim) {}
  		int operator()();
};
#endif /*RANDOMVARIATES_H_*/
