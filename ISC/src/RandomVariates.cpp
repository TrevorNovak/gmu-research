#include "masterheader.h"

RngStream rng("rng");

double Normal(double mean, double variance)
{
    double u1, u2, v1, v2, w, x, y;
    do
    {
        u1 = rng.RandU01();
        u2 = rng.RandU01();
        v1 = 2 * u1 - 1;
        v2 = 2 * u2 - 1;
        w = v1 * v1 + v2 * v2; 
    }
    while (w > 1);
    y = sqrt(-2 * log(w) / w);
    x = v1 * y;
    // x is N(0,1)  
    double normal = mean + sqrt(variance) * x;
    return normal;
}

double Normal(double mean, double variance, RngStream& indrng)
{
    double u1, u2, v1, v2, w, x, y;
    do
    {
        u1 = indrng.RandU01();
        u2 = indrng.RandU01();
        v1 = 2 * u1 - 1;
        v2 = 2 * u2 - 1;
        w = v1 * v1 + v2 * v2; 
    }
    while (w > 1);
    y = sqrt(-2 * log(w) / w);
    x = v1 * y;
    // x is N(0,1)  
    double normal = mean + sqrt(variance) * x;
    return normal;
}

double TruncatedNormal(double mean, double variance, double lowerBound, double upperBound)
{
	double x;
	do
	{
		x = Normal(mean, variance);
	} while (upperBound < x || x < lowerBound);
	return x;
}

double Exponential(double arrivalRate)
{
	double u = rng.RandU01();
	return -log(u) / arrivalRate;
}

// Gumbel  CDF:  exp(-exp(z/(mu + gamma))), where mu is parameter, gamma is Euler's constant,
// gamma = 0.5772. Mean 0, variance mu^2 * pi /6
double Gumbel(double mu)
{
	const double gamma = 0.5772;
	// inverse of CDF is  -(mu + gamma) * log(-log(u))
	double u = rng.RandU01();
	return -(mu + gamma) * log(-log(u));
}

double Gumbel(double mu, RngStream& indrng)
{
	const double gamma = 0.5772;
	// inverse of CDF is  -(mu + gamma) * log(-log(u))
	double u = indrng.RandU01();
	return -(mu + gamma) * log(-log(u));
}


// Extreme value, F(x) = exp(-exp(-x))
double ExtremeValue()
{
	double u = rng.RandU01();
	return -log(-log(u));
}

double ExtremeValue(RngStream& indrng)
{
	double u = indrng.RandU01();
	return -log(-log(u));
}

// Law and Kelton, p478
int Poisson(double lambda)
{
	if (lambda > 100)
	{ 	// Do a normal approximation
		double temp = Normal(lambda, lambda);
		return (temp < 0) ? 0 : int(temp);
	}
	double a = exp(-lambda);
	double b = 1;
	int i = 0;
	while(true)
	{
		b *= rng.RandU01();
		if (b < a)
		{
			return i;
		}
		i++;
	}
}

int Poisson(double lambda, RngStream& indrng)
{
	if (lambda > 100)
	{ 	// Do a normal approximation
		double temp = Normal(lambda, lambda, indrng);
		return (temp < 0) ? 0 : int(temp);
	}
	double a = exp(-lambda);
	double b = 1;
	int i = 0;
	while(true)
	{
		b *= indrng.RandU01();
		if (b < a)
		{
			return i;
		}
		i++;
	}
}

int URandIntGen::operator()()
{
	while(true) 
	{	// Note that RandInt returns integers from 0 to limit inclusively.
		int i = rng.RandInt(0, limit);
	    if(used.find(i) == used.end()) 
	    {
   			used.insert(i);
    	    return i; 
   		}
	}
}

 



