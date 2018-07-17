#ifndef STATALGORITHM_H_
#define STATALGORITHM_H_
#include "myutilities.h"
#include "require.h"
#include <cmath>
#include <limits>
#include <iostream>
/* Algorithm AS66 Applied Statistics (1973) vol.22, no.3 Evaluates the tail area 
 * of the standardised normal curve from x to infinity if upper is .true. or
 * from minus infinity to x if upper is .false.
 */
double alnorm(double x, bool upper);

/* ALGORITHM AS241  APPL. STATIST. (1988) VOL. 37, NO. 3
 * Produces the normal deviate Z corresponding to a given lower tail area of P; 
 * Z is accurate to about 1 part in 10**16.
 * The hash sums below are the sums of the mantissas of the coefficients.   
 * They are included for use in checking transcription.
 */
void  ppnd16(double p, double& normal_dev, int& ifault);

/* Algorithm AS 190.2  Appl. Statist. (1983) Vol.32, No.2
 * Calculates an initial quantile p for a studentized range distribution
 * having v degrees of freedom and r samples for probability p, 0.8 < p < 0.995
 * Uses function ppnd - Algorithm AS 241
 */
double qtrng0(double p, double v, double r);

/* Algorithm AS 190  Appl. Statist. (1983) Vol.32, No. 2
 * Evaluates the probability from 0 to q for a studentized range having v degrees of 
 * freedom and r samples.
 * Uses function ALNORM = algorithm AS66.
 * Arrays vw and qw store transient values used in the quadrature
 * summation.  Node spacing is controlled by step.  pcutj and pcutk control
 * truncation.  Minimum and maximum number of steps are controlled by
 * jmin, jmax, kmin and kmax.  Accuracy can be increased by use of a finer
 * grid - Increase sizes of arrays vw and qw, and jmin, jmax, kmin, kmax and
 * 1/step proportionally.
 */
double prtrng(double q, double v, double r);

/* Algorithm AS 190.1  Appl. Statist. (1983) Vol.32, No. 2
 * Approximates the quantile p for a studentized range distribution
 * having v degrees of freedom and r samples for probability 0.9 < p < 0.99.
 * Uses functions  alnorm, ppnd, prtrng and qtrng0 - Algorithms AS 66, AS 241, AS 190 and AS 190.2
 */
double qtrng(double p, double v, double r);

/* Inverse CDF of student-t distribution. Adapted from A guide to simulation (2nd ed) by Paul Bratley
 * pp. 340
 * Inputs: 	p = probability, df = degree of freedom
 */
double tinv(double p, int df);  

/*Inverse CDF of F
 * Adapted from A guide to simulation (2nd ed) by Paul Bratley
 * pp. 337.
 * Inputs:  dfn= d.f. of numerator, dfd = d.f. of denominator, phi=prob to be inverted
 * Return:  inverse of F cdf evaluated at phi
 */
double Finv(double phi, int dfn, int dfd);

// Inverse CDF of normal. Input: probability p		
double norminv(double p);

/* Inverse CDF of Beta. 
 * Adapted from A guide to simulation (2nd ed) by Paul Bratley
 * pp. 333.
 * Inputs:  phi = prob to be inverted, a, b = 2 shape parmater of the beta
 * Return: the inverse of Beta(a,b) at phi
 */
double Betainv(double phi, double a, double b);

/* t = number of systems under comparison
 * pstar = confidence level, so if user specify tail prob alpha, should do some conversion
 * nu = degree of freedom for the initial sample size, i.e. n0 - 1.
 * Rinott constant
 */
double Rinott(int t, double pstar, int nu);

// This program shows how chi-square pdf function is calculated.
double Chipdf(int N, double c);

// Compute ln(Gamma(xx))
double gammaln(double xx);

//Compute CDF of standard normal
double normcdf(double x);
#endif /*STATALGORITHM_H_*/
