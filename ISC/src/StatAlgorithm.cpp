#include "StatAlgorithm.h"
using namespace std;
/* Algorithm AS66 Applied Statistics (1973) vol.22, no.3 Evaluates the tail area 
 * of the standardised normal curve from x to infinity if upper is .true. or
 * from minus infinity to x if upper is .false.
 */
double alnorm(double x, bool upper)
{
	const double con = 1.28, ltone = 7.0, utzero = 18.66;
	const double p = 0.398942280444, q = 0.39990348504,   \
                        r = 0.398942280385, a1 = 5.75885480458,  \
                        a2 = 2.62433121679, a3 = 5.92885724438,  \
                        b1 = -29.8213557807, b2 = 48.6959930692, \
                        c1 = -3.8052E-8, c2 = 3.98064794E-4,     \
                        c3 = -0.151679116635, c4 = 4.8385912808, \
                        c5 = 0.742380924027, c6 = 3.99019417011, \
                        d1 = 1.00000615302, d2 = 1.98615381364,  \
                        d3 = 5.29330324926, d4 = -15.1508972451, \
                        d5 = 30.789933034;
	double z, y, fn_val;
	bool up = upper;
	z = x;
	if (z < 0)
	{
		up = !up;
		z = -z;
	}
	if (z <= ltone || (up && z <= utzero))
	{
		y = 0.5 * z * z;
		if (con < z)
		{
			fn_val = r * exp(-y) / (z+c1+d1/(z+c2+d2/(z+c3+d3/(z+c4+d4/(z+c5+d5/(z+c6))))));
		}
		else
		{
			fn_val = 0.5 - z * (p - q*y/(y+a1+b1/(y+a2+b2/(y+a3))));
		}
	}
	else
	{
		fn_val = 0;
	}
	if(!up)
	{
		fn_val = 1 - fn_val; 
	}
	return fn_val;
}	


/* ALGORITHM AS241  APPL. STATIST. (1988) VOL. 37, NO. 3
 * Produces the normal deviate Z corresponding to a given lower tail area of P; 
 * Z is accurate to about 1 part in 10**16.
 * The hash sums below are the sums of the mantissas of the coefficients.   
 * They are included for use in checking transcription.
 */
void  ppnd16(double p, double& normal_dev, int& ifault)
{
	const double split1 = 0.425, split2 = 5, const1 = 0.180625, const2 = 1.6;
	double q, r;
	// Coefficients for P close to 0.5
	const double a0 = 3.3871328727963666080e0, \
             a1 = 1.3314166789178437745e+2, \
             a2 = 1.9715909503065514427e+3, \
             a3 = 1.3731693765509461125e+4, \
             a4 = 4.5921953931549871457e+4, \
             a5 = 6.7265770927008700853e+4, \
             a6 = 3.3430575583588128105e+4, \
             a7 = 2.5090809287301226727e+3, \
             b1 = 4.2313330701600911252e+1, \
             b2 = 6.8718700749205790830e+2, \
             b3 = 5.3941960214247511077e+3, \
             b4 = 2.1213794301586595867e+4, \
             b5 = 3.9307895800092710610e+4, \
             b6 = 2.8729085735721942674e+4, \
             b7 = 5.2264952788528545610e+3;
    // HASH SUM AB    55.8831928806149014439
    // Coefficients for P not close to 0, 0.5 or 1.
    const double c0 = 1.42343711074968357734e0, \
             c1 = 4.63033784615654529590e0, \
             c2 = 5.76949722146069140550e0, \
             c3 = 3.64784832476320460504e0, \
             c4 = 1.27045825245236838258e0, \
             c5 = 2.41780725177450611770e-1, \
             c6 = 2.27238449892691845833e-2, \
             c7 = 7.74545014278341407640e-4, \
             d1 = 2.05319162663775882187e0, \
             d2 = 1.67638483018380384940e0, \
             d3 = 6.89767334985100004550e-1, \
             d4 = 1.48103976427480074590e-1, \
             d5 = 1.51986665636164571966e-2, \
             d6 = 5.47593808499534494600e-4, \
             d7 = 1.05075007164441684324e-9;
	// HASH SUM CD    49.33206503301610289036
	// Coefficients for P near 0 or 1.
	const double e0 = 6.65790464350110377720e0, \
             e1 = 5.46378491116411436990e0, \
             e2 = 1.78482653991729133580e0, \
             e3 = 2.96560571828504891230e-1, \
             e4 = 2.65321895265761230930e-2, \
             e5 = 1.24266094738807843860e-3, \
             e6 = 2.71155556874348757815e-5, \
             e7 = 2.01033439929228813265e-7, \
             f1 = 5.99832206555887937690e-1, \
             f2 = 1.36929880922735805310e-1, \
             f3 = 1.48753612908506148525e-2, \
             f4 = 7.86869131145613259100e-4, \
             f5 = 1.84631831751005468180e-5, \
             f6 = 1.42151175831644588870e-7, \
             f7 = 2.04426310338993978564e-15;
	// HASH SUM EF    47.52583317549289671629
	ifault = 0;
	q = p - 0.5;
	if (dabs(q) <= split1)
	{
		r = const1 - q * q;
		normal_dev = q * (((((((a7*r + a6)*r + a5)*r + a4)*r + a3)*r + a2)*r + a1)*r + a0) / \
           (((((((b7*r + b6)*r + b5)*r + b4)*r + b3)*r + b2)*r + b1)*r + 1);
	}
	else
	{
		r = (q < 0) ? p : (1 - p);
		if (r < 0)
		{
			ifault = 1;
			normal_dev = 0;
		}
		else
		{
			r = sqrt(-log(r));
			if (r <= split2)
			{
				r -= const2;
				normal_dev = (((((((c7*r + c6)*r + c5)*r + c4)*r + c3)*r + c2)*r + c1)*r + c0) / \
				             (((((((d7*r + d6)*r + d5)*r + d4)*r + d3)*r + d2)*r + d1)*r + 1);
			}
			else
			{
				r -= split2;
				normal_dev = (((((((e7*r + e6)*r + e5)*r + e4)*r + e3)*r + e2)*r + e1)*r + e0) / \
             				 (((((((f7*r + f6)*r + f5)*r + f4)*r + f3)*r + f2)*r + f1)*r + 1);
			}
			if (q < 0)
			{
				normal_dev = - normal_dev;
			}
		}
	}
}

/* Algorithm AS 190.2  Appl. Statist. (1983) Vol.32, No.2
 * Calculates an initial quantile p for a studentized range distribution
 * having v degrees of freedom and r samples for probability p, 0.8 < p < 0.995
 * Uses function ppnd - Algorithm AS 241
 */
double qtrng0(double p, double v, double r)
{
	const double vmax = 120.0, c1 = 0.8843, c2 = 0.2368, c3 = 1.214, c4 = 1.208, c5 = 1.4142;
	int ifault;
	double q, t;
	// t and ifault will have the returned values from function ppnd16
	ppnd16(0.5 + 0.5 * p, t, ifault);
	if (v < vmax)
	{
		t += (t * t * t + t) / v / 4.0;
	}
	q = c1 - c2 * t;
	if (v < vmax)
	{
		q = q - c3 / v + c4 * t / v;
	}
	return t * (q * log(r - 1.0) + c5);
}

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
double prtrng(double q, double v, double r)
{
	const double pcutj = 0.00003, pcutk = 0.0001, step = 0.45, vmax = 120.0, zero = 0.0, \
			fifth = 0.2, half = 0.5, one = 1.0, two = 2.0, cv1 = 0.193064705, cv2 = 0.293525326, \
			cvmax = 0.39894228;
	const double cv[4] = {0.318309886, -0.268132716e-2, 0.347222222e-2, 0.833333333e-1};
	const int jmin = 3, jmax = 15, kmin = 7, kmax = 15;
	double vw[30], qw[30];
	double g, gmid, r1, c, h, v2, gstep, pk1, pk2, gk, pk;
	double w0, pz, x, hj, ehj, pj;
	int ifault, j, jj, jump, k;
	double prob = zero;
	ifault = 0;
	// check if parametes are invalid, v is d.f. r is # of samples.
	if (v < one || r < two) ifault = 1;
	if (!(q <= zero || ifault == 1))
	{
		// Computing constants, local midpoint, adjusting steps.
		g = step * pow(r, (-fifth));
		gmid = half * log(r);
		r1 = r - one;
		c = log(r * g * cvmax);
		if (v <= vmax)
		{
			h = step * pow(v, (-half));
			v2 = v * half;
			if (v == one) 
			{
				c = cv1;
			}
			else if (v == two) 
			{
				c = cv2;
			}
			else
			{
				c = sqrt(v2) * cv[0] / (one + ((cv[1] / v2 + cv[2]) / v2 + cv[3]) / v2);
			}
			c = log(c * r * g * h);
		}
		/*  Computing integral
		 *  Given a row k, the procedure starts at the midpoint and works
		 *  outward (index j) in calculating the probability at nodes
		 *  symmetric about the midpoint.  The rows (index k) are also
		 *  processed outwards symmetrically about the midpoint.  The centre row is unpaired.
		 */
		gstep = g;
		qw[0] = -one;
		qw[jmax] = -one;
		pk1 = one;
		pk2 = one;
		for (k = 1; k <= kmax; k++)
		{
			gstep = gstep - g;
			do
			{
				gstep = -gstep;
				gk = gmid + gstep;
				pk = zero;
				if (!(pk2 <= pcutk && k > kmin))
				{
					w0 = c - gk * gk * half;
					pz = alnorm(gk, true);
					x = alnorm(gk - q, true) - pz; 
					if (x > zero) pk = exp(w0 + r1 * log(x));
					if (v <= vmax)
					{
						jump = -jmax;
						do
						{
							jump = jump + jmax;
							for (j = 1; j <= jmax; j++)
							{
								jj = j + jump;
								if (qw[jj-1] <= zero)  
								{
									hj = h * j;
									if (j < jmax) qw[jj] = -one;
									ehj = exp(hj);
									qw[jj-1] = q * ehj;
									vw[jj-1] = v * (hj + half - ehj * ehj * half);
								}
								pj = zero;
								x = alnorm(gk - qw[jj-1], true) - pz;
								if (x > zero) pj = exp(w0 + vw[jj-1] + r1 * log(x));
								pk = pk + pj;
								if (pj > pcutj) continue;
								if (jj > jmin || k > kmin) break;
							}
							h = -h;
						} while(h < zero);
					}
				}
				prob += pk;
				if (k > kmin && pk <= pcutk && pk1 <= pcutk) 
				{
					if (ifault != 0) cout << "Error, function prtrng(q,v,r) " << endl;
					return prob;
				}
				pk2 = pk1;
				pk1 = pk;
			}while(gstep > zero);
		}
	}
	if (ifault != 0) cout << "Error, function prtrng(q,v,r) " << endl;
	return prob;	
}			


/* Algorithm AS 190.1  Appl. Statist. (1983) Vol.32, No. 2
 * Approximates the quantile p for a studentized range distribution
 * having v degrees of freedom and r samples for probability 0.9 < p < 0.99.
 * Uses functions  alnorm, ppnd, prtrng and qtrng0 - Algorithms AS 66, AS 241, AS 190 and AS 190.2
 */
double qtrng(double p, double v, double r)
{
	const int jmax = 8;
	const double pcut = 0.001, p75 = 0.75, p80 = 0.80, p90 = 0.9, p99 = 0.99, p995 = 0.995,   \
                         p175 = 1.75, one = 1.0, two = 2.0, five = 5.0, eps = 1.0E-04;
	double quantile;
	double q1, p1, q2, p2, d, e1, e2;
	int ifault, j, nfault;
	ifault = 0;
	nfault = 0;
	// Check if parameters are valid
	if (v < one || r < two) ifault = 1;
	if (p < p90 || p > p99) ifault = 2;
	if (ifault != 0) 
	{
		cout << "Invalid arguments for function qtrng(p,v,r)" << endl;
		return quantile = -1;
	}
	// Obtain initial values
	q1 = qtrng0(p, v, r);
	p1 = prtrng(q1, v, r);
	if (nfault != 0) 
	{
		cout << "Error in function qtrng(p,v,r)" << endl;
		return quantile = -1;
	}
	quantile = q1;
	if (dabs(p1-p) >= pcut)
	{
		if (p1 > p) p1 = p175 * p - p75 * p1;
		if (p1 < p) p2 = p + (p - p1) * (one - p) / (one - p1) * p75;
		if (p2 < p80) p2 = p80;
		if (p2 > p995) p2 = p995;
		q2 = qtrng0(p2, v, r);
		if (nfault != 0) 
		{
			cout << "Error in function qtrng(p,v,r)" << endl;
			return quantile = -1;
		}
		// Refine approximation
		for (j = 2; j <= jmax; j++)
		{
			p2 = prtrng(q2, v, r);
			if (nfault != 0) 
			{
				cout << "Error in function qtrng(p,v,r)" << endl;
				return quantile = -1;
			}
			e1 = p1 - p;
			e2 = p2 - p;
			quantile = (q1 + q2) / two;
			d = e2 - e1;
			if (dabs(d) > eps) quantile = (e2 * q1 - e1 * q2) / d;
			if(dabs(e1) >= dabs(e2)) 
			{
				q1 = q2;
  				p1 = p2;
			}
			if (dabs(p1 - p) < pcut * five) return quantile;
			q2 = quantile;
		}
	}		
	return quantile;
}
		
// Inverse CDF of normal. Input: probability p		
// Accuracy: about 5 decimial points
double norminv(double p)
{
	const double PLIM = 1.0E-18, P0 = -0.322232431088E0, P1 = -1.0, P2 = -0.342242088547E0;
	const double P3 = -0.0204231210245E0, P4 = -0.453642210148E-4, Q0 = 0.099348462606E0, \
		Q1 = 0.588581570495E0, Q2 = 0.531103462366E0, Q3 = 0.10353775285E0, Q4 = 0.38560700634E-2;
	double vtemp;
	require(0 < p && p < 1, "Invalid probability for normal-inverse function");
	if (0.5 < p) p = 1 - p;
	if (PLIM <= p) 
	{
		double y = sqrt(-log(p*p));
		vtemp = y + ((((y*P4 + P3)*y + P2)*y + P1)*y +P0)/((((y*Q4 + Q3)*y + Q2)*y + Q1)*y + Q0);
	}
	else
	{
		// This is as far out in the tails as we go
		vtemp = 8;
	}
	if (0.5 <= p) vtemp = -vtemp;
	return vtemp; 
}
/* Inverse CDF of student-t distribution. Adapted from A guide to simulation (2nd ed) by Paul Bratley
 * pp. 340
 * Inputs: 	p = probability, df = degree of freedom
 * Accuracy in decimal digits:  5 if d.f. >= 8, machine precision if df = 1 or 2, 3 otherwise
 */
double tinv(double p, int df)
{
	const double DWARF = 0.1E-30, PIOVR2 = 1.5707963268E0;
	double t;
	double w[4];
	require(0 < p && p < 1, "Invalid probability for t-inverse function");
	require(df >= 0, "Invalid DF for t-inverse function");
	if (df == 0)
	{ // This may happen with a deterministic simulator, then simply return 0
		return 0;
	}
	// if df = 1, it is a Cauchy, use special method (exact)
	if (df == 1)
	{
		double arg = 3.1415926 * (1 - p);
		double sinarg = (DWARF < sin(arg)) ? sin(arg) : DWARF;
		return cos(arg) / sinarg;
	}
	// if df ==2, use a special method
	if (df == 2)
	{
		t = (p <= 0.5) ? (2 * p) :  (2 * (1 - p));
		if (t <= DWARF) t = DWARF;
		double vtemp = sqrt((2/(t*(2 - t))) - 2);
		if (p <= 0.5) vtemp = -vtemp;
		return vtemp;
	}	
	// Check to see if we are too far out in the tails
	if (DWARF < p && DWARF < (1 -p))
	{// General case
		double px2 = 2 * p;
		double p2tail = (px2 < (2 - px2)) ? px2 : 2 - px2;
		//begin the approximiation. 
		w[1] = 1.0 / (df - 0.5);
		w[0] = 48 / (w[1]*w[1]);
		w[2] = ((20700*w[1]/w[0] - 98)*w[1] -16)*w[1] + 96.36;
		w[3] = ((94.5/(w[0] + w[2]) - 3)/w[0] + 1)*sqrt(w[1]*PIOVR2)*df;
		double x = w[3] * p2tail;
		double c = w[2];
		double y = pow(x, 2.0/double(df));
		if(y <= (0.05 + w[1])) 
		{
			double y2 = (df + 6)/double(df + y) - 0.089*w[3] - 0.822;
			y2 *= (df + 2) * 3;
			y2 = (1/y2 + .5/(df + .4))*y - 1;
			y2 *= double(df + 1)/(df + 2);
			y = y2 + 1/y;
		}
		else
		{
			//asymptotic inverse expansion about normal
			x = norminv(p2tail * .5);
			y = x * x;
			if (df < 5) c += 0.3 * (df - 4.5) * (x + .6);
			c += (((.05*w[3]*x - 5)*x - 7)*x - 2)*x + w[0];
			y = (((((.4*y + 6.3)*y + 36)*y + 94.5)/c - y - 3)/w[0] + 1)*x;
			y = w[1] * y * y;
			if (y <= .002) 
			{
				y = .5 * y * y + y;
			}
			else
			{
				y = exp(y) - 1;
			}
		}
		t = sqrt(df * y);
	}
	else
	{// Tails
		t = 10E30;
	}
	if (p < 0.5) t = -t;
	return t;
}  		


/*Inverse CDF of F
 * Adapted from A guide to simulation (2nd ed) by Paul Bratley
 * pp. 337.
 * Inputs:  dfn= d.f. of numerator, dfd = d.f. of denominator, phi=prob to be inverted
 * Return:  inverse of F cdf evaluated at phi
 * Accuracy: except for dfn or dfd =1, in which case is the same accuracy as tinv, 
 * accuracy is the same as Betainv, which is about 1 decimal point.
 */
double Finv(double phi, int dfn, int dfd)
{
	double dwarf = numeric_limits<double>::epsilon();
	double fv;
	require(0 <= phi && phi <= 1, "Invalid p for Finv");
	require(0 < dfn && 0 < dfd, "Invalid d.f. for Finv");
	if (dfn == 1)
	{
		fv = tinv((1 + phi)/2, dfd);
		return fv * fv;
	}		
	else if (dfd == 1)
	{
		fv = tinv((2 - phi)/2, dfn);
		if (fv < dwarf) fv = dwarf;
		return 1 / (fv*fv);
	}
	else
	{
		fv = 1 - Betainv(phi, dfn/2.0, dfd/2.0);
		if (fv < dwarf) fv = dwarf;
		return double(dfd) / dfn * (1 - fv) / fv;
	}
}

/* Inverse CDF of Beta. 
 * Adapted from A guide to simulation (2nd ed) by Paul Bratley
 * pp. 333.
 * Inputs:  phi = prob to be inverted, a, b = 2 shape parmater of the beta
 * Return: the inverse of Beta(a,b) at phi
 * Accuracy: about 1 decimal point, except a or b = 1, in which case = machine precision
 */
double Betainv(double phi, double a, double b)
{
	require(0 <= phi && phi <= 1, "Invalid p for Betainv");
	double dwarf = 1e-10;
	if (phi < dwarf) return 0;
	if (1 - dwarf < phi) return 1;
	if (b == 1) return pow(phi, 1/a);
	if (a == 1) return 1 - pow(1 - phi, 1/b);
	double yp = -norminv(phi);
	double gl = (yp * yp - 3) / 6;
	double ad = 1 / (2 * a - 1);
	double bd = 1 / (2 * b - 1);
	double h = 2 / (ad + bd);
	double w = (yp * sqrt(h + gl) / h) - (bd - ad) * (gl + 0.83333333 - 0.66666666 / h);
	return a / (a + b * exp(2 * w));
}

/* t = number of systems under comparison
 * pstar = confidence level. so if user specify tail prob alpha, should do some conversion
 * nu = degree of freedom for the initial sample size, i.e. n0 - 1.
 * Rinott constant
 */
double Rinott(int t, double pstar, int nu)
{
	int i, j;
	double temp, ans;
	double h, lowerh, upperh;
	double X[32];
	double w[32];
	double WEX[32];
	X[0] = 0.044489365833267; X[1] = 0.234526109519619; X[2] = 0.576884629301886; X[3] = 1.07244875381782;
    X[4] = 1.72240877644465; X[5] = 2.52833670642579; X[6] = 3.49221327302199; X[7] = 4.61645676974977;
    X[8] = 5.90395850417424; X[9] = 7.35812673318624; X[10] = 8.9829409242126; X[11] = 10.78301863254;
    X[12] = 12.7636979867427; X[13] = 14.9311397555226; X[14] = 17.2924543367153; X[15] = 19.8558609403361;
    X[16] = 22.6308890131968; X[17] = 25.6286360224592; X[18] = 28.8621018163235; X[19] = 32.3466291539647;
    X[20] = 36.100494805752; X[21] = 40.1457197715394; X[22] = 44.5092079957549; X[23] = 49.2243949873086;
    X[24] = 54.3337213333969; X[25] = 59.892509162134; X[26] = 65.975377287935; X[27] = 72.6876280906627;
    X[28] = 80.1874469779135; X[29] = 88.7353404178924; X[30] = 98.829542868284; X[31] = 111.751398097938;
	w[0] = 0.109218341952385;  w[1] = 0.210443107938813; w[2] = 0.235213229669848; w[3]= 0.195903335972881;
    w[4] = 0.129983786286072; w[5] = 7.05786238657174E-02; w[6] = 3.17609125091751E-02; w[7] = 1.19182148348386E-02;
    w[8] = 3.73881629461152E-03; w[9] = 9.80803306614955E-04; w[10] = 2.14864918801364E-04; w[11] = 3.92034196798795E-05;
    w[12] = 5.93454161286863E-06; w[13] = 7.41640457866755E-07; w[14] = 7.60456787912078E-08; w[15] = 6.35060222662581E-09;
    w[16] = 4.28138297104093E-10; w[17] = 2.30589949189134E-11; w[18] = 9.79937928872709E-13; w[19] = 3.23780165772927E-14;
    w[20] = 8.17182344342072E-16; w[21] = 1.54213383339382E-17; w[22] = 2.11979229016362E-19;
    w[23] = 2.05442967378805E-21; w[24] = 1.3469825866374E-23; w[25] = 5.66129413039736E-26;
    w[26] = 1.41856054546304E-28; w[27] = 1.91337549445422E-31; w[28] = 1.19224876009822E-34;
    w[29] = 2.67151121924014E-38; w[30] = 1.33861694210626E-42; w[31] = 4.51053619389897E-48;
    for (i = 0; i < 32; i++)
	{
		WEX[i] = w[i] * exp(X[i]);
	}
	// This code cannot deal with d.f. > 93, in this case, use a closed form for infinite d.f.
	if (60 < nu)
	{
		h = sqrt(2.0) * norminv(pstar);
		return h;
	}
    h = 4;
    lowerh = 0;
    upperh = 20;
    do 
    {
    	ans = 0;
    	for (j = 0; j < 32; j++)
    	{
    		temp = 0;
    		for (i = 0; i < 32; i++)
    		{
    			temp = temp + WEX[i] * normcdf(h / sqrt(nu * (1 / X[i] + 1 / X[j]))) * Chipdf(nu, X[i]);
    		}
    		temp = pow(temp, (t - 1));
    		ans = ans + WEX[j] * temp * Chipdf(nu, X[j]);
    	}
    	if (ans > pstar + 0.000001)
    	{
    		upperh = h;
            h = (lowerh + upperh) / 2;
    	}
    	else if (pstar > ans + 0.000001)
    	{
    		lowerh = h;
    		h = (lowerh + upperh) / 2;
    	}
	} while(0.000001 <= dabs(ans - pstar));
    return h;
}

// This program shows how chi-square pdf function is calculated.
double Chipdf(int N, double c)
{
	double temp;
	temp = -N / 2.0 * log(2.0) - gammaln(N / 2.0) + (N / 2.0 - 1) * log(c) - c / 2;
    return exp(temp);
}

// Compute ln(Gamma(x))
double gammaln(double xx)
{
	require(0 < xx, "Invalid parameter for lngamma");
	double x,y,tmp,ser;
	static double cof[6]={76.18009172947146,-86.50532032941677, 24.01409824083091,-1.231739572450155, 0.1208650973866179e-2,-0.5395239384953e-5};
	int j;
	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<=5;j++) ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
}

//Compute CDF of standard normal
double normcdf(double x)
{ 
	return alnorm(x, false);
}






