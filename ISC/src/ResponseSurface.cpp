#include "masterheader.h"
using namespace std;

double ResponseSurface::F22(const Solution* const solution, bool isNoisy)
{
	int dimension;
	dimension = solution->getDimension();
	require(dimension == 2, "Error of F22 function");
	const int* const x = solution->getDiscreteSolution();
	double f22 = -pow(sin(0.05 * 3.1415926535897932384626 * x[0]), 6) / pow(2, 2 * ((x[0] - 10)/80.0) * ((x[0] - 10)/80.0)) \
			-pow(sin(0.05 * 3.1415926535897932384626 * x[1]), 6) / pow(2, 2 * ((x[1] - 10)/80.0) * ((x[1] - 10)/80.0));
	if (isNoisy) f22 += Normal(0, 0.3 * 0.3);
	return f22;
}

double ResponseSurface::F2n(const Solution* const solution, bool isNoisy)
{
	int dimension;
	dimension = solution->getDimension();
	//require(dimension == 2, "Error of F22 function");
	const int* const x = solution->getDiscreteSolution();
	double f2 = 0;
	for (int d = 0; d < dimension; d++)
	{
		int dmod9 = (d > 9) ? 8 : d+1; 
		f2 -= pow(sin(0.05 / dmod9 * 3.1415926535897932384626 * x[d]), 6) / pow(2, 2 * ((x[d] - 10 * dmod9)/80.0) * ((x[d] - 10 * dmod9)/80.0));
	}
	if (isNoisy) f2 += Normal(0, 0.3 * 0.3);
	return f2;
}

double ResponseSurface::singular(const Solution* const solution, bool isNoisy)
{
	double result = 0;
	int dimension;
	dimension = solution->getDimension();
	require(dimension == 4, "Error of singular function");
	const int* const x = solution->getDiscreteSolution();
	result += 1 + (x[0] + 10*x[1])*(x[0] + 10*x[1]) + 5*(x[2]-x[3])*(x[2]-x[3]) + pow(double(x[1] - 2*x[2]), 4.0) + 10*pow(double(x[0]-x[3]), 4.0);
	double variance = result < 900 ? result : 900;
	//double variance = 900; 
	double noise;
	if (isNoisy)
	{
		noise = Normal(0, variance);
		result += noise;
	}
	return result;
	
}

double ResponseSurface::rosenbrock(const Solution* const solution, bool isNoisy)
{
	double result = 0;
	int dimension;
	dimension = solution->getDimension();
	require(dimension == 2, "Error of singular function");
	const int* const x = solution->getDiscreteSolution();
	result += 100*(static_cast<double>(x[1]) - x[0]*x[0])*(x[1] - x[0]*x[0]) + (1-x[0])*(1-x[0]);
	double variance = result+1;
	double noise;
	if (isNoisy)
	{
		noise = Normal(0, variance);
		result += noise;
	}
	return result;
}

double ResponseSurface::quadratic(const Solution* const solution, bool isNoisy)
{
	// This is a quadratic function with random noise
	double result = 0;
	int dimension;
	dimension = solution->getDimension();
	const int* const x = solution->getDiscreteSolution();
	//RngStream rng("rng1");
	for (int i = 0; i < dimension; i++)
	{
		result += x[i] * x[i];
	} 
	result++;
	double noise;
	if (isNoisy)
	{
		noise = Normal(0, 0.1 * result * 0.1 * result);
		//noise = Normal(0, 30 * 30);
		/*double gx0 = 0;
		for (int i = 0; i < dimension; i++)
		{
			gx0 += 80 * 80;
		}
		gx0++;
		double variance = max(1, gx0 - result);  
		noise = Normal(0, variance);*/
		result += noise;
	}
	return result/1000;
}


double ResponseSurface::highd(const Solution* const solution, bool isNoisy)
{
	int dimension;
	dimension = solution->getDimension();
	const int* const x = solution->getDiscreteSolution();
	int m = 1;
	//double gamma[] = {100, 81, 225};
	//double gamma[] = {5, 6};
	double gamma[] = {1e-3};
	//long xi[] = {36, 168};
	long xi[] = {0};
	//double f = 12500;
	double beta[] = {1e4};
	double g = 0; 
	for (int i = 0; i < m; i++)
	{
		double sum = 0;
		for (int j = 0; j < dimension; j++) 
		{
			//sum += (x[j] - xi[i]) * (x[j] - xi[i]);
			sum += (x[j] - xi[i]) * (x[j] - xi[i]) * gamma[i] * (j+1);
		}
		//g +=  f / gamma[i] * exp(-0.5 / gamma[i] / gamma[i] * sum);
		g +=  beta[i] * exp(-sum);
	} 
	//double variance = 100.0/36.0;
	double variance = g * 0.09;
	variance = mymax(variance, 1.0);
	double noise;
	if (variance < 1e-30)
		noise = 0;
	else
		noise = Normal(0, variance);
	//return -g;
	return -(g + ((isNoisy)? noise : 0));
}

double ResponseSurface::MultiLocalHighd(const Solution* const solution, bool isNoisy)
{
	int dimension;
	dimension = solution->getDimension();
	const int* const x = solution->getDiscreteSolution();
	double gamma[] = {0.001, 0.005};
	size_t m = sizeof(gamma)/sizeof(gamma[0]); //m is the number of normal densities to add up
	//double gamma[] = {5, 6};
	//long xi[] = {36, 168};
	long xi[] = {-38, 56};
	//double f = 12500;
	double beta[] = {300, 500};
	double g = 0; 
	for (int j = 0; j < dimension; j++) 
	{
		for (size_t i = 0; i < m; i++)
		{
			g += beta[i] * exp(-(x[j] - xi[i]) * (x[j] - xi[i]) * gamma[i]);
		}
	} 
	//double variance = 100.0/36.0;
	double variance = 1;
	double noise;
	if (variance < 1e-30)
		noise = 0;
	else
		noise = Normal(0, variance);
	//return -g;
	return -(g + ((isNoisy)? noise : 0));
}

int ResponseSurface::isGlobalOptimum(const Solution& solution)
{
	vector<int> xgopt;
	if (_rsfid == 0)
	{//f22
		xgopt.push_back(10);
		xgopt.push_back(10);
		Solution gopt(xgopt);
		if (solution.isEqual(&gopt))
			return 1;
		else 
			return 0;
	}
	else if (_rsfid == 1)
	{//singular
		for (int i = 0; i < 4; i++) xgopt.push_back(0);
		Solution gopt(xgopt);
		if (solution.isEqual(&gopt))
			return 1;
		else 
			return 0;
	}
	else if (_rsfid == 2)
	{//highd
		for (int i = 0; i < solution.getDimension(); i++) xgopt.push_back(0);
		Solution gopt(xgopt);
		if (solution.isEqual(&gopt))
			return 1;
		else 
			return 0;
	}
	else if (_rsfid == 3)
	{//multilocal highd
		for (int i = 0; i < solution.getDimension(); i++) xgopt.push_back(56);
		Solution gopt(xgopt);
		if (solution.isEqual(&gopt))
			return 1;
		else 
			return 0;
	}
	else if (_rsfid == 6)
	{//quadratic
		for (int i = 0; i < solution.getDimension(); i++) xgopt.push_back(0);
		Solution gopt(xgopt);
		if (solution.isEqual(&gopt))
			return 1;
		else 
			return 0;
	}
	else
		return 0;
} 

double ResponseSurface::rsf(const Solution* const solution, bool isNoisy)
{
	switch (_rsfid)
	{
		case 0 : return F22(solution, isNoisy);
		case 1 : return singular(solution, isNoisy);
		case 2 : return highd(solution, isNoisy);
		case 3 : return MultiLocalHighd(solution, isNoisy);
		case 4 : return rosenbrock(solution, isNoisy);
		case 5 : return F2n(solution, isNoisy);
		case 6 : return quadratic(solution, isNoisy);
		default :
			cout << "Wrong RSF function id: " << _rsfid << endl;
			exit(5);
	}
}

