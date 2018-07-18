#include "myutilities.h"
#include <iostream>
#include <ctime>
using namespace std;

int maximum(double* x, int size)
{
	double maximum = x[0];
	int index = 0;
	for(int i = 1; i < size; i++)
	{
		if (maximum < x[i])
		{
			index = i;
			maximum = x[i];
		}
	}
	return index; 
} 
	
int maximum(int* x, int size)
{
	double maximum = x[0];
	int index = 0;
	for(int i = 1; i < size; i++)
	{
		if (maximum < x[i])
		{
			index = i;
			maximum = x[i];
		}
	}
	return index; 
}

double max(double x, double y)
{
	return (x >= y) ? x : y;
}

int min(int x, int y)
{
	return (x <= y) ? x : y;
}

int min(double* x, int size)
{
	double minimum = x[0];
	int index = 0;
	for(int i = 1; i < size; i++)
	{
		if (minimum > x[i])
		{
			index = i;
			minimum = x[i];
		}
	}
	return index;
}

int minimum(int* x, int size)
{
	int minimum = x[0];
	int index = 0;
	for(int i = 1; i < size; i++)
	{
		if (minimum > x[i])
		{
			index = i;
			minimum = x[i];
		}
	}
	return index;
}


long findMinSolution(const vector<Solution*>& visited, const vector<long>& visitedThisRun)
{
	double minimum = 1e20;
	long index = -1;
	for (long i = 0; i < (long) visitedThisRun.size(); i++)
	{
		if (minimum > visited[visitedThisRun[i]]->getSampleMean())
		{
			minimum = visited[visitedThisRun[i]]->getSampleMean();
			index = i;
		}
	}
	return visitedThisRun[index];
}

long findMinSolution(const vector<Solution*>& visited, long indvIndex, const vector<long>& visitedThisRun)
{
	if (visitedThisRun.size() == 0) return indvIndex;
	double minimum = 1e50;
	long index = -1;
	for (long i = 0; i < (long) visitedThisRun.size(); i++)
	{
		if (minimum > visited[visitedThisRun[i]]->getSampleMean())
		{
			minimum = visited[visitedThisRun[i]]->getSampleMean();
			index = i;
		}
	}
	if (visited[visitedThisRun[index]]->getSampleMean() < visited[indvIndex]->getSampleMean())
		return visitedThisRun[index];
	else
		return indvIndex;
}

long find(std::vector<long> x, long element)
{
	long index = -1;
	for (long i = 0; i < (long) x.size(); i++)
	{
		if (element == x[i])
		{
			index = i;
			break;
		}
	}
	return index;
} 

long findSolution(std::vector<Solution*> visited, Solution* solution)
{
	for (long i = 0; i < (long) visited.size(); i++)
	{
		if (visited[i]->isEqual(solution))
		{
			return i;
		}
	}
	return -1;
}

void copy(int* x, int* y, int dimension)
{
	for (int i = 0; i < dimension; i++)
	{
		y[i] = x[i];
	}
}

double dabs(double x)
{
	return (x > 0) ? x : -x;
}

// Copy the contents of x into y
void vectorCopy(vector<double*> x, vector<double*>& y)
{
	y.clear();
	for (int i = 0; i < (int) x.size(); i++)
	{
		y.push_back(x[i]);
	}	
}


void vectorCopy(vector<double> x, vector<double>& y)
{
	y.clear();
	for (int i = 0; i < (int) x.size(); i++)
	{
		y.push_back(x[i]);
	}
}

void vectorCopy(vector<long> x, vector<long>& y)
{
	y.clear();
	for (int i = 0; i < (int) x.size(); i++)
	{
		y.push_back(x[i]);
	}
}

void vectorCopy(vector<unsigned int> x, vector<unsigned int>& y)
{
	y.clear();
	for (int i = 0; i < (int) x.size(); i++)
	{
		y.push_back(x[i]);
	}
}

void printVector(vector<double> b)
{
	cout << "Vector b size " << b.size() << endl;
	for (long i = 0; i < (long) b.size(); i++)
	{
		cout << b[i] << " ";
	}
	cout << endl;
}

void printVector(vector<int> x)
{
	cout << "Vector x size " << x.size() << endl;
	for (long i = 0; i < (long) x.size(); i++)
	{
		cout << x[i] << " ";
	}
	cout << endl;
}

void printVector(vector<unsigned int> x)
{
	cout << "Vector x size " << x.size() << endl;
	for (long i = 0; i < (long) x.size(); i++)
	{
		cout << x[i] << " ";
	}
	cout << endl;
}

void printVector(std::vector<Record*> r)
{
	cout << "Records size " << r.size() << endl;
	for (long i = 0; i < (long) r.size(); i++)
	{
		r[i]->print();
	}
	cout << endl;
}


void printVector(std::vector<long> W)
{
	for (long i = 0; i < (long) W.size(); i++)
	{
		cout << W[i] << " ";
	}
	cout << endl;
}

void printVector(vector<double*> A, int dimension)
{
	cout << "Matrix A Size " << A.size() << " Dimension " << dimension << endl;
	for (long i = 0; i < (long) A.size(); i++)
	{
		for (int j = 0; j < dimension; j++)
		{
			cout << *(A[i]+j) << " ";
		}
		cout << endl;
	}
}

/*void printVector(double* A, int dimension, const char* name)
{
	cout << name << endl;
	for (int i = 0; i < dimension; i++)
	{
		for (int j = 0; j < dimension; j++)
		{
			cout << A[i*dimension+j] << ", ";
		}
		cout << endl;
	}
}

void printVector(int* A, int dimension, const char* name)
{
	cout << name << endl;
	for (int i = 0; i < dimension; i++)
	{
		for (int j = 0; j < dimension; j++)
		{
			cout << A[i*dimension+j] << ", ";
		}
		cout << endl;
	}
}*/


void printVector(vector<double*> A, int dimension, vector<double> b)
{
	cout << "Matrix A Size " << A.size() << " Dimension " << dimension << endl;
	for (long i = 0; i < (long) A.size(); i++)
	{
		for (int j = 0; j < dimension; j++)
		{
			cout << *(A[i]+j) << " ";
		}
		cout << ">= " << b[i] << endl;
	}
}

void printVector(std::vector<Solution*> X)
{
	cout << "Size of solution vector   " << X.size() << endl;
	for (unsigned int i = 0; i < X.size(); i++)
	{
		X[i]->printObservations();
	}
}

void PrintActiveSolutions(std::vector<Solution*> X, std::vector<long> W)
{
	//cout << "Active solutions" << endl;
	for (unsigned int i = 0; i < W.size(); i++)
	{
		if (W[i] >= long(X.size()))
			cout << "X.size() = " << X.size() << " W[i] = " << W[i] << endl;
		cout << *X[W[i]];
	}
}


void PrintActiveSolutions(std::vector<Solution*> X, std::vector<unsigned int> W)
{
	//cout << "Active solutions" << endl;
	for (unsigned int i = 0; i < W.size(); i++)
	{
		if (W[i] >= X.size())
			cout << "X.size() = " << X.size() << " W[i] = " << W[i] << endl;
		cout << *X[W[i]];
	}
}

void printVectorObservations(std::vector<Solution*> X)
{
	cout << "Size of solution vector   " << X.size() << endl;
	for (unsigned int i = 0; i < X.size(); i++)
	{
		cout << "Solution # " << i << endl;
		X[i]->printObservations();
	}
}

void PrintConstraints(int dimension, std::vector<double*> A, std::vector<double> b, ostream& os, const char* sep)
{
	for (size_t c = 0; c < b.size(); c++)
	{
		for (int d = 0; d < dimension; d++)
		{
			os << A[c][d] << sep;
		}
		os << b[c] << endl;
	}
}

long AddSolution(std::vector<Solution*>& xVector, Solution* x)
{
	for (long i = 0; i < (long) xVector.size(); i++)
	{
		if (x->isEqual(xVector[i]))
		{
			// This genome already exists
			delete x;
			return i;
		}
	}
	xVector.push_back(x);
	return long(xVector.size() - 1);
}


void MySleep(double seconds)
{
	clock_t time1 = clock();
	clock_t time2 = time1 + seconds*CLOCKS_PER_SEC;
	clock_t currtime;
	do
	{
		currtime = clock();
	} while (currtime < time2);
}


// 12/28/11 computes choose k from n
long ChooseKN(long k, long n)
{
	double prod = 1.0;
	if (k > n-k)
	{
		for (long i = 1; i <= n-k; i++)
			prod *= (k+i) / i;
	}
	else
	{
		for (long i = 1; i <= k; i++)
			prod *= (n-k+i) / i;
	}
	return static_cast<long> (prod);
}
