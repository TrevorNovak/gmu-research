#ifndef UTILITIES_H_
#define UTILITIES_H_
#include "Solution.h"
#include "Simulator.h"
#include "require.h"
#include "Statistics.h"
#include "RMD.h"
#include <vector>
#include <list>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <iterator>
#include <math.h>

enum AlgorithmType
{
	compass, 
	revisedcompass	
};

int maximum(double* x, int size);
int maximum(int* x, int size);
double max(double x, double y);

int min(int x, int y);
int min(double* x, int size);
int minimum(int* x, int size);

long findMinSolution(const std::vector<Solution*>& visited, long indvIndex, const std::vector<long>& visitedThisRun);
long findMinSolution(const std::vector<Solution*>& visited, const std::vector<long>& visitedThisRun);
long findSolution(std::vector<Solution*> visited, Solution* solution);
long find(std::vector<long> x, long element); 
double dabs(double x);
// Clear y and copy the contents of x into y
void vectorCopy(std::vector<double*> x, std::vector<double*>& y);
void vectorCopy(std::vector<double> x, std::vector<double>& y); 
void vectorCopy(std::vector<long> x, std::vector<long>& y);
void vectorCopy(std::vector<unsigned int> x, std::vector<unsigned int>& y);
void copy(int* x, int* y, int dimension);
void printVector(std::vector<double> b);
void printVector(std::vector<int> x);
void printVector(std::vector<unsigned int> x);
void printVector(std::vector<long> W);
void printVector(std::vector<double*> A, int dimension);
//void printVector(double* A, int dimension, const char* name);
//void printVector(int* A, int dimension, const char* name);
void printVector(std::vector<double*> A, int dimension, std::vector<double> b);
void printVector(std::vector<Solution*> X);
void printVector(std::vector<Record*> r);
void printVectorObservations(std::vector<Solution*> X);
void PrintActiveSolutions(std::vector<Solution*> X, std::vector<long> W);
void PrintActiveSolutions(std::vector<Solution*> X, std::vector<unsigned int> W);
void PrintConstraints(int dimension, std::vector<double*> A, std::vector<double> b, std::ostream& os = std::cout, const char* sep=",");
void MySleep(double seconds);
// Add a solution to a vector if it is not redundant, return its index in visited
long AddSolution(std::vector<Solution*>& xVector, Solution* x);

// To print an array
template<typename T> 
void printVector(T* A, int dimension, const char* name)
{
	std::cout << name << std::endl;
	for (int i = 0; i < dimension; i++)
	{
		for (int j = 0; j < dimension; j++)
		{
			std::cout << A[i*dimension+j] << ", ";
		}
		std::cout << std::endl;
	}
}

template<typename Iter>
void printcontainer(Iter first, Iter last, const char* nm = "", const char* sep = "\n", std::ostream& os = std::cout) 
{
	if(nm != 0 && *nm != '\0')	os << nm << ": " << "\n";
	typedef typename std::iterator_traits<Iter>::value_type T;
	std::copy(first, last, std::ostream_iterator<T>(os, sep));
  	os << std::endl;
}

template<typename Iter>
void printmap(Iter first, Iter last, const char* nm = "", const char* sep = "\n", std::ostream& os = std::cout) 
{
	if(nm != 0 && *nm != '\0')	os << nm << ": " << "\n";
	Iter it;
	for (it = first; it != last; it++)
	{
		os << (*it).first << ", " << (*it).second << std::endl;
	}
	os << std::endl;
}


template<typename Iter>
void printpointercontainer(Iter first, Iter last, const char* nm = "", const char* sep = "\n", std::ostream& os = std::cout) 
{
	if(nm != 0 && *nm != '\0')
    	os << nm << ": " << sep;
	typedef typename std::iterator_traits<Iter>::value_type T;
	for (; first != last; first++)
		os << *(*first) << sep;
  	os << std::endl;
}

template<class Cont>
void printcontainer(Cont& c, const char* nm = "", const char* sep = "  ", std::ostream& os = std::cout) 
{
	printcontainer(c.begin(), c.end(), nm, sep, os);
}

template<class Cont>
size_t max(const Cont& c)
{
	require(c.size() > 0, "max(c), c container size = 0");
	typedef typename Cont::iterator Iter;
	typedef typename std::iterator_traits<Iter>::value_type T;
	T maximum = c[0];
	size_t index = 0;
	for (size_t i = 1; i < c.size(); i++)
	{
		if (c[i] > maximum)
		{
			maximum = c[i];
			index = i;
		}
	}
	return index;
} 

template<typename Iter>
double summarystat(Iter first, Iter last, double& mean, double& std, double& max, double& min)
{
	double total = 0.0, sumsqr = 0.0, var = 0.0;
	max = *first;
	min = *first;
	size_t num_element = 0; 
	Iter it;
	for (it = first; it != last; it++)
	{
		total += *it;
		num_element++;
		if (*it > max)
			max = *it;
		else if (*it < min)
			min = *it;
	}
	require(num_element > 0, "Take mean of empty container");
	mean = total / num_element;
	for (it = first; it != last; it++)
	{
		var += ((*it)-mean) / num_element * ((*it) - mean) ;
	}
	//var = (sumsqr - num_element * mean * mean) / num_element;
	std = sqrt(var);
	return mean;
}

template<typename T> const T& myabs(const T& a)
{
  return (a < 0) ? -a : a;
}

template<typename T> const T& mymax(const T& a, const T& b)
{
  return (a < b) ? b : a;
}

template<typename T> const T& mymin(const T& a, const T& b)
{
  return (a > b) ? b : a;
}
// 12/28/11 computes choose k from n
long ChooseKN(long k, long n);
#endif /*UTILITIES_H_*/
