#include "masterheader.h"
using namespace std;

Solution::Solution(int n, const int* const x)
{
	_dimension = n;
	_discrete_solution = new int[n];
	for (int i = 0; i < n; i++)
	{
		_discrete_solution[i] = x[i];
	}
	_sumSquared = 0;
	_sum = 0;
	_numOfObservations = 0;
	_mean = 0;
	_sample_variance = 0;
	_sharedFitness = 0;
	_trueObjectiveValue = 0;
}

Solution::Solution(vector<int> x)
{
	_dimension = int(x.size());
	_discrete_solution = new int[_dimension];
	for (int i = 0; i < _dimension; i++)
	{
		_discrete_solution[i] = x[i];
	}
	_sumSquared = 0;
	_sum = 0;
	_numOfObservations = 0;
	_mean = 0;
	_sample_variance = 0;
	_sharedFitness = 0;
	_trueObjectiveValue = 0;
}


Solution::Solution(int n)
{
	_dimension = n;
	_discrete_solution = new int[_dimension];
	_sumSquared = 0;
	_sum = 0;
	_numOfObservations = 0;
	_mean = 0;
	_sample_variance = 0;
	_sharedFitness = 0;
	_trueObjectiveValue = 0;
}

Solution::Solution(const Solution& x)
{
	_dimension = x.getDimension();
	_discrete_solution = new int[_dimension];
	for (int i = 0; i < _dimension; i++)
	{
		_discrete_solution[i] = x.getDiscreteSolution()[i];
	}
	_sumSquared = x.getSumSquared();
	_sum = x.getSum();
	_mean = x.getSampleMean();
	_sample_variance = x.getSampleVariance();
	_numOfObservations = x.getNumOfObservations();
	_sharedFitness = x.getSharedFitness();
	_trueObjectiveValue = x.getTrueObjectiveValue();
}

Solution& Solution::operator=(const Solution& x)
{
	// check self-assignment
	if (this != &x)
	{
		delete []_discrete_solution;
		_dimension = x.getDimension();
		_discrete_solution = new int[_dimension];
		for (int i = 0; i < _dimension; i++)
		{
			_discrete_solution[i] = x.getDiscreteSolution()[i];
		}
		_sumSquared = x.getSumSquared();
		_sum = x.getSum();
		_mean = x.getSampleMean();
		_sample_variance = x.getSampleVariance();
		_numOfObservations = x.getNumOfObservations();
		_sharedFitness = x.getSharedFitness();
		_trueObjectiveValue = x.getTrueObjectiveValue();
	}
	return *this;
}

int Solution::clone(const Solution& x)
{
	for (int i = 0; i < _dimension; i++)
	{
		_discrete_solution[i] = x.getDiscreteSolution()[i];
	}
	_sumSquared = x.getSumSquared();
	_sum = x.getSum();
	_mean = x.getSampleMean();
	_sample_variance = x.getSampleVariance();
	_numOfObservations = x.getNumOfObservations();
	_sharedFitness = x.getSharedFitness();
	_trueObjectiveValue = x.getTrueObjectiveValue();
	return 1;
}

void Solution::clear()
{
	for (int i = 0; i < _dimension; i++)
	{
		_discrete_solution[i] = 0;
	}
	_sumSquared = 0;
	_sum = 0;
	_mean = 0;
	_sample_variance = 0;
	_numOfObservations = 0;
	_sharedFitness = 0;
	_trueObjectiveValue = 0;
}

bool Solution::isFeasible(const std::vector<double*> & A, const std::vector<double> b) const
{
	for (size_t con = 0; con < A.size(); con++)
	{
		double total = 0.0;
		for (int i = 0; i < _dimension; i++)
			total += *(A[con]+i) * *(_discrete_solution+i);
		if (total < b[con])
			return false;
	}
	return true;
}

bool Solution::isEqual(Solution* x) const
{
	for (int i = 0; i < _dimension; i++)
	{
		if (_discrete_solution[i] != x->getDiscreteSolution()[i])
		{
			return false;
		}
	}
	return true;
}

void Solution::setSolution(const int* const x)
{
	for (int i = 0; i < _dimension; i++)
	{
		_discrete_solution[i] = x[i];
	}
}

void Solution::setSolution(const std::vector<int>& x)
{
	for (int i = 0; i < _dimension; i++)
	{
		_discrete_solution[i] = x[i];
	}
}

void Solution::changeOneDimension(int x, int n)
{
	_discrete_solution[n] = x;
}


const int* const Solution::getDiscreteSolution() const
{
	return _discrete_solution;
}


int Solution::getDimension() const
{
	return _dimension;
}

void Solution::setDimension(int dim)
{
	_dimension = dim;
}

void Solution::setNumOfObservations(int num)
{
	_numOfObservations = num;
}

int Solution::getNumOfObservations() const
{
	return _numOfObservations;
}


double Solution::recordObservation(double result)
{
	_objs.push_back(result);
	_sum += result;
	_sumSquared += result * result;
	_numOfObservations++;
	//double xbarn1 = _mean;
	_mean = _sum / _numOfObservations;
	if (_numOfObservations == 1)
	{
		_sample_variance = 0;
	}
	else
	{	// To avoid numerical rounding errors when simulation results are huge, we want to update sample var in the following way
		// S_n^2 = (n-2)*S_{n-1}^2/(n-1) + Xbar_{n-1}^2 + X_n^2/(n-1) - n/(n-1)* Xbar_n^2
		_sample_variance = 0;
		//_sample_variance = (_sumSquared - _numOfObservations * _mean * _mean) / (_numOfObservations - 1);
		//_sample_variance = (_numOfObservations-2)/static_cast<double>(_numOfObservations-1)*sn1 + (xbarn1-_mean)*(xbarn1+_mean) + (result-_mean)*(result+_mean)/(_numOfObservations-1);
		/*double t1 = (_numOfObservations-2)/static_cast<double>(_numOfObservations-1)*sn1;
		double t2 = (xbarn1-_mean)*(xbarn1+_mean);
		double t3 = (result-_mean)*(result+_mean)/(_numOfObservations-1);
		double t5 = t1 + t2 + t3;*/
		for (int i = 0; i < _numOfObservations; i++)
		{
			_sample_variance += (_objs[i]-_mean)*(_objs[i]-_mean)/(_numOfObservations-1);
		}
		if (_sample_variance < 0) 
		{
			double n1 = _numOfObservations * _mean * _mean;
			double n2 = _sumSquared - n1;
			double n3 = n2 / (_numOfObservations - 1);
			cout << this;
			cout << "error, negative sample var" << endl;
			exit(5);
		}
	}
	return _mean;
}


double Solution::getSampleMean() const
{
	return _mean;
}


void Solution::setSampleMean(double mean) 
{
	_mean = mean;
}


double Solution::getSampleVariance() const
{
	return _sample_variance;
}

void Solution::setSampleVariance(double svar)
{
	_sample_variance = svar;
}


double Solution::getSumSquared() const
{
	return _sumSquared;
}

void Solution::setSumSquared(double squared)
{
	_sumSquared = squared;
}


double Solution::getSum() const
{
	return _sum;
}
		
void Solution::setSum(double sum)
{
	_sum = sum;		
}

void Solution::print() const
{
	for (int i = 0; i < _dimension; i++)
	{
		cout << _discrete_solution[i] << " ";
	}
	cout << endl;
}

void Solution::printX(ostream& os) const
{
    for (int i = 0; i < _dimension; i++)
	{
		os << "\t" << _discrete_solution[i];
	}
}

void Solution::printXobs(ostream& os) const
{
	for (int i = 0; i < _dimension; i++)
	{
		os << _discrete_solution[i] << ",";
	}
	os << _mean << "," << _sample_variance << "," << _numOfObservations << endl;
}

void Solution::printObservations() const
{
	cout << "(";
	for (int i = 0; i < _dimension - 1; i++)
	{
		cout << _discrete_solution[i] << " ";
	}
	cout << _discrete_solution[_dimension - 1];
	cout << ")  ";
	//cout << endl;
	cout << "obs " << _numOfObservations;
	cout << "  mean " << _mean;
	cout << "  var " << _sample_variance << endl;
}

const double Solution::getDistanceTo(const Solution* x) const
{
	double distance = 0;
	for (int j = 0; j < x->getDimension(); j++)
	{
		distance += pow(double(_discrete_solution[j] - *(x->getDiscreteSolution() + j)), 2);
	}
	return sqrt(distance);
}

const double Solution::getEuclideanDistanceTo(const Solution& x) const
{
	double distance = 0;
	for (int j = 0; j < x.getDimension(); j++)
	{
		distance += fabs(float(_discrete_solution[j] - *(x.getDiscreteSolution() + j))) * fabs(float(_discrete_solution[j] - *(x.getDiscreteSolution() + j)));
	}
	return sqrt(distance);
}

ostream& operator<<(ostream& os, const Solution& solution)
{
	os << "(";
	for (int i = 0; i < solution._dimension - 1; i++)
	{
		os << solution._discrete_solution[i] << ", ";
	}
	os << solution._discrete_solution[solution._dimension - 1];
	os << ")  ";
	//cout << endl;
	os << "obs " << solution._numOfObservations;
	os << "  mean " << solution._mean;
	os << "  var " << solution._sample_variance << endl;
	return os;
}

ostream& operator<<(ostream& os, const Solution* solution)
{
    os << "(";
	for (int i = 0; i < solution->_dimension - 1; i++)
	{
		os << solution->_discrete_solution[i] << ", ";
	}
	os << solution->_discrete_solution[solution->_dimension - 1];
	os << ")  ";
	//cout << endl;
	os << "obs " << solution->_numOfObservations;
	os << "  mean " << solution->_mean;
	os << "  var " << solution->_sample_variance << endl;
	return os;
}

bool operator<(const Solution& solution1, const Solution& solution2)
{
	return solution1.getSampleMean() < solution2.getSampleMean();
}

bool operator<=(const Solution& solution1, const Solution& solution2)
{
	return solution1.getSampleMean() <= solution2.getSampleMean();
}

bool operator>(const Solution& solution1, const Solution& solution2)
{
	return solution1.getSampleMean() > solution2.getSampleMean();
}

bool operator>=(const Solution& solution1, const Solution& solution2)
{
	return solution1.getSampleMean() >= solution2.getSampleMean();
}

bool operator==(const Solution& solution1, const Solution& solution2)
{
	return solution1.getSampleMean() == solution2.getSampleMean();
}

void Solution::setSharedFitness(double shared)
{
	_sharedFitness = shared;
}

const double Solution::getSharedFitness() const
{
	return _sharedFitness;
}

void Solution::setSharedSampleVar(double shared)
{
	_sharedSampleVar = shared;
}


const double Solution::getSharedSampleVar() const
{
	return _sharedSampleVar;
}


