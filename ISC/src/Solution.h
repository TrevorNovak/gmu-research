#ifndef SOLUTION_H_
#define SOLUTION_H_

#include <vector>
#include <iostream>
const int DISCRETE = 0;
const int CONTINUOUS = 1;
class Solution 
{
	private:
		int _dimension;  // dimension of the solution space
		int* _discrete_solution;  
		double _sumSquared;
		double _sum;
		double _mean;
		double _trueObjectiveValue;
		double _sample_variance;
		int _numOfObservations;
		double _sharedFitness; // This is for niching GA. It's an ugly way to do it...
		double _sharedSampleVar; 
		std::vector<double> _objs; 
	public:
		Solution(int n, const int* const x);
		Solution(int n);
		Solution(const Solution&);
		Solution& operator=(const Solution&);
		bool isFeasible(const std::vector<double*> & A, const std::vector<double> b) const; 
		Solution(std::vector<int> x);
		virtual ~Solution() 
		{
			delete []_discrete_solution;
		}
		int clone(const Solution& x);
		void setSolution(const int* const x);
		void setSolution(const std::vector<int>& x);
		void changeOneDimension(int x, int n);
		const int* const getDiscreteSolution() const;
		int getDimension() const;
		void setDimension(int dim);
		int getNumOfObservations() const;
		void setNumOfObservations(int num);
		double recordObservation(double result);
		double getSampleMean() const;
		void setSampleMean(double mean);
		double getTrueObjectiveValue() const {return _trueObjectiveValue;}
		void setTrueObjectiveValue(double value) {_trueObjectiveValue = value;}  
		double getSampleVariance() const;
		void setSampleVariance(double svar);
		double getSumSquared() const;
		void setSumSquared(double squared);
		double getSum() const;
		void setSum(double sum);
		bool isEqual(Solution* x) const;
		void clear();
		void printObservations() const;
		void print() const;
		void printXobs(std::ostream& os) const; //This function is used to print csv file for a solution used as stochastic kriging input
		void printX(std::ostream& os) const; //This function is used to print txt file for extern simulator to read 
		const double getDistanceTo(const Solution* x) const; 
		const double getEuclideanDistanceTo(const Solution& x) const; 
		void setSharedFitness(double shared);
		const double getSharedFitness() const;
		void setSharedSampleVar(double shared);
		const double getSharedSampleVar() const;
		friend std::ostream& operator<<(std::ostream& os, const Solution& solution);
		friend std::ostream& operator<<(std::ostream& os, const Solution* solution);
		friend bool operator<(const Solution& solution1, const Solution& solution2);
		friend bool operator<=(const Solution& solution1, const Solution& solution2);
		friend bool operator>(const Solution& solution1, const Solution& solution2);
		friend bool operator>=(const Solution& solution1, const Solution& solution2);
		friend bool operator==(const Solution& solution1, const Solution& solution2);
};


#endif /*SOLUTION_H_*/



