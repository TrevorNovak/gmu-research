#ifndef STATISTICS_H_
#define STATISTICS_H_

#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include "Solution.h"
#include "Simulator.h"

class Record
{
	private:
		double sum;
		long number;
		double mean;
		double sumsquared;
		double samplevariance;
	public:
		Record() : sum(0), number(0), mean(0), sumsquared(0), samplevariance(0) {}
		Record(double s) : sum(s), number(1), mean(s) {}
		Record(const Record& r) : sum(r.sum), number(r.number), mean(r.mean) {}
		// Create a new Record object using an existing one and a new observation 
		Record(const Record& r, double s) : sum(r.sum+s), number(r.number+1), mean((r.sum+s)/(r.number+1)) {}
		void writeRecord(double result) 
		{
			sum += result;
			sumsquared += result * result;
			number++;
			mean = sum / number;
			if (number == 1)
			{
				samplevariance = 0;
			}
			else
			{
				samplevariance = (sumsquared - number * mean * mean) / (number - 1);
			}
		}
		void clear()
		{
			sum = mean = sumsquared = samplevariance = 0;
			number = 0;
		}
		double getMean() const {return mean;}
		bool isEmpty() const {return (number == 0);}
		double getVariance() const {return samplevariance;}
		void print()
		{
			std::cout << "mean " << mean << "  #" << number << std::endl;
		}
		friend std::ostream& operator<<(std::ostream& os, const Record& record);
};

void RecordResult(std::vector<Record*>& records, long repCount, double objectiveValue);

void WritePlot(std::ofstream& logFile, std::vector<Record*> records, long budget);
//This function reads the result file written out by ISC and checks if the solution found is a local or global optimum
void CompileRunStat(std::ifstream& resultFile, std::ofstream& summaryfile, int dimension, Simulator& sim, const std::string& caseid, const std::vector<double*> & A, const std::vector<double> b);
#endif /*STATISTICS_H_*/
