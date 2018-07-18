#ifndef MATLABSIM_H_
#define MATLABSIM_H_

#include "masterheader.h"

// This class is used to work with a Matlab simulator through the use of Matlab engine 
using std::vector;

class MatlabSim : public Simulator
{
private:
	//Engine *_ep; //Define matlab engine and test it
	std::string _matlabdir; 
	char* _buffer;
	size_t _maxThreads; 
public:
	MatlabSim(const std::string& matlabdir, size_t maxth);
	~MatlabSim();
	double simulation(const Solution* const solution);
	double simulation(vector<Solution*>& visited, const vector<size_t>& solnid, const vector<size_t> numrep);
	double gettruevalue(const Solution* const solution) {return solution->getSampleMean();}
	int isGlobalOptimum(const Solution& solution);
};

#endif /*MATLABSIM_H_*/
