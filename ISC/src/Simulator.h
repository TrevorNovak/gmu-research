#ifndef SIMULATOR_H_
#define SIMULATOR_H_
#include "Solution.h"
#include "RandomVariates.h"
using std::vector;

double simulate(const Solution* const solution, bool isNoisy = true);

class Simulator
{
public:
	enum simMode { serial, batch };
	virtual ~Simulator() {};
	virtual double simulation(const Solution* const solution) = 0;
	virtual double simulation(vector<Solution*>& visited, const vector<size_t>& solnid, const vector<size_t> numrep) = 0; // This functions considers if the simulator runs in batch mode or not
	virtual double gettruevalue(const Solution* const solution) = 0;
	virtual int isGlobalOptimum(const Solution& solution) = 0;
	virtual int isLocalOptimum(const Solution& solution, const std::vector<double*> & A, const std::vector<double> b);
	void setSimMode(simMode smode) { _simMode = smode; }
	bool isBatchSim() {return (_simMode == batch);}
private: 
	simMode _simMode;
	
};

#endif /*SIMULATOR_H_*/
