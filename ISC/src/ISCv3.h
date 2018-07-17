#ifndef ISCV3_H_
#define ISCV3_H_
#include "masterheader.h"
#include <list>
#include <map>
/* Revised COMPASS algorithm for constrained problem with equal-sample-allocation rule. Linear 
 * constraints take the form Ax>=b. x0 is the starting point, x0=[x01,...x0d]. The objective is a minimization problem.
 * The algorithm takes as input:
 * 
 * epsilon: for stopping test procedure to confirm true local optima
 * alpha: CI level for stopping test
 * 
 * visited: a list of visited solutions, which contains the given starting point if it is the first run. This list is kept throughout 
 * the entire optimization process, i.e., it lives longer than a single run of compass. The motivation is that past simulation output 
 * data should be kept to save effort as well as improve accuracy. 
 * A, b: vectors which specify the original problem constraint. 
 * activeSolutions: candidate active solutions for visited[starIndex]. Need to run MPR to find which are truly are
 * bestSolutions: best solutions found so far by GA
 * starIndex: the index in visited as the starting point as well as the local optimum solution returned by one run.
 * budget: the number of available simulation replications, when the algorithm finishes, budget becomes the remaining computing power
 * reocrds: to record evolution of COMPASS
 * aType: original COMPASS or revised COMPASS
 * logFile: reference to the file where statistics of the algorithm run are kept. 
 * 
 * initialNumReps: initial number of reps assigned to a new solution
 * numCandidates: size of S_k, i.e. # of samplings within MPR (repetion ok)
 * alphac: CI-level for generalized COMPASS constraint placement procedure
 * alphae: CI-level for stopping test in Ver 2 of revised COMPASS to ensure that a optimal is indeed optimal before proceeding with COMPASS steps
 * deltae: indifference zone parameter for stopping test in Ver. 2 of revised COMPASS 
 * Return value: index of xStar.
 */

/*long ISCv3(Simulator* simulator, bool stocsim, double deltam, double alpham, std::vector<Solution*>& visited, std::vector<double*> A, \
		std::vector<double> b, std::vector<long>& activeSolutions, \
		long& budget, long starIndex, std::ofstream& logFile, std::vector<Record*>& records, \
		AlgorithmType aType, int initialNumReps, int numCandidates, \
		double alphac, double alphae, double deltae, int pruning_freq, int mpamode,int rmdmode, long lastRep, \
		double lastObjectiveValue, double bestValue, double trueBestObjectiveValue, \
		bool backTrackingTest, bool ocba, bool firstRunOfCOMPASS, bool& optimum, int metamodel);*/
long ISCv3(Simulator* simulator, std::vector<Solution*>& visited, std::vector<double*> A, std::vector<double> b, ISCParameters* iscparams,\
		std::vector<long>& activeSolutions, long& budget, long starIndex, std::ofstream& logFile, std::vector<Record*>& records, \
		long lastRep, double lastObjectiveValue, double bestValue, double trueBestObjectiveValue, bool firstRunOfCOMPASS, bool& optimum);

const int COMPASSMPA = 0, BOXMPA = 1;
// key is the coordinate position, which is an int. value is the index of the solution in the visited vector, a size_t
typedef std::multimap<int, size_t>  CornerMap;
typedef CornerMap::iterator CMIter;
typedef std::multimap<double, size_t>  ObjMap;
typedef ObjMap::iterator OMIter;

#endif /*ISCV3_H_*/
