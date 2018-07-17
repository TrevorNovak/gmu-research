#ifndef ISCCONTROL_H_
#define ISCCONTROL_H_
#include "masterheader.h"

#include "Solution.h"
#include <vector>
#include "ISCv3.h"

/*long ISCControl(Simulator* simulator, long& N, double P, int G, double globaldelta, double globalalpha, double alphadt, \
	double alphac, double deltae, double alphae, double deltam, double alpham, double deltacu, double alphacu, \
	int initialNumReps, int numCandidates, int dimension, std::vector<Solution*>& visited, const std::vector<double*>& A, \
	const std::vector<double>& b, std::ofstream& repFile, std::ofstream& resultFile, std::vector<Record*>& records, bool backTrackingTest, bool ocba, \
	bool dntest, bool docleanup, bool stocsim, bool elitism, int pruning_freq, int mpamode,int rmdmode, int metamodel);*/
long ISCControl(Simulator* simulator, long& N, int dimension, std::vector<Solution*>& visited, const std::vector<double*>& A, \
	const std::vector<double>& b, std::ofstream& repFile, std::ofstream& resultFile, std::vector<Record*>& records, ISCParameters* iscparams);
	
/* simulator: simulator used for this problem
 * N: budget, # of reps
 * P: percentage of budget spent on GA
 * G: numer of generations without improvement for NGA to quit
 * globaldelta: indifference parameter used in grouping procedure of the NGA
 * globalalpha: confidence parameter used in grouping procedure of the NGA 
 * alphadt: dominant niche test alpha parameter of the NGA
 * alphac: CI-level for generalized COMPASS constraint placement procedure
 * deltae, alphae: backtracking test parameters 
 * deltam, alpham:  for optimality stopping test
 * deltacu:  indif parameter for clean-up
 * alphacu: CI level for clean up
 * initialNumReps: initial number of reps assigned to a new solution
 * numCandidates: size of S_k, i.e. # of samplings within MPR and NGA (repetion ok)
 * dimension: dimension of solution space
 * visited: vector of visited solutions. It has the starting point passed from the program calling ISCControl.
 * A, b: problem constraints
 * repFile: write sample path of this rep to this file
 * records: accumulated statistics for averaged sample path after all reps
 * backTrackingTest: true if backTracking test is used 
 * ocba: true if ocba replications application rule is used
 * dntest: true if dominant niche rule is used
 * docleanup : true if cleanup is desired
 */
 
double EvaluateObjectiveValue(const Solution* const solution, bool truevalue);

#endif /*ISCCONTROL_H_*/
