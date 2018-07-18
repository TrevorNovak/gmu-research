#ifndef RMD_H_
#define RMD_H_

#include <vector>
#include "Simulator.h"
#include "Solution.h"
//using namespace std;

const double EPSILON_STABILITY = 1e-20;
const double UPPER_BOUND = 1e31;
const double LARGE_NUMBER = 1e32;

// candidates has the newly sampled solutions returned from RMD. They may have repetitive entries, and so 
// this function checks and removes those redundnant copies. Finally "unique" has all unique new solutions.
size_t UniqueCandidates(std::vector<Solution*> candidates, std::vector<Solution*>& unique, Solution* x0);
size_t UniqueCandidates(std::vector<Solution*>& candidates, Solution* x0, bool isunique);
size_t UniqueCandidates(std::vector<Solution*> candidates, std::vector<Solution*>& unique);
void RMD(const std::vector<double*>& A, const std::vector<double>& b, const Solution* const X0, int n, int T, std::vector<Solution*>& candidates);
void SimpleRMD(const std::vector<double*>& A, const std::vector<double>& b, const Solution* const X0, int n, int T, std::vector<Solution*>& candidates);
// void RMD(Solution X0);
// Ax>=b;
// x0-- starting solution
// n -- number of candidate solution
// T -- warm-up period

#endif /*RMD_H_*/
