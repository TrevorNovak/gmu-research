#ifndef MPR_H_
#define MPR_H_

#include "masterheader.h"
#include <map>

// visited contains all visited solutions
// starIndex is the index of the current best solution in the vector "visited"
// Ak is the coefficients of the linear constraints, including MPR constraints and original problem formulation constraints
// bk is the RHS of constraints
// A is original problem formulation constraints, only linear constraints are considered
// b is the RHS of the original problem formulation constraints. 
// Constraints are Ak * x >= bk
// Ak, bk, W are the returned values by this function.
// activeSolutions has the indices of those solutions in "visited" that define MPR 

int ConstructMPR(std::vector<Solution*> visited, long starIndex, long oldStarIndex, std::vector<double*>& Ak, \
			std::vector<double>& bk, const std::vector<double*>& A, const std::vector<double>& b, \
			std::vector<long>& activeSolutions, const std::vector<long>& bestSolutions);

void MPRv2(bool stocsim, std::vector<Solution*> visited, long starIndex, long oldStarIndex, std::vector<double*>& Ak, \
			std::vector<double>& bk, const std::vector<double*>& A, const std::vector<double>& b, \
			std::vector<long>& activeSolutions, const std::vector<long>& bestSolutions, double alphac, \
			bool constraintPruning = true);

//12/27/11 adds a function to compute the volume of an MPA specified in Ax>=b constraints, right now only hyperbox supported
double calcVolMPA(const std::vector<double*>& A, const std::vector<double>& b, unsigned int dimension);

//12/28/11 computes the alignment prob of a blind picking rule with N total solutions, top G solutions, select S and with at least k among top G
double AlgnProb(double N, double G, double S, double k, int method);

//12/28/11 determines sample size for a blind picking OO 
double SampleSizeOO(double N, double G, double k);


#endif /*MPR_H_*/
