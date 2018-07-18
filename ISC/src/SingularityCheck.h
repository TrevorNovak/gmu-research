#ifndef SINGULARITYCHECK_H_
#define SINGULARITYCHECK_H_

#include "Solution.h"
#include <vector>

// For COMPASS, singular means there is only one solution in the MPA.
/* For hyperbox algorithm, singularity means there is only one interior solution in the MPA, although there may be boundary
   solutions inside this MPA. If that happens, this function will not only declare singularity, but make sure all of the 
   interior solution's neighbors have been visited. If not, it will visit these solutions one by one and add them to 
   list of visited solutions and update coordinate position map
   */

// A and b are MPA constraints. A0 and b0 are original problem constraints
bool SigularityCheck(std::vector<double*> A, std::vector<double> b, std::vector<double*> A0, std::vector<double> b0, Solution* X0, int mpamode, \
					 std::vector<Solution*>& visited, std::vector<long>& activeSolutions, long& budget, \
					 long& numOfTotalSamples, std::vector<long> visitedThisRun, Simulator* simulator, std::vector<CornerMap*> corners);

// Check if solution x is within the region 
bool CheckFeasibility(const Solution& x, std::vector<double*> A, std::vector<double> b);

#endif /*SINGULARITYCHECK_H_*/
