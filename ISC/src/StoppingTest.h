#ifndef STOPPINGTEST_H_
#define STOPPINGTEST_H_

#include <vector>
#include "Solution.h"
#include "Simulator.h"
/* The StoppingTest is given the current sample best solution indexed by starIndex in the list of 
 * visited solutions. The function then identifies the 1-neighborhood of the current sample best solution. 
 * The function then declares the current sample best as a local optimal solution with probability > 1-alpha
 * if g(xStar) <= min g(neighborhood of xStar), and thus stop the search; or declare xStar is not a local optimal
 * solution and notifies Compass that search should continue.
 * 
 * epsilon: parameters for stopping test, as used in epsilon-local optimizer definition.
 * alpha: type-one error and power of the stopping test.
 */ 
long StoppingTest(Simulator* simulator, bool stocsim, double epsilon, double alpha, std::vector<Solution*>& visited, long starIndex, long& budget, int mpamode);
long StoppingTest(Simulator* simulator, double epsilon, double alpha, std::vector<Solution*>& visited, const std::vector<long>& activeSolutions, long oldStarIndex, long starIndex, long& budget);
long StoppingTest(Simulator* simulator, bool stocsim, double epsilon, double alpha, std::vector<Solution*>& visited, long starIndex, const std::vector<long>& activeSolutions, long& budget, int mpamode);
#endif /*STOPPINGTEST_H_*/
