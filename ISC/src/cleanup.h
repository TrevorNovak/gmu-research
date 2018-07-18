#ifndef CLEANUP_H_
#define CLEANUP_H_
#include "masterheader.h"

/* Given the list of visited solutions (visited), a list of local minimum (localmin), a constant delta, a 
 * constant alpha, cleanup() performs a NSGC procedure and report the best minimum.
 * Return value is the index of the best in visited.
 * Note that this cleanup procedure is meant for maximization problems instead. To use it for minimization
 * problem, in the begining of the procedure, need to negate the sample means and work on 
 * negated sample means instead.
 */
long cleanup(Simulator* simulator, std::vector<Solution*>& visited, std::vector<long> M, double delta, double alpha);

#endif /*CLEANUP_H_*/
