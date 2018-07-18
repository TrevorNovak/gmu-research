#ifndef RESPONSESURFACE_H_
#define RESPONSESURFACE_H_

#include "masterheader.h"

using std::vector;
 
class ResponseSurface : public Simulator
{
	private :
		int _rsfid; //indicate which function to use
		double F22(const Solution* const solution, bool isNoisy);
		double F2n(const Solution* const solution, bool isNoisy);
		double rosenbrock(const Solution* const solution, bool isNoisy);
		double singular(const Solution* const solution, bool isNoisy);
		double quadratic(const Solution* const solution, bool isNoisy);
		double highd(const Solution* const solution, bool isNoisy);
		double MultiLocalHighd(const Solution* const solution, bool isNoisy);
		double rsf(const Solution* const solution, bool isNoisy);
	public :
		ResponseSurface(int id) : _rsfid(id) {setSimMode(serial);}
		~ResponseSurface() {}
		double simulation(const Solution* const solution)	{return rsf(solution, true);}
		double gettruevalue(const Solution* const solution)	{return rsf(solution, false);}
		int isGlobalOptimum(const Solution& solution);
}; 

#endif /*RESPONSESURFACE_H_*/
