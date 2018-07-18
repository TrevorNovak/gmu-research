#ifndef SUBSTITUTIONINVENTORY_H_
#define SUBSTITUTIONINVENTORY_H_
// Operations research vol 49 no 3 334-351
// Determine inventory level given know prices with dynamic substitution.
#include "DES.h"
#include "RandomVariates.h"
#include "RngStream.h"
#include "require.h"
#include "myutilities.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
// int numOfVariants: number of product variants

class SubstitutionInventory : public Simulator
{
	private :
		// To use common random numbers, need two random number generators
		RngStream _rngpoisson;
		RngStream _rngutility;
		double substitutioninventory(const Solution* const solution);
		double ContinuousSubstitutionInventory(const double* const initialInventory, int numOfVariants);		
	public :
		SubstitutionInventory() : _rngpoisson("rngpoisson"), _rngutility("rngutility") {setSimMode(serial);}
		double simulation(const Solution* const solution);
		double gettruevalue(const Solution* const solution);
};
 
#endif /*SUBSTITUTIONINVENTORY_H_*/
