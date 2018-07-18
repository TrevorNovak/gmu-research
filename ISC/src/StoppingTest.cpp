#include "StoppingTest.h"
#include <cmath>
#include "myutilities.h"
#include "Simulator.h"
using namespace std;
/* The StoppingTest is given the current sample best solution indexed by starIndex in the list of 
 * visited solutions. The function then identifies the 1-neighborhood of the current sample best solution. 
 * The function then declares the current sample best as a local optimal solution with probability > 1-alpha
 * if g(xStar) <= min g(neighborhood of xStar), and thus stop the search; or declare xStar is not a local optimal
 * solution and notifies Compass that search should continue.
 * 
 * epsilon: parameters for stopping test, as used in epsilon-local optimizer definition.
 * alpha: type-one error and power of the stopping test.
 * 
 * Return value:  the index of the optimal solution in the vector visited. If it is the same as starIndex, it means the current
 * sample best is locally optimal.
 */

long StoppingTestProcedure(Simulator* simulator, double epsilon, double alpha, vector<Solution*>& visited, const vector<long>& solutionsToCompare, long& budget, bool stocsim = true ); 
long StoppingTest(Simulator* simulator, bool stocsim, double epsilon, double alpha, vector<Solution*>& visited, long starIndex, long& budget, int mpamode)
{
	// Setup: identifies the neighborhood of the current sample best solution.
	if (budget < 0)
	{
		throw long(-1); 
	}
	vector<long> solutionsToCompare;
	int dimension = visited[starIndex]->getDimension();
	solutionsToCompare.push_back(starIndex);
	for (long i = 0; i < (long) visited.size(); i++)
	{
		if (i != starIndex)
		{
			double distance = visited[starIndex]->getDistanceTo(visited[i]);
			if (distance <= 1)
			{
				// This solution is within the neighborhood of xStar
				solutionsToCompare.push_back(i);
			}
		}
	}
	// Debug use.
	if (solutionsToCompare.size() != 2*dimension + 1)
	{
		cout << "Xstar" << endl << *visited[starIndex];
		PrintActiveSolutions(visited, solutionsToCompare);
	}
	return StoppingTestProcedure(simulator, epsilon, alpha, visited, solutionsToCompare, budget, stocsim);
}

long StoppingTest(Simulator* simulator, bool stocsim, double epsilon, double alpha, vector<Solution*>& visited, long starIndex, const vector<long>& activeSolutions, long& budget, int mpamode)
{
	// Setup: identifies the neighborhood of the current sample best solution.
	if (budget < 0)
	{
		throw long(-1);
	}
	vector<long> solutionsToCompare;
	solutionsToCompare.push_back(starIndex);
	int dimension = visited[starIndex]->getDimension();
	for (long i = 0; i < (long) activeSolutions.size(); i++)
	{
		if (i != starIndex)	solutionsToCompare.push_back(activeSolutions[i]);
	}
	// Debug use.
	/*if (solutionsToCompare.size() < 2*dimension + 1)
	{
		PrintActiveSolutions(visited, solutionsToCompare);
	}*/
	return StoppingTestProcedure(simulator, epsilon, alpha, visited, solutionsToCompare, budget, stocsim);
}

// An overloaded version for step 3 in LCRS.
/*long StoppingTest(Simulator* simulator, double epsilon, double alpha, vector<Solution*>& visited, const vector<long>& activeSolutions, long oldStarIndex, long starIndex, long& budget)
{
	// Setup: identifies the neighborhood of the current sample best solution.
	if (budget < 0)
	{
		throw long(-1);
	}
	vector<long> solutionsToCompare;
	// starIndex is the standard, as indicated by x0 in the procedure. so it takes the first position
	solutionsToCompare.push_back(starIndex);
	for (long count = 0; count < (long) activeSolutions.size(); count++)
	{
		solutionsToCompare.push_back(activeSolutions[count]);
	}
	solutionsToCompare.push_back(oldStarIndex);
	return StoppingTestProcedure(simulator, epsilon, alpha, visited, solutionsToCompare, budget);
}*/


long StoppingTestProcedure(Simulator* simulator, double epsilon, double alpha, vector<Solution*>& visited, const vector<long>& solutionsToCompare, long& budget, bool stocsim)
{
	if (!stocsim)
	{ // Deterministic case, simply find the minimum
		return findMinSolution(visited, solutionsToCompare);
	}
	// solutionsToCompare[0] is the standard
	// Compute parameters
	// A CATCH: here p is actually p+1 in the Coordinate Sampling paper by Jeff L Hong
	long p = long(solutionsToCompare.size());
	//cout << "Compare these solutions\n";
	//PrintActiveSolutions(visited, solutionsToCompare);
	long numOfSamples = 0;
	double* S2 = new double[p*p];
	for (long i = 0; i < p; i++)
	{
		for (long j = 0; j < p; j++)
		{ 
			S2[i*p + j] = visited[solutionsToCompare[i]]->getSampleVariance() + visited[solutionsToCompare[j]]->getSampleVariance();
		}
	}
	//printVector(S2, p, "S2");
	int* f = new int[p*p];
	for (long i = 0; i < p; i++)
	{
		for (long j = 0; j < p; j++)
		{
			f[i*p + j] = min(visited[solutionsToCompare[i]]->getNumOfObservations(), visited[solutionsToCompare[j]]->getNumOfObservations()) - 1;
		}
	}
	//printVector(f, p, "f");
	double* epsilonij = new double[p*p];
	for (long i = 0; i < p; i++)
	{
		for (long j = 0; j < p; j++)
		{
			if (i == 0 || j == 0)
				epsilonij[i*p + j] = epsilon / 2;
			else
				epsilonij[i*p + j] = epsilon;
		}
	}
	//printVector(epsilonij, p, "epsilon");
	double* eta = new double[p*p];
	for (long i = 0; i < p; i++)
	{
		for (long j = 0; j < p; j++)
		{
			// Note p here is the size of all Iset to compare, but p in the paper
			// of this stopping test procedure means the systems to be compared with the best
			// so p-1 should be used here
			eta[i*p + j] = pow((p-1)/alpha/2.0, 2.0/double(f[i*p + j])) - 1;
		}
	}
	//printVector(eta,p, "eta");
	double* lambda = new double[p*p];
	for (long i = 0; i < p; i++)
	{
		for (long j = 0; j < p; j++)
		{
			lambda[i*p + j] = epsilonij[i*p + j] / 2;
		}
	}
	//printVector(lambda,p, "lambda");
	double* a = new double[p*p];
	for (long i = 0; i < p; i++)
	{
		for (long j = 0; j < p; j++)
		{
			a[i*p + j] =  eta[i*p + j] * f[i*p + j] * S2[i*p + j] / epsilonij[i*p + j] / 2;
		}
	}
	//printVector(a,p, "a");
	/*
	 *  Initialization
	 */
	// Let Iset be the set of solutions still in contention, which is denoted as I in the CSA paper. 
	// Iset has the indexes of those solutions with respect to the list solutionsToCompare. So Iset should always
	// have values between 0 and p. To actually index the solution, use visited[solutionsToCompare[Iset[i]]].
	vector<long> Iset;
	vector<long> oldIset;
	for (long i = 0; i < p; i++)
	{
		Iset.push_back(i);
	}
	// The step counter r = min(n_00, n_10, ... n_p0)
	int stepCounter = visited[solutionsToCompare[Iset[0]]]->getNumOfObservations();
	for(long i = 1; i < p; i++)
	{
		if (stepCounter > visited[solutionsToCompare[Iset[i]]]->getNumOfObservations())
		{
			stepCounter = visited[solutionsToCompare[Iset[i]]]->getNumOfObservations();
		}
	}
	// observation counter
	vector<long> observationCounter; 
	long returnValue = solutionsToCompare[0];
	long iteration = 0;
	/* ***********************************************************************************
	 * Screening procedure
	 * ***********************************************************************************
	 */
	// I try to avoid any possibility of a dead loop. That's why I use here a maximum number of iterations.
	while(0 < budget)
	{
		// Calculate those "H" values as in Jeff's procedure
		iteration++;
		if (iteration % 5000 == 0)
		{ 
			cout << "Iteration " << iteration << endl;
			for (size_t i = 0; i < Iset.size (); i++)
			{
				cout << *visited[solutionsToCompare[Iset[i]]];
			}
			cout << endl;
		}
		double* H = new double[Iset.size()];
		// Note here the index for i actually is between 0 and p, that is, the indexes for the remaining solutions
		// in the list solutionsToCompare. 
		for (long i = 0; i < (long) Iset.size(); i++)
		{
			if (Iset[i] == 0)
			{
				H[0] = -visited[solutionsToCompare[Iset[0]]]->getSampleMean() + epsilon / 2;
			}
			else
			{
				H[i] = -visited[solutionsToCompare[Iset[i]]]->getSampleMean();
			}
		}
		// Iset is I, oldIset is I^{old}
		// Debug use
		/*if (iteration > 100)
		{
			PrintActiveSolutions(visited, Iset);
			printcontainer(H, H+Iset.size());
			for (long i = 0; i < (long) Iset.size(); i++)
			{
				for (long j = 0; j < (long) Iset.size(); j++)
				{
					cout << stepCounter * (H[i] - H[j]) << ", " << a[i*p + j] - stepCounter * lambda[i*p + j] << endl;
				}
			}
		}*/
		vectorCopy(Iset, oldIset);
		Iset.clear();
		observationCounter.clear();
		// Identify the new "I" set
		for (long i = 0; i < (long) oldIset.size(); i++)
		{
			bool flag = true;
			for (long j = 0; j < (long) oldIset.size(); j++)
			{	// For the ith and jth solutions in the list oldIset, their indexes in visited are solutionsToCompare[oldIset[i]]
				// Their indexes in the parameter matrices eta, a..., are given by oldIset[i], as a value between 0 and p
				// We should not use i and j directly to access solution or parameters
				double temp1 = stepCounter * (H[i] - H[j]);
				double temp2 = a[oldIset[i]*p + oldIset[j]];
				double temp3 = stepCounter * lambda[oldIset[i]*p + oldIset[j]];
				if ((i != j) && (stepCounter * (H[i] - H[j]) < -max(0, a[oldIset[i]*p + oldIset[j]] - stepCounter * lambda[oldIset[i]*p + oldIset[j]])))
				{
					flag = false;
					break;
				}
			}
			if (flag)
			{
				Iset.push_back(oldIset[i]);
				observationCounter.push_back(visited[solutionsToCompare[oldIset[i]]]->getNumOfObservations());
			}
		}
		/* ******************************************************************************************
		 * Stopping rule
		 * *****************************************************************************************
		 */
		double minimum = visited[solutionsToCompare[Iset[0]]]->getSampleMean();
		long minX = 0;
		for (long i = 1; i < (long) Iset.size(); i++)
		{
			if (minimum > visited[solutionsToCompare[Iset[i]]]->getSampleMean())
			{
				minimum = visited[solutionsToCompare[Iset[i]]]->getSampleMean();
				minX = i;
			}
		}
		if (minimum < visited[solutionsToCompare[0]]->getSampleMean())
		{
			returnValue = solutionsToCompare[Iset[minX]];
			double oldstarobj = visited[solutionsToCompare[0]]->getSampleMean();
			// The current sample best is actually worse than a neighbor, return the index of that neighbor so that search can continue there.
			break;
		}
		else if (Iset.size() == 1)
		{
			// The current sample best declared as the local optimum.
			returnValue = solutionsToCompare[0];
			double oldstarobj = visited[solutionsToCompare[0]]->getSampleMean();
			break;
		}
		else
		{	// Take additional reps for remaining solutions indexed by Iset, and go back to screening
			for (long i = 0; i < (long) Iset.size(); i++)
			{
				if (observationCounter[i] == stepCounter)
				{
					visited[solutionsToCompare[Iset[i]]]->recordObservation(simulator->simulation(visited[solutionsToCompare[Iset[i]]]));
					++numOfSamples;
					if (--budget <= 0)
					{
						returnValue = solutionsToCompare[0];
						throw numOfSamples;
					}
				}
			}
			stepCounter++;	
		}
		delete []H;
	}
	delete []S2;
	delete []f;
	delete []epsilonij;
	delete []eta;
	delete []lambda;
	delete []a;
	//cout << "Returned as the winner of the stopping test" << endl;
	//cout << *visited[returnValue];
	//cout << endl << "Now the solutions become " << endl;
	//PrintActiveSolutions(visited, solutionsToCompare);
	return returnValue;
}



