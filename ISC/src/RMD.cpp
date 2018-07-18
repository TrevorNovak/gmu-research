#include "RMD.h"
#include "RngStream.h"
#include <math.h>
#include "myutilities.h"
#include <iostream>

using namespace std;
// void RMD(Solution X0);
// Ax>=b;
// x0-- starting solution
// n -- number of candidate solution
// T -- warm-up period
extern RngStream rng;

void RMD(const vector<double*>& A, const vector<double>& b, const Solution* const X0, int n, int T, vector<Solution*>& candidates)
{
	int dimension = X0->getDimension();
	int numOfConstraints = (int) A.size();
	// Debug use
	//printVector(A, dimension);
	//copy(b.begin(), b.end(), ostream_iterator<double>(cout, ", "));
	// Start from the current sample best.
	vector<Solution*> newSolutions;
	Solution* newSolution0 = new Solution(dimension);
	newSolutions.push_back(newSolution0);
	newSolutions[0]->setSolution(X0->getDiscreteSolution());
	for (int i = 0; i < n; i++)
	{	// Whenever a new solution is pushed into the vector candidate, a new solution is created and 
		// initially it has the same content as the solution just pushed into the vector.
		if (i > 0)
		{
			Solution* newSolutioni = new Solution(dimension);
			newSolutions.push_back(newSolutioni);
			newSolutions[i]->setSolution(newSolutions[i-1]->getDiscreteSolution());
		}
		for (int j = 0; j < T; j++)
		{	// To warm up
			// Randomly pick up a dimension to move along
			int directionToMove = rng.RandInt(0, dimension-1);
			double* b1 = new double[numOfConstraints];
			for (int k = 0; k < numOfConstraints; k++)	
			{// For each constraint
				double sum = 0;
				for (int l = 0; l < dimension; l++)
				{	// Do a matrix multiplication
					if (l != directionToMove)
					{
						sum += A[k][l] * newSolutions[i]->getDiscreteSolution()[l];
					}
				}
				b1[k] = b[k] - sum;
			}
			// Now check which constraint is tight
			double upper = UPPER_BOUND, lower = UPPER_BOUND;
			for (int k = 0; k < numOfConstraints; k++)
			{
				double temp = 0;
				// temp is the temporary value of x_i to make the jth constraint tight
				if (dabs(A[k][directionToMove]) > EPSILON_STABILITY)
				{
					temp = b1[k] / A[k][directionToMove];
				}
				else
				{
					temp = LARGE_NUMBER;
				}
				if (temp > newSolutions[i]->getDiscreteSolution()[directionToMove] + EPSILON_STABILITY)
				{
					//If the value to make the constraint tight is greater than the value of the current point,
					//it means that there is space "above" the current point, and the upper bound could be shrinked, until
					//the upper bound becomes the current point itself or cannot be smaller than 1.
					if (temp - newSolutions[i]->getDiscreteSolution()[directionToMove] < upper)
					{
						upper= temp - newSolutions[i]->getDiscreteSolution()[directionToMove];
					}
				}
				else if (temp < newSolutions[i]->getDiscreteSolution()[directionToMove] - EPSILON_STABILITY)
				{
					if (newSolutions[i]->getDiscreteSolution()[directionToMove] - temp < lower)
					{
						lower = newSolutions[i]->getDiscreteSolution()[directionToMove] - temp;
					}
				}
				else
				{
					// The constraint is already tight at current value, i.e., the point is now on the boundary. !!!!!!!!!!!
					// If the coefficient is positive, then increasing, i.e., moving "up" will reenter feasible region because the 
					// inequalitys are Ax>=b
					if (A[k][directionToMove] > 0)
					{
						lower = 0;
					}
					else
					{
						// Don't need to worry about A[k][directionToMove] = 0, because in that case temp will be a large number
						upper = 0;
					}
				}
			}
			int int_maxX_directionToMove = (int) floor(upper) + newSolutions[i]->getDiscreteSolution()[directionToMove];
			int int_minX_directionToMove = newSolutions[i]->getDiscreteSolution()[directionToMove] - (int) floor(lower);
			int length = int_maxX_directionToMove - int_minX_directionToMove;
			int step = rng.RandInt(0, length);
			newSolutions[i]->changeOneDimension(int_minX_directionToMove + step, directionToMove);
			delete []b1;
		}	
		// Warmup over, take this solution as a candidate
		candidates.push_back(newSolutions[i]);
	}
}


// Simple RMD samples along the coordinates of X0. 

void SimpleRMD(const vector<double*>& A, const vector<double>& b, const Solution* const X0, int n, int T, vector<Solution*>& candidates)
{
	int dimension = X0->getDimension();
	int numOfConstraints = (int) A.size();
	// Start from the current sample best.
	vector<Solution*> newSolutions;
	for (int i = 0; i < n; i++)
	{	// Whenever a new solution is pushed into the vector candidate, a new solution is created and 
		// initially it has the same content as the solution just pushed into the vector.
		Solution* newSolutioni = new Solution(dimension);
		newSolutions.push_back(newSolutioni);
		newSolutions[i]->setSolution(X0->getDiscreteSolution());
		// Randomly pick up a dimension to move along
		int directionToMove = rng.RandInt(0, dimension-1);
		double* b1 = new double[numOfConstraints];
		for (int k = 0; k < numOfConstraints; k++)	
		{// For each constraint
			double sum = 0;
			for (int l = 0; l < dimension; l++)
			{	// Do a matrix multiplication
				if (l != directionToMove)
				{
					sum += A[k][l] * newSolutions[i]->getDiscreteSolution()[l];
				}
			}
			b1[k] = b[k] - sum;
		}
		// Now check which constraint is tight
		double upper = UPPER_BOUND, lower = UPPER_BOUND;
		for (int k = 0; k < numOfConstraints; k++)
		{
			double temp = 0;
			// temp is the temporary value of x_i to make the jth constraint tight
			if (dabs(A[k][directionToMove]) > EPSILON_STABILITY)
			{
				temp = b1[k] / A[k][directionToMove];
			}
			else
			{
				temp = LARGE_NUMBER;
			}
			if (temp > newSolutions[i]->getDiscreteSolution()[directionToMove] + EPSILON_STABILITY)
			{
				//If the value to make the constraint tight is greater than the value of the current point,
				//it means that there is space "above" the current point, and the upper bound could be shrinked, until
				//the upper bound becomes the current point itself or cannot be smaller than 1.
				if (temp - newSolutions[i]->getDiscreteSolution()[directionToMove] < upper)
				{
					upper= temp - newSolutions[i]->getDiscreteSolution()[directionToMove];
				}
			}
			else if (temp < newSolutions[i]->getDiscreteSolution()[directionToMove] - EPSILON_STABILITY)
			{
				if (newSolutions[i]->getDiscreteSolution()[directionToMove] - temp < lower)
				{
					lower = newSolutions[i]->getDiscreteSolution()[directionToMove] - temp;
				}
			}
			else
			{
				// The constraint is already tight at current value, i.e., the point is now on the boundary. !!!!!!!!!!!
				// If the coefficient is positive, then increasing, i.e., moving "up" will reenter feasible region because the 
				// inequalitys are Ax>=b
				if (A[k][directionToMove] > 0)
				{
					lower = 0;
				}
				else
				{
					// Don't need to worry about A[k][directionToMove] = 0, because in that case temp will be a large number
					upper = 0;
				}
			}
		}
		int int_maxX_directionToMove = (int) floor(upper) + newSolutions[i]->getDiscreteSolution()[directionToMove];
		int int_minX_directionToMove = newSolutions[i]->getDiscreteSolution()[directionToMove] - (int) floor(lower);
		int length = int_maxX_directionToMove - int_minX_directionToMove;
		int step = rng.RandInt(0, length);
		newSolutions[i]->changeOneDimension(int_minX_directionToMove + step, directionToMove);
		delete []b1;
		// Warmup over, take this solution as a candidate
		candidates.push_back(newSolutions[i]);
	}
}


// candidates has the newly sampled solutions returned from RMD. They may have repetitive entries, and so 
// this function checks and removes those redundnant copies. Finally "unique" has all unique new solutions.
size_t UniqueCandidates(vector<Solution*> candidates, vector<Solution*>& unique, Solution* x0)
{
	for (int i = 0; i < (int) candidates.size(); i++)
	{
		bool isUnique = true;
		for (int j = 0; j < (int) unique.size(); j++)
		{
			// check if the solution being inspected already exists in the unique list
			if (candidates[i]->isEqual(unique[j]))
			{
				isUnique = false;
				break;
			}
		}
		if (candidates[i]->isEqual(x0))
		{
			isUnique = false;
		}
		if (isUnique)
		{	// Keep this solution
			unique.push_back(candidates[i]);
		}
		else
		{	// This entry is redundant, delete it.
			delete candidates[i];
		}
	}
	return unique.size();
}

size_t UniqueCandidates(vector<Solution*>& candidates, Solution* x0, bool isunique) // This one excludes solution x0 from the list of candidates, which can be unique already if "isunique"=true
{
	for (size_t i = 0; i < candidates.size(); i++)
	{
		if (candidates[i]->isEqual(x0))
		{
			delete candidates[i];
			candidates.erase(candidates.begin()+i);
			if (isunique) break;
		}
	}
	return candidates.size();
}


size_t UniqueCandidates(vector<Solution*> candidates, vector<Solution*>& unique)
{
	for (int i = 0; i < (int) candidates.size(); i++)
	{
		bool isUnique = true;
		for (int j = 0; j < (int) unique.size(); j++)
		{
			// check if the solution being inspected already exists in the unique list
			if (candidates[i]->isEqual(unique[j]))
			{
				isUnique = false;
				break;
			}
		}
		if (isUnique)
		{	// Keep this solution
			unique.push_back(candidates[i]);
		}
		else
		{	// This entry is redundant, delete it.
			delete candidates[i];
		}
	}
	return unique.size();
}



