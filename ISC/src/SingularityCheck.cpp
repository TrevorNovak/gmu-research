#include "masterheader.h"
using namespace std;

// This function determines whether x0 is the only integer point in the region defined by constraints Ax>=b; 
// It returns 1 if it is, 0 o/w 
// mpamode determines how to check singularity. Right now I include original COMPASS and box constraint. 
// budget is the remaining simulation reps, numOfTotalSamples is the ISC consumed budget so far
extern RngStream rng;

bool CheckFeasibility(const Solution& x, vector<double*> A, vector<double> b)
{
	int dimension = x.getDimension();
	int numOfConstraints = (int) A.size();
	for (int k = 0; k < numOfConstraints; k++)	
	{	// For each constraint
		double sum = 0;
		for (int l = 0; l < dimension; l++)
		{	// Do a matrix multiplication
			sum += A[k][l] * x.getDiscreteSolution()[l];
		}
		if (b[k] > sum - EPSILON_STABILITY)
		{
			return false;
		}
	}
	return true;
}

bool SigularityCheck(vector<double*> A, vector<double> b, vector<double*> A0, vector<double> b0, Solution* X0, int mpamode, vector<Solution*>& visited,\
					 vector<long>& activeSolutions, long& budget, long& numOfTotalSamples, vector<long> visitedThisRun, Simulator* simulator, vector<CornerMap*> corners)
{
	vector<int> dimToSample;
	int dimension = X0->getDimension();
	int observations = 10;
	int numOfConstraints = (int) A.size();
	// flag indicates whether it is a singleton
	bool flag = true;
	for (int i = 0; i < dimension; i++)
	{
		double* b1 = new double[numOfConstraints];
		for (int k = 0; k < numOfConstraints; k++)	
		{	// For each constraint
			double sum = 0;
			for (int l = 0; l < dimension; l++)
			{	// Do a matrix multiplication
				if (l != i)
				{
					sum += A[k][l] * X0->getDiscreteSolution()[l];
				}
			}
			b1[k] = b[k] - sum;
		}
		// b1 has the values b - (A minus the ith column for xi)' * (x0 minus the ith element) 
		double upper = UPPER_BOUND, lower = UPPER_BOUND;
		for (int j = 0; j < numOfConstraints; j++)
		{	//for each constraint, see which one bounds the movement of the solution along this dimension
			double temp = 0;
			// temp is the temporary value of x_i to make the jth constraint tight
			if (dabs(A[j][i]) > EPSILON_STABILITY)
			{
				temp = b1[j] / A[j][i];
			}
			else
			{	// If the coefficient of x[i] for constraint A[j] is 0, then this variable has no effect on the
				// feasibility of solution, it can be adjusted arbitrarily. Set temp to a large number
				temp = 1 / EPSILON_STABILITY;
			}
			if (mpamode == BOXMPA)
			{	// For box constraints, singularity means the x0 is the only interior solution and all of its neighbors have been visited. 
				// If x0 is less than or equal to 1 away from the boundary, I consider x0 as the only interior solution. But to make it singular,
				// all of its neighbors need to be visited.
				if (dabs(temp - X0->getDiscreteSolution()[i]) > 1 + EPSILON_STABILITY)
				{
					if (temp > X0->getDiscreteSolution()[i])
					{
						//If the difference between the value to make the constraint tight and the value of the current point
						//is greater than "upper", it means that there is space "above" the current point, and the upper bound could be shrinked, until
						//the upper bound becomes the current point itself or cannot be smaller than 1.
						if (temp - X0->getDiscreteSolution()[i] < upper)
						{
							upper= temp - X0->getDiscreteSolution()[i];
						}
					}
					else
					{
						if (X0->getDiscreteSolution()[i] - temp < lower)
						{
							lower = X0->getDiscreteSolution()[i] - temp;
						}
					}
				}
				else
				{
					// If x0 is less than or equal to 1 from boundary along this coordinate direction, we consider it on the boundary. 
					// Now look at the sign of the coefficient, if it is positive, since its Ax>=b, can move upward as much as
					// we want. So set lower = 0, leave upper alone; else, reverse it. 
					// Doesn't have to worry about A[j][i] == 0 because if that's the case, then temp = large number, and will not
					// come to this place.
					if (A[j][i] > 0)
						lower = 0;
					else
						upper = 0;
				}
			}
			else // Default is COMPASS mode
			{
				if (dabs(temp - X0->getDiscreteSolution()[i]) > EPSILON_STABILITY) 					
				{	//First make sure that this constraint is not tight at X0->getDiscreteSolution()[i]. If it is tight, it means x0 is on the boundary.
					//Moving along both directions is still inside feasible region if the current solution is
					//in the interior. But a solution on the boundary is different, only one direction can be followed.
					if (temp > X0->getDiscreteSolution()[i])
					{
						//If the difference between the value to make the constraint tight and the value of the current point
						//is greater than "upper", it means that there is space "above" the current point, and the upper bound could be shrinked, until
						//the upper bound becomes the current point itself or cannot be smaller than 1.
						if (temp - X0->getDiscreteSolution()[i] < upper)
						{
							upper= temp - X0->getDiscreteSolution()[i];
						}
					}
					else
					{
						if (X0->getDiscreteSolution()[i] - temp < lower)
						{
							lower = X0->getDiscreteSolution()[i] - temp;
						}
					}
				}
				else
				{
					// temp and X0->getDiscreteSolution()[i] are equal, so the current solution is right on this constraint
					// Now look at the sign of the coefficient, if it is positive, since its Ax>=b, can move upward as much as
					// we want. So set lower = 0, leave upper alone; else, reverse it. 
					// Doesn't have to worry about A[j][i] == 0 because if that's the case, then temp = large number, and will not
					// come to this place.
					if (A[j][i] > 0)
						lower = 0;
					else
						upper = 0;
				}
			}
		}
		// if the upper or the lower boundary is greater than 1 from the current point along this coordinate position, it is not a singleton yet. 
		if (upper >= 1 || lower >= 1)
        {
        	flag = false;
			break;
			delete []b1;
        }
        delete []b1;
	}
	// Now need to make sure all of x0's neighbors have been visited for the box mode. It won't hurt to do so for COMPASS, too.
	if (mpamode == BOXMPA && flag)
	{
		/*cout << *X0;
		printVector(A, dimension);
		printVector(b);
		cout << endl;*/
		// Note that total neighbors is not 2*dimension, they need to be feasible
		int neighborvisited = 0, totalneighbors = 0;
		dimToSample.clear();
		for (int i = 0; i < dimension; i++)
		{
			bool tovisit = false;
			Solution* x0nbr1 = new Solution(dimension, X0->getDiscreteSolution());
			x0nbr1->changeOneDimension(*(X0->getDiscreteSolution()+i)+1, i);
			// Check the feasibility of this neighbor
			if (CheckFeasibility(*x0nbr1, A0, b0))
			{
				tovisit = true;
				totalneighbors++;
				long index = findSolution(visited, x0nbr1);
				if (index < 0)
				{	
					delete x0nbr1;
				}
				else
				{
					neighborvisited++;
					delete x0nbr1;
				}
			}
			Solution* x0nbr2 = new Solution(dimension, X0->getDiscreteSolution());
			x0nbr2->changeOneDimension(*(X0->getDiscreteSolution()+i)-1, i);
			if (CheckFeasibility(*x0nbr2, A0, b0))
			{
				tovisit = true;
				totalneighbors++;
				long index = findSolution(visited, x0nbr2);
				if (index < 0)
				{	
					delete x0nbr2;
				}
				else
				{
					delete x0nbr2;
					neighborvisited++;
				}
			}
			if (tovisit)
			{
				dimToSample.push_back(i);
			}
		}
		// Now we know how many feasible neighbors have been visited. We invoke stopping test when all feasible
		// neighbors have been visited. If not, we randomly select a direction that has unsampled neighbors
		if (neighborvisited < totalneighbors)
		{
			flag = false;
			int directionToMove = dimToSample[rng.RandInt(0, dimToSample.size()-1)];
			//for (int i = 0; i < dimension; i++)
			{
				Solution* x0nbr1 = new Solution(dimension, X0->getDiscreteSolution());
				//x0nbr1->changeOneDimension(*(X0->getDiscreteSolution()+i)+1, i);
				x0nbr1->changeOneDimension(*(X0->getDiscreteSolution()+directionToMove)+1, directionToMove);
				// Check the feasibility of this neighbor
				if (CheckFeasibility(*x0nbr1, A0, b0))
				{
					long index = findSolution(visited, x0nbr1);
					if (index < 0)
					{	// If one of its neighbors has not been visited, add it to visited
						// To avoid unnecessary complication with stopping test procedure, require the new solution has at least one rep. 
						for (int newrepi = 0; newrepi < observations; newrepi++)
						{
							x0nbr1->recordObservation(simulator->simulation(x0nbr1));
						}
						numOfTotalSamples+=observations;
						budget-=observations;
						visited.push_back(x0nbr1);	
						visitedThisRun.push_back(long(visited.size() - 1));
						activeSolutions.push_back(long(visited.size() - 1));
						// Also need to update the coordinate position map
						for (int j = 0; j < dimension; j++)
						{	// The index in the visited vector is the content. The key is its coordinate position
							corners[j]->insert(make_pair(*(visited[long(visited.size() - 1)]->getDiscreteSolution() + j), long(visited.size() - 1)));
						}
					}
					else
					{
						delete x0nbr1;
					}
				}
				Solution* x0nbr2 = new Solution(dimension, X0->getDiscreteSolution());
				x0nbr2->changeOneDimension(*(X0->getDiscreteSolution()+directionToMove)-1, directionToMove);
				if (CheckFeasibility(*x0nbr2, A0, b0))
				{
					long index = findSolution(visited, x0nbr2);
					if (index < 0)
					{	// If one of its neighbors has not been visited, add it to visited
						// To avoid unnecessary complication with stopping test procedure, require the new solution has at least one rep. 
						for (int newrepi = 0; newrepi < observations; newrepi++)
						{	
							x0nbr2->recordObservation(simulator->simulation(x0nbr2));
						}
						numOfTotalSamples+=observations;
						budget-=observations;
						visited.push_back(x0nbr2);	
						visitedThisRun.push_back(long(visited.size() - 1));
						activeSolutions.push_back(long(visited.size() - 1));
						for (int j = 0; j < dimension; j++)
						{	// The index in the visited vector is the content. The key is its coordinate position
							corners[j]->insert(make_pair(*(visited[long(visited.size() - 1)]->getDiscreteSolution() + j), long(visited.size() - 1)));
						}
					}
					else
					{
						delete x0nbr2;
					}
				}
			}
		}
		else
		{
			flag = true;
		}
	}
	return flag;
}


