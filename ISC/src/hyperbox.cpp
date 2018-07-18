#include "masterheader.h"
#include "hyperbox.h"
#include <algorithm>

using namespace std;
/* 
stocsim: whether it is a stochastic simulation and thus should keep on adding replications to solutions defining the box
visited: the vector of all visited solutions
starIndex: the current sample best, the box will encompass it
oldStarIndex: may be useful to reduce the number of solutions that need to be simulated further in future work. Now it is not useful
Ak: The constraint coefficients defining the box. This is the output
bk: the rhs for those constraints. This is output
A: original problem constraints
b: original constraints rhs
activeSolutions: vector of solutions defining the hyperbox. This is the output
corners: the multimap containing coordinate positions for all visited solutions
*/

void hyperbox(bool stocsim, const std::vector<Solution*>& visited, long starIndex, long oldStarIndex, std::vector<double*>& Ak, \
			std::vector<double>& bk, const std::vector<double*>& A, const std::vector<double>& b, \
			std::vector<long>& activeSolutions, const std::vector<CornerMap*>& corners, bool boxalone)
{
	// Each time we construct a hyperbox, the previous one, defined by Ak, bk are no longer useful. 
	// Release memory allocated to arrays pointed by pointers in Ak
	// But be very, very careful, Ak also contains the pointers in A. If they are deleted, original problem
	// constraints are lost and memory problems occur!!!!!
	//static long iterationcount = 0;
	//clock_t beginning = clock();
	//static ofstream lpsizefile("../lpsize.dat");
	unsigned int dimension = visited[0]->getDimension();
	pair<CMIter, CMIter> ptt;
	if (boxalone)
	{
		for (size_t i = A.size(); i < Ak.size(); i++)
		{
			delete []Ak[i];
		}
		activeSolutions.clear();
		// First the original problem formulation constraints
		// copy the contents of A into Ak
		vectorCopy(A, Ak);
		vectorCopy(b, bk);
	}
	// In the future, I may use oldactive solutions to help reduce the amount of working in search for the new hyperbox. It is not used for now.
	vector<long> oldActiveSolutions;
	// for each coordinate d, find xup_d and xlow_d that are closest to xstar_d.If one or both of xup_d and xlow_d do not exist. It is still ok because
	// we have original problem constraints to bound the search region.
	for (size_t i = 0; i < dimension; i++)
	{	// First find the key, which is the coordinate position of xstar
		int key = *(visited[starIndex]->getDiscreteSolution() + i);
		/* lower_bound() returns the iterator to the smallest element whose key is 
		   equal to or BIGGER than the argument "key". Decreasing it will give the largest element 
		   with a key smaller than the argument "key", if such an element exists. 
		*/
		CMIter cmit = corners[i]->lower_bound(key);
		if (cmit == corners[i]->end())
		{
			cout << "Error. Cannot find coord " << i << " in the coordinate position map\n" << *visited[starIndex];
			for (cmit = corners[i]->begin(); cmit != corners[i]->end(); cmit++)
			{
				cout << *visited[(*cmit).second];
			}
			exit(2);
		}
		if ((*cmit).first == key)
		{	// This should always be the case because at least xstar has this key.
			// But if there is no element with a smaller key, that's fine, we simply don't set up a constraint along this direction.
			if (cmit != corners[i]->begin())
			{	// So we fetch the previous item and use its key to do a search for all solutions with this key, if there is such a key smaller than the key of xstar
				cmit--;
				/*ptt = corners[i]->equal_range((*cmit).first);
				//printmap(ptt.first, ptt.second);
				// Write Ak, bk and activesolutions
				for (CMIter activeit = ptt.first; activeit != ptt.second; activeit++)
				{	// Now I treat every solution with the same key as an active solution. This is an overkill. Needs to be refined in the future.
					activeSolutions.push_back((*activeit).second);
				}*/

				/*	Aug 24 2009. We do not need to make solutions acroos the edges of the hyperbox active solutions.
					They are only needed to define the hyperbox. What we need to simulate are the new solutions and 
					the neighbors of current sample best that have been active solutions before. So experiment with not
					adding these edge solutions as active solutions and add neighbors that have been visited as active solutions.
				*/
				/* Aug 31, 2009, do not add neighbors that have been visited as active solutions. This step may not
				   be necessary for convergence and slows down AHA too much
				   */
				/*for (int dim = 0; dim < dimension; dim++)
				{
					Solution* x0nbr1 = new Solution(dimension, visited[starIndex]->getDiscreteSolution());
					x0nbr1->changeOneDimension(*(visited[starIndex]->getDiscreteSolution()+dim)+1, dim);
					// Check the feasibility of this neighbor
					if (CheckFeasibility(*x0nbr1, A, b))
					{	// If a feasible neighbor has been visited, it must have been in the estimation set 
						// before. So add it to the active solution set.
						long index = findSolution(visited, x0nbr1);
						if (index < 0)
						{
							delete x0nbr1;
						}
						else
						{
							delete x0nbr1;
							activeSolutions.push_back(index);
						}
					}
					Solution* x0nbr2 = new Solution(dimension, visited[starIndex]->getDiscreteSolution());
					x0nbr2->changeOneDimension(*(visited[starIndex]->getDiscreteSolution()+dim)-1, dim);
					if (CheckFeasibility(*x0nbr2, A, b))
					{
						long index = findSolution(visited, x0nbr2);
						if (index < 0)
						{	
							delete x0nbr2;
						}
						else
						{
							delete x0nbr2;
							activeSolutions.push_back(index);
						}
					}
				}*/
				activeSolutions.push_back((*cmit).second);
				// Debug use
				/*if (i == 0) 
				{
					cout << visited[(*cmit).second] << endl;
				}*/
				double* Atemp = new double[dimension];
				fill(Atemp, &Atemp[dimension], 0);
				Atemp[i] = 1;
				// The key is the coordinate position. So it is the rhs of the constraint
				//double btemp = (key+(*cmit).first)/2.0;
				double btemp = (*cmit).first;
				Ak.push_back(Atemp);
				bk.push_back(btemp);
			}
		}
		else
		{
			// This should not happen. It is a problem
			cout << "Problem in constructing hyperbox" << endl;
			printmap(corners[i]->begin(), corners[i]->end());
			exit(2);
		}
		// upper_bound returns the smallest element whose key is bigger than (excluding equal to) "key", if such an element exists
		cmit = corners[i]->upper_bound(key);
		if (cmit != corners[i]->end())
		{
			/*ptt = corners[i]->equal_range((*cmit).first);
			//printmap(ptt.first, ptt.second);
			for (CMIter activeit = ptt.first; activeit != ptt.second; activeit++)
			{	// Now I treat every solution with the same key as an active solution. This is an overkill. Needs to be refined in the future.
				activeSolutions.push_back((*activeit).second);
			}*/
			/*	Aug 24 2009. We do not need to make solutions acroos the edges of the hyperbox active solutions.
					They are only needed to define the hyperbox. What we need to simulate are the new solutions and 
					the neighbors of current sample best that have been active solutions before. So experiment with not
					adding these edge solutions as active solutions, but need to add already visited neighbors.
			*/
			/* Aug 31, 2009, do not add neighbors that have been visited as active solutions. This step may not
			   be necessary for convergence and slows down AHA too much
			   */
			/*for (int dim = 0; dim < dimension; dim++)
			{
				Solution* x0nbr1 = new Solution(dimension, visited[starIndex]->getDiscreteSolution());
				x0nbr1->changeOneDimension(*(visited[starIndex]->getDiscreteSolution()+dim)+1, dim);
				// Check the feasibility of this neighbor
				if (CheckFeasibility(*x0nbr1, A, b))
				{	// If a feasible neighbor has been visited, it must have been in the estimation set 
					// before. So add it to the active solution set.
					long index = findSolution(visited, x0nbr1);
					if (index < 0)
					{
						delete x0nbr1;
					}
					else
					{
						delete x0nbr1;
						activeSolutions.push_back(index);
					}
				}
				Solution* x0nbr2 = new Solution(dimension, visited[starIndex]->getDiscreteSolution());
				x0nbr2->changeOneDimension(*(visited[starIndex]->getDiscreteSolution()+dim)-1, dim);
				if (CheckFeasibility(*x0nbr2, A, b))
				{
					long index = findSolution(visited, x0nbr2);
					if (index < 0)
					{	
						delete x0nbr2;
					}
					else
					{
						delete x0nbr2;
						activeSolutions.push_back(index);
					}
				}
			}*/
			activeSolutions.push_back((*cmit).second);
			// Debug use
			/*if (i == 0) 
			{
				cout << visited[(*cmit).second];
				cout << endl;
			}*/
			double* Atemp = new double[dimension];
			fill(Atemp, &Atemp[dimension], 0);
			Atemp[i] = -1;
			// The key is the coordinate position. So it is the rhs of the constraint
			//double btemp = -((*cmit).first+key)/2.0;
			double btemp = -(*cmit).first;
			Ak.push_back(Atemp);
			bk.push_back(btemp);
		}
	}
}