#include "masterheader.h"
using namespace std;


// Ak, bk, activeSolutions are the returned values by this function.
/* When this function is called, activeConstraints has the indices of those solutions in "visited" 
   that define MPR in the previous iteration. When the function call finishes, activeSolutions
   contains active constraints for this iteration, which will be used when RMD is called to sample
   new solutions, and it doesn't include current sample best inside it(activeSolutions) */ 
int ConstructMPR(vector<Solution*> visited, long starIndex, long oldStarIndex, vector<double*>& Ak, vector<double>& bk, \
			const vector<double*>& A, const vector<double>& b, vector<long>& activeSolutions, const vector<long>& bestSolutions)
{
	// Once RMD is done, Ak, bk are no longer useful. Release memory allocated to arrays pointed by pointers in Ak
	// But be very, very careful, Ak also contains the pointers in A. If they are deleted, original problem
	// constraints are lost and memory problems occur!!!!!
	for (size_t i = A.size(); i < Ak.size(); i++)
	{
		delete []Ak[i];
	}
	// First the original problem formulation constraints
	// copy the contents of A into Ak
	vectorCopy(A, Ak);
	vectorCopy(b, bk);
	vector<long> oldActiveSolutions;
	/* I am NOT using all other visited solutions when writing down the LP formulations for determining 
	 * new active solutions because doing that may make the size of LP too large as the algorithm goes on. 
	 * Instead, the center is indexed by starIndex, and possible candidates for activeSolutions are
	 * 1. Those "old" active solutions which are held in activeSolutions at the entry of this function. 
	 *    If the current sample best happens to be one of those old active solutions, should remove it from
	 *    activeSolutions.
	 * 2. The "old" best solution if it is no longer the current sample best.
	 * 3. 
	 */
	// First keep all old active solutions except possibly the one now being current sample best
	for (long i = 0; i < (long) activeSolutions.size(); i++)
	{
		if (activeSolutions[i] != starIndex) oldActiveSolutions.push_back(activeSolutions[i]);
	}
	// Second, add oldStarindex if it is not the same as starIndex not already inside oldActiveSolutions
	bool found = false;
	if (oldStarIndex != starIndex)
	{
		for (long i = 0; i < (long) oldActiveSolutions.size(); i++)
		{
			if (oldActiveSolutions[i] == oldStarIndex)
			{
				found = true;
				break;
			}
		}	
		if (!found)	oldActiveSolutions.push_back(oldStarIndex);
	}
	// Third, add best solutions so that in case of revisedcompass restarting, compass won't end up with
	// too few active solutions and thus make the MPR unnecessarily large.
	for (long i = 0; i < (long) bestSolutions.size(); i++)
	{
		found = false;
		if (bestSolutions[i] != starIndex)
		{
			for(long j = 0; j < (long) oldActiveSolutions.size(); i++)
			{
				if (oldActiveSolutions[j] == bestSolutions[i])
				{
					found = true;
					break;
				}
			}
			if (!found)	oldActiveSolutions.push_back(bestSolutions[i]);
		}
	}
	activeSolutions.clear();
	lprec *lp;
	int dimension = visited[starIndex]->getDimension();
	int* colno = new int[dimension];
	for (int i = 0; i < dimension; i++)	colno[i] = i + 1;
	/* MPR constraints are of the form:  
	   (xStar - xi)' * (x - (xStar+xi)/2) >= 0
	   So 
	   Ak = xStar - xi
	   bk = (xStar - xi)' * ((xStar+xi)/2)
	*/
	// Implement Michael Trick's idea about removing redundant constraints in RevisedCompass when MPR is constructed. See research notes p3. 
	/* For each active visited solution, check if it is redundant in terms of defining MPR.
	 * Note that when solving the LP, only those solutions that are still active plus new solutions enter
	 * the LP constraints formulation because otherwise the size of LP may be unnecessarily too large.
	 */
	/* Constraints for the MPR come from two sets of solutions: active solutions and new solutions.
	 * Form constraints from those two sets and determine which ones might be redundant. First consider
	 * active solutions to see if they are redundant.
	 */
	for (long i = 0; i < (long) oldActiveSolutions.size(); i++)
	{
		/* This is an active solution forming the MPR constraints. If this visited solution 
		 * is redundant, its objective will be positive. Form a LP for this active solution indexed by
		 * visited[oldActiveSolutions[i]]
		 */
		lp = make_lp(0, dimension);
		set_verbose(lp, CRITICAL);
		set_add_rowmode(lp, TRUE);
		double* Atemp = new double[dimension];
		for (long j = 0; j < (long) oldActiveSolutions.size(); j++)
		{
			if (i != j)
			{
				double btemp = 0;
				for (int n = 0; n < dimension; n++)
				{
					Atemp[n] = (visited[starIndex]->getDiscreteSolution())[n] - \
								(visited[oldActiveSolutions[j]]->getDiscreteSolution())[n];
					btemp += Atemp[n] * ((visited[starIndex]->getDiscreteSolution())[n] + \
							(visited[oldActiveSolutions[j]]->getDiscreteSolution())[n]) / 2;
				}
				add_constraintex(lp, dimension, Atemp, colno, GE, btemp);
			}
		}
		// Now I finish adding constraints from other active solutions. Don't forget original problem formulation constraints!
		for (unsigned int k = 0; k < A.size(); k++)
		{
			for (int n = 0; n < dimension; n++)
			{
				Atemp[n] = *(A[k]+n);
			}
			add_constraintex(lp, dimension, Atemp, colno, GE, b[k]);
		}
		set_add_rowmode(lp, FALSE);
		for (int k = 0; k < dimension; k++)
		{
			Atemp[k] = 0;
		}			
		double objective = 0;
		for (int n = 0; n < dimension; n++)
		{
			Atemp[n] = (visited[starIndex]->getDiscreteSolution())[n] - \
						(visited[oldActiveSolutions[i]]->getDiscreteSolution())[n];
			objective += Atemp[n] * ((visited[starIndex]->getDiscreteSolution())[n] + \
						(visited[oldActiveSolutions[i]]->getDiscreteSolution())[n]) / 2;
		}
		set_obj_fnex(lp, dimension, Atemp, colno);
		for (int n = 0; n < dimension; n++)
		{
			set_bounds(lp, n + 1, -get_infinite(lp), get_infinite(lp));        
		}
		set_minim(lp);
	#ifdef MYDEBUG	
		write_LP(lp, stdout);
		set_verbose(lp, IMPORTANT);
	#endif
		int ret = solve(lp);
		switch (ret)
		{
			case OPTIMAL: 
				// If the objective value is less than btemp, then the constraint is not redundant
				if (get_objective(lp) < objective)
				{
					Ak.push_back(Atemp);
					bk.push_back(objective);
					activeSolutions.push_back(oldActiveSolutions[i]);
				}
				break;
			case SUBOPTIMAL: 
				cout << " !!!!!!!!!!!!!!!!!!!  Lpsolve cannot solve to optimum" << endl;
				delete []Atemp;
				break;		
			case INFEASIBLE: 
				cout << " !!!!!!!!!!!!!!!!!!!  Lpsolve encounters an infeasible problem" << endl;
				Ak.push_back(Atemp);
				bk.push_back(objective);
				activeSolutions.push_back(oldActiveSolutions[i]);
				//PrintActiveSolutions(visited, oldActiveSolutions);
				//cout << "Xstar =  " << *(visited[starIndex]);
				write_lp(lp, "lp.txt");
				ret = solve(lp);
				/*if (ret == OPTIMAL && !constraintPruning || get_objective(lp) < objective)
				{
					Ak.push_back(Atemp);
					bk.push_back(objective);
					activeSolutions.push_back(oldActiveSolutions[i]);
				}*/
				cout << get_objective(lp) << endl;
				//delete []Atemp;
				break;   
			case UNBOUNDED: 		
    			cout << " Unbounded !!!!!!!!,  keep this solution as an active constraint" << endl;
    			//cout << visited.size() << " visited" << endl;
				//cout << oldActiveSolutions.size() << " old active constraints" << endl;
				Ak.push_back(Atemp);
				bk.push_back(objective);
				activeSolutions.push_back(oldActiveSolutions[i]);
				//PrintActiveSolutions(visited, oldActiveSolutions);
				//cout << "Xstar =  " << *(visited[starIndex]);
				write_lp(lp, "lp.txt");
				ret = solve(lp);
				/*if (ret == OPTIMAL && !constraintPruning || get_objective(lp) < objective)
				{
					Ak.push_back(Atemp);
					bk.push_back(objective);
					activeSolutions.push_back(oldActiveSolutions[i]);
				}*/
				cout << get_objective(lp) << endl;
				//cout << "XStar" << endl;
				//delete []Atemp;
				break;
			default:  
				cout << "LpSolve encounters problems, code " << ret << endl; 
				Ak.push_back(Atemp);
				bk.push_back(objective);
				activeSolutions.push_back(oldActiveSolutions[i]);
				//PrintActiveSolutions(visited, oldActiveSolutions);
				//cout << "Xstar =  " << *(visited[starIndex]);
				write_lp(lp, "lp.txt");
				ret = solve(lp);
				/*if (ret == OPTIMAL && !constraintPruning || get_objective(lp) < objective)
				{
					Ak.push_back(Atemp);
					bk.push_back(objective);
					activeSolutions.push_back(oldActiveSolutions[i]);
				}*/
				cout << get_objective(lp) << endl;
				//delete []Atemp;
    	}
    	delete_lp(lp);
	}
	//cout << activeSolutions.size() << " new active constraints" << endl;
	//PrintActiveSolutions(visited, activeSolutions);
	delete []colno;
	return 1;
}

// Ak, bk, activeSolutions are the returned values by this function.
/* When this function is called, activeConstraints has the indices of those solutions in "visited" 
   that define MPR in the previous iteration. When the function call finishes, activeSolutions
   contains active constraints for this iteration, which will be used when RMD is called to sample
   new solutions, and it doesn't include current sample best inside it(activeSolutions) */ 
void MPRv2(bool stocsim, vector<Solution*> visited, long starIndex, long oldStarIndex, vector<double*>& Ak, vector<double>& bk, \
			const vector<double*>& A, const vector<double>& b, vector<long>& activeSolutions, \
			const vector<long>& bestSolutions, double alphac, bool constraintPruning) 
{
	// Once RMD is done, Ak, bk are no longer useful. Release memory allocated to arrays pointed by pointers in Ak
	// But be very, very careful, Ak also contains the pointers in A. If they are deleted, original problem
	// constraints are lost and memory problems occur!!!!!
	//static long iterationcount = 0;
	//clock_t beginning = clock();
	//static ofstream lpsizefile("../lpsize.dat");
	for (size_t i = A.size(); i < Ak.size(); i++)
	{
		delete []Ak[i];
	}
	// First the original problem formulation constraints
	// copy the contents of A into Ak
	vectorCopy(A, Ak);
	vectorCopy(b, bk);
	vector<long> oldActiveSolutions;
	/* I am NOT using all other visited solutions when writing down the LP formulations for determining 
	 * new active solutions because doing that may make the size of LP too large as the algorithm goes on. 
	 * Instead, the center is indexed by starIndex, and possible candidates for activeSolutions are
	 * 1. Those "old" active solutions which are held in activeSolutions at the entry of this function. 
	 *    If the current sample best happens to be one of those old active solutions, should remove it from
	 *    activeSolutions.
	 * 2. The "old" best solution if it is no longer the current sample best.
	 * 3. 
	 */
	// First keep all old active solutions except possibly the one now being current sample best
	for (long i = 0; i < (long) activeSolutions.size(); i++)
	{
		if (activeSolutions[i] != starIndex) oldActiveSolutions.push_back(activeSolutions[i]);
	}
	// Second, add oldStarindex if it is not the same as starIndex and not already inside oldActiveSolutions
	// bool found = false;
	if (oldStarIndex != starIndex)
	{
		if (find(oldActiveSolutions.begin(), oldActiveSolutions.end(), oldStarIndex) == oldActiveSolutions.end())
		{
			oldActiveSolutions.push_back(oldStarIndex);
		}
	}
	// Third, add best solutions so that in case of revisedcompass restarting, compass won't end up with
	// too few active solutions and thus make the MPR unnecessarily large.
	for (size_t i = 0; i < bestSolutions.size(); i++)
	{
		//found = false;
		if (bestSolutions[i] != starIndex)
		{
			if (find(oldActiveSolutions.begin(), oldActiveSolutions.end(), bestSolutions[i]) == oldActiveSolutions.end())
			{
				oldActiveSolutions.push_back(bestSolutions[i]);
			}
		}
	}
	activeSolutions.clear();
	lprec *lp;
	int dimension = visited[starIndex]->getDimension();
	int* colno = new int[dimension];
	for (int i = 0; i < dimension; i++)
	{
		colno[i] = i + 1;
	}
	bool* redundant = new bool[oldActiveSolutions.size()];
	for(size_t i = 0; i < oldActiveSolutions.size(); i++)
	{
		redundant[i] = false;
	}
	/* MPR constraints are of the form:  
	   (xStar - xi)' * (x - (beta*xStar + (1-beta)*xi)) >= 0,  0<=beta<=0.5
	   So 
	   Ak = xStar - xi
	   bk = (xStar - xi)' * (beta*xStar + (1-beta)*xi)) = Ak * (beta*xStar + (1-beta)*xi))
	*/
	/* For each active visited solution, check if it is redundant in terms of defining MPR.
	 * Note that when solving the LP, only those solutions that are still active plus new solutions enter
	 * the LP constraints formulation because otherwise the size of LP may be unnecessarily too large.
	 */

	/* Constraints for the MPR come from two sets of solutions: active solutions and new solutions.
	 * Form constraints from those two sets and determine which ones might be redundant. First consider
	 * active solutions to see if they are redundant.
	 */
	/* To incorporate the modified COMPASS constraint placement scheme, need to compute beta for every 
	 * candidate of active solution, which is for every element in oldActiveSolutions. 
	 */
	vector<double> beta(oldActiveSolutions.size());
	double S2star = visited[starIndex]->getSampleVariance();
	int nstar = visited[starIndex]->getNumOfObservations();
	for (long i = 0; i < (long) oldActiveSolutions.size(); i++)
	{
		int ni = visited[oldActiveSolutions[i]]->getNumOfObservations();
		// Step 0
		require(0 < ni && 0 < nstar, "Unexpected error: A solution has 0 replication in MPR");
		if (!stocsim || visited[oldActiveSolutions[i]]->getDistanceTo(visited[starIndex]) < sqrt(double(dimension)))
		{	// If this is a deterministic simulator, simply set beta = 0.5
			beta[i] = 0.5;
			continue;
		}
		// Step 1, calculate Sd and nu
		double S2i = visited[oldActiveSolutions[i]]->getSampleVariance();
		double Sd = sqrt(S2star/nstar + S2i/ni);
		// Need to be careful. nStar and ni cannot be 1 here. But if some solution only has one observation taken. 
		// cannot do a divison by 0
		double denomstar, denomi;
		denomstar = (nstar == 1) ? 0 : S2star/nstar * S2star/nstar / (nstar - 1);
		denomi = (ni == 1) ? 0 : S2i/ni * S2i/ni / (ni - 1);
		// If both points only have one observation, treat it as deterministic and don't do test
		if (denomstar < numeric_limits<double>::epsilon() && denomi < numeric_limits<double>::epsilon())
		{
			beta[i] = 0.5;
		}
		else
		{	//else at least one of them has more than 1 observation, 
			int	nu = int(ceil((S2star/nstar + S2i/ni) * (S2star/nstar + S2i/ni) / (denomstar + denomi)));
			double t = tinv(1 - alphac, nu);
			if (dabs(visited[starIndex]->getSampleMean() - visited[oldActiveSolutions[i]]->getSampleMean()) < t * Sd)
			{   // Set beta to a minimum value rather than 0 so that a singleton can be built.
				beta[i] = pow(0.05, dimension);
			}
			else
			{
				double tempbeta = 1- t * Sd / dabs(visited[starIndex]->getSampleMean() - visited[oldActiveSolutions[i]]->getSampleMean());
				beta[i] = (tempbeta < 0.5) ? tempbeta : 0.5;
			}
		}
	}	
	vector<double*> ConstraintTable(oldActiveSolutions.size() + A.size());
	for (size_t i = 0; i < ConstraintTable.size(); i++) ConstraintTable[i] = new double[dimension];
	// First add constraints from old active solutions
	for (long i = 0; i < (long) oldActiveSolutions.size(); i++)
	{
		for (int n = 0; n < dimension; n++)
		{
			*(ConstraintTable[i]+n) = (visited[starIndex]->getDiscreteSolution())[n] - (visited[oldActiveSolutions[i]]->getDiscreteSolution())[n];
		}
	}
	// Then add constraints from original problem formulation
	for (unsigned int k = 0; k < A.size(); k++)
	{
		for (int n = 0; n < dimension; n++)
		{
			*(ConstraintTable[oldActiveSolutions.size()+k] + n) = *(A[k]+n);
		}
	}
	for (long i = 0; i < (long) oldActiveSolutions.size(); i++)
	{
		/* This is an active solution forming the MPR constraints. If this visited solution 
		 * is redundant, its objective will be positive. Form a LP for this active solution indexed by
		 * visited[oldActiveSolutions[i]]
		 */
		double* Atemp = new double[dimension];
		double objective = 0;
		if (!constraintPruning)
		{ // If no constraint pruning is needed, simply keep this constraint
			for (int n = 0; n < dimension; n++)
			{
				Atemp[n] = (visited[starIndex]->getDiscreteSolution())[n] - \
							(visited[oldActiveSolutions[i]]->getDiscreteSolution())[n];
				objective += Atemp[n] * (beta[i] * (*(visited[starIndex]->getDiscreteSolution() + n)) + \
							(1 - beta[i]) * (*(visited[oldActiveSolutions[i]]->getDiscreteSolution() + n)));
			}
			Ak.push_back(Atemp);
			bk.push_back(objective);
			activeSolutions.push_back(oldActiveSolutions[i]);
		}
		else
		{// Construct LP to determine if this constraint should be pruned
			lp = make_lp(0, dimension);
			set_verbose(lp, CRITICAL);
			set_add_rowmode(lp, TRUE);
			for (long j = 0; j < (long) oldActiveSolutions.size(); j++)
			{
				if (!redundant[j] && i != j)
				{
					double btemp = 0;
					for (int n = 0; n < dimension; n++)
					{
						btemp += *(ConstraintTable[j] + n) * (beta[j] * (*(visited[starIndex]->getDiscreteSolution() + n)) + (1 - beta[j]) * (*(visited[oldActiveSolutions[j]]->getDiscreteSolution() + n)));
					}
					add_constraintex(lp, dimension, ConstraintTable[j], colno, GE, btemp);
				}
			}
			// Now I finish adding constraints from other active solutions. Don't forget original problem formulation constraints!
			for (unsigned int k = 0; k < A.size(); k++)
			{
				add_constraintex(lp, dimension, ConstraintTable[oldActiveSolutions.size()+k], colno, GE, b[k]);
			}
			set_add_rowmode(lp, FALSE);
			for (int n = 0; n < dimension; n++)
			{
				Atemp[n] = (visited[starIndex]->getDiscreteSolution())[n] - \
						(visited[oldActiveSolutions[i]]->getDiscreteSolution())[n];
				objective += Atemp[n] * (beta[i] * (*(visited[starIndex]->getDiscreteSolution() + n)) + \
						(1 - beta[i]) * (*(visited[oldActiveSolutions[i]]->getDiscreteSolution() + n)));
			}
			set_obj_fnex(lp, dimension, Atemp, colno);
			for (int n = 0; n < dimension; n++)
			{
				set_bounds(lp, n + 1, -get_infinite(lp), get_infinite(lp));        
			}
			set_minim(lp);
			// write_LP(lp, stdout);
			// set_verbose(lp, IMPORTANT);
			int ret = solve(lp);
			switch (ret)
			{
				case OPTIMAL: 
					// If the objective value is less than btemp, then the constraint is not redundant
					// or if it doesn't require pruning, simply keep it.
					if (!constraintPruning || get_objective(lp) < objective)
					{
						Ak.push_back(Atemp);
						bk.push_back(objective);
						activeSolutions.push_back(oldActiveSolutions[i]);
					}
					else
					{// This constraint is redundant, remove it from old active solutions.
						redundant[i] = true;
					}
					break;
				case SUBOPTIMAL: 
					cout << " !!!!!!!!!!!!!!!!!!!  Lpsolve cannot solve to optimum" << endl;
					Ak.push_back(Atemp);
					bk.push_back(objective);
					activeSolutions.push_back(oldActiveSolutions[i]);
					//PrintActiveSolutions(visited, oldActiveSolutions);
					//cout << "Xstar =  " << *(visited[starIndex]);
					write_lp(lp, "lp.txt");
					ret = solve(lp);
					/*if (ret == OPTIMAL && !constraintPruning || get_objective(lp) < objective)
					{
						Ak.push_back(Atemp);
						bk.push_back(objective);
						activeSolutions.push_back(oldActiveSolutions[i]);
					}*/
					cout << get_objective(lp) << endl;
					//delete []Atemp;
    				break;		
				case INFEASIBLE: 
					cout << " !!!!!!!!!!!!!!!!!!!  Lpsolve encounters an infeasible problem" << endl;
					Ak.push_back(Atemp);
					bk.push_back(objective);
					activeSolutions.push_back(oldActiveSolutions[i]);
					//PrintActiveSolutions(visited, oldActiveSolutions);
					//cout << "Xstar =  " << *(visited[starIndex]);
					write_lp(lp, "lp.txt");
					ret = solve(lp);
					/*if (ret == OPTIMAL && !constraintPruning || get_objective(lp) < objective)
					{
						Ak.push_back(Atemp);
						bk.push_back(objective);
						activeSolutions.push_back(oldActiveSolutions[i]);
					}*/
					cout << get_objective(lp) << endl;
					//delete []Atemp;
    				break;
				case UNBOUNDED: 		
    				cout << " Unbounded !!!!!!!!!!!!!!!" << endl;
    				//cout << visited.size() << " visited" << endl;
					//cout << oldActiveSolutions.size() << " old active constraints" << endl;
					Ak.push_back(Atemp);
					bk.push_back(objective);
					activeSolutions.push_back(oldActiveSolutions[i]);
					//PrintActiveSolutions(visited, oldActiveSolutions);
					//cout << "Xstar =  " << *(visited[starIndex]);
					write_lp(lp, "lp.txt");
					ret = solve(lp);
					/*if (ret == OPTIMAL && !constraintPruning || get_objective(lp) < objective)
					{
						Ak.push_back(Atemp);
						bk.push_back(objective);
						activeSolutions.push_back(oldActiveSolutions[i]);
					}*/
					cout << get_objective(lp) << endl;
					//cout << "XStar" << endl;
					//delete []Atemp;
					break;
			default:
					cout << "LpSolve encounters problems!!!!" << endl; 
					Ak.push_back(Atemp);
					bk.push_back(objective);
					activeSolutions.push_back(oldActiveSolutions[i]);
					//PrintActiveSolutions(visited, oldActiveSolutions);
					//cout << "Xstar =  " << *(visited[starIndex]);
					write_lp(lp, "lp.txt");
					ret = solve(lp);
					/*if (ret == OPTIMAL && !constraintPruning || get_objective(lp) < objective)
					{
						Ak.push_back(Atemp);
						bk.push_back(objective);
						activeSolutions.push_back(oldActiveSolutions[i]);
					}*/
					cout << get_objective(lp) << endl;
					//delete []Atemp;
    		}
    		delete_lp(lp);
		}
	}
	delete []colno;
	for (size_t i = 0; i < ConstraintTable.size(); i++) delete ConstraintTable[i];
}

//12/27/11 adds a function to compute the volume of an MPA specified in Ax>=b constraints, right now only hyperbox supported
double calcVolMPA(const std::vector<double*>& A, const std::vector<double>& b, unsigned int dimension)
{
	vector<set<long>> mpaub, mpalb;
	for (size_t i = 0; i < dimension; i++)
	{
		set<long> dimi;
		mpaub.push_back(dimi);
		mpalb.push_back(dimi);
	}
	for (size_t i = 0; i < A.size(); i++)
	{
		double asum = 0;
		for (size_t d = 0; d < dimension; d++)
		{
			if (myabs(A[i][d]) > 0 && myabs(A[i][d]) != 1)
			{
				cout << "Wrong constraints\n";
				exit(8);
			}
			else if (myabs(A[i][d]) > 0)
			{
				if (A[i][d] > 0) 
					mpalb[d].insert(b[i]);
				else
					mpaub[d].insert(-b[i]);
				asum += myabs(A[i][d]);
			}
		}
		if (asum > 1) 
		{
			cout << "Not a hyperbox constraint \n";
			exit(8);
		}
	}
	double mpavol = 1;
	for (size_t d = 0; d < dimension; d++)
	{
		double dimlength = *(mpaub[d].begin()) - (*mpalb[d].rbegin());
		mpavol *= *(mpaub[d].begin()) - (*mpalb[d].rbegin());
	}
	return mpavol;
}

//12/28/11 computes the alignment prob of a blind picking rule with N total solutions, top G solutions, select S and with at least k among top G
//Blind picking top G from a total of N with a sample size of S, is 1 - (N-g choose s) / (N choose s) (c.f. OO, p27, (2.13)
// Stirling's approx gives n! ~ sqrt(2*pi*n) (n/e)^n
double AlgnProb(double N, double G, double S, double k, int method)
{ // The current implementation only supports k = 1
  // c.f. OO book p27
	double prob = 1.0;
	if (method == 0)
	{ // Exact computation when N is small
		for (double i = 0; i < S; i++)
		{	// To compute the blink picking alignment probability, need to cancel terms in the formula to avoid extremly large numbers
			prob *= (N-G-i) / (N-i);
		}
	}
	else
	{ // Use Stirling's formula to approximate binomial coefficient n choose k, c.f. Wikipedia, sqrt(n) mn choose n = pow(m, m(n-1)+1)/pow(m-1, (m-1)(n-1))
		double n1 = (N-G)/(N-S-G), n2 = (N-G)/N, n3 = (N-S)/N;
		double sqrtn = sqrt(n1*n3);
		double temp1 = (N-S) * log(n3);
		double temp2 = S * log(n2);
		double temp3 = (N-G-S) * log(n1);
		prob = sqrtn * exp(temp1 + temp2 + temp3);
	}
	return 1-prob;
}

//12/28/11 determines sample size for a blind picking OO 
double SampleSizeOO(double N, double G, double k)
{ // Now only supports k = 1
	require(k==1, "OO only supports k == 1");
	double Sl = 1, Su = N, S; // hard code predication max number to 1M
	double ap;
	//double tp1 = AlgnProb(900, 5, 5, 1, 1);
	//double tp0 = AlgnProb(900, 5, 5, 1, 0);
	while (Sl + 1 < Su) 
	{
		S = floor((Sl + Su)/2);
		// The last arg of AlgnProb indicates whether to use exact factorial (0) or stirling's formula (1)
		if (N > 100000)
			ap = AlgnProb(N, G, S, 1, 1);
		else
			//ap = AlgnProb(N, G, S, 1, 0);
			ap = AlgnProb(N, G, S, 1, 0);
		if (ap < 0.9) // 1/2/12 hardcode AP guarantee to 0.9 
			Sl = S;
		else
			Su = S;
	}
	return S;
}