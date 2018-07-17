#include "masterheader.h"
#include "NichingGAWTR.h"
#include "NichePopulationWTR.h"

using namespace std;
/* simulator: simulator used for this problem
 * N: budget, # of reps
 * P: percentage of budget spent on GA
 * G: numer of generations without improvement for NGA to quit
 * globaldelta: indifference parameter used in grouping procedure of the NGA
 * globalalpha: confidence parameter used in grouping procedure of the NGA 
 * alphadt: dominant niche test alpha parameter of the NGA
 * alphac: CI-level for generalized COMPASS constraint placement procedure
 * deltae, alphae: backtracking test parameters 
 * deltam, alpham:  for optimality stopping test
 * deltacu:  indif parameter for clean-up
 * alphacu: CI level for clean up
 * initialNumReps: initial number of reps assigned to a new solution
 * numCandidates: size of S_k, i.e. # of samplings within MPR and NGA (repetion ok)
 * dimension: dimension of solution space
 * visited: vector of visited solutions. It has the starting point passed from the program calling ISCControl.
 * A, b: problem constraints
 * repFile: write sample path of this rep to this file
 * records: accumulated statistics for averaged sample path after all reps
 * backTrackingTest: true if backTracking test is used 
 * ocba: true if ocba replications application rule is used
 * dntest: true if dominant niche rule is used
 * docleanup : true if cleanup is desired
 */
const double EPSILON = 1e-6;

/*long ISCControl(Simulator* simulator, long& N, double P, int G, double globaldelta, double globalalpha, double alphadt, \
	double alphac, double deltae, double alphae, double deltam, double alpham, double deltacu, double alphacu, \
	int initialNumReps, int numCandidates, int dimension, vector<Solution*>& visited, const vector<double*>& A, \
	const vector<double>& b, ofstream& repFile, ofstream& resultFile, vector<Record*>& records, bool backTrackingTest, bool ocba, \
	bool dntest, bool docleanup, bool stocsim, bool elitism, int pruning_freq, int mpamode,int rmdmode, int metamodel)*/
long ISCControl(Simulator* simulator, long& N, int dimension, vector<Solution*>& visited, const vector<double*>& A, \
	const vector<double>& b, ofstream& repFile, ofstream& resultFile, vector<Record*>& records, ISCParameters* iscparams)
{
	vector<long> M;
	long budgetforGA = iscparams->get_gabudget() < EPSILON ? 0 : long(ceil(iscparams->get_gabudget() * iscparams->get_budget())); 
	long remainingBudget = iscparams->get_budget(); 
	double bestObjectiveValue, trueBestObjectiveValue;
	size_t bestSolutionIndex;
	bool optimum;
	long min;
	if (budgetforGA > 1-EPSILON)
	{
		/*NichingGAWTR<Solution> ga(budgetforGA, 50, iscparams->get_gagen(), alphadt, 50, 1.0, 0.328, dimension, 2, 1.5, 10, globaldelta, globalalpha, \
			3, initialNumReps, visited, A, b, repFile, records, dntest, simulator, stocsim, elitism);*/
		NichingGAWTR<Solution> ga(budgetforGA, 50, 20, 1.0, 0.328, dimension, 2, 1.5, 10, 3, visited, A, b, repFile, records, simulator, iscparams);
		long usedreps = ga.solve();
		remainingBudget = iscparams->get_budget() - usedreps;
		cout << "GA evolved " << ga.getGeneration() << " generations and used " << usedreps << " reps" << endl;
		resultFile << "GA evolved " << ga.getGeneration() << " generations and used " << usedreps << " reps" << endl;
		// This gives us the index of best solution found by GA in visited
		bestSolutionIndex = ga.getBestGenome();
		trueBestObjectiveValue = simulator->simulation(visited[bestSolutionIndex]);  
		bestObjectiveValue = visited[bestSolutionIndex]->getSampleMean();
		vector<NichePopulationWTR<Solution>::Niche*> niches = ga.getNiches();
		vector<uint> nonPeakNiches = ga.getNonPeakNiches();
		// If there is a dominant niche, starts from it, otherwise, just use the order in niches, 
		// which is from the best niche center to the worst as a result of the greedy niche identification procedure.
		vector<size_t> order;
		for (size_t i = 0; i < niches.size(); i++) order.push_back(i);
		if (0 <= ga.getDominantNiche())  
		{
			order[0] = ga.getDominantNiche();
			for (int i = 1; i <= ga.getDominantNiche(); i++)
				order[i] = i - 1;
		}
		for (size_t i = 0; i < order.size(); i++)
		{
			cout << "Niche size "  << niches[order[i]]->getNicheCount() << ",  ";
			cout << "center: "  << *(visited[niches[order[i]]->getCenter()]);
			//const vector<uint> residents = niches[order[i]]->getResidents();
			//std::cout << "Residents " << std::endl;
			//for (size_t j = 0; j < residents.size(); j++)	cout << *(visited[residents[j]]);
		}	
		/*cout << "Nonpeak niceh " << endl;
		for (size_t i = 0; i < nonPeakNiches.size(); i++)
		{
			cout << *(visited[nonPeakNiches[i]]);
		}
		cout << endl;*/
		// For each niche center, as long as budget allows, run COMPASS starting from it.
		for (size_t i = 0; 0 < remainingBudget && i < order.size(); i++)
		{
			// Those two variables lastObjectiveVAlue, lastRepCount, are used to record results
			double lastObjectiveValue = simulator->simulation(visited[niches[order[i]]->getCenter()]);
			// iscparams->get_budget() is the overall budget, so iscparams->get_budget() - remainingBudget is the used # of reps
			long lastRepCount = iscparams->get_budget() - remainingBudget;
			vector<long> activeSolutions;
			const vector<uint> residents = niches[order[i]]->getResidents();
			cout << endl << "Now starting from niche center " << *visited[niches[order[i]]->getCenter()];
			for (size_t j = 0; j < residents.size(); j++)
			{
				if (find(activeSolutions.begin(), activeSolutions.end(), long(residents[j])) == \
						activeSolutions.end() && residents[j] != niches[i]->getCenter()) 
					activeSolutions.push_back(residents[j]);
			}
			// If there are fewer then 2*dimension many other residents in this niche, then 
			// look at the entire history and use closet ones as active solution candidates to make the size
			// 2*dimension.
			if (int(activeSolutions.size()) < 2 * dimension)
			{
				double* distanceToCenter = new double[visited.size()];
				for (size_t disti = 0; disti < visited.size(); disti++) 
					distanceToCenter[disti] = visited[disti]->getDistanceTo(visited[niches[order[i]]->getCenter()]);
				size_t* rank = new size_t[visited.size()];
				for (size_t ranki = 0; ranki < visited.size(); ranki++) rank[ranki] = ranki;
				for (size_t dest = visited.size() - 1; 0 < dest; dest--)
				{
					for (size_t src = 0; src < dest; src++)
					{
						if (distanceToCenter[rank[src+1]] < distanceToCenter[rank[src]])
						{
							size_t temp = rank[src+1];
							rank[src+1] = rank[src];
							rank[src] = temp;
						}
					}
				}
				require(rank[0] == niches[order[i]]->getCenter(), "Error in distance ranking in ISCControl.cpp");			 
				for (size_t ai = 1; int(activeSolutions.size()) < 2*dimension && ai < visited.size(); ai++)
				{
					if (find(activeSolutions.begin(), activeSolutions.end(), long(rank[ai])) == \
						activeSolutions.end() && rank[ai] != niches[i]->getCenter()) 
						activeSolutions.push_back(long(rank[ai]));
				}
				delete []distanceToCenter;
				delete []rank;
			}
			/*long xstar = ISCv3(simulator, stocsim, deltam, alpham, visited, A, b, activeSolutions, remainingBudget, niches[order[i]]->getCenter(), \
				repFile, records, revisedcompass, initialNumReps, numCandidates, alphac, alphae, deltae, pruning_freq, mpamode, rmdmode,\
				lastRepCount, lastObjectiveValue, bestObjectiveValue, trueBestObjectiveValue, \
				backTrackingTest, ocba, i == 0, optimum, metamodel);*/
			long xstar = ISCv3(simulator, visited, A, b, iscparams, activeSolutions, remainingBudget, niches[order[i]]->getCenter(), repFile, records, \
				lastRepCount, lastObjectiveValue, bestObjectiveValue, trueBestObjectiveValue, i == 0, optimum);
			// xstar is the index of the best found by this COMPASS run in visited. 	
			if (find(M.begin(), M.end(), xstar) == M.end()) M.push_back(xstar);
			cout << "Local search returns " << *visited[xstar] << "after " << iscparams->get_budget() - remainingBudget - lastRepCount << "  replications" << endl;
			if (optimum) 
				resultFile << "local search returns a local minimum " << *visited[xstar];
			if (i > 0)
			{ // If this is not the first run of Compass, use best value before to record the result.
		  // Only up to the replications this COMPASS finishes, can we use a possibly new best
		  // objective value.
				for (long repcount = lastRepCount + 1; repcount < iscparams->get_budget() - remainingBudget; repcount++)
				{
					RecordResult(records, repcount, trueBestObjectiveValue);
				}
			}
			// Update best objective value found so far
			double tempObjectiveValue = visited[xstar]->getSampleMean();
			if (i == 0 || tempObjectiveValue < bestObjectiveValue)
			{ // Only results returned by COMPASS are reliable. So if this is the first run, bestObjectiveValue
		  // is returned by GA, not reliable. And thus should not compare it with COMPASS result.			
				bestObjectiveValue = tempObjectiveValue;
				bestSolutionIndex = xstar;
				trueBestObjectiveValue = simulator->gettruevalue(visited[bestSolutionIndex]); 
			}
		}
		if (remainingBudget <= 0) 
		{
			cout << endl << "All " << iscparams->get_budget() << " reps used without finding a local min: " << endl;
		}
		else
		{
			cout << endl << "Found local min. Used " << iscparams->get_budget() - remainingBudget << " reps" << endl;
			// N records how many replications have been used, excluding clean-up.
			N = iscparams->get_budget() - remainingBudget;
		}
		// Run clean-up on M
		cout << endl << "# of local min " << M.size() << endl;
		PrintActiveSolutions(visited, M);
		if (iscparams->get_cleanup() && iscparams->get_stocsim())
		{//cleanup(Simulator* simulator, vector<Solution*>& visited, vector<long> M, double delta, double alpha)
			min = cleanup(simulator, visited, M, iscparams->get_cleanupdelta(), iscparams->get_cleanupalpha());
			if (min < 0) 
			{
				cout << "No local min found. No clean-up performed" << endl;
			}
			else
			{
				cout << "With confidence level " << 1-iscparams->get_cleanupalpha() << ", we obtain the best in the cleanup." << endl;
				cout << "The CI = ";
				cout << "[" << visited[min]->getSampleMean() - iscparams->get_cleanupdelta() << ", " << visited[min]->getSampleMean() + iscparams->get_cleanupdelta() << "]" << endl;
			}
			resultFile << "With confidence level " << 1-iscparams->get_cleanupalpha() << ", we obtain the best in the cleanup." << endl;
			resultFile << "The CI = ";
			resultFile << "[" << visited[min]->getSampleMean() - iscparams->get_cleanupdelta() << ", " << visited[min]->getSampleMean() + iscparams->get_cleanupdelta() << "]" << endl;
		}
		else
		{	// No cleanup, simply report the current sample best.
			double best = visited[M[0]]->getSampleMean();
			min = M[0];
			for (size_t i = 1; i < M.size(); i++)
			{
				if (visited[M[i]]->getSampleMean() < best)
				{
					best = visited[M[i]]->getSampleMean();
					min = M[i];
				}
			}	 
			cout << "Best among local min (no clean up):  " << *visited[min] << endl; 
		}
	}
	else
	{ // Just run the local search algorithm when the budget for GA is 0
		cout << "Local search only... " << endl;
		vector<long> activeSolutions;
		visited[0]->recordObservation(simulator->simulation(visited[0]));
		remainingBudget = iscparams->get_budget()-1;
		double lastObjectiveValue = simulator->gettruevalue(visited[0]);
		long lastRepCount = iscparams->get_budget() - remainingBudget;
		trueBestObjectiveValue = simulator->gettruevalue(visited[0]);  
		//Output the first replication objective value
		repFile << "1 , " << trueBestObjectiveValue;
		bestObjectiveValue = visited[0]->getSampleMean();
		/*long xstar = ISCv3(simulator, stocsim, deltam, alpham, visited, A, b, activeSolutions, remainingBudget, 0, \
				repFile, records, revisedcompass, initialNumReps, numCandidates, alphac, alphae, deltae,\
				pruning_freq, mpamode, rmdmode, lastRepCount, lastObjectiveValue, bestObjectiveValue, \
				trueBestObjectiveValue, backTrackingTest, ocba, true, optimum, metamodel);*/
		long xstar = ISCv3(simulator, visited, A, b, iscparams, activeSolutions, remainingBudget, 0, repFile, records, \
				lastRepCount, lastObjectiveValue, bestObjectiveValue, trueBestObjectiveValue, true, optimum);
		// xstar is the index of the best found by this COMPASS run in visited. 	
		if (find(M.begin(), M.end(), xstar) == M.end()) M.push_back(xstar);
		cout << "Local search returns " << *visited[xstar] << "after " << iscparams->get_budget() - remainingBudget - lastRepCount << "  replications" << endl;
		if (optimum) 
			resultFile << "Local search returns a local minimum " << *visited[xstar];
		if (remainingBudget <= 0) 
		{
			cout << endl << "All " << iscparams->get_budget() << " reps used without finding a local min: " << endl;
		}
		else
		{
			cout << endl << "Found local min. Used " << iscparams->get_budget() - remainingBudget << " reps" << endl;
			// N records how many replications have been used, excluding clean-up.
			N = iscparams->get_budget() - remainingBudget;
		}
		// Run clean-up on M
		cout << endl << "# of local min " << M.size() << endl;
		PrintActiveSolutions(visited, M);
		if (iscparams->get_cleanup() && iscparams->get_stocsim())
		{//cleanup(Simulator* simulator, vector<Solution*>& visited, vector<long> M, double delta, double alpha)
			min = cleanup(simulator, visited, M, iscparams->get_cleanupdelta(), iscparams->get_cleanupalpha());
			if (min < 0) 
			{
				cout << "No local min found. No clean-up performed" << endl;
			}
			else
			{
				cout << "With confidence level " << 1-iscparams->get_cleanupalpha() << ", we obtain the best in the cleanup." << endl;
				cout << "The CI = ";
				cout << "[" << visited[min]->getSampleMean() - iscparams->get_cleanupdelta() << ", " << visited[min]->getSampleMean() + iscparams->get_cleanupdelta() << "]" << endl;
			}
			resultFile << "With confidence level " << 1-iscparams->get_cleanupalpha() << ", we obtain the best in the cleanup." << endl;
			resultFile << "The CI = ";
			resultFile << "[" << visited[min]->getSampleMean() - iscparams->get_cleanupdelta() << ", " << visited[min]->getSampleMean() + iscparams->get_cleanupdelta() << "]" << endl;
		}
		else
		{	// No cleanup, simply report the current sample best.
			double best = visited[M[0]]->getSampleMean();
			min = M[0];
			for (size_t i = 1; i < M.size(); i++)
			{
				if (visited[M[i]]->getSampleMean() < best)
				{
					best = visited[M[i]]->getSampleMean();
					min = M[i];
				}
			}	 
			cout << "Best among local min (no clean up):  " << *visited[min] << endl; 
		}
	}
	return min;
}

