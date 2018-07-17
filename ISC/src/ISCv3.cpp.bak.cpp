#include "masterheader.h"
#include <ctime>
#include <map>

using namespace std; 
const int WARMUP = 10;


/* Revised COMPASS algorithm for constrained problem with equal-sample-allocation rule. Linear 
 * constraints take the form Ax>=b. x0 is the starting point, x0=[x01,...x0d]. The objective is a minimization problem.
 * The algorithm takes as input:
 *  
 * epsilon: for stopping test procedure to confirm true local optima
 * alpha: CI level for stopping test
 * 
 * visited: a list of visited solutions, which contains the given starting point if it is the first run. This list is kept throughout 
 * the entire optimization process, i.e., it lives longer than a single run of compass. The motivation is that past simulation output 
 * data should be kept to save effort as well as improve accuracy. 
 * A, b: vectors which specify the original problem constraint. 
 * activeSolutions: candidate active solutions for visited[starIndex]. Need to run MPR to find which are truly are
 * bestSolutions: best solutions found so far by GA
 * starIndex: the index in visited as the starting point as well as the local optimum solution returned by one run.
 * budget: the number of available simulation replications, when the algorithm finishes, budget becomes the remaining computing power
 * reocrds: to record evolution of COMPASS
 * recordInterval: determine the interval to record result
 * aType: original COMPASS or revised COMPASS
 * logFile: reference to the file where statistics of the algorithm run are kept. 
 * 
 * initialNumReps: initial number of reps assigned to a new solution
 * numCandidates: size of S_k, i.e. # of samplings within MPR (repetion ok)
 * alphac: CI-level for generalized COMPASS constraint placement procedure
 * alphae: CI-level for stopping test in Ver 2 of revised COMPASS to ensure that a optimal is indeed optimal before proceeding with COMPASS steps
 * deltae: indifference zone parameter for stopping test in Ver. 2 of revised COMPASS 
 * Return value: index of xStar.
 */

// Used in GetBestSolutions to sort visited solutions (pointers) 
bool ltsolution(Solution* x, Solution* y)
{
	return x->getSampleMean() < y->getSampleMean();
}
 
// bestSolutions has index in visited.
void GetBestSolutions(vector<long>& bestSolutions, const vector<Solution*>& visited, const vector<long> visitedThisRun)
{
	size_t size = visitedThisRun.size();
	unsigned int dimension = visited[0]->getDimension();
	list<Solution*> solutions;
	for (unsigned int i = 0; i < size; i++) solutions.push_back(visited[visitedThisRun[i]]);
	solutions.sort(ltsolution);
    bestSolutions.clear();
    list<Solution*>::iterator it = solutions.begin();
    for (unsigned int i = 0; i < size && i < dimension; i++)
    {
    	for (unsigned int j = 0; j < size; j++)
    	{
    		if (*it == visited[visitedThisRun[j]])
    		{
    			bestSolutions.push_back(visitedThisRun[j]);
    			break;
    		}
    	}
    	++it;
    }
}

long ISCv3(Simulator* simulator, bool stocsim, double deltam, double alpham, vector<Solution*>& visited, vector<double*> A, \
		vector<double> b, vector<long>& activeSolutions, long& budget, long starIndex, \
		ofstream& logFile, vector<Record*>& records, AlgorithmType aType, int initialNumReps, \
		int numCandidates, double alphac, double alphae, double deltae, int pruning_freq,int mpamode,int rmdmode,\
		long lastRep, double lastObjValue, double bestObjectiveValue, double trueBestObjectiveValue, \
		bool backTrackingTest, bool ocba, bool firstRunOfCompass, bool& optimum)
{
	double lastObjectiveValue = lastObjValue;
	double currentBestObjectiveValue = bestObjectiveValue;
	double trueCurrentBestObjectiveValue = trueBestObjectiveValue;
	long lastRepCount = 0;
	long baselineRepCount = lastRep;
	vector<long> Ek;
	vector<long> Sk;
	vector<double*> Ak;
	vector<double> bk;
	vector<long> bestSolutions;
	vector<long> visitedThisRun;
	vector<CornerMap*> corners;
	ObjMap objvalues; //a map<double, size_t) for predicated objective values of candidate solutions
	bool singularityFlag = false, optimumFound = false;
	long numOfUniqueCandidates = 1,	numOfTotalSamples = 0;
	long oldStarIndex = starIndex;
	long Nk = stocsim ? initialNumReps : 1;
	logFile << "\n";
	// Construct multimap for each coordinate
	unsigned int dimension = visited[0]->getDimension();
	for (size_t i = 0; i < dimension; i++)
	{
		corners.push_back(new CornerMap);
	}
	// For each visited solution, record its coordinate positions in the multimap
	for (size_t i = 0; i < visited.size(); i++)
	{
		for (size_t j = 0; j < dimension; j++)
		{	// The index in the visited vector is the content. The key is its coordinate position
			corners[j]->insert(make_pair(*(visited[i]->getDiscreteSolution() + j), i));
		}
	}
	optimum = false;
	// Step 0 of COMPASS
	/* activeSolutions contains the set of indices of visited solutions that actually 
	 * form MRP. Unlike in algorithm description, it is not exactly Ak because it doent's include
	 * original problem constraints, which will be added inside MPRv1(). 
	 * Its usage is as follows:
	 * 1. Before each call to ConstructMpr to construct a new MPR, activeSolutions contain the indexes
	 *    of those solutions that define MPR as judged by previous iteration (but excluding current sample best
	 * 	  plus newly sampled solutions returned by RMD();
	 * 2. After a call to ConstructMpr, it contains those that actually form MPR, excluding current sample best. 
	 * */
	/* Algorithm 11 (COMPASS in the integrated GA + COMPASS framework
	 * Step 1.1: 
	 * Ak * X >= bk is the constraint, including original problem formulation and MPR
	 * At the start of COMPASS, a set of active solutions and best solutions are given, so need to
	 * run MPR first to get the MPR.
	 */
	if (mpamode == 1)
	{
		hyperbox(stocsim, visited, starIndex, oldStarIndex, Ak, bk, A, b, activeSolutions, corners, true);
	}
	else
	{
		MPRv2(stocsim, visited, starIndex, oldStarIndex, Ak, bk, A, b, activeSolutions, bestSolutions, alphac);
	}
	// Record the starting point and activesolutions as those visited during this COMPASS run
	visitedThisRun.push_back(starIndex);
	for (size_t i = 0; i < activeSolutions.size(); i++)
	{
		if (find(visitedThisRun.begin(), visitedThisRun.end(), activeSolutions[i]) == visitedThisRun.end())
			visitedThisRun.push_back(activeSolutions[i]);
	}
	// Keep a record of best solutions	
	GetBestSolutions(bestSolutions, visited, visitedThisRun);
	// End of step 1.1
	clock_t time1 = clock();
	double seconds = time1/CLOCKS_PER_SEC;
	for (long iteration = 1; !optimumFound && 0 < budget; iteration++)
	{
		vector<Solution*> candidates;
		vector<Solution*> uniqueCandidates;
		vector<Solution*> candidatesToSim;
		/* STEP 2: Use RMD to sample S_k from MPR defined by activeSolutions. 
		 * Then Apply SAR to x^*_{k-1} and A_{k-1}
		 * Nk is the SAR. Apply SAR to xStar and activeSolutions. 
		 */
		Sk.clear();
		if (rmdmode == 1)
		{	// This is the coordinate sampling approach
			SimpleRMD(Ak, bk, visited[starIndex], numCandidates, WARMUP, candidates);
		}
		else
		{	// This is the default uniform sampling distribution
			RMD(Ak, bk, visited[starIndex], numCandidates, WARMUP, candidates);
		}
		// For debug only, look at newly sampled solutions.
		//	if (iteration > 1500 && iteration % 10 == 0) printVector(candidates);
		// Print out visited solution information for stochastic kriging input
		// Use Stochastic Kriging to predict the sample solution's quality and then only simulate selected ones
		ofstream iscmatlab("iscoutput.csv");
		require(iscmatlab.is_open(), "iscoutput.csv is not open");
		for (size_t ipr = 0; ipr < visited.size(); ipr++)
		{
			visited[ipr]->printX(iscmatlab);
		}
		iscmatlab.close();
		numOfUniqueCandidates = UniqueCandidates(candidates, uniqueCandidates, visited[starIndex]);
		// Use stoch kriging to determine what candidates to actually simulate 
		set<long> candidatesToKeep;
		ofstream predfile("isc_pred.csv");
		for (long count = 0; count < numOfUniqueCandidates; count++)
		{
			require(predfile.is_open(), "isc_pred.csv is not open");
			uniqueCandidates[count]->printX(predfile);
			candidatesToKeep.insert(count);
		}
		predfile.close();
		if (visited.size() >= 3) 
		{
			int ti=system ("runmatlab.bat");
			cout << "\nThe value returned was " << ti << endl;
			ifstream matlabout;
			do
			{
				matlabout.open("C:\\Users\\jxu13.SITE\\Documents\\MATLAB\\matlabout.csv");
			} while(!matlabout.is_open());
			MySleep(0.2);
			double tempy;
			char delimiter;
			objvalues.clear();
			for (size_t count = 0; count < numOfUniqueCandidates; count++) 
			{
				matlabout >> tempy; 
				objvalues.insert(make_pair(tempy, count));
			}
			matlabout.close();
			//predMean has the predicated mean of the solutions and we want to pick, say the top 30% for actual simulation 
			int sim_soln_cnt = 0;
			candidatesToKeep.clear();
			for (OMIter omit = objvalues.begin(); omit != objvalues.end() && sim_soln_cnt < numCandidates; omit++)
			{
				//cout << *uniqueCandidates[omit->second] << endl;
				candidatesToKeep.insert(omit->second);
				sim_soln_cnt++;
			}
		}
		if (numOfUniqueCandidates >= 1)
		{
			//Now only keep a portion of the solutions in uniqueCandidates to simulate
			//Need to release the memory taken by those we do not simulate
			for (long count = 0; count < numOfUniqueCandidates; count++)
			{
				if (candidatesToKeep.find(count) != candidatesToKeep.end()
				{

				}
			}
			for (long count = 0; count < numOfUniqueCandidates; count++)
			{
				// Check if this solution has already been visited. 
				long index = findSolution(visited, uniqueCandidates[count]);
				if (index < 0)
				{
					// New solutions need to be put into the global vector of visited solutions, and one for the local.
					visited.push_back(uniqueCandidates[count]);	
					visitedThisRun.push_back(long(visited.size() - 1));
					index = long(visited.size() - 1);
					for (size_t j = 0; j < dimension; j++)
					{	// The index in the visited vector is the content. The key is its coordinate position
						corners[j]->insert(make_pair(*(visited[index]->getDiscreteSolution() + j), index));
					}
					
				}
				else
				{
					// If this solution has already been visited, need to find out whether 
					// this solution was visited in this run or in previous runs
					long indexThisRun = find(visitedThisRun, index);
					if (indexThisRun == -1)	visitedThisRun.push_back(index);
					for (size_t j = 0; j < dimension; j++)
					{	// The index in the visited vector is the content. The key is its coordinate position
						corners[j]->insert(make_pair(*(visited[index]->getDiscreteSolution() + j), index));
					}
				}
				// Record the index into Sk
				Sk.push_back(index);
				// Debug use
				//cout << "New " << visited[index] << endl;
			}
		}
		// To avoid unnecessary complication with stopping test procedure, require that Sk has at least 2 reps. 
		for (size_t i = 0; i < Sk.size(); i++) 
		{
			int addrep = 2 - visited[Sk[i]]->getNumOfObservations();
			for (int repi = 0; repi < addrep; repi++)
			{
				visited[Sk[i]]->recordObservation(simulator->simulation(visited[Sk[i]]));
				numOfTotalSamples++;
				budget--;
			}
		}
		if (ocba && stocsim)
		// ocba only makes sense when it is a stoc simulation 
		{	// allocate some reps to active solutions and sample best according to an ocba like heuristic if ocba option is turned on.
			// There are deltan more reps to allocate, where deltan = sizeof(activesolutions).
			int deltan = int(activeSolutions.size());
			// Always add two more reps to current sample best
			for (size_t i = 0; i < 2; i++) visited[starIndex]->recordObservation(simulator->simulation(visited[starIndex]));	
			numOfTotalSamples += 2;
			budget -= 2;
			deltan -= 2;
			if (deltan > 0)
			{	// Get R
				double R = 0;
				for (size_t i = 0; i < activeSolutions.size(); i++) 
				{
					R += visited[activeSolutions[i]]->getSampleVariance() / max(1e-10, dabs(visited[activeSolutions[i]]->getSampleMean() - visited[starIndex]->getSampleMean()));
				}
				for (size_t i = 0; i < activeSolutions.size(); i++) 
				{
					double deltani = visited[activeSolutions[i]]->getSampleVariance() / max(1e-10, dabs(visited[activeSolutions[i]]->getSampleMean() - visited[starIndex]->getSampleMean())) / R * deltan;
					if (deltani >= 1)
					{
						long rounded_deltani = static_cast<long>(deltani);
						// Debug use
						if (rounded_deltani < 0)
						{
							cout << *visited[activeSolutions[i]];
							cout << *visited[starIndex];
							cout << rounded_deltani;
						}
						for (int j = 0; j < rounded_deltani; j++) 
						{
							visited[activeSolutions[i]]->recordObservation(simulator->simulation(visited[activeSolutions[i]]));
						}
						numOfTotalSamples += rounded_deltani;
						budget -= rounded_deltani;
					}
				}  
			}
		} 
		// If it is a deterministic simulation, only one rep
		long newReps = (long) ceil(initialNumReps * max(1,pow(log((double) iteration/2), 1.01)));
		if (stocsim)
		{
			Nk = (Nk >= newReps) ? Nk : newReps;
		}
		else
		{
			Nk = 1;
		}		
		// Now do the simulations for activesolutions
		for (size_t count = 0; count < activeSolutions.size(); count++)
		{
			if (visited[activeSolutions[count]]->getNumOfObservations() < Nk)
			{
				int newrep = Nk - visited[activeSolutions[count]]->getNumOfObservations();
				for (long rep = 0; rep < newrep; rep++)
				{
					visited[activeSolutions[count]]->recordObservation(simulator->simulation(visited[activeSolutions[count]]));	
				}
				numOfTotalSamples += newrep;
				budget -= newrep;
			}
		}
		// Simulate current sample best
		if (visited[starIndex]->getNumOfObservations() < Nk)
		{
			int newrep = Nk - visited[starIndex]->getNumOfObservations();
			for (long rep = 0; rep < newrep; rep++)
			{
				visited[starIndex]->recordObservation(simulator->simulation(visited[starIndex]));	
			}
			numOfTotalSamples += newrep;
			budget -= newrep;
		}
		/* STEP 2 ends here. 
		 * The vector uniqueCandidates contains new solutions, 
		 * and activeSolutions contain visited solutions that define MPR
		 */ 
		/* STEP 3: 
		 * Determine if x^*_{k-1} is still the current sample best compared with activeSolutions A_{k-1}
		 */ 
		// findMinSolution defined in myutilities.cpp. It compares solutions in activeSolutions with starIndex.
		long xtilde = findMinSolution(visited, starIndex, activeSolutions);
		if (!backTrackingTest || xtilde == starIndex)
		{	/* If x^*_{k-1} is still the current sample best, Step 5 begins.
			 * Apply SAR to S_k. 
			 */
			for (size_t i = 0; i < Sk.size(); i++)
			{	// Simulate this solution sampled within MPR, as indexed by "index"
				if (visited[Sk[i]]->getNumOfObservations() < Nk)
				{
					int newrep = Nk - visited[Sk[i]]->getNumOfObservations();
					for (int rep = 0; rep < newrep; rep++)
					{
						visited[Sk[i]]->recordObservation(simulator->simulation(visited[Sk[i]]));	
					}
					numOfTotalSamples += newrep;
					budget -= newrep;
				}
			} 
			// End of Step 5
			/* Step 6: 
			 * Determine new current sample best from E_k.
			 */
			oldStarIndex = starIndex;
			Ek.clear();
			for (size_t i = 0; i < Sk.size(); i++) Ek.push_back(Sk[i]);
			for (size_t i = 0; i < activeSolutions.size(); i++) Ek.push_back(activeSolutions[i]);
			Ek.push_back(starIndex);
			starIndex = findMinSolution(visited, Ek);
			// End of Step 6. Proceed to step 7.
		}
		else
		{
			/* If x^*_{k-1} is no longer the current sample best, Step 4 begins.
			 * Run stopping test on x^*_{k-1} (starIndex), S_k (uniqueCandidates), A_{k-1} (activeSolutions)
			 * with xtilde as the standard. Let x0 be the result returned by stopping test
			 */
			Ek.clear();
			for (size_t i = 0; i < activeSolutions.size(); i++) 
			{
				if (activeSolutions[i] != xtilde) Ek.push_back(activeSolutions[i]);
			}
			Ek.push_back(starIndex);
			for (size_t i = 0; i < Sk.size(); i++)
			{
				if (find(Ek.begin(), Ek.end(), Sk[i]) == Ek.end()) Ek.push_back(Sk[i]);
			}
			long testbudget = budget;
			// testbudget upon calling StoppingTest() gives the budget, and upon return, gives remaining budget
			// Step 4:  
			long x0;
			try
			{	
				x0 = StoppingTest(simulator, stocsim, deltae, alphae, visited, xtilde, Ek, testbudget, mpamode);
				numOfTotalSamples += budget - testbudget;
				budget = testbudget;
				// If budget exhausted, terminate iterations
				if (testbudget <= 0)	break;
			}
			catch(long)
			{
				numOfTotalSamples += (budget > 0) ? budget : 0;
				budget = 0;
				break;
			}
			// Record # of reps spent on stopping test
			oldStarIndex = starIndex;
			// Update xstar
			starIndex = x0;
			// Update A_{k-1} 
			for (size_t i = 0; i < activeSolutions.size(); i++) 
			{
				if (activeSolutions[i] == xtilde) activeSolutions[i] = oldStarIndex;
			}
			// End of Step 4.
		}
		/* Step 7
		 * Update B_k (bestSolutions).
		 * Then construct A_k
		 */
		// No need to use best solutions to construct MPR unless false optimum declared.
		bestSolutions.clear();
		// Right now activeSolutions only has A_{k-1}, need to combine S_k
		for (size_t i = 0; i < Sk.size(); i++) activeSolutions.push_back(Sk[i]);
		// Ak and bk will have all the constraints, including original problem formulation and MPR
		// A b are original problem formulation constraints
		// activeSolutions will then have the indices for those solutions already visited and define MPR
		// excluding current best solution
		// Do constraint pruning, i.e., determine redundant solutions in activeSolutions, every 5 iterations 
		if (mpamode == 1)
		{
			hyperbox(stocsim, visited, starIndex, oldStarIndex, Ak, bk, A, b, activeSolutions, corners, true);
		}
		else
		{
			MPRv2(stocsim, visited, starIndex, oldStarIndex, Ak, bk, A, b, activeSolutions, bestSolutions, alphac);
		}
		//MPRticks = clock() - current;
		// Record the optimization search process. First write to the repfile to record sample path
		for (long i = lastRepCount + 1; i < numOfTotalSamples; i++)
		{   // first write sample path
			logFile << baselineRepCount + i << " , " << lastObjectiveValue << "\n";
			// Next write records for averaged sample paths
			if (firstRunOfCompass) RecordResult(records, i + baselineRepCount, trueCurrentBestObjectiveValue);	
		}
		// Now update best value found so far and repcount
		lastObjectiveValue = simulator->gettruevalue(visited[starIndex]);
		// If at this stage of optimization, the current sample best seems to be better than
		// what we have kept so far, we'll think the current sample best changes. But when
		// recording the performance curve, we should use the true value of that solution.
		if (visited[starIndex]->getSampleMean() < currentBestObjectiveValue) 
		{	
			currentBestObjectiveValue = visited[starIndex]->getSampleMean();
			trueCurrentBestObjectiveValue = simulator->gettruevalue(visited[starIndex]);
		}
		lastRepCount = numOfTotalSamples;
		logFile << baselineRepCount + numOfTotalSamples << " , " << lastObjectiveValue << "\n";
		// Debug use only begin
		if (iteration % 50 == 0)
		{
			cout << "Iteration " << iteration << "  " << *visited[starIndex] << "   Rep #:  " << baselineRepCount + lastRepCount << endl << endl;
		}
		// Debug use only end
		// COMPASS may be run multiple times because NGA may identify several niches. If so, we should 
		if (firstRunOfCompass) RecordResult(records, baselineRepCount + lastRepCount, trueCurrentBestObjectiveValue);
		/* Step 8:
		 * Apply singleton test. If so, apply stopping test. If it is true optimal, exit and return the found
		 * solution; if it is declared false optimum, restart from the point the stopping procedure returns by updating
		 * xStar and activeSolutions and this while loop continues.
		 */
		// Debug use only begin
		//PrintActiveSolutions(visited, activeSolutions);
		// Debug use only end
		singularityFlag = SigularityCheck(Ak, bk, A, b, visited[starIndex], mpamode, visited, activeSolutions, budget, numOfTotalSamples,visitedThisRun, simulator, corners);
		while (singularityFlag)
		{
			// This solution has a singleton neighborhood. Now call statistical testing procedure to determine whether it is true local optimal
			// within the given statistical guarantee.
			//cout << "Singleton" << endl;
			//cout << *visited[starIndex];
			try 
			{
				long testbudget = budget;
				// Before the call, testBudget specifies the number of available samples the test procedure can sample
				// After the call, testBudget contains the # of samples taken.
				if (testbudget > 0)
				{
					cout << "Test local optimality of " << *visited[starIndex] << endl;
					long localOptimalIndex = StoppingTest(simulator, stocsim, deltam, alpham, visited, starIndex, testbudget, mpamode);
					cout << "The optimal after the test is " << *visited[localOptimalIndex] << endl;
					numOfTotalSamples += budget - testbudget;
					// Record the optimization search process after stopping test. First write to the repfile to record sample path
					for (long i = lastRepCount + 1; i < numOfTotalSamples; i++)
					{   // first write sample path
						logFile << baselineRepCount + i << " , " << lastObjectiveValue << "\n";
						// Next write records for averaged sample paths
						if (firstRunOfCompass) RecordResult(records, i + baselineRepCount, trueCurrentBestObjectiveValue);	
					}
					// Now update best value found so far and repcount
					lastObjectiveValue = simulator->gettruevalue(visited[localOptimalIndex]);
					if (visited[starIndex]->getSampleMean() < currentBestObjectiveValue) 
					{	
						currentBestObjectiveValue = visited[starIndex]->getSampleMean();
						trueCurrentBestObjectiveValue = simulator->gettruevalue(visited[starIndex]);
					}
					lastRepCount = numOfTotalSamples;
					logFile << baselineRepCount + numOfTotalSamples << " , " << lastObjectiveValue << "\n";
					if (firstRunOfCompass) RecordResult(records, baselineRepCount + lastRepCount, trueCurrentBestObjectiveValue);
					budget = testbudget;
					if (testbudget <= 0)	break;
					if (localOptimalIndex == starIndex)
					{ // The current sample best is believed to be the true local optimal
						optimumFound = true;
						optimum = true;
						break;
					}
					else
					{ 	// Otherwise, search continues with the new xStar as determined by the stopping test procedure
				  		// After stopping test, objective values change and should update Bk
				  		//GetBestSolutions(bestSolutions, visited, visitedThisRun);
						// Change xStar, then active solution set also needs be changed for revised compass. Switch old xstar and new xstar.
						// This is handled by mpr.cpp 
						oldStarIndex = starIndex;
						starIndex = localOptimalIndex;
						// Construct new MPR
						//current = clock();
						if (mpamode == 1)
						{
							hyperbox(stocsim, visited, starIndex, oldStarIndex, Ak, bk, A, b, activeSolutions, corners, true);
						}
						else
						{
							MPRv2(stocsim, visited, starIndex, oldStarIndex, Ak, bk, A, b, activeSolutions, bestSolutions, alphac);
						}
						// There is a MPR in the main loop, so need to add it up to MPRticks
						// Check if MPR is now again a singleton, if so, loop continues,
						// else, go back to step 1 (exit this while loop)
						singularityFlag = SigularityCheck(Ak, bk, A, b, visited[starIndex], mpamode, visited, activeSolutions, budget, numOfTotalSamples, visitedThisRun, simulator, corners);
					}
				}
				else
				{	
					// Record the optimization search process. First write to the repfile to record sample path
					for (long i = lastRepCount + 1; i < numOfTotalSamples; i++)
					{   // first write sample path
						logFile << baselineRepCount + i << " , " << lastObjectiveValue << "\n";
						// Next write records for averaged sample paths
						if (firstRunOfCompass) RecordResult(records, i + baselineRepCount, trueCurrentBestObjectiveValue);	
					}
					// Now update best value found so far and repcount
					lastObjectiveValue = simulator->gettruevalue(visited[starIndex]);
					if (visited[starIndex]->getSampleMean() < currentBestObjectiveValue) 
					{	
						currentBestObjectiveValue = visited[starIndex]->getSampleMean();
						trueCurrentBestObjectiveValue = simulator->gettruevalue(visited[starIndex]);
					}
					lastRepCount = numOfTotalSamples;
					logFile << baselineRepCount + numOfTotalSamples << " , " << lastObjectiveValue << "\n";
					if (firstRunOfCompass) RecordResult(records, baselineRepCount + lastRepCount, trueCurrentBestObjectiveValue);
					break;
				}
			}
			catch(long)
			{
				numOfTotalSamples +=  (budget > 0) ? budget : 0;
				budget = 0;
				// Record the optimization search process. First write to the repfile to record sample path
				for (long i = lastRepCount + 1; i < numOfTotalSamples; i++)
				{   // first write sample path
					logFile << baselineRepCount + i << " , " << lastObjectiveValue << "\n";
					// Next write records for averaged sample paths
					if (firstRunOfCompass) RecordResult(records, i + baselineRepCount, trueCurrentBestObjectiveValue);	
				}
				// Now update best value found so far and repcount
				lastObjectiveValue = simulator->gettruevalue(visited[starIndex]);
				if (visited[starIndex]->getSampleMean() < currentBestObjectiveValue) 
				{	
					currentBestObjectiveValue = visited[starIndex]->getSampleMean();
					trueCurrentBestObjectiveValue = simulator->gettruevalue(visited[starIndex]);
				}
				lastRepCount = numOfTotalSamples;
				logFile << baselineRepCount + numOfTotalSamples << " , " << lastObjectiveValue << "\n";
				if (firstRunOfCompass) RecordResult(records, baselineRepCount + lastRepCount, trueCurrentBestObjectiveValue);
				break;
			}	
		}
	}
	return starIndex;
}













