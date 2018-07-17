#include "masterheader.h"
#include <ctime>
#include <map>

using namespace std; 
const int WARMUP = 10;

//extern Engine *ep;
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

// Used in GetTrueBestSolutions to sort a given set of solutions (pointers) 
bool ltsolutionTrueObj(Solution* x, Solution* y)
{
	return x->getTrueObjectiveValue() < y->getTrueObjectiveValue();
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

// return the pointer to best solutions in the subset of solutions in solns, indexed by [is, ie-1]
Solution* GetTrueBestSolutions(Simulator* simulator, vector<Solution*>& solns, size_t is, size_t ie)
{
	size_t size = solns.size();
	unsigned int dimension = solns[0]->getDimension();
	list<Solution*> solutions;
	for (size_t i = is; i < ie; i++) 
	{
		double trueobj = simulator->gettruevalue(solns[i]);
		solns[i]->setTrueObjectiveValue(trueobj);
		solutions.push_back(solns[i]);
	}
	solutions.sort(ltsolutionTrueObj);
    list<Solution*>::iterator it = solutions.begin();
    return *it;
}

/*long ISCv3(Simulator* simulator, bool stocsim, double deltam, double alpham, vector<Solution*>& visited, vector<double*> A, \
		vector<double> b, vector<long>& activeSolutions, long& budget, long starIndex, \
		ofstream& logFile, vector<Record*>& records, AlgorithmType aType, int initialNumReps, \
		int numCandidates, double alphac, double alphae, double deltae, int pruning_freq,int mpamode,int rmdmode,\
		long lastRep, double lastObjValue, double bestObjectiveValue, double trueBestObjectiveValue, \
		bool backTrackingTest, bool ocba, bool firstRunOfCompass, bool& optimum, int metamodel)*/
long ISCv3(Simulator* simulator, vector<Solution*>& visited, vector<double*> A, vector<double> b, ISCParameters* iscparams, \
	vector<long>& activeSolutions, long& budget, long starIndex, ofstream& logFile, vector<Record*>& records, long lastRep, \
	double lastObjValue, double bestObjectiveValue, double trueBestObjectiveValue, bool firstRunOfCOMPASS, bool& optimum)

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
	ObjMap eivalues; //a map<double, size_t) for predicated expected improvement of candidate solutions
	vector<double> SKPredMean; // Mean of SK predictions, indexed by solution index as in SK output file
	vector<double> SKPredStd; // Std of SK predictions, indexed by solution index as in SK output file
	vector<size_t> PredRanking; //The indexes of prediction points ranked by sampling from the MVN conditional distribution of prediction points
	bool singularityFlag = false, optimumFound = false;
	size_t numOfUniqueCandidates = 1;	
	long numOfTotalSamples = 0;
	long oldStarIndex = starIndex;
	long Nk = iscparams->get_stocsim() ? iscparams->get_initialNumReps() : 1;
	logFile << "\n";
	/*Engine *ep; //Define matlab engine and test it
	if (!(ep=engOpen(NULL)))  // Nov 11,2012 
	{
		cout <<"Can't start Matlab engine!" <<endl;
		//cin.get();
		exit(1);
	}
	engSetVisible(ep, 0); // Set matlab engine command window invisible*/
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
	if (iscparams->get_mpamode() == 1)
	{
		hyperbox(iscparams->get_stocsim(), visited, starIndex, oldStarIndex, Ak, bk, A, b, activeSolutions, corners, true);
	}
	else
	{
		MPRv2(iscparams->get_stocsim(), visited, starIndex, oldStarIndex, Ak, bk, A, b, activeSolutions, bestSolutions, iscparams->get_alphac());
	}
	// In case of local search only, simulate the initial solution
	cout << visited[starIndex];
	for (int repi = visited[starIndex]->getNumOfObservations(); repi < iscparams->get_initialNumReps(); repi++)
	{
		double obj = simulator->simulation(visited[starIndex]);
		visited[starIndex]->recordObservation(obj);
		numOfTotalSamples++;
		budget--;
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
	ofstream skdebug("skdebug.csv");
	require(skdebug.is_open(), "skdebug not open");
	ofstream skinfo("skinfo.csv");
	require(skinfo.is_open(), "skinfo not open");
	skdebug.close();		
	skinfo.close();
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
		int numSampleSolns = (iscparams->get_metamodel() > 0) ? iscparams->get_metaPredNum() : iscparams->get_numOfCandidates();
		double mpavol;
		int lhsnum;
		//For OO option, the number of sampled solution comes from the alignment prob calculation
		mpavol = calcVolMPA(Ak, bk, dimension);
		// We turn Metamodel off if mpavol is too small, else turn it on
		if (mpavol < 100)
		{
			iscparams->set_metamodel(0);
		}
		else
		{
			iscparams->set_metamodel(iscparams->get_metamodelbak());
		}
		if (iscparams->get_metamodel() > 2 && mpavol > 3 * iscparams->get_numOfCandidates())
		{	// Use SK design only when MPA is three times larger than number of solutions to sample w/o SK 
			// Hardcode the goal of OO G to be 0.01% of all solutions, if it is smaller than 1, make it one
			double G = mymax(mpavol*0.0001, 1.0); 
			double OOSize = SampleSizeOO(mpavol, G, 1);
			// Note alignment level k is hardcoded to 1
			numSampleSolns = static_cast<int>(mymin(OOSize, 5000.0)); // 1/2/2012 hardcode the max number of solutions to simulate
			numSampleSolns = static_cast<int>(mymin(static_cast<double>(numSampleSolns), mpavol*0.1));
		}
		if (iscparams->get_rmdmode() == 1)
		{	// This is the coordinate sampling approach
			SimpleRMD(Ak, bk, visited[starIndex], numSampleSolns, WARMUP, candidates);
		}
		else
		{	// This is the default uniform sampling distribution
			RMD(Ak, bk, visited[starIndex], numSampleSolns, WARMUP, candidates);
		}
		numOfUniqueCandidates = UniqueCandidates(candidates, uniqueCandidates, visited[starIndex]);
		// For debug only, look at newly sampled solutions.
		//	if (iteration > 1500 && iteration % 10 == 0) printVector(candidates);
		// Print out visited solution information for stochastic kriging input
		// Use Stochastic Kriging to predict the sample solution's quality and then only simulate selected ones
		//ofstream skdebug("skdebug.csv", ios::app);
		ofstream skdebug("skdebug.csv");
		require(skdebug.is_open(), "skdebug not open");
		//ofstream skinfo("skdinfo.csv", ios::app);
		ofstream skinfo("skdinfo.csv");
		require(skinfo.is_open(), "skinfo not open");
		// JX 11/8/11 adds one more metamodel option: 1 for using prediction mean to rank, 2 for using expected improvement to rank	
		if (iscparams->get_metamodel() == 1 || iscparams->get_metamodel() == 2)
		{	// Below is for EI or predicted value to rank solutions 
			/*ofstream iscmatlab("iscoutput.csv");
			require(iscmatlab.is_open(), "iscoutput.csv is not open");
			skdebug << endl << "\nIteration " << iteration << endl;
			for (size_t ipr = 0; ipr < visited.size(); ipr++)
			{
				visited[ipr]->printX(iscmatlab);
			}
			iscmatlab.close();
			// Use stoch kriging to determine what candidates to actually simulate 
			set<long> candidatesToKeep;
			ofstream predfile("isc_pred.csv");
			for (long count = 0; count < numOfUniqueCandidates; count++)
			{
				require(predfile.is_open(), "isc_pred.csv is not open");
				uniqueCandidates[count]->printX(predfile);
				uniqueCandidates[count]->printX(skdebug);
			}
			for (long count = 0; count < numOfUniqueCandidates && count < iscparams->get_numOfCandidates(); count++)
			{
				candidatesToKeep.insert(count); //In case SK does not work, we randomly pick some sampled solutions for simulation
			}
			predfile.close();
			if (visited.size() >= 3) 
			{
				int ti=system ("runmatlab.bat");
				cout << "\nThe value returned was " << ti << endl;
				ifstream matlabout;
				// Wait for Matlab for a maximum of 20 seconds in case there is a problem and it does not finish
				clock_t time1 = clock();
				clock_t time2 = time1 + 20*CLOCKS_PER_SEC;
				clock_t currtime;
				do
				{
					matlabout.open("C:\\Users\\jxu13.SITE\\Documents\\MATLAB\\matlabout.csv");
					currtime = clock();
				} while(!matlabout.is_open() && currtime < time2);
				MySleep(0.2);
				if (matlabout.is_open())
				{
					double tempy, tempsd;
					objvalues.clear();
					ifstream skstd("C:\\Users\\jxu13.SITE\\Documents\\MATLAB\\std.csv");
					require(skstd.is_open(), "SK std file not open");
					eivalues.clear();
					for (int count = 0; count < numOfUniqueCandidates; count++) 
					{
						matlabout >> tempy; 
						skstd >> tempsd;
						objvalues.insert(make_pair(tempy, count));
						eivalues.insert(make_pair(tempy-tempsd, count));
					}
					matlabout.close();
					skstd.close();
					//predMean has the predicated mean of the solutions and we want to pick, say the top 30% for actual simulation 
					int sim_soln_cnt = 0;
					candidatesToKeep.clear();
					if (iscparams->get_metamodel() == 1)
					{
						for (OMIter omit = objvalues.begin(); omit != objvalues.end() && sim_soln_cnt < iscparams->get_numOfCandidates(); omit++)
						{
							cout << omit->first  << "  pred,  " << *uniqueCandidates[omit->second] << endl;
							skdebug << omit->first  << "  pred  " << *uniqueCandidates[omit->second] << endl;
							candidatesToKeep.insert(omit->second);
							sim_soln_cnt++;
						}
					}
					else if (iscparams->get_metamodel() == 2)
					{
						for (OMIter omit = eivalues.begin(); omit != eivalues.end() && sim_soln_cnt < iscparams->get_numOfCandidates(); omit++)
						{
							cout << omit->first  << "  EI,  " << *uniqueCandidates[omit->second] << endl;
							skdebug << omit->first  << "  EI  " << *uniqueCandidates[omit->second] << endl;
							candidatesToKeep.insert(omit->second);
							sim_soln_cnt++;
						}
					}
					else
					{ 
						cout << "Error, Metamodel parameter input not supported\n";
						exit(5);
					}
				}
				else 
				{
					cout << "Matlab does not finish in 20 seconds\n";
				}
			}
			if (numOfUniqueCandidates >= 1)
			{
				//Now only keep a portion of the solutions in uniqueCandidates to simulate
				//Need to release the memory taken by those we do not simulate
				candidatesToSim.clear();
				for (long count = 0; count < numOfUniqueCandidates; count++)
				{
					if (candidatesToKeep.find(count) != candidatesToKeep.end())
					{
						candidatesToSim.push_back(uniqueCandidates[count]);
					}
					else
					{
					delete uniqueCandidates[count];
					}
				}
				uniqueCandidates.clear();
				uniqueCandidates = candidatesToSim;
				numOfUniqueCandidates = uniqueCandidates.size();
			}*/
		}
		//else if ((iscparams->get_metamodel() == 3 || iscparams->get_metamodel() == 5) && iscparams->get_mpamode()==1) 
		else if ((iscparams->get_metamodel() == 3 || iscparams->get_metamodel() == 4 || iscparams->get_metamodel() == 5) && iscparams->get_mpamode()==1) 
		{	// Below is for OO + SK, right now it is only for AHA because I need to compute the volume of the MPA, 12/26/2011
			// Use SK only when MPA is at least 3 times larger than number of solutions to sample w/o SK 
			// First we generate design points within the current MPA 
			// We use AHA w/o global phase with latin hypercube design. Before we enter AHA, use a latin hypercube to come up with design points 
			// Write latin hypercube design space (only for hyperbox geometry), defined by Ak (lhs) and bk (rhs) 
			//ofstream lhs_spec("isclhs.csv");
			ofstream lhs_spec("C:\\Users\\jxu13\\Documents\\MATLAB\\isclhs.csv");
			require(lhs_spec.is_open(), "Error: isclhs.csv not open");
			// Require lhs design points to be given by get_numOfCandidates(), but at most 5% of the MPA volume
			if (iscparams->get_numOfCandidates() > mpavol*0.05)
				lhsnum = static_cast<int>(mpavol*0.05); 
			else
				lhsnum = iscparams->get_numOfCandidates();
			lhs_spec << lhsnum << "," << visited[0]->getDimension() << endl;
			// Print Ax >= b constraints in csv format
			PrintConstraints(visited[0]->getDimension(), Ak, bk, lhs_spec);
			lhs_spec.close();
			//int tilhs=system ("runmatlablhs.bat");
			//cout << "\nThe value returned was " << tilhs << endl;
			/*int engRt = engEvalString(ep, "sk4isc_lhs");   // Nov 11,2012 
			if (engRt == 1)
			{
				cout << "Matlab engine error" << endl;
				engRt = engEvalString(ep, "sk4isc_lhs");
				//engClose(ep); // Close matlab engin,Nov 11,2012 
				return -1;
			}*/
			//engClose(ep); // Close matlab engin,Nov 11,2012 
			ifstream matlabout_lhs;
			// Wait for Matlab for a maximum of 60 seconds in case there is a problem and it does not finish
			clock_t time1_lhs = clock();
			clock_t time2_lhs = time1_lhs + 300*CLOCKS_PER_SEC;
			clock_t currtime_lhs;
			do
			{
				matlabout_lhs.open("C:\\Users\\jxu13\\Documents\\MATLAB\\lhsdesign.txt");
				currtime_lhs = clock();
			} while(!matlabout_lhs.is_open() && currtime_lhs < time2_lhs);
			//MySleep(10);
			if (matlabout_lhs.is_open())
			{
				for (int numdes = 0; numdes < lhsnum; numdes++)
				{
					vector<int> xd(visited[0]->getDimension(), 0);
					for (int dd = 0; dd < visited[0]->getDimension(); dd++)
					{
						matlabout_lhs >> xd[dd];
					}
					Solution* newdesign = new Solution(xd);
					// Check if this design point is also selected as prediction point. If so, remove it from prediction point list
					numOfUniqueCandidates = UniqueCandidates(uniqueCandidates, newdesign, true);
					// Check if this solution has already been visited. 
					long index = findSolution(visited, newdesign);
					if (index < 0)
					{
						// New solutions need to be put into the global vector of visited solutions, and one for the local.
						visited.push_back(newdesign);	
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
					Sk.push_back(index); // Sk now contains LHS design points
				}
				// We simulate LHS design points to make sure they have at least an initial # of reps
				for (size_t i = 0; i < Sk.size(); i++) 
				{
					for (int repi = visited[Sk[i]]->getNumOfObservations(); repi < iscparams->get_initialNumReps(); repi++)
					{
						double obj = simulator->simulation(visited[Sk[i]]);
						visited[Sk[i]]->recordObservation(obj);
						numOfTotalSamples++;
						budget--;
					}
				}
			}
			else
			{
				cout << "Error, Matlab doesn't finish LHS" << endl;
				exit(6);
			}
			// The following file contains all visited solutions info as SK design points 
			ofstream iscmatlab("C:\\Users\\jxu13\\Documents\\MATLAB\\iscoutput.csv");
			require(iscmatlab.is_open(), "iscoutput.csv is not open");
			skdebug << endl << "\nIteration " << iteration << endl;
			for (size_t ipr = 0; ipr < visited.size(); ipr++)
			{
				visited[ipr]->printX(iscmatlab);
			}
			iscmatlab.close();
			// Use stoch kriging to determine what candidates to actually simulate 
			set<long> candidatesToKeep;
			// The following file contains points sampled uniformly from MPA to be predicted by SK
			ofstream predfile("C:\\Users\\jxu13\\Documents\\MATLAB\\isc_pred.csv");
			require(predfile.is_open(), "isc_pred.csv is not open");
			for (long count = 0; count < numOfUniqueCandidates; count++)
			{
				require(predfile.is_open(), "isc_pred.csv is not open");
				uniqueCandidates[count]->printX(predfile);
				uniqueCandidates[count]->printX(skdebug);
			}
			// The following file contains points sampled uniformly from MPA to be predicted by SK
			ofstream prednum("C:\\Users\\jxu13\\Documents\\MATLAB\\isc_prednum.csv");
			require(prednum.is_open(), "isc_prednum.csv is not open");
			prednum << iscparams->get_numOfCandidates();
			prednum.close();
			for (long count = 0; count < numOfUniqueCandidates && count < iscparams->get_numOfCandidates(); count++)
			{
				candidatesToKeep.insert(count); //In case SK does not work, we randomly pick some sampled solutions for simulation
			}
			predfile.close();
			SKPredMean.clear();
			SKPredStd.clear();
			PredRanking.clear();
			// Print out the current sample best true obj, the true best of the newly sampled w/o SK, the true best obj of all predicted x using SK, its SK predicted value, SK 
			// predicted best mean and its true value, SK predicted most expected improvement and its true value 
			//skinfo << endl << iteration << "," << simulator->gettruevalue(visited[starIndex]) << ",";
			//cout << endl << iteration << "," << simulator->gettruevalue(visited[starIndex]) << ",";
			// Find out the best of the newly sampled if we use the old AHA/COMPASS way
			//skinfo << GetTrueBestSolutions(simulator, uniqueCandidates, 0, iscparams->get_numOfCandidates())->getTrueObjectiveValue() << ",";
			//cout << GetTrueBestSolutions(simulator, uniqueCandidates, 0, iscparams->get_numOfCandidates())->getTrueObjectiveValue() << ",";
			//skinfo << GetTrueBestSolutions(simulator, uniqueCandidates, 0, uniqueCandidates.size())->getTrueObjectiveValue() << ",";
			//cout << GetTrueBestSolutions(simulator, uniqueCandidates, 0, uniqueCandidates.size())->getTrueObjectiveValue() << ",";
			//int ti=system ("runmatlab.bat");
			//cout << "\nThe value returned was " << ti << endl;
			
			/*engRt = engEvalString(ep, "sk4isc1");   // Nov 11,2012 
			if (engRt == 1)
			{
				cout << "Matlabe eng sk4isc error" << endl;
				exit(6);
			}*/
			ifstream matlabout, skstd, skranking;
			// Wait for Matlab for a maximum of 20 seconds in case there is a problem and it does not finish
			clock_t time1 = clock();
			clock_t time2 = time1 + 120*CLOCKS_PER_SEC;
			clock_t currtime;
			do
			{
				matlabout.open("C:\\Users\\jxu13\\Documents\\MATLAB\\matlabout.csv");
				skstd.open("C:\\Users\\jxu13\\Documents\\MATLAB\\std.csv");
				skranking.open("C:\\Users\\jxu13\\Documents\\MATLAB\\ranking.csv");
				currtime = clock();
			} while((!matlabout.is_open() || !skstd.is_open() || !skranking.is_open()) && currtime < time2);
			//MySleep(30);
			matlabout.close(); skstd.close(); skranking.close();
			matlabout.open("C:\\Users\\jxu13\\Documents\\MATLAB\\matlabout.csv");
			skstd.open("C:\\Users\\jxu13\\Documents\\MATLAB\\std.csv");
			skranking.open("C:\\Users\\jxu13\\Documents\\MATLAB\\ranking.csv");
			if (matlabout.is_open() && skranking.is_open())
			{
				double tempy, tempsd;
				size_t predindex, numtopxpred;
				objvalues.clear();
				eivalues.clear();
				for (int count = 0; count < numOfUniqueCandidates; count++) 
				{
					matlabout >> tempy; 
					skstd >> tempsd;
					objvalues.insert(make_pair(tempy, count));
					eivalues.insert(make_pair(tempy-tempsd, count));
					SKPredMean.push_back(tempy);
					SKPredStd.push_back(tempsd);
				}
				matlabout.close();
				skstd.close();
				skranking >> numtopxpred; // The first number in the ranking file is the number of indexes to read 
				for (size_t ip = 0; ip < numtopxpred; ip++)
				{
					skranking >> predindex;
					PredRanking.push_back(predindex-1); //The ranking is computed by Matlab, so need to offset index by 1
				}
				skranking.close();
				//predMean has the predicated mean of the solutions and we want to pick, say the top 30% for actual simulation 
				int sim_soln_cnt = 0;
				candidatesToKeep.clear();
				if (iscparams->get_metamodel() == 1 || iscparams->get_metamodel() == 3)
				{
					for (OMIter omit = objvalues.begin(); omit != objvalues.end() && sim_soln_cnt < iscparams->get_numOfCandidates(); omit++)
					{
						cout << omit->first  << "  pred,  " << *uniqueCandidates[omit->second] << " actual, " << simulator->gettruevalue(uniqueCandidates[omit->second]) << endl;
						skdebug << omit->first  << "  pred  " << *uniqueCandidates[omit->second] << " actual, " << simulator->gettruevalue(uniqueCandidates[omit->second]) << endl;
						candidatesToKeep.insert(omit->second);
						sim_soln_cnt++;
					}
					//skinfo << objvalues.begin()->first << "," << simulator->gettruevalue(uniqueCandidates[objvalues.begin()->second]) << "," << SKPredStd[objvalues.begin()->second] << ",";
					//cout << objvalues.begin()->first << "," << simulator->gettruevalue(uniqueCandidates[objvalues.begin()->second]) << "," << SKPredStd[objvalues.begin()->second] << ",";
				}
				else if (iscparams->get_metamodel() == 2 || iscparams->get_metamodel() == 4)
				{
					for (OMIter omit = eivalues.begin(); omit != eivalues.end() && sim_soln_cnt < iscparams->get_numOfCandidates(); omit++)
					{
						cout << omit->first  << "  EI,  " << *uniqueCandidates[omit->second] << " actual, " << simulator->gettruevalue(uniqueCandidates[omit->second]) << endl;
						skdebug << omit->first  << "  EI  " << *uniqueCandidates[omit->second] << " actual, " << simulator->gettruevalue(uniqueCandidates[omit->second]) << endl;
						candidatesToKeep.insert(omit->second);
						sim_soln_cnt++;
					}
					// skinfo records the ei/obj value of the best solution picked by SK
					//skinfo << eivalues.begin()->first << "," << simulator->gettruevalue(uniqueCandidates[objvalues.begin()->second]) << "," << SKPredStd[eivalues.begin()->second] << ",";
					//cout << eivalues.begin()->first << "," << simulator->gettruevalue(uniqueCandidates[objvalues.begin()->second]) << "," << SKPredStd[eivalues.begin()->second] << ",";
				}
				else if (iscparams->get_metamodel() == 5)
				{	// Another option is to use the simulation ranking of prediction points 
					for (size_t iprk = 0; iprk < PredRanking.size(); iprk++)
					{
						candidatesToKeep.insert(PredRanking[iprk]);
					}
				}
				else
				{ 
					cout << "Error, Metamodel parameter input not supported\n";
					exit(5);
				}
			} //end-if for determining which prediction points to include in the sampling set to be simulated
			else 
			{
				cout << "Matlab does not finish in 20 seconds\n";
				skinfo << ",,,," << endl;
				cout << ",,,," << endl;
				// If SK does not work, most likely it will not work any further, disable it
				//iscparams->set_metamodel(0);
			}
			if (numOfUniqueCandidates >= 1)
			{
				//Now only keep a portion of the solutions in uniqueCandidates to simulate
				//Need to release the memory taken by those we do not simulate
				candidatesToSim.clear();
				for (long count = 0; count < numOfUniqueCandidates; count++)
				{
					if (candidatesToKeep.find(count) != candidatesToKeep.end())
					{
						candidatesToSim.push_back(uniqueCandidates[count]);
					}
					else
					{
						delete uniqueCandidates[count];
					}
				}
				uniqueCandidates.clear();
				uniqueCandidates = candidatesToSim;
				numOfUniqueCandidates = uniqueCandidates.size();
			}
			// Now print the obj of the best solution in the newly sampled solution Sk
			//skinfo << GetTrueBestSolutions(simulator, uniqueCandidates, 0, uniqueCandidates.size())->getTrueObjectiveValue() << endl;
			//cout << GetTrueBestSolutions(simulator, uniqueCandidates, 0, uniqueCandidates.size())->getTrueObjectiveValue() << endl;
		} // end else-if using OO with SK to select points to simulate
		// Now uniqueCandidates has the solutions we sampled within MPA to simualte. If we don't use SK, it's simply the # of solutions we sample from MPA. 
		// If we use SK, they are the "top" solutions from a lot of prediction points selected via SK predictions, note Sk now already has LHS design points inside the MPA 
		// This is the end of codes added for SK. Below are COMPASS/AHA codes 
		if (numOfUniqueCandidates >= 1)
		{
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
				// Record the index into Sk, note if we use SK, LHS design within MPA are already in Sk
				Sk.push_back(index);
				// Debug use
				//cout << "Sk " << visited[index] << endl;
				//skdebug << "Sk " << visited[index] << endl;
			}
		}
		skdebug.close();
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
		if (iscparams->get_ocba() && iscparams->get_stocsim())
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
		long newReps = (long) ceil(iscparams->get_initialNumReps() * max(1,pow(log((double) iteration/2), 1.01)));
		if (iscparams->get_stocsim())
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
		if (!iscparams->get_backtracking() || xtilde == starIndex)
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
				x0 = StoppingTest(simulator, iscparams->get_stocsim(), iscparams->get_backdelta(), iscparams->get_backalpha(), \
					visited, xtilde, Ek, testbudget, iscparams->get_mpamode());
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
		if (iscparams->get_mpamode() == 1)
		{
			hyperbox(iscparams->get_stocsim(), visited, starIndex, oldStarIndex, Ak, bk, A, b, activeSolutions, corners, true);
		}
		else
		{
			MPRv2(iscparams->get_stocsim(), visited, starIndex, oldStarIndex, Ak, bk, A, b, activeSolutions, bestSolutions, iscparams->get_alphac());
		}
		//MPRticks = clock() - current;
		// Record the optimization search process. First write to the repfile to record sample path
		for (long i = lastRepCount + 1; i < numOfTotalSamples; i++)
		{   // first write sample path
			logFile << baselineRepCount + i << " , " << lastObjectiveValue << "\n";
			// Next write records for averaged sample paths
			if (firstRunOfCOMPASS) RecordResult(records, i + baselineRepCount, trueCurrentBestObjectiveValue);	
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
		if (firstRunOfCOMPASS) RecordResult(records, baselineRepCount + lastRepCount, trueCurrentBestObjectiveValue);
		/* Step 8:
		 * Apply singleton test. If so, apply stopping test. If it is true optimal, exit and return the found
		 * solution; if it is declared false optimum, restart from the point the stopping procedure returns by updating
		 * xStar and activeSolutions and this while loop continues.
		 */
		// Debug use only begin
		//PrintActiveSolutions(visited, activeSolutions);
		// Debug use only end
		singularityFlag = SigularityCheck(Ak, bk, A, b, visited[starIndex], iscparams->get_mpamode(), visited, activeSolutions, budget, \
			numOfTotalSamples,visitedThisRun, simulator, corners);
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
					long localOptimalIndex = StoppingTest(simulator, iscparams->get_stocsim(), iscparams->get_localoptdelta(), iscparams->get_localoptalpha(), \
						visited, starIndex, testbudget, iscparams->get_mpamode());
					cout << "The optimal after the test is " << *visited[localOptimalIndex] << endl;
					numOfTotalSamples += budget - testbudget;
					// Record the optimization search process after stopping test. First write to the repfile to record sample path
					for (long i = lastRepCount + 1; i < numOfTotalSamples; i++)
					{   // first write sample path
						logFile << baselineRepCount + i << " , " << lastObjectiveValue << "\n";
						// Next write records for averaged sample paths
						if (firstRunOfCOMPASS) RecordResult(records, i + baselineRepCount, trueCurrentBestObjectiveValue);	
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
					if (firstRunOfCOMPASS) RecordResult(records, baselineRepCount + lastRepCount, trueCurrentBestObjectiveValue);
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
						if (iscparams->get_mpamode() == 1)
						{
							hyperbox(iscparams->get_stocsim(), visited, starIndex, oldStarIndex, Ak, bk, A, b, activeSolutions, corners, true);
						}
						else
						{
							MPRv2(iscparams->get_stocsim(), visited, starIndex, oldStarIndex, Ak, bk, A, b, activeSolutions, bestSolutions, iscparams->get_alphac());
						}
						// There is a MPR in the main loop, so need to add it up to MPRticks
						// Check if MPR is now again a singleton, if so, loop continues,
						// else, go back to step 1 (exit this while loop)
						singularityFlag = SigularityCheck(Ak, bk, A, b, visited[starIndex], iscparams->get_mpamode(), visited, activeSolutions, budget, \
							numOfTotalSamples, visitedThisRun, simulator, corners);
					}
				}
				else
				{	
					// Record the optimization search process. First write to the repfile to record sample path
					for (long i = lastRepCount + 1; i < numOfTotalSamples; i++)
					{   // first write sample path
						logFile << baselineRepCount + i << " , " << lastObjectiveValue << "\n";
						// Next write records for averaged sample paths
						if (firstRunOfCOMPASS) RecordResult(records, i + baselineRepCount, trueCurrentBestObjectiveValue);	
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
					if (firstRunOfCOMPASS) RecordResult(records, baselineRepCount + lastRepCount, trueCurrentBestObjectiveValue);
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
					if (firstRunOfCOMPASS) RecordResult(records, i + baselineRepCount, trueCurrentBestObjectiveValue);	
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
				if (firstRunOfCOMPASS) RecordResult(records, baselineRepCount + lastRepCount, trueCurrentBestObjectiveValue);
				break;
			}	
		}
	}
	//engClose(ep); // Close matlab engin,Nov 11,2012 
	return starIndex;
}













