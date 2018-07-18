#include "masterheader.h"

using namespace std;

ISCParameters* ReadInput(std::ifstream& inputFile, string& output, std::vector<int>& x0, std::vector<double*>& A, std::vector<double>& b)
{
	// Constraints are of the form Ax>=b. The first line is dimension and size of constraints
	// The second line is the starting point
	// Then follow constraints per line.
	int dimension, size;
	long budget;
	int reps, backtracking, ocba, gagen, dominantniche, cleanup, stocsim, initialNumReps, numOfCandidates, elitism, pruning_freq;
	double  globaldelta, backalpha, backdelta, localoptalpha, localoptdelta, cleanupalpha, cleanupdelta, gabudget; 
	int mpamode, rmdmode, metamodel, metaprednum;
	if (!inputFile.is_open())
	{
		cout << "Error: invalid input file" << endl;
		return NULL;
	}
	string temp;
	inputFile >> temp >> output;
	inputFile >> temp >> budget;
	inputFile >> temp >> reps;
	inputFile >> temp >> backtracking;
	inputFile >> temp >> ocba;
	inputFile >> temp >> gagen;
	inputFile >> temp >> dominantniche;
	inputFile >> temp >> cleanup;
	inputFile >> temp >> stocsim;
	inputFile >> temp >> globaldelta;
	inputFile >> temp >> backalpha;
	inputFile >> temp >> backdelta;
	inputFile >> temp >> localoptalpha;
	inputFile >> temp >> localoptdelta;
	inputFile >> temp >> cleanupalpha;
	inputFile >> temp >> cleanupdelta;
	inputFile >> temp >> initialNumReps;
	inputFile >> temp >> numOfCandidates;
	inputFile >> temp >> elitism;
	inputFile >> temp >> pruning_freq;
	inputFile >> temp >> mpamode;
	inputFile >> temp >> rmdmode;
	inputFile >> temp >> gabudget;
	inputFile >> temp >> metamodel;
	inputFile >> temp >> metaprednum;
	require(cleanupdelta >= 0, "Indifference parameter for clean-up must be positive");
	ISCParameters* iscparam = new ISCParameters(cleanupdelta,  cleanupalpha,  localoptdelta,  localoptalpha,  backdelta,  backalpha, \
		 globaldelta, budget, gabudget, reps, gagen, initialNumReps, numOfCandidates, pruning_freq, \
		  mpamode, rmdmode, metamodel, metaprednum, backtracking, ocba, dominantniche, cleanup, stocsim, elitism);
	/*budget = (budget < 0) ? BUDGET : budget;
	reps = (reps < 0) ? REPS : reps;
	backtracking = (backtemp < 0) ? BACKTRACKING : static_cast<bool>(backtemp);
	ocba = (ocbatemp < 0) ? OCBA : static_cast<bool>(ocbatemp);
	gagen = (gagen < 0) ? GAGEN : gagen;
	dominantniche = (dntemp < 0) ? DOMINANTNICHE : static_cast<bool>(dntemp);
	cleanup = (cleanuptemp < 0) ? CLEANUP : static_cast<bool>(cleanuptemp);
	stocsim = (stocsimtemp == 0) ? false : true;
	globaldelta = (globaldelta < 0) ? cleanupdelta : globaldelta;
	backalpha = (backalpha < 0) ? BACKALPHA : backalpha;
	backdelta = (backdelta < 0) ? cleanupdelta : backdelta;
	localoptalpha = (localoptalpha < 0) ? LOCALOPTALPHA : localoptalpha;
	localoptdelta = (localoptdelta < 0) ? cleanupdelta : localoptdelta;
	cleanupalpha = (cleanupalpha < 0) ? CLEANUPALPHA : cleanupalpha;
	initialNumReps = (initialNumReps < 0) ? INITIALNUMREPS : initialNumReps;
	// If it is a stochastic sim problem, require at least 2 initial reps
	if (initialNumReps < 2 && stocsim) initialNumReps = 2;
	numOfCandidates = (numOfCandidates < 0) ? NUMOFCANDIDATES : numOfCandidates;
	elitism = (elitismtemp < 0) ? ELITISM : static_cast<bool>(elitismtemp);
	pruning_freq = (pruning_freq < 0) ? PRUNINGFREQ : pruning_freq;
	mpamode = (mpatemp < 0) ? MPA : mpatemp;
	rmdmode = (rmdtemp < 0) ? RMDCHOICE : rmdtemp;
	gabudget = (gabudget < 0 || gabudget > 1) ? GABUDGETPROP : gabudget;
	metamodel = (metamodel == 1) ? metamodel : METAMODEL; */
	inputFile >> temp >> dimension;
	inputFile >> temp >> size;
	for (int i = 0; i < dimension; i++)
	{
		// First read starting point
		int input;
		inputFile >> input;
		x0.push_back(input);
	}
	for (int i = 0; i < size; i++)
	{
		double* Ai = new double[dimension];
		for (int j = 0; j < dimension; j++)
		{
			double input;
			inputFile >> input;
			Ai[j] = input;
		}
		A.push_back(Ai);
		double input;
		inputFile >> input;
		b.push_back(input);
	}
	return iscparam;
}

void ISCParameters::set_backtracking(int svar) {_backtracking = (svar < 0) ? BACKTRACKING : static_cast<bool>(svar);}
	
void ISCParameters::set_ocba(int svar) {_ocba = (svar < 0) ? OCBA : static_cast<bool>(svar);}

void ISCParameters::set_dominantniche(int svar) {_dominantniche = (svar < 0) ? DOMINANTNICHE : static_cast<bool>(svar);}

void ISCParameters::set_cleanup(int svar) {_cleanup = (svar < 0) ? CLEANUP : static_cast<bool>(svar);}

void ISCParameters::set_stocsim(int svar) {_stocsim = (svar < 0) ? STOCSIM : static_cast<bool>(svar);}

void ISCParameters::set_elitism(int svar) {_elitism = (svar < 0) ? ELITISM : static_cast<bool>(svar);}

ISCParameters::ISCParameters(double cleanupdelta, double cleanupalpha, double localoptdelta, double localoptalpha, double backdelta, double backalpha, \
		double globaldelta, long budget, double gabudget, int reps,	int gagen, int initialNumReps, int numOfCandidates, int pruning_freq, \
		int mpamode, int rmdmode, int metamodel, int metaprednum, int backtracking, int ocba, int dominantniche, int cleanup, int stocsim, int elitism) \
		: _cleanupdelta(cleanupdelta)
	{
		if (cleanupdelta <= 0)
		{
			std::cout << "Error, cleanup delta is negative\n";
			exit(5);
		}
		set_cleanupalpha(cleanupalpha);
		set_localoptdelta(localoptdelta);
		set_localoptalpha(localoptalpha);
		set_backdelta(backdelta);
		set_backalpha(backalpha);
		set_globaldelta(globaldelta);
		set_gabudget(gabudget);
		set_budget(budget);
		set_reps(reps);
		set_gagen(gagen);
		// If it is a stochastic sim problem, require at least 2 initial reps
		if (initialNumReps < 2 && stocsim) 
			set_initialNumReps(2);
		else
			set_initialNumReps(initialNumReps);
		set_numOfCandidates(numOfCandidates);
		set_pruning_freq(pruning_freq);
		set_mpamode(mpamode);
		set_rmdmode(rmdmode);
		set_metamodel(metamodel);
		_metamodelbak = metamodel;
		set_metaPredNum(metaprednum);

		set_backtracking(backtracking);
		set_ocba(ocba);
		set_dominantniche(dominantniche);
		set_cleanup(cleanup);
		set_stocsim(stocsim);
		set_elitism(elitism);
		
		set_globalalpha(-1); //USE default value
		set_alphadt(-1);
		set_alphac(-1);
	}
