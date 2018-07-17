#ifndef READINPUT_H_
#define READINPUT_H_

#include "masterheader.h"

/*bool ReadInput(std::ifstream& inputFile, std::string& output, std::vector<int>& x0, std::vector<double*>& A, std::vector<double>& b, \
		long& budget, int& reps, bool& backtracking, bool& ocba, int& gagen, bool& dominantniche, bool& cleanup, \
		bool& stocsim, double& globaldelta, double& backalpha, double& backdelta, double& localoptapha, \
		double& localoptdelta, double& cleanupalpha, double& cleanupdelta, int& initialNumReps, \
		int& numOfCandidates, bool& elitism, int& pruning_freq, int& mpamode, int& rmdmode, double& gabudget,  int& metamodel);*/

class ISCParameters
{
private:
	double _cleanupdelta;  // This parameter must be specified by user, no default, clean-up indifference zone parameter 
	double _cleanupalpha; // clean up test type I error prob
	double _localoptdelta; // local optimality test indifference parameter, default to cleanup delta
	double _localoptalpha; // local optimality test type I error prob
	double _backdelta; // Back tracking test in COMPASS indifference zone parameter
	double _backalpha; // Back tracking test in COMPASS type I error prob
	double _globaldelta; // The indifference zone parameter in the global NGA phase, default to cleanup delta
	double _globalalpha; // confidence parameter used in grouping procedure of the NGA 
	double _alphadt; // dominant niche test alpha parameter of the NGA
	double _alphac; // CI-level for generalized COMPASS constraint placement procedure
	long _budget; // Total simulation budget in # of reps
	double _gabudget;  // percentage of budget spent on GA
	int _reps; // # of macro reps
	int _gagen; // numer of generations without improvement for NGA to quit
	int _initialNumReps; // initial number of reps assigned to a new solution
	int _numOfCandidates;  // size of S_k, i.e. # of samplings within MPR and NGA (repetion ok)
	int _pruning_freq; // COMPASS constraint pruning frequency
	int _mpamode; // MPA gemoetry, 0 for COMPASS and 1 for AHA
	int _rmdmode; // sampling distribution, 0 for uniform and 1 for coordinate
	int _metamodel; // 1 to use metamodel prediction to rank, 2 to use expected improvement to rank, 3 to use OO, and 0 to switch off
	int _metamodelbak; // Back up original model choice
	// For OO, _numOfCandidates is the top good solutions acceptable, ie, if it is 5, it means we are happy with inding top 5 solutions
	int _metaPredNum; // Number of prediction points for a metamodel, for OO, it is the alignment level, i.e., the number of top solutions to include

	bool _backtracking; // whether to perform the backtracking test for COMPASS
	bool _ocba; // whether to use OCBA in COMPASS 
	bool _dominantniche; // the dominant niche transition test for NGA
	bool _cleanup; // to perform final cleanup of multiple local optima found
	bool _stocsim; // indicate if the problem is stochastic or not
	bool _elitism; // use elitism in NGA step or not
	
	static const int CLEANUPALPHA = 5; //static const cannot be doulbe and thus we use the percents here to initialize these 4 parameters
	static const int LOCALOPTALPHA = 5; 
	static const int BACKALPHA = 10;
	static const int GLOBALALPHA = 10;
	static const int ALPHADT = 10;
	static const int ALPHAC = 10;
	static const int GABUDGETPROP = 18;

	static const long BUDGET = 100000;
	static const int REPS = 5; 
	static const int GAGEN = 3; 
	static const int INITIALNUMREPS = 5; 
	static const int NUMOFCANDIDATES = 5; 
	static const int PRUNINGFREQ = 5; 
	static const int MPA =0; 
	static const int RMDCHOICE = 0; 
	static const int METAMODEL = 0;
	static const int METAPREDNUM = 10;
	static const int OOAPCNT = 95; // Alignment probability as percents

	static const bool BACKTRACKING = false; 
	static const bool OCBA = true; 
	static const bool DOMINANTNICHE = false; 
	static const bool CLEANUP = true; 
	static const bool STOCSIM = true; 
	static const bool ELITISM = true;
	
public:
	ISCParameters(double cleanupdelta, double cleanupalpha, double localoptdelta, double localoptalpha, double backdelta, double backalpha, \
		double globaldelta, long budget, double gabudget, int reps,	int gagen, int initialNumReps, int numOfCandidates, int pruning_freq, \
		int mpamode, int rmdmode, int metamodel, int metaprednum, int backtracking, int ocba, int dominantniche, int cleanup, int stocsim, int elitism);
	/*ISCParameters(double cleanupdelta, double cleanupalpha, double localoptdelta, double localoptalpha, double backdelta, double backalpha, \
		double globaldelta, long budget, double gabudget, int reps,	int gagen, int initialNumReps, int numOfCandidates, int pruning_freq, \
		int mpamode, int rmdmode, int metamodel, int metaprednum, bool backtracking, bool ocba, bool dominantniche, bool cleanup, bool stocsim, bool elitism) \
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
	}*/
	
	void set_cleanupdelta(double svar) {_cleanupdelta = (svar < 0) ? _cleanupdelta : svar;}
	double get_cleanupdelta() const {return _cleanupdelta;}
	void set_cleanupalpha(double svar) {_cleanupalpha = (svar < 0) ? static_cast<double>(CLEANUPALPHA)/100.0 : svar;}
	double get_cleanupalpha() const {return _cleanupalpha;}
	void set_localoptdelta(double svar) {_localoptdelta = (svar < 0) ? _cleanupdelta : svar;}
	double get_localoptdelta() const {return _localoptdelta;}
	void set_localoptalpha(double svar) {_localoptalpha = (svar < 0) ? static_cast<double>(LOCALOPTALPHA)/100.0 : svar;}
	double get_localoptalpha() const {return _localoptalpha;}
	void set_backdelta(double svar) {_backdelta = (svar < 0) ? _cleanupdelta : svar;}
	double get_backdelta() const {return _backdelta;}
	void set_backalpha(double svar) {_backalpha = (svar < 0) ? static_cast<double>(BACKALPHA)/100.0 : svar;}
	double get_backalpha() const {return _backalpha;}
	void set_globaldelta(double svar) {_globaldelta = (svar < 0) ? _cleanupdelta : svar;}
	double get_globaldelta() const {return _globaldelta;}
	void set_gabudget(double svar) {_gabudget = (svar < 0) ? static_cast<double>(GABUDGETPROP/100.0) : svar;}
	double get_gabudget() const {return _gabudget;}

	void set_globalalpha(double svar) {_globalalpha = (svar < 0) ? static_cast<double>(GLOBALALPHA)/100.0 : svar;}
	double get_globalalpha() const {return _globalalpha;}
	void set_alphadt(double svar) {_alphadt = (svar < 0) ? static_cast<double>(ALPHADT)/100.0 : svar;}
	double get_alphadt() const {return _alphadt;}
	void set_alphac(double svar) {_alphac = (svar < 0) ? static_cast<double>(ALPHAC)/100.0 : svar;}
	double get_alphac() const {return _alphac;}

	void set_budget(long svar) {_budget = (svar < 0) ? BUDGET : svar;}
	long get_budget() const {return _budget;}
	void set_reps(int svar) {_reps = (svar < 0) ? REPS : svar;}
	int get_reps() const {return _reps;}
	void set_gagen(int svar) {_gagen = (svar < 0) ? GAGEN : svar;}
	int get_gagen() const {return _gagen;}
	void set_initialNumReps(int svar) {_initialNumReps = (svar < 0) ? INITIALNUMREPS : svar;}
	int get_initialNumReps() const {return _initialNumReps;}
	void set_numOfCandidates(int svar) {_numOfCandidates = (svar < 0) ? NUMOFCANDIDATES : svar;}
	int get_numOfCandidates() const {return _numOfCandidates;}
	void set_pruning_freq(int svar) {_pruning_freq = (svar < 0) ? PRUNINGFREQ : svar;}
	int get_pruning_freq() const {return _pruning_freq;}
	void set_mpamode(int svar) {_mpamode = (svar < 0) ? MPA : svar;}
	int get_mpamode() const {return _mpamode;}
	void set_rmdmode(int svar) {_rmdmode = (svar < 0) ? RMDCHOICE : svar;}
	int get_rmdmode() const {return _rmdmode;}
	void set_metamodel(int svar) {_metamodel = (svar < 0) ? METAMODEL : svar;}
	int get_metamodel() const {return _metamodel;}
	int get_metamodelbak() const {return _metamodelbak;}
	void set_metaPredNum(int svar) {_metaPredNum = (svar < 0) ? METAPREDNUM : svar;}
	int get_metaPredNum() const {return _metaPredNum;}

	void set_backtracking(int svar);
	bool get_backtracking() const {return _backtracking;}
	void set_ocba(int svar);
	bool get_ocba() const {return _ocba;}
	void set_dominantniche(int svar);
	bool get_dominantniche() const {return _dominantniche;}
	void set_cleanup(int svar);
	bool get_cleanup() const {return _cleanup;}
	void set_stocsim(int svar);
	bool get_stocsim() const {return _stocsim;}
	void set_elitism(int svar);
	bool get_elitism() const {return _elitism;}
};

ISCParameters* ReadInput(std::ifstream& inputFile, std::string& output, std::vector<int>& x0, std::vector<double*>& A, std::vector<double>& b);

#endif /*READINPUT_H_*/
