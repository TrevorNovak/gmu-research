#ifndef NICHINGGAWTR_H_
#define NICHINGGAWTR_H_

// When GA is constructed, make sure that the vector of visited solution already contains the initial point from which
// RMD will sample other points as the initial population.
#include "NichePopulationWTR.h"
#include "masterheader.h"

template<class Genome>
class NichingGAWTR 
{
    private :
    	NichePopulationWTR<Genome>* _currentPopulation;
    	std::vector<Genome*>& _visited;
    	long _budget;
    	int _maxGeneration;
    	int _maxGenerationNoImprovement;
		unsigned int _lastBestGenome;
		unsigned int _bestGenome;
    	int _noOfGenerationNoImprovement;
    	int _populationSize, _genomeSize, _tournamentSize;
    	double _crossoverProbability, _mutationProbability;
    	double _alphadt;
        std::ofstream& _out; 
        std::vector<Record*>& _records;
        int _noOfEvaluations;
    	double _highestFitness, _oldHighestFitness;
    	double _nicheRadius;
    	// Dominant niche center
    	int _dominantNiche;
    	// Last replication count when an objective value is recorded.
    	long _lastRepCount;
    	double _lastObjectiveValue;
    	bool _dntest;
    	Simulator* _simulator;
    	bool _stocsim; 
    	bool _elitism;
    	std::vector<typename NichePopulationWTR<Genome>::Niche*> _niches; 
    	// Perfrom one iteration of the evolution process
    	void evolve();
    	bool terminate();
    	bool nichebasedTR();
		void writeStatistics();
		double evaluateTrueValue(Genome* genome)
		{
			return _simulator->simulation(genome);
		}
		void start(); 
    public :
    	// Constructor initializes maximum number of gnerations, maximum number of generations without improvement
	    // population size, mutation probability, the length of genome string, tournament size, and the GA log file.
	    // alphadat: dominant niche test
	    // alphag, deltag:  for grouping individuals into groups
    	/*NichingGAWTR(long budget, int maxG, int maxGNI, double alphadt, int pSize, double cProb, double mProb, int gSize, int tSize, \
   			double linearRankingCoeff, int inbreedingMF, double deltag, double alphag, \
	   		int minNumGroups, int initialNumRep, std::vector<Genome*>& visited, const std::vector<double*>& A, \
   			const std::vector<double>& b, std::ofstream& logFile, std::vector<Record*>& records, \
   			bool dntest, Simulator* simulator, bool stocsim, bool elitism) :  */
   		NichingGAWTR(long budget, int maxG, int pSize, double cProb, double mProb, int gSize, int tSize, double linearRankingCoeff, int inbreedingMF, \
			int minNumGroups, std::vector<Genome*>& visited, const std::vector<double*>& A, const std::vector<double>& b, \
			std::ofstream& logFile, std::vector<Record*>& records, Simulator* simulator, ISCParameters* iscparams) :  _visited(visited), _budget(iscparams->get_budget()), _maxGeneration(maxG), _maxGenerationNoImprovement(iscparams->get_gagen()), \
	    	_populationSize(pSize), _genomeSize(gSize), _tournamentSize(tSize), _lastBestGenome(0), \
		    _crossoverProbability(cProb), _mutationProbability(mProb), _alphadt(iscparams->get_alphadt()), _out(logFile), \
		    _records(records), _noOfEvaluations(0), _highestFitness(0), _oldHighestFitness(0), _simulator(simulator)
		{
			/* NichePopulationWTR(int pSize, double cProb, double mProb, int gSize, int tSize, int maxG, \
	    		std::vector<Genome*>& visited, const std::vector<double*>& A, \
		    	const std::vector<double>& b, const Genome& x0, const int crossoverType, const int mutationType, \
		    	const int selectionType, double lrc, int inbreedingMF, double delta, double alpha_g, \
	    		double delta_m, double alpha_m, int minNumGroups, int initialNumRep) 
		    */	
		    _dominantNiche = -1;
		    _stocsim = iscparams->get_stocsim();
		    _elitism = iscparams->get_elitism();
		    _lastRepCount = 0;
			_dntest = iscparams->get_dominantniche();
		    _currentPopulation = new NichePopulationWTR<Solution>(pSize, cProb, mProb, gSize, tSize, maxG, visited, A, b, _simulator, \
				1, 2, 2, linearRankingCoeff, inbreedingMF, iscparams->get_globaldelta(), iscparams->get_globalalpha(), \
				minNumGroups, iscparams->get_initialNumReps(), iscparams->get_stocsim(), iscparams->get_elitism());
			_lastObjectiveValue = _simulator->gettruevalue(_visited[_currentPopulation->getBestGenome()]);
			writeStatistics();
		}
		~NichingGAWTR() {delete _currentPopulation;}
		// Start the GA iteration process 
    	long solve();
	    // Collect statistics during GA run. 
    	unsigned int getBestGenome() const
	    {
    		return _currentPopulation->getBestGenome();
	    }
	    const std::vector<typename NichePopulationWTR<Genome>::Niche*>& getNiches() const {return _currentPopulation->getNiches();}
	    const std::vector<uint>& getNonPeakNiches() const {return _currentPopulation->getNonPeakNiches();}
	    const int getGeneration() const {return _currentPopulation->getGeneration();}
	    const int getDominantNiche() const {return _dominantNiche;}
	    double getLastObjectiveValue() const {return _lastObjectiveValue;}
	    long getLastRepCount() const {return _lastRepCount;} 
};

template<class Genome> 
void NichingGAWTR<Genome>::start()  
{
	do{
		// Record the fitness of the best genome in the parent generation 
	    _oldHighestFitness = _highestFitness;
    	// Determine niche radius
	    _currentPopulation->idNicheRadius();
		// Construct dynamic niches by a greedy algorithm similar to Miller 1995 and Li 2002
	   	_currentPopulation->identifyNiches();
   		// Apply transit rules based on niche construction here, so that unnecessary further selection and crossover operators won't happen
	   	
		// ********************************************************************************************************
		// Disable the niche based transition rule for experiments with binary variables, experiment with binary variables 06/04/2018 
		// if (nichebasedTR()) break;
		// ********************************************************************************************************

   		// Apply fitness sharing procedure, in particular, use the number of individuals in a niche as the discounting factor
		_currentPopulation->fitnessSharing();   	
	   	// Select parents
   		std::vector<unsigned int> parentPool;
		_currentPopulation->selectParentPool(parentPool);
		std::vector<unsigned int> offspringIndex;
   		std::vector<Genome*> offspring;
	   	/* Crossover P(t) with dynamic inbreeding restrictions. i.e., the first parent is selected 
   		 * using SUS + linear ranking, and then randomly sample M individuals from the parent solution, 
	   	 * and select the best individual that belongs to the same niche as the first parent as the 
   		 * second parent; if there is no such individual, simply select the one that is closest to the 
	   	 * first parent as the second parent
   		 */
	   	_currentPopulation->crossover(parentPool, offspring);
		// Mutate
	   	_currentPopulation->mutate(offspring, offspringIndex);
	   	// Replace parents with offsprings and evaluate new generation
	    // This is the only place new reps are taken
    	_currentPopulation->getNewGeneration(offspringIndex);
	    _highestFitness = _visited[_currentPopulation->getBestGenome()]->getSampleMean();

		// Debug use only
		std::cout << "Generation " << getGeneration() << " " << *_visited[_currentPopulation->getBestGenome()] << std::endl;
		
		writeStatistics();
	} while (!terminate());  
}

template<class Genome> 
bool NichingGAWTR<Genome>::nichebasedTR()
{
	// Rule 1: If there is only one niche, switch to COMPASS
    _niches = _currentPopulation->getNiches();
    if (_niches.size() == 1) 
    {
		std::cout << "NGA only detects one niche, quit" << std::endl;
    	return true;
    }
    if (!_dntest)	return false;
    // Rule 4: Dominance test
    size_t N = _niches.size();
    double beta = pow(1 - _alphadt, 1 / double(N -1));
    double* SSE = new double[N];
    double* SST = new double[N];
    double* t = new double[N];
    int* n = new int[N];
    size_t* M = new size_t[N];
    int* nu = new int[N];
    double* muhat = new double[N];
    // muhat should be negated because the subset selection procedure is for maximization. So I use minusmuhat
    double* minusmuhat = new double[N];
    double* ybar = new double[N];
    double* w = new double[N*N];
    for(size_t i = 0; i < N; i++) M[i] = _niches[i]->getNicheCount();
   	for(size_t i = 0; i < N; i++)
   	{
   		SSE[i] = SST[i] = muhat[i] = n[i] = 0;
   		const std::vector<uint>& residents = _niches[i]->getResidents();
   		for (size_t j = 0; j < M[i]; j++) n[i] += _visited[residents[j]]->getNumOfObservations();
   		double sum = 0;
   		for (size_t j = 0; j < M[i]; j++) sum += _visited[residents[j]]->getSampleMean();
   		ybar[i] = sum / M[i];
   		for (size_t j = 0; j < M[i]; j++) muhat[i] += _visited[residents[j]]->getNumOfObservations() * _visited[residents[j]]->getSampleMean();
   		muhat[i] /= n[i];
   		minusmuhat[i] = -muhat[i];
   		nu[i] = int(n[i] - M[i]);
   		t[i] = tinv(beta, nu[i]);
   		for (size_t j = 0; j < M[i]; j++) SSE[i] += _visited[residents[j]]->getSumSquared() - \
   			_visited[residents[j]]->getNumOfObservations() * _visited[residents[j]]->getSampleMean() * _visited[residents[j]]->getSampleMean();
   		for (size_t j = 0; j < M[i]; j++) SST[i] += _visited[residents[j]]->getNumOfObservations() * _visited[residents[j]]->getSampleMean() \
   			* _visited[residents[j]]->getSampleMean();
   		SST[i] -= n[i] * ybar[i] * ybar[i];
   	}
   	for(size_t i = 0; i < N; i++)
   	{
   		for (size_t j = 0; j < N; j++)
   		{
   			w[i*N + j] = sqrt(t[i]*t[i]*SSE[i]/nu[i]/n[i] + t[j]*t[j]*SSE[j]/nu[j]/n[j]);
   		}
   	}
   	// Now do a subset selection, note this procedure is for maximization, so negate muhat.
   	std::vector<size_t> I;
   	for(size_t i = 0; i < N; i++)
   	{
   		bool addi = true;
   		for (size_t j = 0; j < N; j++)
   		{
   			if (i != j && minusmuhat[i] < minusmuhat[j] - w[i*N + j])
   			{
   				addi = false;
   				break;
   			}
   		}
   		if (addi) I.push_back(i);
   	}
   	delete []ybar;
    delete []SSE;
    delete []SST;
    delete []n;
    delete []nu;
    delete []muhat;
    delete []minusmuhat;
    delete []w;
    delete []M;
    delete []t;
    if (I.size() == 1) 
    {
    	_dominantNiche = I[0];
		std::cout << "Dominant niche found, quit NGA" << std::endl;
    	return true;
    }
    // If none of those rules stops GA, continue
    return false;
}
template<class Genome> 
bool NichingGAWTR<Genome>::terminate()
{
	// Rule 0: budget based, GA must exit when exceeding budget
    if (_budget < _currentPopulation->getNumOfReps()) 
    {
    	cout << "GA Budget exceeded" << endl;
    	return true;
    }
    // Rule 3: If no GA improvement in _maxGenerationNoImprovement generations, swithc to COMPASS
	// This is stochastic simulation, so we should check if it is still the same solution rather than looking
	// at the objective function. 
	_bestGenome = getBestGenome();
    if (_lastBestGenome != _bestGenome)
    {	// Remeber this is a minimization problem.
        // There is improvement, reset the counter
        _noOfGenerationNoImprovement = 0;
		_lastBestGenome = _bestGenome;
    }
    if (++_noOfGenerationNoImprovement >= _maxGenerationNoImprovement)
    {
        std::cout << _maxGenerationNoImprovement << " generations without improvement, quit" << std::endl;
        return true;
    }
    return false;	
}

template<class Genome> 
long NichingGAWTR<Genome>::solve()
{
	_currentPopulation->initialize();
	writeStatistics();
	_highestFitness = _currentPopulation->getBestFitness();
    _oldHighestFitness = 0;
    _noOfGenerationNoImprovement = 0;
    start();
    // Now _currentPopulation is the final population. 
    return _currentPopulation->getNumOfReps();
}

template<class Genome> 
void NichingGAWTR<Genome>::writeStatistics()
{
	// When I get a new replication count, use the last replication count and last objective value
	// to fill in records up to the new replication count. 
	for (int i = _lastRepCount + 1; i < _currentPopulation->getNumOfReps(); i++)
	{
		// Two records to keep, one is sample path of this rep, written into the logfile. The
		// other is for average of all reps, kept in records.
		_out << i << " ,  " << _lastObjectiveValue << "\n";
		// Note that records index = repcount - 1 since the index starts at 0 while rep count starts at 1
		// If records size is less than repcount, add new records, otherwise, add data to exiting
		// records.
		RecordResult(_records, i, _lastObjectiveValue);	
	}
	// After filling records and out with previous data, now process the result of this generation
	_lastObjectiveValue = _simulator->gettruevalue(_visited[_currentPopulation->getBestGenome()]);
	_out << _currentPopulation->getNumOfReps() << " , " << _lastObjectiveValue << "\n";
	_lastRepCount =  _currentPopulation->getNumOfReps();
	RecordResult(_records,  _lastRepCount, _lastObjectiveValue);		
}
#endif /*NICHINGGAWTR_H_*/
