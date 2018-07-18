#ifndef NICHEPOPULATIONWTR_H_
#define NICHEPOPULATIONWTR_H_
#include "masterheader.h"

extern RngStream rng;
typedef unsigned int uint;
typedef std::vector<uint>:: iterator vuit;
using std::cout;
using std::endl;
template<class Genome> 

class NichePopulationWTR
{
    private :
    	// const for simple crossover operator
    	//static const double RO = 10;
    	// const for heurstic crossover operator
    	static const int W = 10;
    	// const for non-uniform mutation
    	//static const double B = 1.5;  
    	int _populationSize;
    	double _crossoverProbability;
    	double _mutationProbability;
    	int _genomeSize;
    	int _tournamentSize;
    	int _generation;
    	int _maxGeneration;
    	int _numberOfDistinctGenomes; 
    	// There are 3 types of crossover, 1-arithmetical, 2-simple, 3-heuristic
    	int _crossoverType;
    	// There are 2 types of mutation, 1-uniform, 2-non-uniform
    	int _mutationType;
    	int _selectionType;
    	double _nicheRadius;
    	double _linearRankingCoefficient;
    	int _inbreedingMF;
    	int _maxNumNiches;
    	// The indifference zone parameter delta
    	double _delta_g;
    	// Alpha quantile for grouping
    	double _alpha_g;
    	// Minimum number of groups to get in a generation for GA in the presence of noise
    	int _minNumGroups;
    	// Initial number of replications assigned to a newly visited individual
    	int _initialNumRep;
    	// Number of replications spent so far
    	long _numOfReps;
    	// Note that a global database of visited solution is passed into the population object.
    	// _genomes[i] is the index of the ith genome in _visited. All newly visited solutions are added to "_visited"
    	std::vector<uint> _genomes;
    	// _bestGenome has the index of the best genome in _visited.
    	uint _bestGenome;
    	// _distinctFlag indicates whether a genome is distinct. If the ith genome is distinct,
    	// _distinctFlag[i] is i, otherwise, it will be the index of the first occurence 
    	// of this genome in the population.
    	uint* _distinctFlag;
    	std::vector<Genome*>& _visited;
    	// Constraints are Ax >= b
    	const std::vector<double*>& _A;
    	const std::vector<double>& _b;
    	Simulator* _simulator; 
    	bool _stocsim;
    	bool _elitism;
    	// This method scans the population and finds out how many distinct genomes there
    	// are in the population. It rewrites the array _distintFlag to mark out which genome is unique and which
	    // is the same as another genome. 
    	void getDistinctGenomes();
    	// Evaluate the population and returns the index of the best 
    	uint evaluatePopulation();
    	// Evaluate each individual without noise for response surface problems.
    	void evaluateNoiseFree();
    	uint tournamentSelection();
    	// This method assigns each individual a selection prob according to linear ranking scheme. prob[i] has
		// the selection prob of the ith individual in _genomes. prob must be initialized to all 0.
    	void linearRanking(std::vector<double>& prob);
    	/* grouping procedure takes as input the index(as in _genomes, so they are from 0 to
 		 * _populationSize, NOT indexes in _visited!!!) of solutions to be grouped. Also as an input
		 * is the minimum number of groups to get from. toBeGrouped initially has all individuals inside it.
		 * Then each invocation of grouping() partition individuals in toBeGrouped into groups, and if one group needs
		 * to be breaked up further, it recursively calls itself, indicating which group needs to be divided further
		 * by input parameter index, i.e., toBeGrouped[index] is the group to be divided further.
		 */
    	void grouping(std::vector<std::vector<uint> >& toBeGrouped, int index, int minNumGroup);
    	/* This method assigns each individual a selection probability adjusted for a stochastic environment based on 
    	 * its deterministic linear ranking prob and grouping result. std::vector<double>& prob first holds linear ranking
    	 * prob as input, then it holds grouping ranking prob as output. 
    	 */
    	void groupRanking(std::vector<double>& prob, std::vector<std::vector<uint> >& groups);
    	void SUSLinearRanking(std::vector<uint>& parentPool); // stochastic uniform sampling
		// Use a mating restriction procedure to find the other mate for the given parent1. 
		// At this stage, I implement inbreeding. 
    	uint matingRestriction(uint parent1);
    	uint inbreeding(uint parent1);
    	uint linebreeding(uint parent1);
    	void arithmeticalCrossover(uint parent1, uint parent2, Genome*& child1, Genome*& child2); 
    	void simpleCrossover(uint parent1, uint parent2, Genome*& child1, Genome*& child2);
    	void heuristicCrossover(uint parent1, uint parent2, Genome*& child);
    	void uniformMutate(Genome* offspring);
    	void nonuniformMutate(Genome* offspring);
    	// Below are methods to get new generation
    	void noElitism(const std::vector<uint>& offspring);
    	void conserveAllNicheCenters(const std::vector<uint>& offspring);
    	void selectiveElitism(const std::vector<uint>& offspring);
    	// This is the list of dominant individuals in a niche, or "species seeds"
    	std::vector<uint> _nonPeakNiche;
    	int delta(int generation, int maxGeneration, int xk, int bound)
    	{
    		double factor;
    		if (maxGeneration <= generation) 
    			factor = 0.005;
    		else
    			factor = pow(1 - generation / maxGeneration, 1.5); //B);
    		double temp = abs(xk - bound) * rng.RandU01() * factor;
    		int int_temp = int(ceil(temp));
    		if (bound <= xk && xk - int_temp < bound)
    			int_temp = xk - bound;
    		else if (xk <= bound && bound < xk + int_temp)
    			int_temp = bound - xk;
    		return int_temp;
    	}
    	bool isFeasible(Genome* genome)
    	{	// Check if the given genome is feasible or not, Remember the constraints are Ax >= b
    		for (uint i = 0; i < _A.size(); i++)
    		{
    			double sum = 0;
    			for (int j = 0; j < _genomeSize; j++)
    			{
    				sum += *(genome->getDiscreteSolution() + j) * (*(_A[i] + j));
    			}
    			if (sum < _b[i])
    				return false;
    		}
    		return true;
    	}
    	uint findMin(double* fitnessScores, size_t size)
    	{
    		double min = fitnessScores[0];
    		uint index = 0;
    		for (uint i = 0; i < size; i++)
    		{
    			if (fitnessScores[i] < min)
    			{
    				min = fitnessScores[i];
    				index = i;
    			}
    		}
    		return index;
    	}
    	void fathom(const std::vector<double*>& A, const std::vector<double>& b, const Genome* const genome, int allele, int& int_lower, int& int_upper);
    	double distance(uint i, uint j) const;
    public : 
    	// Initialize a population with populationSize many genomes'
	    NichePopulationWTR(int pSize, double cProb, double mProb, int gSize, int tSize, int maxG, \
	    	std::vector<Genome*>& visited, const std::vector<double*>& A, const std::vector<double>& b, \
	    	Simulator* simulator, const int crossoverType, const int mutationType, \
	    	const int selectionType, double lrc, int inbreedingMF, double delta_g, double alpha_g, \
	    	int minNumGroups, int initialNumRep, bool stocsim, bool elitism) : \
	    	
	    	_populationSize(pSize), _crossoverProbability(cProb), _mutationProbability(mProb), \
			_genomeSize(gSize), _tournamentSize(tSize), _generation(1), _maxGeneration(maxG), \
			_numberOfDistinctGenomes(0), _crossoverType(crossoverType), _mutationType(mutationType), \
			_selectionType(selectionType), _linearRankingCoefficient(lrc), _inbreedingMF(inbreedingMF), \
			_delta_g(delta_g), _alpha_g(alpha_g), _minNumGroups(minNumGroups), \
			_bestGenome(0), _visited(visited), _A(A), _b(b), _simulator(simulator), _stocsim(stocsim),
			_elitism(elitism)
	    {
	    	_numOfReps = 0;
	    	_initialNumRep = (_stocsim) ? initialNumRep : 1;
	    	_distinctFlag = new uint[_populationSize];
	    	if ( _tournamentSize % 2 == 1)
	    	{ // tournament size must not be odd
	    		_tournamentSize++;
	    	}
	    	else if (_tournamentSize == 0)
	    	{
	    		_tournamentSize = 2;
	    	}
	    	// Evaluate the starting point so that I can use it as the best objective value
			// when writing statistics in NichingGAWTR.h
			int additionalrep = _initialNumRep - _visited[0]->getNumOfObservations();
			// call simulator in a batch mode 
			vector<size_t> numrep(1,additionalrep);
			vector<size_t> solnid(1,0);
			_simulator->simulation(_visited, solnid, numrep);
			/*for (int rep = 0; rep < additionalrep; rep++)
			{
				_visited[0]->recordObservation(_simulator->simulation(_visited[0]));
			}
			*/
			_numOfReps += (0 < additionalrep) ? additionalrep : 0;
    		_bestGenome = 0;
    		// If this is a deterministic simulator, use only one rep
    	}
    	
	    ~NichePopulationWTR() 
	    {
	    	delete[] _distinctFlag;
	    	for (uint i = 0; i < _niches.size(); i++)
			{	
				delete _niches[i];
			}
	    	_distinctFlag = 0;
	    }
	    class Niche
		{
			private :
				// center is the index of the center in the vector _visited
				uint center;
				// residents has the indexes of individuals in this niche in _visited
				std::vector<uint> residents;
				uint nicheCount;
			public :
			 	Niche(uint center) : center(center), nicheCount(1) {residents.push_back(center);}
			 	~Niche() {residents.clear();}
			 	void add(uint index) {residents.push_back(index); nicheCount++;}
			 	const uint getCenter() const {return center;}
			 	const int getNicheCount() const {return nicheCount;}
			 	const std::vector<uint>& getResidents() const {return residents;}
			 	// Note that r1 and r2 here are the indexes in _visited!!!
			 	const bool sameNiche(uint r1, uint r2) const
			 	{
			 		bool inside1 = false, inside2 = false;
			 		for (uint i = 0; i < residents.size(); i++)
			 		{
			 			if (residents[i] == r1)
			 			{
			 				inside1 = true;
			 				break;
			 			}
			 		}
			 		for (uint i = 0; i < residents.size(); i++)
			 		{
			 			if (residents[i] == r2)
			 			{
			 				inside2 = true;
			 				break;
			 			}
			 		}
			 		return (inside1 && inside2);
			 		//return true;
			 	}
		};
		std::vector<Niche*> _niches;
   		// Randomly generate an intial population and evaluate it
   		void initialize();
   		// Use nearest neighbord info to identify local mins, 
   		void idNicheRadius();
   		// Identfiy niches/species seeds
   		void identifyNiches();
   		void dynamicNiche();
   		// Apply fitness sharing 
   		void fitnessSharing();
   		void dynamicNicheFitnessSharing();
	    // Select a candidate pool for mating. Tournament selection, size given by _tournamentSize.
	    // candidatePool is given as an empty vector when this method is called, and upon returning, 
	    // it has the indexes of genomes to perform crossover inside.
    	void selectParentPool(std::vector<uint>& parentPool);
        // Perform crossover over the selected candidate pool and return new generation.
        // offspring is an empty vector when this method is called, it then has the offspring upon finishing.
    	void crossover(const std::vector<uint>& parentPool, std::vector<Genome*>& offspring); 
	    void mutate(std::vector<Genome*>& offspring, std::vector<uint>& offspringIndex);
        // Replace the old generation with new offsprings
    	void getNewGeneration(const std::vector<uint>& offspring);
    	int getGeneration() const {return _generation;}
        int getNumberOfDistinctGenomes() const {return _numberOfDistinctGenomes;}
        uint getBestGenome() const {return _bestGenome;}
        double getBestFitness() const {return _visited[_bestGenome]->getSampleMean();}
        long getNumOfReps() const {return _numOfReps;}
        const std::vector<Niche*>& getNiches() const {return _niches;}
        const std::vector<uint>& getNonPeakNiches() const {return _nonPeakNiche;}
        const uint* const getDistinctFlags() const {return _distinctFlag;}
        void print() const;
        void print(const std::vector<uint>& genomes) const;
        void print(const std::vector<uint>& genomes, bool) const;
};

template<class Genome> 
void NichePopulationWTR<Genome>::getDistinctGenomes()
{
    std::vector<uint> distinctGenomes;
    for (uint i = 0; i < _genomes.size(); i++)
    {
    	bool isDistinct = true;
        for (uint j = 0; j < distinctGenomes.size(); j++)
        {
        	if (_genomes[i] == _genomes[distinctGenomes[j]])
            {
                isDistinct = false;
                _distinctFlag[i] = distinctGenomes[j];
                break;
            }
        }
    	if (isDistinct)
        {
            distinctGenomes.push_back(i);
            _distinctFlag[i] = i;
        }
    }
    _numberOfDistinctGenomes = int(distinctGenomes.size());
}

template<class Genome> 
void NichePopulationWTR<Genome>::initialize()
{
	// Randomly generate initial population. This is not a trivial job. Need to consider constraints.
	// Call RMD for compass to get this done.
	std::vector<Genome*> candidates;
	RMD(_A, _b, _visited[0], _populationSize-1, 1000, candidates);
	// In visited, the starting point x0 is already pushed in and has index 0. So _genomes has 0 in it.
	_genomes.push_back(0);
	for (uint i = 0; i < candidates.size(); i++)
	{
		_genomes.push_back(AddSolution(_visited, candidates[i]));
	}
	// Evaluate the initial population
   	_bestGenome = evaluatePopulation();
}

template<class Genome>     
uint NichePopulationWTR<Genome>::evaluatePopulation()
{
	// When evaluting a genome, first make sure it is distinct. Otherwise,
    // simply fetch the fitness score from the same genome already evaluated
    getDistinctGenomes();
    double* fitnessScores = new double[_genomes.size()];
	vector<size_t> numrep;
	vector<size_t> solnid;
	for (uint i = 0; i < _genomes.size(); i++)
    {
    	if (_distinctFlag[i] == i)
        {
       		// Here this is an arbitrary choice of initial number of replications
	        int additionalrep = _initialNumRep - _visited[_genomes[i]]->getNumOfObservations();
	        additionalrep = (additionalrep > 0) ? additionalrep : 0;
	        // Implement a log10(iteration count) rule to determine additional rep. So we maintain
	        // global convergence property of NGA without spending too many iterations there.
	        if (_stocsim) additionalrep += static_cast<int>(log10(static_cast<double>(_generation)));
			numrep.push_back(additionalrep);
			solnid.push_back(_genomes[i]);
	        /*for (int rep = 0; rep < additionalrep; rep++)
    		{
    			_visited[_genomes[i]]->recordObservation(_simulator->simulation(_visited[_genomes[i]]));	
    		}*/
    		_numOfReps += (0 < additionalrep) ? additionalrep : 0;
    		fitnessScores[i] = _visited[_genomes[i]]->getSampleMean();
        }
        else
        {
            // The fitness is already known. However, the fitness of the genome still needs be set.
            fitnessScores[i] = fitnessScores[_distinctFlag[i]];
        }
    }
	_simulator->simulation(_visited, solnid, numrep);
	for (uint i = 0; i < _genomes.size(); i++)
	{
		if (_distinctFlag[i] == i)
		{
			fitnessScores[i] = _visited[_genomes[i]]->getSampleMean();
		}
		else
		{
			// The fitness is already known. However, the fitness of the genome still needs be set.
			fitnessScores[i] = fitnessScores[_distinctFlag[i]];
		}
	}
    size_t min = findMin(fitnessScores, _genomes.size());
    delete []fitnessScores;
    return _genomes[min];
}


template<class Genome>   
void NichePopulationWTR<Genome>::idNicheRadius()
{
	// First do a sort in decreasing order of fitness. Note this is a minimization problem.
	// So individual with largest objective value is at the end
	uint* rank = new uint[_populationSize];
    for (int i = 0; i < _populationSize; i++)
    	rank[i] = i;
    for (int i = _populationSize - 1; 0 < i; i--)
    {
    	for (int j = 0; j < i; j++)
    	{
    		if (_visited[_genomes[rank[j]]]->getSampleMean() > _visited[_genomes[rank[j + 1]]]->getSampleMean())
    		{
    			unsigned int temp = rank[j];
    			rank[j] = rank[j + 1];
    			rank[j + 1] = temp;
    		}
    	}
    }	
    //List I and vector L have the indexes as in _visited
   	std::list<uint> I;
   	std::vector<uint> L;
   	for (int i = 0; i < _populationSize; i++)
    {
    	// If there are duplicate copies of a genome in the population, igore those duplicates.
   		if (_distinctFlag[rank[i]] == rank[i])
	    	I.push_back(_genomes[rank[i]]);
    }
    std::list<uint> oldI(I.begin(), I.end());
	int counter = 0;
	while(0 < I.size())
    { 
    	// Ak, bk, activeConstraints are the returned values by this function.
		/* When this function is called, activeConstraints has the indices of those solutions in "visited" 
   		 * that define MPR in the previous iteration. When the function call finishes, activeConstraints
	   	 * contains active constraints for this iteration, which will be used when RMD is called to sample
   		 * new solutions 
   		*/ 
	   	/* To apply this function here, activeSolutions contains all other points, starIndex and oldStarIndex
   		 * both set to the index of the solution of interest, A, b are original problem formulation constraints
	   	 * bestSolutions, Ak, bk are empty vectors. Note that activeSolutions has indexes in _visited
   		 */ 
    	std::vector<double*> Ak;
   	 	Ak.clear();
   	 	std::vector<double> bk;
   	 	bk.clear();
   	 	std::vector<long> bestSolutions;
   	 	bestSolutions.clear();
   	 	std::vector<long> activeSolutions;
   	 	for (uint i = 0; i < (uint) _populationSize; i++)
		{
			if (_genomes[i] != I.front() && _distinctFlag[i] == i)
			{	// Note activeSolutions has indexes in _visited
				activeSolutions.push_back(_genomes[i]);
			}
		}		
		ConstructMPR(_visited, I.front(), I.front(), Ak, bk, _A, _b, activeSolutions, bestSolutions);
		// Debug use only
		//std::cout << _visited[*I.begin()];
		//PrintActiveSolutions(_visited, activeSolutions);
		// Debug use only

		// Now activeSolutions has the set Aj, those that forms active COMPASS constraints for Xi
		bool isLocalMin = true;
		for (uint i = 0; i < activeSolutions.size(); i++)
		{	// The stopping test procedure might be used here to get some statistical gurantee for this comparison
			if (_visited[*I.begin()]->getSampleMean() >= _visited[activeSolutions[i]]->getSampleMean())
			{
				isLocalMin = false;
				break;
			}  
		}
		// The first individual in list I has the best fitness, though there may be a tie. Make it a local min automatically. 
		if (counter++ < 1) isLocalMin = true;
		if (isLocalMin)
		{	// L has the indexes of local mins in _visited
			L.push_back(I.front());
			// Now exclude the neighbors of this localmin from I, which are listed in activesolutions. Because
			// they cannot be local mins
			for (uint i = 0; i < activeSolutions.size(); i++)
			{
				for (std::list<uint>::iterator it = I.begin(); it != I.end(); it++)
				{
					if (long(*it) == activeSolutions[i])
					{
						I.erase(it);
						break;
					}
				}
			}
		}
		I.erase(I.begin());
    } 
	// Determine Niche Radius
    double R = 1e50; //std::numeric_limits<double>::max();
	for (uint i = 0; i < L.size() - 1; i++)
	{
    	for (uint j = i + 1; j < L.size(); j++)
    	{
    		if (distance(L[i], L[j]) < R)
    		{
    			R = distance(L[i], L[j]);
    		}
    	}
	}
	_nicheRadius = R / 2;
    _maxNumNiches = int(L.size()); 
    delete []rank;
}

/* This method scans the population and identifies niche peaks in a greedy manner. It clears the
 * vector _niches and put new niche centers into it 
 */
template<class Genome>   
void NichePopulationWTR<Genome>::identifyNiches()
{
	dynamicNiche();
}

template<class Genome>   
void NichePopulationWTR<Genome>::dynamicNiche()
{
	// First do a sort in decreasing order of fitness. Note this is a minimization problem.
	uint* rank = new uint[_populationSize];
    for (int i = 0; i < _populationSize; i++)
    	rank[i] = i;
    for (int i = _populationSize - 1; 0 < i; i--)
    	for (int j = 0; j < i; j++)
    		if (_visited[_genomes[rank[j]]]->getSampleMean() > _visited[_genomes[rank[j + 1]]]->getSampleMean())
    		{
    			uint temp = rank[j];
    			rank[j] = rank[j + 1];
    			rank[j + 1] = temp;
    		}
    for (uint i = 0; i < _niches.size(); i++)
		delete _niches[i];
	_niches.clear();
	_nonPeakNiche.clear();
	int count = 0, numPeaks = 0;
	while(numPeaks < _maxNumNiches && count < _populationSize)
	{
		bool nicheFound = false;
		for (uint i = 0; i < _niches.size(); i++)
		{	// it's the index in _visited that is recorded in _niches. Need be consistent
			double dist = distance(_niches[i]->getCenter(), _genomes[rank[count]]);
			if (dist <= _nicheRadius)
			{
				nicheFound = true;
				// It's the index in _visited that is recorded in _niches. Need be consistent
				_niches[i]->add(_genomes[rank[count]]);
				break;
			}
		}
		if (!nicheFound)
		{	// Note that I push the index in _visited into _niches for the niche center
			_niches.push_back(new Niche(_genomes[rank[count]]));
			numPeaks++;
		}
		count++; 
	}
	// Add remaining individuals to niches, if they don't belong to any niche, leave them alone
	for (; count < _populationSize; count++)
	{
		bool nicheFound = false;
		for (uint i = 0; i < (uint) _maxNumNiches; i++)
		{	// it's the index in _visited that is recorded in _niches. Need be consistent
			double dist = distance(_niches[i]->getCenter(), _genomes[rank[count]]);
			if (dist <= _nicheRadius)
			{
				nicheFound = true;
				// It's the index in _visited that is recorded in _niches. Need be consistent
				_niches[i]->add(_genomes[rank[count]]);
				break;
			}
		}
		if (!nicheFound)
			_nonPeakNiche.push_back(_genomes[rank[count]]);
	}
	delete []rank;
}

// Apply fitness sharing procedure, in particular, use the number of individuals in a niche as the discounting factor
template<class Genome>   
void NichePopulationWTR<Genome>::fitnessSharing()
{
	dynamicNicheFitnessSharing();
}

template<class Genome>   
void NichePopulationWTR<Genome>::dynamicNicheFitnessSharing()
{
	// If an individual is within a niche, discount it by the size of the niche
	// otherwise, use standard sharing function (triangular)
	for (size_t i = 0; i < _niches.size(); i++)
	{
		const std::vector<uint> residents = _niches[i]->getResidents();
		// Remember that residents[i] is the index in _visited
		for (size_t j = 0; j < residents.size(); j++)
		{
			_visited[residents[j]]->setSharedFitness(_visited[residents[j]]->getSampleMean() / _niches[i]->getNicheCount());
			_visited[residents[j]]->setSharedSampleVar(_visited[residents[j]]->getSampleVariance() / _niches[i]->getNicheCount() / _niches[i]->getNicheCount());			
		}
		/* Do not change the fitness of the niche peak. It might be advisable to leave the fitness 
		 * value of a niche peak unchanged because it may happen that there are many singleton 
		 * niches that are located in unpromising regions, and thus those niche peaks will be 
		 * given artificially too much a chance to take part in recombination than those true 
		 * niche peaks around local optima. 
		 */
	}
	for (size_t i = 0; i < _nonPeakNiche.size(); i++)
	{
		double mi = 1;
		for (size_t j = 0; j < _nonPeakNiche.size(); j++)
			if (distance(_nonPeakNiche[i], _nonPeakNiche[j]) <= _nicheRadius)
				mi += 1 - distance(_nonPeakNiche[i], _nonPeakNiche[j]) / _nicheRadius;
		_visited[_nonPeakNiche[i]]->setSharedFitness(_visited[ _nonPeakNiche[i]]->getSampleMean() / mi);
		_visited[_nonPeakNiche[i]]->setSharedSampleVar(_visited[ _nonPeakNiche[i]]->getSampleVariance() / mi / mi);
	}
}

template<class Genome>   
uint NichePopulationWTR<Genome>::tournamentSelection()
{	// This is tournament selection with continuous updating
	uint* rivals = new uint[_tournamentSize];
    double* fitnessScores = new double[_tournamentSize];
    for (int j = 0; j < _tournamentSize; j++)
    {
    	// Randomly pick up a genome
        rivals[j] = rng.RandInt(0, _populationSize - 1);
        // May need to ensure rivals are unique
        fitnessScores[j] = _visited[_genomes[rivals[j]]]->getSharedFitness();
    } 
    // Select the winner of this tournament 
    uint winner = findMin(fitnessScores, _tournamentSize);
    delete []rivals;
    delete []fitnessScores;
	return _genomes[rivals[winner]];
}

// This method assigns each individual a selection prob according to linear ranking scheme. prob[i] has
// the selection prob of the ith individual in _genomes. prob must be initialized to all 0.
template<class Genome>   
void NichePopulationWTR<Genome>::linearRanking(std::vector<double>& prob)
{
	// First, rank genomes. After the procedure is done, rank has the indexes of _genomes[]
    // such that rank[0] gives the index of the best individual in _genomes[] and rank[_populationSize-1] the worst
    // So if rank[j] = i, it means that the ith genome in _genomes has a rank of (j+1) in the population
    // Note that the index starts from 0.
    uint* rank = new uint[_populationSize];
    for (int i = 0; i < _populationSize; i++)
    	rank[i] = i;
    for (int i = _populationSize - 1; 0 < i; i--)
    {
    	for (int j = 0; j < i; j++)
    	{	// remember that we are working on a minimization problem
    		if (_visited[_genomes[rank[j + 1]]]->getSharedFitness() < _visited[_genomes[rank[j]]]->getSharedFitness())
    		{
    			uint temp = rank[j];
    			rank[j] = rank[j + 1];
    			rank[j + 1] = temp;
    		}
    	}
    }
    for (int j = 0; j < _populationSize; j++)
    {	// Assign selection probabilites of the jth ranked genome
		prob[rank[j]] = (_linearRankingCoefficient - 2*(_linearRankingCoefficient - 1)* j/(_populationSize - 1)) / _populationSize;
    }
	delete []rank;   
}

/* This method assigns each individual a selection probability adjusted for a stochastic environment based on 
 * its deterministic linear ranking prob and grouping result. std::vector<double>& prob first holds linear ranking
 * prob as input, then it holds grouping ranking prob as output. Notice that groups has indexes as in _genomes, not
 * _visited!!!
 */
template<class Genome>  
void NichePopulationWTR<Genome>::groupRanking(std::vector<double>& prob, std::vector<std::vector<uint> >& groups)
{
	size_t groupsize = groups.size();
	double* gprob = new double[groupsize];
	// First calculate each group's average selection probs
	for (size_t i = 0; i < groupsize; i++)
	{
		double sum = 0;
		for (size_t j = 0; j < groups[i].size(); j++)
		{
			sum += prob[groups[i].at(j)];
		}
		gprob[i] = sum / groups[i].size();
	}
	for (size_t i = 0; i < groupsize; i++)
	{
		for (size_t j = 0; j < groups[i].size(); j++)
		{
			prob[groups[i].at(j)] = gprob[i];
		}
	}
	delete []gprob;
}
    	
/* grouping procedure takes as input the index(as in _genomes, so they are from 0 to
 * _populationSize, NOT indexes in _visited!!!) of solutions to be grouped. Also as an input
 * is the minimum number of groups to get from. toBeGrouped initially has all individuals inside it.
 * Then each invocation of grouping() partition individuals in toBeGrouped into groups, and if one group needs
 * to be breaked up further, it recursively calls itself, indicating which group needs to be divided further
 * by input parameter index, i.e., toBeGrouped[index] is the group to be divided further.
 */
template<class Genome>   
void NichePopulationWTR<Genome>::grouping(std::vector<std::vector<uint> >& groups, int index, int minNumGroup)
{
	// 1. Given current individuals to group in toBeGrouped, with sample means, variances, and # of observations
	// 2. Sort so that individual 1's fitness > individual 2's fitness, etc.. 
	uint m = uint(groups[index].size());
	// Make a local copy of the vector that has the indexes of individuals in _genomes that needs to be grouped
	std::vector<uint> toBeGrouped(groups[index]);
	groups[index].clear();
	// Sort individuals in toBeGrouped in the decreasing order of objective values, so inferior systems are at the front
	for (uint i = m - 1; 0 < i; i--)
    {
    	for (uint j = 0; j < i; j++)
    	{	// remember that we are working on a minimization problem
    		// Note that toBeGrouped has indexes in _genomes
    		if (_visited[_genomes[toBeGrouped[j + 1]]]->getSharedFitness() > _visited[_genomes[toBeGrouped[j]]]->getSharedFitness())
    		{
    			uint temp = toBeGrouped[j];
    			toBeGrouped[j] = toBeGrouped[j + 1];
    			toBeGrouped[j + 1] = temp;
    		}
    	}
    }
    // 3. Compute average variance and distance to divide group R
    double s2 = 0;
    int nbar = 0, nsum = 0;
    for (uint i = 0; i < m; i++)
    {
    	s2 += _visited[_genomes[toBeGrouped[i]]]->getSharedSampleVar();
    	nsum += _visited[_genomes[toBeGrouped[i]]]->getNumOfObservations();
    }
    s2 /= m;
    nbar = nsum / m;
    // qtrng(p, v, r) Approximates the quantile p for a studentized range distribution
 	// having v degrees of freedom and r samples for probability 0.9 < p < 0.99.
    double R = qtrng(1 - _alpha_g, nsum - m, m) * sqrt(double(s2)) / sqrt(double(nbar));
    // 4. do the grouping
    // g is the index of that group in the vector groups. So the first subgroup divided from group[index] will
    // be held at the same place, and subsequent groups will be pushed to the back end of vector groups.
    uint g = index;
    for(uint i = 0; i < m;)
    {
    	// If this is not the first group that uses the existing vector, need to do initilization
    	if (g >= groups.size())
    	{
    		std::vector<uint> newgroup;
    		// Note that standard containers hold objects by value, so it's safe to define a local vector
    		// newgroup here and push it into groups, and then simply work on the newly added vector in groups
    		groups.push_back(newgroup);
    		g = uint(groups.size() - 1);
    	}
    	groups[g].push_back(toBeGrouped[i]);
    	uint bottom = i++;
    	while (i < m && (_visited[_genomes[toBeGrouped[bottom]]]->getSharedFitness() - _visited[_genomes[toBeGrouped[i]]]->getSharedFitness() < R))
    	{
    		if (_visited[_genomes[toBeGrouped[i]]]->getSharedFitness() > _visited[_genomes[toBeGrouped[bottom]]]->getSharedFitness()) 
    		{
    			cout << _visited[_genomes[toBeGrouped[i]]];	
    			cout << _visited[_genomes[toBeGrouped[i]]]->getSharedFitness() << endl;
    			cout << _visited[_genomes[toBeGrouped[bottom]]];
    			cout << _visited[_genomes[toBeGrouped[bottom]]]->getSharedFitness() << endl;
    			
    		}
    		require(_visited[_genomes[toBeGrouped[i]]]->getSharedFitness() <= _visited[_genomes[toBeGrouped[bottom]]]->getSharedFitness() \
    			, "The group being grouped is not sorted in decreasing order of fitness values");
    		groups[g].push_back(toBeGrouped[i++]);
    	}
    	// Change the value of g so that the next iteration a new group will be added
   		g = uint(groups.size());
    }
    // 5. We now have g groups group[0], ... group[g-1]. See if g is large enough 
    if (groups.size() < (uint) minNumGroup)
    {	// If there are too few groups, do grouping again on a selected group with largest gap
		double maxRHat = dabs(_visited[_genomes[*(groups[0].begin())]]->getSharedFitness() - _visited[_genomes[*(--(groups[0].end()))]]->getSharedFitness());	
    	uint lhat = 0;
    	for (uint j = 1; j < groups.size(); j++)
    	{
    		if (maxRHat < dabs(_visited[_genomes[*(groups[j].begin())]]->getSharedFitness() - _visited[_genomes[*(--(groups[j].end()))]]->getSharedFitness()))
    		{
    			maxRHat = dabs(_visited[_genomes[*(groups[j].begin())]]->getSharedFitness() - _visited[_genomes[*(--(groups[j].end()))]]->getSharedFitness());
    			lhat = j;
    		}
    	}
    	if (_delta_g < maxRHat) 
    	{	//If the difference within a group is large enough to be detected, continue
    		// Compute Q for this subgroup with max Rhat, So need to know # of systems in that group and 
	    	// the totalnumber of observations
    		uint mhatlhat = uint(groups[lhat].size());
    		double s2lhat = 0;
	    	int nbarlhat = 0, nsumlhat = 0;
    		for (uint j = 0; j < mhatlhat; j++)
    		{
    			nsumlhat += _visited[_genomes[groups[lhat].at(j)]]->getNumOfObservations();
	    		s2lhat += _visited[_genomes[groups[lhat].at(j)]]->getSharedSampleVar();
    		}
    		s2lhat /= mhatlhat;
	    	nbarlhat = nsumlhat / mhatlhat;
    		double Q = qtrng(1 - _alpha_g, nsumlhat - mhatlhat, mhatlhat);
    		double nhat = ceil(Q * Q * s2lhat / maxRHat / maxRHat);
	    	// obtain max(nhat - ni more observations from i in groups[lhat]
			vector<size_t> numrep;
			vector<size_t> solnid;
			int additionalRep;
    		for (uint j = 0; j < mhatlhat; j++)
    		{
				if (_distinctFlag[groups[lhat].at(j)] == groups[lhat].at(j))
				{ 
    				additionalRep = int(nhat) - _visited[_genomes[groups[lhat].at(j)]]->getNumOfObservations();
					if (additionalRep > 0)
					{
						numrep.push_back(additionalRep);
						solnid.push_back(_genomes[groups[lhat].at(j)]);
						_numOfReps += additionalRep;
					}
				}
				/*	for (int rep = 0; rep < additionalRep; rep++)
    			{
    				_visited[_genomes[groups[lhat].at(j)]]->recordObservation(_simulator->simulation(_visited[_genomes[groups[lhat].at(j)]]));	
    			}*/
    		}
			_simulator->simulation(_visited, solnid, numrep);
			// Now fitness values may have changed, need to do fitness sharing again.
    		fitnessSharing();
	    	// Now repeat step 1-4 on this group
    		grouping(groups, lhat, minNumGroup);
    	}
    }
}
    	
// The vector parentPool upon finish contains the indices of parents in _visited, not in _genomes!
// There are _populationSize/2 many individual selected 
template<class Genome>   
void NichePopulationWTR<Genome>::SUSLinearRanking(std::vector<uint>& parentPool)
{
    // prob is the selection probability of individual indexed by _genomes[i]
    std::vector<double> prob(_populationSize, 0);
    // prob is calculated by invoking a specific routine, it could be linear ranking, or the grouping procedure
    // First treat everything as in the deterministic case. If this is indeed a deterministic problem, then skip
    // subsequent grouping steps. If not, this prob vector will be used in groupranking method.
    linearRanking(prob);
    if (_stocsim)
    {	// In noisy environment, first do a grouping, then do group ranking to assign modified selection probabilities
    	std::vector<uint> firstgroup;
    	for (uint i = 0; i < (uint) _populationSize; i++) firstgroup.push_back(i);
    	std::vector<std::vector<uint> > groups;
    	groups.push_back(firstgroup);
    	grouping(groups, 0, _populationSize / 6);
    	groupRanking(prob, groups);
    }
    double* chart = new double[_populationSize];
    chart[0] = prob[0];
    for (int j = 1; j < _populationSize; j++)
    {
		chart[j] = chart[j-1] + prob[j];   
    }
	// Do only one spin of the roulettewheel
	double pointer = rng.RandU01();
	// Pointer spacing determines how many samples we take. Because of mating restriction, I only
	// select half population size. So the for-loop only does _populationSize/2 loop
	double pointerSpacing = 1.0 / _populationSize * 2; 
	for (int j = 0; j < _populationSize/2; j++)
    {
    	int counter = 0;
    	for (counter = 0; counter < _populationSize; counter++)
    	{
    		if (pointer < chart[counter])
    		{
    			break;
    		}
    	}
    	// counter is where the pointer points to on the roulette wheel after the spin
    	// Note that counter is actually an index in _genomes, while parentPool holds index in _visited
    	parentPool.push_back(_genomes[counter]);
    	// Now move one-pointer interval to next reading
    	pointer = (1 < pointer + pointerSpacing) ? (pointer + pointerSpacing - 1) : (pointer + pointerSpacing);
    }	
    delete []chart;
}
 
template<class Genome>   
void NichePopulationWTR<Genome>::selectParentPool(std::vector<uint>& parents)
{	// parentPool has the indexes of those genomes in the vector _visited.
	std::vector<uint> parentPool;
	switch (_selectionType)
	{
		case 1 : // tournament selection
		    for (int i = 0; i < _populationSize; i++)
    		{
		    	parentPool.push_back(tournamentSelection());
		    }
		    break;
		case 2 : // SUS + linear ranking
			SUSLinearRanking(parentPool);
			break;
		default :
			break;
	}
    // parentPool has the indexes of those genomes in the vector _visited. 
	for (uint i = 0; i < parentPool.size(); i++)
    { // Select a parent according to crossover probability
		if (rng.RandU01() < _crossoverProbability)
		{
			parents.push_back(parentPool[i]);
		}
    }
}  

// Mating restriction returns the second parent for crossover in niching GA. The returned value is the
// index of the second parent in _visited.
template<class Genome>   
uint NichePopulationWTR<Genome>::matingRestriction(uint parent1)
{
	return inbreeding(parent1);
}

/* Mating restriction technique. Use inbreeding right now.
 * dynamic inbreeding restrictions. i.e., the first parent is selected 
 * using SUS + linear ranking, and then randomly sample M individuals from the parent solution, 
 * and select the best individual that belongs to the same niche as the first parent as the 
 * second parent; if there is no such individual, simply select the one that is closest to the 
 * first parent as the second parent. The returned value is the index in _visited.
 */
template<class Genome>   
uint NichePopulationWTR<Genome>::inbreeding(uint parent1)
{
	// randomly sample M individuals without replacement from the parent solution
	std::vector<int> MFparents(_inbreedingMF);
	// This function is declared in RandomVariates.h. It generates unique random integers from 0 to the argument inclusively
	URandIntGen urg(_populationSize - 1);
	// Note that _inbreedingMF has the indexes in _genomes, not _visited!!
 	std::generate_n(MFparents.begin(), _inbreedingMF, urg);
 	// select the best individual that belongs to the same niche as the first parent as the second parent
	// First see if there is an individual that is within the same niche as parent1
	int bestIndex = -1;
	double bestFitness = 1e50; //std::numeric_limits<double>::max();
	for (int i = 0; i <  _inbreedingMF; i++)
	{
		for (uint j = 0; j < _niches.size(); j++)
		{	//Note it's the index in _visited that is recorded in _niches. And sameNiche() also expects
			//indices of _visited. So needs be consistent
			if (_niches[j]->sameNiche(parent1, _genomes[MFparents[i]]))
			{
				break;
			}
		} 
		if (_visited[_genomes[MFparents[i]]]->getSharedFitness() < bestFitness)
		{
			bestIndex = i;
			bestFitness = _visited[_genomes[MFparents[i]]]->getSharedFitness();
		}
	}
	if (0 <= bestIndex)
	{
		return uint(_genomes[MFparents[bestIndex]]);
	}
	else
	{
		// simply select the one that is closest to the first parent as the second parent
		int closestIndex = -1;
		double closestDistance = 1e50; //std::numeric_limits<double>::max();
		for (int i = 0; i <  _inbreedingMF; i++)
		{
			if (distance(parent1, _genomes[MFparents[i]]) < closestDistance)
			{
				closestIndex = i;
				closestDistance = distance(parent1, _genomes[MFparents[i]]);
			}
		}
		assert(0 <= closestIndex);
		return uint(_genomes[MFparents[closestIndex]]);	
	}
}


template<class Genome>   
uint NichePopulationWTR<Genome>::linebreeding(uint parent1)
{
	return 1;
}
/* Crossover P(t) with mating restriction
 */
template<class Genome>   
void NichePopulationWTR<Genome>::crossover(const std::vector<uint>& parents, std::vector<Genome*>& offspring)
{	//Make sure there are even number of parents
    uint numOfParents = uint(parents.size());
    Genome* child1;
    Genome* child2;
    for (uint i = 0; i < numOfParents; i++)
    { 
		// The first parent is selected as usual according to its shared fitness.
		// Note that parents holds indexes in_visited, not in _genomes!
		uint parent1 = parents[i];
		// The second parent is found by matingRestriction
		// matingRestriction also needs to return the index in _visited.
		uint parent2 = matingRestriction(parent1);
		child1 = 0;
		child2 = 0;
		switch (_crossoverType)
		{
			case 1 : //Perform arithmetical crossover.
				arithmeticalCrossover(parent1, parent2, child1, child2);
				break;
			case 2 : //Perform simple crossover.
				simpleCrossover(parent1, parent2, child1, child2);
				break;
			case 3 : //Perform heuristic crossover.
				heuristicCrossover(parent1, parent2, child1);
				child2 = 0;
				break;
			default :
				break;
		}
	// Add those genomes to _visited, and record their indexes in _visited in "offspring" if the pointers are not empty
		if (child1) 
			offspring.push_back(child1);
		if (child2)
			offspring.push_back(child2);
    }
}


template<class Genome>  
void  NichePopulationWTR<Genome>::arithmeticalCrossover(uint parent1, uint parent2, Genome*& child1, Genome*& child2)
{
	double alpha = rng.RandU01();
	int* newGenome1 = new int[_genomeSize]; 
	int* newGenome2 = new int[_genomeSize]; 
	for (int dim = 0; dim < _genomeSize; dim++)
	{ 
		double temp = *(_visited[parent1]->getDiscreteSolution() + dim) * alpha \
			+ *(_visited[parent2]->getDiscreteSolution() + dim) * (1 - alpha);
		if (rng.RandU01() < 0.5)
		{ // Round fractions up and down with equal probability
			newGenome1[dim] = (int) ceil(temp);
		}
		else
		{
			newGenome1[dim] = (int) floor(temp);
		}
		temp = *(_visited[parent1]->getDiscreteSolution() + dim) * (1 - alpha) \
			+ *(_visited[parent2]->getDiscreteSolution() + dim) * alpha;
		if (rng.RandU01() < 0.5)
		{ // Round fractions up and down with equal probability
			newGenome2[dim] = (int) ceil(temp);
		}
		else
		{
			newGenome2[dim] = (int) floor(temp);
		}
	}
	// Now construct new genomes
	child1 = new Genome(_genomeSize, newGenome1);
	child2 = new Genome(_genomeSize, newGenome2);
	if (!isFeasible(child1)) 
	{
		delete child1;
		child1 = new Genome(_genomeSize, _visited[parent1]->getDiscreteSolution());
	}
	if (!isFeasible(child2))
	{
		delete child2;
		child2 = new Genome(_genomeSize, _visited[parent2]->getDiscreteSolution());
	} 
	delete []newGenome1;
    delete []newGenome2;
}

template<class Genome>  
void NichePopulationWTR<Genome>::simpleCrossover(uint parent1, uint parent2, Genome*& child1, Genome*& child2)
{
	int crossoverPoint = rng.RandInt(0, _genomeSize - 1);
	int* newGenome1 = new int[_genomeSize]; 
	int* newGenome2 = new int[_genomeSize]; 
	for (int dim = 0; dim <= crossoverPoint; dim++)
	{ 
		newGenome1[dim] = *(_visited[parent1]->getDiscreteSolution() + dim);
		newGenome2[dim] = *(_visited[parent2]->getDiscreteSolution() + dim);
	}
	for (double alpha = 1; 0 < alpha; alpha -= 1 / 10.0)
	{
		for (int dim = crossoverPoint + 1; dim < _genomeSize; dim++)
		{
			double temp = *(_visited[parent1]->getDiscreteSolution() + dim) * (1 - alpha) \
				+ *(_visited[parent2]->getDiscreteSolution() + dim) * alpha;
			if (rng.RandU01() < 0.5)
			{ // Round fractions up and down with equal probability
				newGenome1[dim] = (int) ceil(temp);
			}
			else
			{
				newGenome1[dim] = (int) floor(temp);
			}
			temp = *(_visited[parent1]->getDiscreteSolution() + dim) * alpha \
				+ *(_visited[parent2]->getDiscreteSolution() + dim) * (1 - alpha);
			if (rng.RandU01() < 0.5)
			{ // Round fractions up and down with equal probability
				newGenome2[dim] = (int) ceil(temp);
			}
			else
			{
				newGenome2[dim] = (int) floor(temp);
			}
		}
		// Need to check if I obtain 2 feasible solutions. If not, decrease alpha by 1/RO, and repeat.
		// If no feasible solution can be obtained, return the parents. 
		child1 = new Genome(_genomeSize, newGenome1);
		child2 = new Genome(_genomeSize, newGenome2);	
		if (isFeasible(child1) && isFeasible(child2))
		{
			break;
		}
		else
		{
			delete child1;
			child1 = 0;
			delete child2;
			child2 = 0;
		}
	}
	// If no  new genomes is constructed, simply return the parents
	if (child1 == 0 || child2 == 0)
	{
		for (int dim = 0; dim <= _genomeSize; dim++)
		{	 
			newGenome1[dim] = *(_visited[parent1]->getDiscreteSolution() + dim);
			newGenome2[dim] = *(_visited[parent2]->getDiscreteSolution() + dim);
		}
		child1 = new Genome(_genomeSize, newGenome1);
		child2 = new Genome(_genomeSize, newGenome2);	
	}
	delete []newGenome1;
    delete []newGenome2;
}

template<class Genome>  
void NichePopulationWTR<Genome>::heuristicCrossover(uint parent1, uint parent2, Genome*& child)
{
	int* newGenome = new int[_genomeSize]; 
	// Parent2 is no worse than parent1. Note this is a minimization problem
	if (_visited[parent1]->getSampleMean() < _visited[parent2]->getSampleMean())
	{
		uint temp = parent1;
		parent1 = parent2;
		parent2 = temp;
	}
	for (int w = 0; w < W; w++)
	{	// Give W trials for heuristic crossover
		for (int dim = 0; dim < _genomeSize; dim++)
		{
			double temp = (*(_visited[parent2]->getDiscreteSolution() + dim) - \
				*(_visited[parent1]->getDiscreteSolution() + dim)) * rng.RandU01() + 
				*(_visited[parent2]->getDiscreteSolution() + dim);
			if (rng.RandU01() < 0.5)
			{ // Round fractions up and down with equal probability
				newGenome[dim] = (int) ceil(temp);
			}
			else
			{
				newGenome[dim] = (int) floor(temp);
			}
		}
		// Need to check if I obtain a feasible solutions. If not, try again.
		// If no feasible solution can be obtained, no offspring is produced. 
		child = new Genome(_genomeSize, newGenome);
		if (isFeasible(child))
		{
			break;
		}
		else
		{
			delete child;
			child = 0;
		}
	}
	delete []newGenome;
}

template<class Genome>   
void NichePopulationWTR<Genome>::mutate(std::vector<Genome*>& offspring, std::vector<uint>& offspringIndex)
{	// Uniform mutation is used here
	/* First determine if this genome needs be mutated. Note that here I only operate on new genomes.
	 * Be very careful here. I put all "new" genomes produced by crossover into the global list _visited. And offspring
	 * has the indexes of those "new" genomes. The are "new", not new, because crossover may very likely produce genomes already 
	 * existing.  Therefore, if I inadvertently directly modify that genome indexed by an element of offspring, I will change the genes
	 * of an existing genome, which is also indexed by many other individuals. So here is a problem of ownership. There are many vectors,
	 * such as offspring, _genomes, that have elements that point to the same element in _visited. So if one mutation is to carry out, needs to
	 * get a new one and insert it into _visited, rather than changing it directly.
	 */
	for (uint pop = 0; pop < offspring.size(); pop++)
	{	
		if (rng.RandU01() < _mutationProbability)
		{	// Then randomly select a copmonent of the genome
			switch (_mutationType)
			{
				case 1 :  // Uniform mutation
					// Then insert this mutated gene into _visited. It may already exist, in that case, AddSolution will delete this pointer
					uniformMutate(offspring[pop]);
					break;
				case 2 :  // Non-Uniform mutation
					nonuniformMutate(offspring[pop]);
					break;
				default :
					break;
			}
		}
		// Then insert this new genome, wehter mutated or not, into _visited. It may already exist, in that case, AddSolution will delete this pointer
		offspringIndex.push_back(AddSolution(_visited, offspring[pop]));
		// This is a very ugly thing to do. Ideally, anything inserted into _visited should already
		// be evaluated. But some of those offspring may be replaced by niche centers during niche 
		// conservation procedure. During that procedure, the fitness values of individuals are 
		// examined to determine whether they should be replaced or not. So there is no way to avoid this 
		// waste. The evaluation is done at the beginning of getNewGeneration().
	}	
}

template<class Genome>   
void NichePopulationWTR<Genome>::uniformMutate(Genome* offspring)
{
	// Randomly select an allele of the genome
	int mutatedGene = rng.RandInt(0, _genomeSize - 1);
	int int_lower = 0, int_upper = 0;
	fathom(_A, _b, offspring, mutatedGene, int_lower, int_upper);
	int length = int_upper - int_lower;
	int step = rng.RandInt(0, length);
	// To mutate this gene
	offspring->changeOneDimension(int_lower + step, mutatedGene);
}

template<class Genome>   
void NichePopulationWTR<Genome>::nonuniformMutate(Genome* offspring)
{
	// Randomly select an allele of the genome
	int mutatedGene = rng.RandInt(0, _genomeSize - 1);
	int int_lower = 0, int_upper = 0;
	fathom(_A, _b, offspring, mutatedGene, int_lower, int_upper);
	// xk is the position of the to be mutated gene. 
	int xk =  *(offspring->getDiscreteSolution() + mutatedGene);
	int xkprime;
	int maxgen = (_generation < _maxGeneration) ? _maxGeneration : _generation;
	if (rng.RandU01() < 0.5)
	{
		xkprime = xk + delta(_generation, maxgen, xk, int_upper);
	}
	else
	{
		xkprime = xk - delta(_generation, maxgen, xk, int_lower);
	}
	// To mutate this gene
	offspring->changeOneDimension(xkprime, mutatedGene);
}


template<class Genome>   
void NichePopulationWTR<Genome>::getNewGeneration(const std::vector<uint>& offspring)
{
	// It is possible to implement various replacement schemes here. I simply replace the worst individuals in the old generation
    // As a result, elitism is also implemented here.
    // Note that here it is a minimization problem.
    // rank has the ordered _genomes, with worst in the front and best in the end. The content of rank is the index of the
    // genome as in the vector _genomes.
    // First do a sort in decreasing order of fitness
    uint* rank = new uint[_populationSize];
    for (int i = 0; i < _populationSize; i++)
    	rank[i] = i;
    if (_elitism)
    {
    	for (int i = _populationSize - 1; 0 < i; i--)
    	{
    		for (int j = 0; j < i; j++)
    		{
    			if (_visited[_genomes[rank[j]]]->getSampleMean() < _visited[_genomes[rank[j + 1]]]->getSampleMean())
    			{
    				uint temp = rank[j];
    				rank[j] = rank[j + 1];
    				rank[j + 1] = temp;
    			}
    		}
    	}
    }
    // First do elitism, replace worst individuals in the previous generation with new offspring
    for (uint i = 0; i < offspring.size(); i++)
    {
    	_genomes[rank[i]] = offspring[i];    
    }
    delete []rank;
    _generation++;
    // Need to evaluate the new population before it is subject to niche conservation.
    _bestGenome = evaluatePopulation();
    if (_elitism)
    {
    	conserveAllNicheCenters(offspring);
    	// After conservation, it is somewhat redundant to evaluate it again. Just do it to be safe
    	// and consistent with the logic.
    	_bestGenome = evaluatePopulation();
    }
}

/* This method combines the new offspring generated by crossover and mutation and the old 
 * population to generate the new generation. There are many ways to do it, for example, whether to 
 * keep old niche centers or not. This method should take the vector containing new offsprings as input
 * parameter and replace _genomes with the new generation upon finishing.
 */
template<class Genome>   
void NichePopulationWTR<Genome>::conserveAllNicheCenters(const std::vector<uint>& offspring)
{
	// Now conserve niche centers in the parent population
    bool* unprocessed = new bool[_populationSize];
    for (int i = 0; i < _populationSize; i++)
    	unprocessed[i] = true;
    for (uint i = 0; i < _niches.size(); i++)
    {
    	// Select the worst unprocessed individual from the new population that lies within the niche of that niche center
    	// If such an individual exists, replace it with nichecenter.
    	// worstWithinNiche & worst have the index of the worst individual within a niche as in the vector _genomes, not _visited!
    	int worstWithinNiche = -1;
    	int worst = -1;
    	// Rember I am dealing with minimization problem. So worst fitness should be bigger
    	double worstFitnessWithinNiche = -1e50; //std::numeric_limits<double>::max();
    	double worstFitness = -1e50; //std::numeric_limits<double>::max();
    	for (int j = 0; j < _populationSize; j++)
    	{
    		if (unprocessed[j] && distance(_genomes[j], _niches[i]->getCenter()) <= _nicheRadius \
    			&& worstFitnessWithinNiche < _visited[_genomes[j]]->getSampleMean())
    		{
    			worstWithinNiche = j;
    			worstFitnessWithinNiche = _visited[_genomes[j]]->getSampleMean();
    		}
    		if (unprocessed[j] && worstFitness < _visited[_genomes[j]]->getSampleMean())
    		{
    			worst = j;
    			worstFitness = _visited[_genomes[j]]->getSampleMean();
    		}
    	}
    	if (-1 < worstWithinNiche)
    	{
    		// Replace worstWithinNiche with old niche center if it is worse than old niche center
    		// Otherwise, discard the old niche center because a better solution is found in its vicinity.
    		if (_visited[_niches[i]->getCenter()]->getSampleMean() < _visited[_genomes[worstWithinNiche]]->getSampleMean())
	    	{
	    		_genomes[worstWithinNiche] = _niches[i]->getCenter();
	    	}
    		// Mark this individual as processed
    		unprocessed[worstWithinNiche] = false;
    	}
    	else
    	{
    		// Else, select the worst unprocessed individual from the new population and replace it with the niche center.
    		_genomes[worst] = _niches[i]->getCenter();
    		// Mark this individual as processed
    		unprocessed[worst] = false;
    	}
    }
  	delete []unprocessed;
}

template<class Genome>   
void NichePopulationWTR<Genome>::fathom(const std::vector<double*>& A, const std::vector<double>& b, const Genome* const genome, int allele, int& int_lower, int& int_upper)
{
	// Now needs to find the range <left, right>. The idea is the same as in RMD
	uint numOfConstraints = uint(A.size());
	double* b1 = new double[numOfConstraints];
	for (uint k = 0; k < numOfConstraints; k++)	
	{// For each constraint
		double sum = 0;
		for (int l = 0; l < _genomeSize; l++)
		{	// Do a matrix multiplication
			if (l != allele)
			{
				sum += A[k][l] * genome->getDiscreteSolution()[l];
			}
		}
		b1[k] = b[k] - sum;
	}
	double upper = UPPER_BOUND, lower = UPPER_BOUND;
	for (uint k = 0; k < numOfConstraints; k++)
	{
		double temp = 0;
		// temp is the temporary value of x_i to make the jth constraint tight
		if (dabs(A[k][allele]) > EPSILON_STABILITY)
		{
			temp = b1[k] / A[k][allele];
		}
		else
		{
			temp = LARGE_NUMBER;
		}
		if (temp > genome->getDiscreteSolution()[allele] + EPSILON_STABILITY)
		{
			//If the value to make the constraint tight is greater than the value of the current point,
			//it means that there is space "above" the current point, and the upper bound could be shrinked, until
			//the upper bound becomes the current point itself or cannot be smaller than 1.
			if (temp - genome->getDiscreteSolution()[allele] < upper)
			{
				upper= temp - genome->getDiscreteSolution()[allele];
			}
		}
		else if (temp < genome->getDiscreteSolution()[allele] - EPSILON_STABILITY)
		{
			if ( genome->getDiscreteSolution()[allele] - temp < lower)
			{
				lower =  genome->getDiscreteSolution()[allele] - temp;
			}
		}
		else
		{
			// The constraint is already tight at current value, i.e., the point is now on the boundary. !!!!!!!!!!!
			// If the coefficient is positive, then increasing, i.e., moving "up" will reenter feasible region because the 
			// inequalitys are Ax>=b
			if (A[k][allele] > 0)
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
	int_upper = (int) floor(upper) + genome->getDiscreteSolution()[allele];
	int_lower = genome->getDiscreteSolution()[allele] - (int) floor(lower);
	delete []b1;
}

template<class Genome>   
double NichePopulationWTR<Genome>::distance(uint i, uint j) const
{
	return _visited[i]->getEuclideanDistanceTo(*_visited[j]);
}


template<class Genome>   
void NichePopulationWTR<Genome>::print() const
{
	for (uint i = 0; i < _genomes.size(); i++)
		_visited[_genomes[i]]->printObservations();
}

template<class Genome>   
void NichePopulationWTR<Genome>::print(const std::vector<uint>& genomes) const
{
	printVector(genomes);
	for (uint i = 0; i < genomes.size(); i++)
		_visited[genomes[i]]->printObservations();
	std::cout << std::endl;
}

// This method prints out the individuals indexed twice, i.e., its position is _genomes[second_level_indexes[i]]
template<class Genome>   
void NichePopulationWTR<Genome>::print(const std::vector<uint>& second_level_indexes, bool) const
{
	for (uint i = 0; i < second_level_indexes.size(); i++)
		_visited[_genomes[second_level_indexes[i]]]->printObservations();
}

#endif /*NICHEPOPULATIONWTR_H_*/
