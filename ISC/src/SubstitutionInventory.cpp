#include "SubstitutionInventory.h"

using namespace std;

double SubstitutionInventory::simulation(const Solution* const solution)
{
	return -substitutioninventory(solution); 
}

/*double SubstitutionInventory::sim_batch(vector<Solution*>& visited, vector<size_t>& solnid, vector<size_t> numrep)
{
	return 0;
}*/

double SubstitutionInventory::gettruevalue(const Solution* const solution)
{
	return solution->getSampleMean();
}

double SubstitutionInventory::substitutioninventory(const Solution* const solution)
{
	// The number of replications a solution has received determines the index of the substream of  
	// random number generators used in this class.
	int substreamindex = solution->getNumOfObservations();
	// First reset random number generators to its starting point
	_rngpoisson.ResetStartStream();  
	_rngutility.ResetStartStream();
	for (int i = 0; i < substreamindex; i++)
	{	// Then advance it to the proper substream for this replication. So that common random numbers
		// are used for one rep for different systems. 
		_rngpoisson.ResetNextSubstream();
		_rngutility.ResetNextSubstream();
	}
	int numOfVariants = solution->getDimension();
	int* initialInventory = new int[numOfVariants];
	for (int i = 0; i < numOfVariants; i++) initialInventory[i] = *(solution->getDiscreteSolution() + i);
	// require(numOfVariants == 10, "Wrong number of product variants");
	// mu determines variance of error term
	const double mu = 1.5;
	// price of products
	const double pp = 8;
	// customer arrival rate
	const double CustomerArrivalRate = 1000;
	// mean order quantity by each customer
	const double Q = 1; 
	// quality of each variant and no buy option
	double* a = new double[numOfVariants + 1];
	a[0] = 0; //4.0;
	for (int i = 1; i < numOfVariants + 1; i++)	a[i] = 12.25 - 0.5 * (i - 1);
	// cost of each variant
	double* c = new double[numOfVariants+1];
	c[0] = 0;
	for (int i = 1; i < numOfVariants + 1; i++)	c[i] = 3;
	// price of each variant 
	double* p = new double[numOfVariants+1];
	p[0] = 0;
	for (int i = 1; i < numOfVariants + 1; i++)	p[i] = pp; 
	// Compute procurement cost for this initial inventory setup
	double procurementcost = 0;
	// initialInventory only has numOfVariants elements, which correspond to stock levels for these items
	for (int i = 1; i <= numOfVariants; i++)	procurementcost += c[i] * initialInventory[i - 1];
	// Set up inventory
	int* stock = new int[numOfVariants];
	for (size_t i = 0; int(i) < numOfVariants; i++)	stock[i] = initialInventory[i];  
	double* utilities = new double[numOfVariants + 1];
	double revenue = 0;
	// Start simulation
	// First determine number of customers
	unsigned int T = Poisson(CustomerArrivalRate, _rngpoisson);
	//unsigned int T = Poisson(CustomerArrivalRate);
	for (unsigned int i = 0; i < T; i++)
	{	// for each customer, first generate order quantity
		unsigned int quantity = Poisson(Q);
		if (quantity == 0) continue;
		// for each customer, generate utility vector
		for(int j = 0; j <= numOfVariants; j++)	utilities[j] = a[j] - p[j] + Gumbel(mu, _rngutility);
		//for(int j = 0; j <= numOfVariants; j++)	utilities[j] = a[j] - p[j] + Gumbel(mu);
		// for(int j = 0; j <= numOfVariants; j++)	utilities[j] = a[j] - p[j] + ExtremeValue();
		for (unsigned int q = 0; q < quantity; q++)
		{
			size_t maxindex = 0;
			double maxutility = utilities[0];
			for (int j = 1; j <= numOfVariants; j++)
			{
				if (utilities[j] > utilities[maxindex] && stock[j - 1] > 0)
				{
					maxindex = j;
					maxutility = utilities[j];
				}
			}
			if (maxindex > 0)
			{	// customer chose to buy
		  		// Note that utitlities are indexed from 0 to numOfVariants + 1, whereas stock is indexed
		  		// from 0 to numOfVariants - 1. So the index of variant j in stock should be j-1
				require(stock[maxindex - 1] > 0, "Cannot be an empty stock");
				--stock[maxindex - 1];
			}
			revenue += p[maxindex];
		}		
	}
	double profit = revenue - procurementcost;
	delete []initialInventory;
	delete []stock;
	delete []utilities;
	delete []p;
	delete []a;
	delete []c;
	return profit;
}


double SubstitutionInventory::ContinuousSubstitutionInventory(const double* const initialInventory, int numOfVariants)
{
	require(numOfVariants == 10, "Wrong number of product variants");
	// mu determines variance of error term
	const double mu = 1.5;
	// price of products
	const double pp = 8;
	// customer arrival rate
	const double CustomerArrivalRate = 30;
	// mean order quantity by each customer
	const double Q = 1; 
	// quality of each variant and no buy option
	double* a = new double[numOfVariants + 1];
	a[0] = 4.0;
	for (int i = 1; i < numOfVariants + 1; i++)	a[i] = 12.25 - 0.5 * (i - 1);
	// cost of each variant
	double* c = new double[numOfVariants+1];
	c[0] = 0;
	for (int i = 1; i < numOfVariants + 1; i++)	c[i] = 3;
	// price of each variant 
	double* p = new double[numOfVariants+1];
	p[0] = 0;
	for (int i = 1; i < numOfVariants + 1; i++)	p[i] = pp; 
	// Compute procurement cost for this initial inventory setup
	double procurementcost = 0;
	for (int i = 1; i <= numOfVariants; i++)	procurementcost += c[i] * initialInventory[i - 1];
	// Set up inventory
	double* stock = new double[numOfVariants];
	for (size_t i = 0; int(i) < numOfVariants; i++)	stock[i] = initialInventory[i];  
	double* utilities = new double[numOfVariants + 1];
	utilities[0] = a[0] - p[0];
	double revenue = 0;
	// Start simulation
	// First determine number of customers
	unsigned int T = Poisson(CustomerArrivalRate);
	for (unsigned int i = 0; i < T; i++)
	{	// for each customer, first generate order quantity
		double quantity = Exponential(Q);
		// for each customer, generate utility vector
		for(int j = 1; j <= numOfVariants; j++)	utilities[j] = a[j] - p[j] + Gumbel(mu);
		while(quantity > 0)
		{
			double fulfilledQuantity = 0;
			size_t maxindex = 0;
			double maxutility = utilities[0];
			for (size_t j = 1; int(j) <= numOfVariants; j++)
			{
				if (utilities[j] > utilities[maxindex] && stock[j - 1] > 0)
				{
					maxindex = j;
					maxutility = utilities[j];
				}
			}
			if (maxindex > 0)
			{	// customer chose to buy
		  		// Note that utitlities are indexed from 0 to numOfVariants + 1, whereas stock is indexed
		  		// from 0 to numOfVariants - 1. So the index of variant j in stock should be j-1
				require(stock[maxindex - 1] > 0, "Cannot be an empty stock");
				if (stock[maxindex - 1] > quantity)
				{
					stock[maxindex - 1] -= quantity;
					fulfilledQuantity = quantity;
					quantity = 0;
				}
				else
				{
					quantity -= stock[maxindex - 1];
					fulfilledQuantity = stock[maxindex - 1];
					stock[maxindex - 1] = 0;
				}
			}
			else
			{
				quantity = 0;
			}
			revenue += p[maxindex] * fulfilledQuantity;
		}		
	}
	double profit = revenue - procurementcost;
	delete []stock;
	delete []utilities;
	delete []p;
	delete []a;
	delete []c;
	return profit;
}




