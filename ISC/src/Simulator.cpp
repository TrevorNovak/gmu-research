#include "masterheader.h"
using namespace std;

int Simulator::isLocalOptimum(const Solution& solution, const vector<double*> & A, const vector<double> b)
{
	int dimension = solution.getDimension();
	double nv, sv = gettruevalue(&solution);;
	for (int i = 0; i < dimension; i++)
	{
		Solution neighbor(solution);
		neighbor.changeOneDimension(*(solution.getDiscreteSolution()+i)+1, i);
		bool fs = neighbor.isFeasible(A, b);
		if (fs)
		{
			nv = gettruevalue(&neighbor);
		}
		
		if (neighbor.isFeasible(A, b) && gettruevalue(&neighbor) < gettruevalue(&solution))
			return 0;
		neighbor.changeOneDimension(*(solution.getDiscreteSolution()+i)-1, i);
		fs = neighbor.isFeasible(A, b);
		if (fs)
			nv = gettruevalue(&neighbor);
		if (neighbor.isFeasible(A,b) && gettruevalue(&neighbor) < gettruevalue(&solution))
			return 0;
	}
	return 1;
};