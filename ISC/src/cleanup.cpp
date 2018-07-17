#include "masterheader.h"
using namespace std;
/* Given the list of visited solutions (visited), a list of local minimum (localmin), a constant delta, a 
 * constant alpha, cleanup() performs a NSGC procedure and report the best minimum.
 * Return value is the index of the best in visited.
 * Note that this cleanup procedure is meant for maximization problems instead. To use it for minimization
 * problem, in the begining of the procedure, need to negate the sample means and work on 
 * negated sample means instead.
 */
long cleanup(Simulator* simulator, vector<Solution*>& visited, vector<long> M, double delta, double alpha)
{
	// Input
	cout << "Clean up...     ";
	size_t k = M.size(); 
	if (k == 0)
	{
		cout << "Empty local minimum set, error" << endl;
		return -1;
	}
	int* num = new int[k];
	double *x = new double[k];
	double *S2 = new double[k];
	double *t = new double[k];
	double *w = new double[k * k];
	for (size_t i = 0; i < k; i++)
	{
		num[i] = visited[M[i]]->getNumOfObservations();
		// Note that NSGS is for maximization, so need to negate objective values here!!!
		x[i] = -visited[M[i]]->getSampleMean();
		S2[i] = visited[M[i]]->getSampleVariance();
	}
	// Step 1
	int n = num[minimum(num, int(k))];
	double plevel = pow(1 - alpha/2, 1/max(1, (double)(k-1)));
	//cout << "Rinott computation, 1-alpha/2 = " << 1-alpha/2 << " pow " << 1/max(1, (double)(k-1)) << endl;
	//cout << "n= " << n << ", " << "p= " << plevel << " alpha= " << alpha << endl;
	//Rinott(systems under comparison, CI level, d.f.)
	double h = Rinott(2, plevel, n);
	//cout << "Rinott " << h << endl;
	for (size_t i = 0; i < k; i++)
	{
		// Note that there should be at least 2 observation, if there is only 1, then variance must be 0
		// and this t[i] value has no significance. But tinv() is not well defined for 0 as the 2nd arg.
		// So set it as below
		t[i] = tinv(plevel, (1 <= num[i] - 1) ? num[i] - 1 : 1);
	}
	for (size_t i = 0; i < k; i++)
	{
		for (size_t l = 0; l < k; l++)
		{
			if (i == l)
			{
				w[i*k + l] = 0;
			}
			else
			{
				w[i*k + l] = sqrt(t[i]*t[i]*S2[i]/num[i] + t[l]*t[l]*S2[l]/num[l]);
			}
		}
	}
	// Step 2: Get set I.
	// Note that I has indexes in M, not visited!!!
	vector<long> I;
	for (size_t i = 0; i < k; i++)
	{
		bool recordi = true;
		for (size_t l = 0; l < k ; l++)
		{
			if (i != l && x[l] - w[i*k + l] > x[i])
			{
				recordi = false; 
				break;
			}
		}
		if (recordi) I.push_back(long(i));
	}
	// Step 3: for all i in I, compute new number of samples and collect more observations
	size_t kI = I.size();
	double* xI = new double[kI];
	if (kI == 0) 
	{
		cout << "Clean-up error: Empty I" << endl;
		return -1;
	}
	double addreps = 0;
	for (size_t i = 0; i < kI; i++)
	{
		int newni = (int) ceil(h*h*S2[I[i]]/delta/delta);
		// I[i] is the index in M.
		int N = (num[I[i]] < newni) ? newni : num[I[i]];
		addreps += N - num[I[i]];
		//cout << "Add " << N - num[I[i]] << " more reps for " << *visited[M[I[i]]] << endl;
		for (int rep = 0; rep < N - num[I[i]]; rep++)
		{
			visited[M[I[i]]]->recordObservation(simulator->simulation(visited[M[I[i]]]));
		}
		// Note this clean up procedure is for maximization problem. So use -sampleMean here!!!!!
		xI[i] = -visited[M[I[i]]]->getSampleMean();
	}
	cout << "Cleanup used " << addreps << " reps" << endl; 
	// Step 4: Find the best. Note bestIndex is the index in vector I, not in M! not in visited!
	int bestindex = maximum(xI, int(kI));
	delete []num;
	delete []x;
	delete []S2;
	delete []t;
	delete []w;
	delete []xI;
	return M[I[bestindex]];
}



