#include "masterheader.h"
#include <map>
#include <iomanip>
#include <string>
#include <utility>
#include "FlowLine.h"
#include "SubstitutionInventory.h"
#include <Python.h>


using namespace std; 

#pragma comment( lib, "libeng.lib" )   // Nov 11,2012 
#pragma comment( lib, "libmx.lib" )    // Nov 11,2012 
#pragma comment( lib, "libmat.lib" )   // Nov 11,2012 
#pragma comment( lib, "lpsolve55.lib" )   

//Engine *ep; //Define matlab engine and test it

int main()   
{   
//	Py_Initialize();
	/*ofstream iscmatlab("iscoutput.csv");
	iscmatlab << 1 << ", " << 2 << endl;
	iscmatlab.close();
	int ti=system ("runmatlab.bat");
	printf ("\nThe value returned was: %d.\n",ti);
	  
	ifstream matlabout;
	do
	{
		matlabout.open("C:\\Users\\jxu13.SITE\\Documents\\MATLAB\\matlabout.csv");
	} while(!matlabout.is_open());
	MySleep(0.2);
	int tempa, tempb;
	char delimiter;
	matlabout >> tempa >> delimiter >> tempb; 
	cout << tempa << delimiter << tempb << endl; 
	return 1;*/
	/**************************************************/
	// Make change here for your own input paramter file 
	//const char* inputFiles = {"./inputmd10.txt"}; 
	const char* inputFiles = {"./inputcascade.txt"}; // changed by  Nov 6,2012
	vector<double*> A; 
	vector<double> b;
	ifstream inputFile(inputFiles);
	require(inputFile.is_open(), "Error: invalid input file"); 
	vector<int> xx;
	//long iscbudget;
	//int reps, initialNumReps, numOfCandidates;
	string output; 
	//bool backtracking, ocba, dominantniche, cleanup, stocsim, elitism;
	//int gagen, pruning_freq, mpamode, rmdmode, metamodel;
	double globaldelta, backalpha, backdelta, localoptalpha, localoptdelta, cleanupalpha, cleanupdelta, gabudget;
	/*ReadInput(inputFile, output, xx, A, b, iscbudget, reps, backtracking, \
		ocba, gagen, dominantniche, cleanup, stocsim, globaldelta, backalpha, backdelta, \
		localoptalpha, localoptdelta, cleanupalpha, cleanupdelta, initialNumReps, numOfCandidates, \
		elitism, pruning_freq, mpamode, rmdmode, gabudget, metamodel);*/
	ISCParameters* iscparams = ReadInput(inputFile, output, xx, A, b);
	
	
	int dimension = (int) xx.size(); 
	Solution* x0 = new Solution(xx);
	vector<Solution*> visited;
	vector<Record*> records; 
	/**************************************************/
	// Make change here to use your own simulator
	// set up the simulator, the argument to the ResponseSurface constructor tells which function to use, see its constructor code for detail
	// 0-f22, 1-singular, 2-highd 3-multihighd 4-rosenbrock 5-f2n, 6-quadratic
	//ResponseSurface sim(2);
	//Flowline sim;
	//SubstitutionInventory sim;
	//MatlabSim sim("C:\\Users\\jxu13\\Documents\\MATLAB\\ANL_par\\CascadeSim");
	MatlabSim sim("/home/cferris/ISC/simulation", 20);
	// Set up ISC by putting the starting point into the vector visited.
	visited.push_back(x0); 
	string outputfilename = output + ".dat";
	string resultfilename = output + "result.dat"; 
	ofstream plotFile(outputfilename.c_str());
	require(plotFile.is_open(), "Error: invalid plot file");
	ofstream resultFile(resultfilename.c_str());
	require(resultFile.is_open(), "Error: invalid plot file");
	long globalmins = 0;
	for (int i = 0; i < iscparams->get_reps(); i++) 
	{
		ostringstream logfilename;
		// logfilename is for individual sample path, distinguished by i
		logfilename << "./" << output << i << ".dat";
		ofstream repfile((logfilename.str()).c_str());
		require(repfile.is_open(), "Error: invalid rep file"); 
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
		resultFile << endl << endl << endl << "Macro rep #:"  << i+1 << endl << "ISC starts from" << *visited[0] << endl << "g(x0) = " << sim.gettruevalue(visited[0]) << endl << endl;
		long N; //N will record the # of reps used by ISC once ISCControl finishes
		clock_t isctime1 = clock();
		/*globalmins = ISCControl(&sim, N, gabudget, gagen, globaldelta, 0.1, 0.1, 0.1, \
			backdelta, backalpha, localoptdelta, localoptalpha, cleanupdelta, cleanupalpha, \
			initialNumReps, numOfCandidates, dimension, visited, A, b, repfile, resultFile, records, \
			backtracking, ocba, dominantniche, cleanup, stocsim, elitism, pruning_freq, mpamode, rmdmode, metamodel);*/
		globalmins = ISCControl(&sim, N, dimension, visited, A, b, repfile, resultFile, records, iscparams);
		clock_t isctime2 = clock();
		double seconds = (isctime2-isctime1)/CLOCKS_PER_SEC;
		resultFile << "ISC used " << seconds << " seconds for global and local phase, excluding cleanup\n" << endl;
		cout << "ISC used " << seconds << " seconds for global and local phase, excluding cleanup\n" << endl;
		cout << "Rep " << i << " " << visited[globalmins] << endl;
		resultFile << "ISC used   " << N << " reps, " << visited[globalmins] << endl;
		// resultFile << "Macrorep   Rep   Iter   x*   g(x*)   N(x*)   var(x*)   x0  g(x0)" << endl;
		for (size_t j = 0; j < visited.size(); j++) delete visited[j];
		visited.clear();
  		x0 = new Solution(xx);
		visited.push_back(x0);
		repfile.close();
	}	
	WritePlot(plotFile, records, iscparams->get_budget());
	plotFile.close();
	resultFile.close();
	for (size_t j = 0; j < records.size(); j++) delete records[j];
	records.clear();
	
	ifstream resultFileIn(resultfilename.c_str());
	require(resultFileIn.is_open(), "Error: invalid result file");
	string summaryfilename = output + "_summary.txt";
	ofstream summaryFile(summaryfilename.c_str());
	require(summaryFile.is_open(), "Error: invalid summary file");
	size_t casepos = output.find_last_of("/");
	string caseid = output.substr(casepos+1, output.length()-casepos-1);
	CompileRunStat(resultFileIn, summaryFile, dimension, sim, caseid, A, b);
	
	for (unsigned int i = 0; i < A.size(); i++)	delete [](A[i]);
	A.clear();
	//engClose(ep); // Close matlab engin,Nov 11,2012 
//	Py_Finalize();
	return 0;

} ///:~

