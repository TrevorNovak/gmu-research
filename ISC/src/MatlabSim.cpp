#include "masterheader.h"
#include <Python.h>
#include <omp.h>
#include "MatlabEngine.hpp"
//Engine *ep;

using namespace std;
int MatlabSim::isGlobalOptimum(const Solution& solution)
{
	return 0; 
} 

double MatlabSim::simulation(const Solution* const solution)
{
	string iscfilename = _matlabdir + "/iscoutput.txt";
	ofstream iscoutFile(iscfilename.c_str());
	require(iscoutFile.is_open(), "Error: invalid iscoutput file");
	iscoutFile << -1 << "\t" << 1 << "\t" << solution->getNumOfObservations();
	solution->printX(iscoutFile);
	iscoutFile << endl;
	iscoutFile.close();
	ifstream matlabout;
	string matfilename = _matlabdir + "./SimulatorOutput.txt";
	clock_t time1_mat = clock();
	clock_t time2_mat = time1_mat + 3600 * CLOCKS_PER_SEC;
	clock_t currtime_mat;
	double simresult = 0;
	// Executes the simulation using parfor with input argument 1 

	/*int engRt = engEvalString(ep, "fm_cas_parfor_enginetest(1,8);");
	if (engRt == 1)
	{
		cout << "Matlab engine error" << endl;
		engClose(ep);
		exit(7);
	}
	do
	{
		matlabout.open(matfilename);
		currtime_mat = clock();
	} while (!matlabout.is_open() && currtime_mat < time2_mat);
	//MySleep(10);
	size_t solutionid, repcount, taskid;
	if (matlabout.is_open())
	{
		// 1st col: id of the solution in the vector of visited solution
		// 2nd col: number of replications 
		// 3rd col and after: simulation results
		// This function is for individual simulation call, and thus repcount will always be 1 and only one simulation result would be available
		matlabout >> taskid >> solutionid >> repcount >> simresult;
	}*/
	return simresult;
}

double MatlabSim::simulation(vector<Solution*>& visited, const vector<size_t>& solnid, const vector<size_t> numrep)
{
	size_t solncnt = 0;
	if (isBatchSim())
	{	// Write iscoput.txt for the Matlab simulator to use
		string iscfilename = _matlabdir + "/iscoutput.txt";
		ofstream iscoutFile(iscfilename.c_str());
		require(iscoutFile.is_open(), "Error: invalid iscoutput file"); 
		for (size_t si = 0; si < solnid.size(); si++)
		{
			// 1st col: id of the solution in the vector of visited solution
			// 2nd col: number of replications to run 
			// 3rd col: number of replications already spent on this solution 
			// 4th col and after: solution vector 
			if (numrep[si] > 0)
			{
				iscoutFile << solnid[si] << "\t" << numrep[si] << "\t" << visited[solnid[si]]->getNumOfObservations();
				visited[solnid[si]]->printX(iscoutFile);
				iscoutFile << endl;
				solncnt += numrep[si]; // keeps track of total reps to run 
			}
		}
		iscoutFile.close(); 
		if (solncnt == 0) return 0; 
		//int taskID; // the current task ID
		std::stringstream command; // The matlab function to call
		int numberOfTasks; // the total number of threads created
		omp_set_dynamic(0); // disables dynamic teams
		omp_set_num_threads(min(_maxThreads,solncnt)); // sets the number of threads to use
		#pragma omp parallel // defines a parallel section
		{
			int taskID = omp_get_thread_num(); // Gives the threads its task number
			numberOfTasks = omp_get_num_threads(); //Tells the thread how many threads there are
			taskID++; // Increments the TaskID since Matlab is expecting the taskID to start at 1 and not 0
			 // Opens a matlab engine
			matlab::engine::FutureResult<std::unique_ptr<matlab::engine::MATLABEngine>> futureMatlab = matlab::engine::startMATLABAsync();
			std::unique_ptr<matlab::engine::MATLABEngine> matlabPtr = futureMatlab.get();
			cout << "Starting Thread " << taskID << " out of " << numberOfTasks << "threads\n";
			// Moves the matlab engine to the correct working directory
			matlabPtr->eval(matlab::engine::convertUTF8StringToUTF16String("cd /home/cferris/ISC/simulation;"));
			//command << "fm_cas_parfor_enginetest2(" << taskID << "," << numberOfTasks << ");"; // Defines the matlab command to run
			// Calls the simulation
			std::vector<matlab::data::Array> args;
			matlab::data::ArrayFactory factory;
			args.push_back(factory.createScalar<double>(taskID));
			args.push_back(factory.createScalar<double>(numberOfTasks));
			matlabPtr->feval(matlab::engine::convertUTF8StringToUTF16String("fm_cas_parfor_enginetest2"), 0, args);
		}
		cout << "Done Running Threads\n";

		// Wait for Matlab simulator to simulate all solutions in iscoutput.txt
		// Wait for a maximum of 60 minutes in case there is a problem and it does not finish
		/*clock_t time1_mat = clock();
		clock_t time2_mat = time1_mat + 3600 * CLOCKS_PER_SEC;
		clock_t currtime_mat;*/
		// Executes the simulation using parfor with input argument 1 
		/*int engRt = engEvalString(ep, "fm_cas_parfor_enginetest(1,18);"); 
		if (engRt == 1)
		{
			cout << "Matlab engine error" << endl;
			engClose(ep);  
			exit(7);
		}*/

		Py_Initialize();

		// Allocates memory for the commandline arguements of the python scripts
		wchar_t** wargv;
		wargv = (wchar_t**)malloc(1 * sizeof(wchar_t*));
		*wargv = (wchar_t*)malloc(14 * sizeof(wchar_t));
		**wargv = L'aggregator.py';  
		FILE* file = _Py_fopen("./log_analyzer/aggregator.py", "r+"); // opens the python script you want to run
		if (file != NULL) // checks the file was able to be opened
		{
			//Py_SetProgramName((wchar_t*)*wargv);
			//PySys_SetArgv(1, (wchar_t**)wargv);
			PyRun_SimpleFile(file, "./log_analyzer/aggregator.py"); // Tells python to run the file as "__main__"
		}
		file = _Py_fopen("./log_analyzer/main.py", "r+"); // opens the python script you want to run
		if (file != NULL) // checks the file was able to be opened
		{
			free(*wargv);
			*wargv = (wchar_t*)malloc(8 * sizeof(wchar_t));
			**wargv = L'main.py';
			//Py_SetProgramName((wchar_t*)*wargv);
			PySys_SetArgv(1, (wchar_t**)wargv);
			PyRun_SimpleFile(file, "./log_analyzer/main.py"); // Tells python to run the file as "__main__"
		}
		free(*wargv);
		free(wargv);
		
		Py_Finalize();

		string matfilename = "./SimulateOut.txt", stringinput;
		ifstream matlabout(matfilename.c_str());
		require(matlabout.is_open(), "Error: invalid iscoutput file");
		/*do
		{ 
			matlabout.open(matfilename);
			currtime_mat = clock();
		} while (!matlabout.is_open() && currtime_mat < time2_mat);
		MySleep(10);
		if (matlabout.is_open())
		{ */
		size_t solutionid, repcount, lastrep, taskid;
		double simresult;
		//for (size_t si = 0; si < solncnt; si++)
		while(getline(matlabout,stringinput))
		{
			istringstream sistrm(stringinput);
			//matlabout >> taskid >> solutionid >> lastrep >> repcount;
			sistrm >> taskid >> solutionid >> lastrep >> repcount;
			for (size_t ri = 0; ri < repcount; ++ri)
			{	// 1st col: id of the solution in the vector of visited solution
				// 2nd col: number of replications 
				// 3rd col and after: simulation results
				sistrm >> simresult;
				visited[solutionid]->recordObservation(simresult);
			}
		}
		//}
	}
	else
	{
		for (size_t si = 0; si < solncnt; si++)
		{
			for (size_t ri = 0; ri < numrep[si]; ++ri)
			{
				visited[solnid[si]]->recordObservation(simulation(visited[solnid[si]]));
			}
		}
	}
	return 0;
}

MatlabSim::MatlabSim(const string& matlabdir, size_t maxth) : _maxThreads(maxth)
{
	_matlabdir = matlabdir; 
	string cdmatlabdir = "cd " + matlabdir;
	setSimMode(batch);  
	/*if (!(ep = engOpen(NULL)))  
	{
		cout << "Can't start Matlab engine!" << endl;
		//cin.get();
		exit(1);
	}
	//engSetVisible(ep, 0); // Set matlab engine command window invisible
	//_buffer = static_cast<char*>(malloc(12800000 * sizeof(char))); // Creates a buffer to store the stdout
	_buffer = static_cast<char*>(malloc(1024 * sizeof(char))); // Creates a buffer to store the stdout
	//engOutputBuffer(ep, _buffer, 12800000); // Starts recording the standard output from the Matlab engine
	engOutputBuffer(ep, _buffer, 1024); // Starts recording the standard output from the Matlab engine
	int engRt = engEvalString(ep, cdmatlabdir.c_str()); // cd Matlab working directory
	if (engRt == 1) // checks if the 
	{
		fprintf(stderr, "Matlab error!");
		exit(8);
	}
	engRt = engEvalString(ep, "pwd"); 
	if (engRt == 1) // checks if the 
	{
		fprintf(stderr, "Matlab error!");
		exit(8);
	}
	// Prints out the output buffer for pwd, showing current Matlab directory
	if (_buffer[0] != '?') {
		fprintf(stdout, "%s\n", _buffer);
	}
	else {
		fprintf(stderr, "%s\n", _buffer);
	}
	// create a persistent parpool using the input argument 0 
	engRt = engEvalString(ep, "fm_cas_parfor_enginetest(0,18);"); 
	if (engRt == 1) // checks if the 
	{
		fprintf(stderr, "Matlab error!");
		exit(7);
		//return -1;
	}
	if (_buffer[0] != '?') {
		fprintf(stdout, "%s\n", _buffer);
	}
	else {
		fprintf(stderr, "%s\n", _buffer);
	}
	free(_buffer);*/
}

MatlabSim::~MatlabSim()
{
	_buffer = static_cast<char*>(malloc(1024 * sizeof(char))); // Creates a buffer to store the stdout
	/*engOutputBuffer(ep, _buffer, 1024); // Starts recording the standard output from the Matlab engine
	int engRt = engEvalString(ep, "fm_cas_parfor_enginetest(2,18);");
	if (engRt == 1) // checks if the 
	{
		fprintf(stderr, "Matlab error!");
		exit(7);
		//return -1;
	}
	if (_buffer[0] != '?') {
		fprintf(stdout, "%s\n", _buffer);
	}
	engClose(ep); */
	free(_buffer);
}
