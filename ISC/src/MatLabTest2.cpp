/*
	A test file to test the functionality of the the Matlab and Python libraries
*/

#include "MatlabEngine.hpp"
#include <stdio.h>
#include <Python.h>
#include <string>
#include <sstream>
#include <fstream>
#include <omp.h>
#include <vector>

#pragma comment( lib, "libeng.lib" )   // Nov 11,2012 
#pragma comment( lib, "libmx.lib" )    // Nov 11,2012 
#pragma comment( lib, "libmat.lib" )   // Nov 11,2012 
#pragma comment( lib, "lpsolve55.lib" )   

int main()
{
/*	std::stringstream command; // The matlab function to call
	omp_set_dynamic(0); // disables dynamic teams
	omp_set_num_threads(8); // sets the number of threads to use
	
#pragma omp parallel // defines a parallel section
	{
		int taskID = omp_get_thread_num(); // Gives the threads its task number
		int numberOfTasks = omp_get_num_threads(); //Tells the thread how many threads there are
		taskID++; // Increments the TaskID since Matlab is expecting the taskID to start at 1 and not 0

		// Opens a matlab engine
		matlab::engine::FutureResult<std::unique_ptr<matlab::engine::MATLABEngine>> futureMatlab = matlab::engine::startMATLABAsync();
		std::unique_ptr<matlab::engine::MATLABEngine> matlabPtr = futureMatlab.get();

		printf("Starting Thread %d/%d\n", taskID, numberOfTasks);

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
	printf("Done Running Threads\n");
*/
	Py_Initialize();
	wchar_t** wargv;
	wargv = (wchar_t**)malloc(1 * sizeof(wchar_t*));
	*wargv = (wchar_t*)malloc(14 * sizeof(wchar_t));
	**wargv = L'aggregator.py'; 
	for(int i = 0; i < 10; i++)
	{
	//Py_Initialize(); 

	FILE* file = _Py_fopen("./log_analyzer/aggregator.py", "r+"); // opens the python script you want to run
	if (file != NULL) // checks the file was able to be opened
	{
		//char* p_argv[1]; // creates an arguement list. This is like the argv used by the c++ main function
		//p_argv[0] = "aggregator.py"; // Sets the first parameter of the argument list to the name of the program
		//Py_SetProgramName((wchar_t*)p_argv[0]); // Sets the program name to the name put in the arguement list
		//PySys_SetArgv(1, (wchar_t**)p_argv); // gives python the arguement list
		//Py_SetProgramName((wchar_t*)*wargv);
		//PySys_SetArgv(1, (wchar_t**)wargv);
		PyRun_SimpleFile(file, "./log_analyzer/aggregator.py"); // Tells python to run the file as "__main__"
	}

	file = _Py_fopen("./log_analyzer/main.py", "r+"); // opens the python script you want to run
	if (file != NULL) // checks the file was able to be opened
	{
		//char* p_argv[1]; // creates an arguement list. This is like the argv used by the c++ main function
		//p_argv[0] = "main.py"; // Sets the first parameter of the argument list to the name of the program
		//Py_SetProgramName((wchar_t*)p_argv[0]); // Sets the program name to the name put in the arguement list
		//PySys_SetArgv(1, (wchar_t**)p_argv); // gives python the arguement list
		free(*wargv);
		*wargv = (wchar_t*)malloc(8 * sizeof(wchar_t));
		**wargv = L'main.py';
		//Py_SetProgramName((wchar_t*)*wargv);
		PySys_SetArgv(1, (wchar_t**)wargv);
		PyRun_SimpleFile(file, "./log_analyzer/main.py"); // Tells python to run the file as "__main__"
		//free(*wargv);
	}

	//Py_Finalize();

	//free(*wargv);
	//free(wargv);
	}
	free(wargv);
	Py_Finalize();
	//system("PAUSE");
	return 0;
}
