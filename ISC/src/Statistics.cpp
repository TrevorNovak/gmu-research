#include "masterheader.h"

using namespace std;

ostream& operator<<(ostream& os, const Record& record)
{
	return os << "mean " << record.mean << " var " << record.samplevariance << endl;
}


void RecordResult(vector<Record*>& records, long repCount, double objectiveValue)
{
	if ((long)records.size() <  repCount)
	{
		Record* newRecord = new Record(objectiveValue);		
		records.push_back(newRecord);
	}
	else
	{
		records[repCount - 1]->writeRecord(objectiveValue);
	}		
}


void WritePlot(std::ofstream& logFile, std::vector<Record*> records, long budget)
{
	for (long i = 0; i < (long) records.size(); i++)
	{
		if (i+1 <= budget) logFile << 1 + i << ", " << records[i]->getMean() << "\n";
	}
}
//This function reads the result file written out by ISC and checks if the solution found is a local or global optimum
void CompileRunStat(std::ifstream& resultFile, std::ofstream& summaryfile, int dimension, Simulator& sim, const string& caseid, const vector<double*> & A, const vector<double> b)
{
	string oneline;
	int run = 0;
	vector<long> reps;
	vector<int> islocal;
	vector<int> isglobal;
	vector<double> objvalue;
	vector<int> x(dimension, 0);
	while(getline(resultFile, oneline))
	{
		if (oneline.find("ISC used") != string::npos && oneline.find("reps") != string::npos)
		{	
			run++;
			size_t found = oneline.find("used");
			//replace all ',' with blank
			replace(oneline.begin(), oneline.end(), ',', ' ');
			string temp = oneline.substr(found+4, oneline.length()-found-4);
			istringstream tmpstrm(temp);
			long tmprep;
			tmpstrm >> tmprep;
			reps.push_back(tmprep);
			size_t prpos = temp.find("(");
			string temp1 = temp.substr(prpos+1, temp.length()-prpos-1);
			istringstream temp1strm(temp1);
			for (int i = 0; i < dimension; i++) 
				//temp1strm >> x[i] >> temp1; //use this one if solution is printed with comma
				temp1strm >> x[i]; //use this one if printed with space
			Solution xstar(x);
			int islocalopt = sim.isLocalOptimum(xstar, A, b);
			int isglobalopt = sim.isGlobalOptimum(xstar);
			double finalobj = sim.gettruevalue(&xstar);
			if (islocalopt != isglobalopt)
			{
				cout << xstar;
				cout << "local optimum and global optimum status do not match" << endl;
			}
			islocal.push_back(islocalopt);
			isglobal.push_back(isglobalopt);
			objvalue.push_back(finalobj);
		}
	}
	summaryfile << "Row\tN\tMean\tStd\tMax\tMin\tCase" << endl;
	double mean, std, max, min;
	summarystat(reps.begin(), reps.end(), mean, std, max, min);
	summaryfile << "NumReps\t" << reps.size() << "\t" << mean  << "\t" << std << "\t" << max << "\t" << min << "\t" << caseid << endl;
	summarystat(islocal.begin(), islocal.end(), mean, std, max, min);
	summaryfile << "IsLocalOpt\t" << reps.size() << "\t" << mean  << "\t" << std << "\t" << max << "\t" << min << "\t" << caseid << endl;
	summarystat(isglobal.begin(), isglobal.end(), mean, std, max, min);
	summaryfile << "IsGlobalOpt\t" << reps.size() << "\t" << mean  << "\t" << std << "\t" << max << "\t" << min << "\t" << caseid << endl;
	summarystat(objvalue.begin(), objvalue.end(), mean, std, max, min);
	summaryfile << "ObjValue\t" << reps.size() << "\t" << mean  << "\t" << std << "\t" << max << "\t" << min << "\t" << caseid;
}