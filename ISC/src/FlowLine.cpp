#include "FlowLine.h"
using namespace std;

const int JOBARRIVAL = 1;
const int JOBCOMPLETION = 2;
const int UNBLOCK = 3;
// Service rate is in the unit of  # of items/time unit. TOTALTIME and WARMUP both use the
// same time unit.

JobArrivalEvent::JobArrivalEvent(int station, double time) : Event(JOBARRIVAL, time), _station(station) 
{}

JobCompletionEvent::JobCompletionEvent(int station, double time) : Event(JOBCOMPLETION, time), _station(station)
{}
 
UnblockEvent::UnblockEvent(int station, double time) : Event(UNBLOCK, time), _station(station) {}

typedef vector<Configuration> vectorconfig;
typedef vector<Configuration>::iterator vectorcit;

double Flowline::simulation(const Solution* const solution)
{
	int config[5];
	for (int i = 0; i < 4; i++) config[i] = *(solution->getDiscreteSolution() + i);
	config[4] = 20 - config[3];
	// The last parameter indicates a simulation is required.
	// Remember we work on minimization problem, while FlowLine returns the throughput we want to maximize
	return -flowline(3, config, 50, 1000, false);;
}


double Flowline::gettruevalue(const Solution* const solution)
{
	int config[5];
	for (int i = 0; i < 4; i++) config[i] = *(solution->getDiscreteSolution() + i);
	config[4] = 20 - config[3];
	// The last parameter indicates a true value is required.
	return -flowline(3, config, 50, 1000, true);;
}

// linelength: number of stations on this line
// config: solution array, the first linelength many elements are processing rate of station[0] 
// to station[linelength-1], followed by buffer sizes for station[1] to station[linelength]
double Flowline::flowline(int linelength, const int* const config, const double WARMUP, const double TOTAL, bool truevalue)
{
	static bool firstRun = true;
	static vectorconfig trueValueTable;
	if (firstRun && truevalue) 
	{	// If this is the first run of FlowLine requiring truevale, read truevalue
		firstRun = false;
		ifstream inputfile("./truevalue.txt");
		assure(inputfile, "./truevalue.txt");
		string s;
		while(getline(inputfile, s))
		{
			istringstream iss(s);
			int temp[5];
			double value;
			for (int i = 0; i < 5; i++) iss >> temp[i];
			vector<int> thisconfig(temp, temp+5);
			iss >> value;
			Configuration newconfig(thisconfig, value);
			//cout << "Insert " << newconfig << endl;
			/*if (newconfig == newconfig) 
				cout << "==" << endl;
			else
				cout << "!=" << endl;*/
			trueValueTable.push_back(newconfig); 
		}
		/*for (vectorcit it = trueValueTable.begin(); it != trueValueTable.end(); it++)
			cout << *it << endl;*/
	}
	if (truevalue)
	{
		require(linelength == 3, "Wrong line length");
		vector<int> thisconfig(config, config+5);
		Configuration tobefound(thisconfig, 0);
		vectorcit it = find(trueValueTable.begin(), trueValueTable.end(), tobefound);
		if (it == trueValueTable.end())
		{
			cout << "To be found " << tobefound << endl;
		}
		require(it != trueValueTable.end(), "Error, cannot find this config in true value table");
		return (*it).getValue();
	}
	require(1 < linelength, "Error, less than 2 machines in the line");
	int id;
	//cout << endl << endl << "Config" << endl;
	//for (int i = 0; i < 2*linelength - 1; i++) cout << config[i] << ", ";
	//cout << endl;
	EventList eventList;
	// schedule the end event
	//eventList.schedule(new Event(END, TOTAL)); 
	int rawthroughput = 0;
	int totalthroughput = 0;
	double clock = 0;
	int* servicerate = new int[linelength];
	for (int i = 0; i < linelength; i++) servicerate[i] = config[i];
	int *buffersize = new int[linelength];
	buffersize[0] = 0;
	for (int j = linelength; j < 2*linelength - 1; j++) buffersize[j - linelength + 1] = config[j];
	// Now set up machines
	vector<Station*> stations;
	for (int i = 0; i < linelength; i++) stations.push_back(new Station(i, servicerate[i], buffersize[i]));
	// Calculate the first job release time and schedule the event
	double eventtime = Exponential(servicerate[0]);
	eventList.schedule(new JobCompletionEvent(0, clock + eventtime));
	//printpointercontainer(stations.begin(), stations.end(), "Stations");
	//for (int i = 0; i < linelength; i++) stations[i]->print();
	//while (rawthroughput < TOTAL)
	while (clock < TOTAL)
	{
		// event is the type of event that is gonna happen next, clock is its time
		unsigned int event = eventList.nextEvent();
		clock = eventList[event]->getTime();
		int eventType = eventList[event]->getType();
		switch (eventType)
		{	
			case JOBARRIVAL :
				// This is a job arrival event. 
				// Note that in a communication blocking protocol, a job can begin only when
				// the immediate downstream machine has at least one buffer space available. So
				// when this job is completed, it can at least enter the queue. But there is no
				// guarantee that the next stage machine is not locked, so need to check this
				// If the next stage station is not locked and not busy, begin 
				// processing the job and set the station busy, and schedule
				// the job completion event; else, put enque the job
				id = static_cast<JobArrivalEvent*>(eventList[event])->getStation();
				//cout << "Job arrival at # " << id << endl;
				require(stations[id]->bufferAvailable() || !stations[id]->isBusy(), "Error: no buffer space or station is not idle upon job arrival");
				require(0 < id, "Error, should not be a job arrival event for station 0"); 
				if (stations[id]->isBufferEmpty())
				{ // If the buffer is empty, two possibilities
					if (stations[id]->isBlocked() || stations[id]->isBusy())
					{ // If the machine is blocked, enque the job
						stations[id]->enqueue();
						// If the buffer is full now, block the immediate predecessor station.
						if (stations[id]->isBufferFull()) stations[id-1]->block();
					}
					else
					{ // Else start processing this job and schedule a job completion event
						stations[id]->setBusy();
						eventList.schedule(new JobCompletionEvent(id, clock + Exponential(servicerate[id])));
					}
				}
				else
				{	// The buffer is not empty and thus needs to enqueue
					stations[id]->enqueue();
				}
				// Check if the buffer is full or not; if not, unblock predecessor.
				// Note its predecessor blocks itself upon completion of this job, so no need
				// to block it explicitly here.
				if (!stations[id]->isBufferFull())
				{
					eventList.schedule(new UnblockEvent(id - 1, clock));
				}
				break;
			case JOBCOMPLETION :
				id = static_cast<JobArrivalEvent*>(eventList[event])->getStation();
				//cout << "Job completion at # " << id << endl;
				stations[id]->setIdle();
				if (id == linelength - 1) 
				{ // Last machine, count throughput after warmup. And schedule next job 
				  // completion event.
				  	rawthroughput++;
					//if (WARMUP < rawthroughput) ++totalthroughput;
					if (WARMUP < clock) ++totalthroughput;
					// If there is a job waiting in the buffer, do it.
					if (!stations[id]->isBufferEmpty()) 
					{
						stations[id]->dequeue();
						stations[id]->setBusy();
						eventList.schedule(new JobCompletionEvent(id, clock + Exponential(servicerate[id])));
					}
				 // Unblock predecessor if it is blocked since at least one new buffer space becomes available now.
				 	if (stations[id - 1]->isBlocked())
					 	eventList.schedule(new UnblockEvent(id - 1, clock));
				}
				else
				{ // schedule arrival to next station
					eventList.schedule(new JobArrivalEvent(id + 1, clock));
					// Block the station because after sending this job to next stage, the buffer
					// there might be full and thus blocking happens. So I choose to implement in 
					// this way: always block a station by itself and wait for unblocking from its
					// immediate successor, which is scheduled as soon as job arrival is processed.
					stations[id]->block();
				}
				break;
			case UNBLOCK :
				id = static_cast<UnblockEvent*>(eventList[event])->getStation();
				//cout << "Unblock # " << id << endl;
				require(stations[id]->isBlocked(), "Error, when unblcok happens, station is not blocked");
				stations[id]->unblock();
				require(!stations[id]->isBusy(), "Error, when unblcok happens, station is busy");
				// If this is not the first machine.
				// Check out if there is job waiting in the queue, and the station is not blocked,
				// , if so, continue working on the next job, and unblock its predecessor because now
				// there is at least one vacant buffer. If the station is blocked instead, do 
				// nothing and wait for an unblcok event.
				if (0 < id) 
				{
					if (!stations[id]->isBufferEmpty())
					{	// if the buffer is empty, do nth; else, do the job
						stations[id]->setBusy();
						stations[id]->dequeue();
						eventList.schedule(new JobCompletionEvent(id, clock + Exponential(servicerate[id])));
						// If predecessor is blocked, schedule an unblock event since this is not station 0, 
						// and now there is buffer avaialbe at this stage, and the unblcok event will happen right away
					}
					// No matter there is a job to do or not, there must be at least one buffer space at this
					// station that becomes available now, so unblock predecessor.
					if (stations[id - 1]->isBlocked())
						eventList.schedule(new UnblockEvent(id - 1, clock));
				}
				else
				{ // This is station 0, simply schedule another job completion if it is not blocked
					stations[id]->setBusy();
					eventList.schedule(new JobCompletionEvent(id, clock + Exponential(servicerate[id])));
				}			
				// If this station is unblocked and has job to do, do it
				break;
			default :	
				break;
		}
		eventList.removeEvent(event);
		//for (int i = 0; i < linelength; i++) stations[i]->print();
	}
	delete []servicerate;
	delete []buffersize;
	for (size_t i = 0; i < stations.size(); i++) delete stations[i];
	//cout << "Throughput " << totalthroughput / (TOTALTIME - WARMUP) << endl;
	//return totalthroughput / clock;
	return totalthroughput / (TOTAL - WARMUP);
}

int Flowline::isGlobalOptimum(const Solution& solution)
{
	vector<int> xgopt1;
	vector<int> xgopt2;
	xgopt1.push_back(6);
	xgopt1.push_back(7);
	xgopt1.push_back(7);
	xgopt1.push_back(12);
	xgopt2.push_back(7);
	xgopt2.push_back(7);
	xgopt2.push_back(6);
	xgopt2.push_back(8);
	Solution gopt1(xgopt1);
	Solution gopt2(xgopt2);
	if (solution.isEqual(&gopt1) || solution.isEqual(&gopt2))
		return 1;
	else 
		return 0;
}
