#ifndef FLOWLINE_H_
#define FLOWLINE_H_
/* The flow line example from "Stochastic Models of Manufacturing Systems"
 * by John A. Buzacott, p 189.
 * Communication blocking protocol used
 */
#include "DES.h"
#include "RandomVariates.h"
#include "RngStream.h"
#include "require.h"
#include "myutilities.h"
#include "Simulator.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>

using std::vector;

// linelength: number of stations on this line
// config: solution array, the first linelength many elements are processing rate of station[0] 
// to station[linelength-1], followed by buffer sizes for station[1] to station[linelength]
// There are altogether TOTAL many units produced, and the first WARMUP many units discarded in
// calculating throughput.
class Flowline : public Simulator
{
	private :
		double flowline(int linelength, const int* const config, const double WARMUP, const double TOTAL, bool truevalue);
	public :
		Flowline() { setSimMode(serial); }
		double simulation(const Solution* const solution);
		double gettruevalue(const Solution* const solution);
		int isGlobalOptimum(const Solution& solution);
};


class JobArrivalEvent : public Event
{
	private :
		int _station;
	public :
		JobArrivalEvent(int station, double time);
		const int getStation() const {return _station;}
};
	
class JobCompletionEvent : public Event
{
	int _station;
	public:
		JobCompletionEvent(int station, double time);
		const int getStation() const {return _station;}
};

class UnblockEvent : public Event
{
	int _station;
	public:
		UnblockEvent(int station, double time);
		const int getStation() const {return _station;}
};

class Configuration
{
	private :
		std::vector<int> _config;
		double _value;
	public :
		Configuration(std::vector<int>& data, double value)
		{
			_config.reserve(data.size());
			_config = data;
			_value = value;
		}
		/*Configuration(std::vector<int>& data)
		{
			_config.reserve(data.size());
			_config = data;
		}*/
		const std::vector<int>& getConfig() const {return _config;}
		double getValue() const {return _value;}
		friend bool operator==(const Configuration& config1, const Configuration& config2)
		{
			const std::vector<int> data1 = config1.getConfig();
			const std::vector<int> data2 = config2.getConfig();
			if (data1.size() != data2.size()) return false;
			for (size_t i = 0; i < data1.size(); i++)
			{
				if (data1[i] != data2[i]) return false;
			}
			return true;
		}
		friend bool operator<(const Configuration& config1, const Configuration& config2)
		{
			return config1.getValue() < config2.getValue();
			//return true;
		}
		friend std::ostream& operator<<(std::ostream& os, const Configuration& config)
		{
			for (size_t i = 0; i < (config.getConfig()).size(); i++) 
				os << (config.getConfig())[i] << " ";
			return os << " " << config.getValue();
		}
};
	
class Station
{
	private :
		int _id;
		const int _servicerate;
		const int _buffersize;
		bool _blocked;
		bool _busy;
		int _numOfJobsInBuffer;
	public :
		Station(int id, int servicerate, int buffersize) : _id(id), _servicerate(servicerate), \
			_buffersize(buffersize), _blocked(false), _busy(false), _numOfJobsInBuffer(0) {}
		const bool isBlocked() const {return _blocked;}
		void block() {_blocked = true;}
		void unblock() {_blocked = false;}	
		const bool isBusy() const {return _busy;}
		void setBusy() {_busy = true;}
		void setIdle() {_busy = false;}
		const bool bufferAvailable() const {return _numOfJobsInBuffer < _buffersize;} 
		const bool isBufferEmpty() const {return _numOfJobsInBuffer == 0;}
		const bool isBufferFull() const {return _numOfJobsInBuffer == _buffersize;}
		void enqueue() 
		{
			require(_numOfJobsInBuffer < _buffersize, "Error, buffer full while enqueuing");
			++_numOfJobsInBuffer;
		}
		void dequeue() 
		{
			if (0 < _numOfJobsInBuffer)
			{
				_busy = true;
				--_numOfJobsInBuffer;
			}
		}
		friend std::ostream& operator<<(std::ostream& os, const Station* station);
		void print()
		{
			if (_busy) 
				std::cout << "Station " << _id << " busy";
			else
				std::cout << "Station " << _id << " idle";
			std::cout << "# jobs in buffer " << _numOfJobsInBuffer << std::endl;
		}
};

struct StationGen
{
	Station operator()(int id, int servicerate, int buffersize) {return Station(id, servicerate, buffersize);}
};
#endif /*FLOWLINE_H_*/
