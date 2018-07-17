#ifndef DESS_H_
#define DESS_H_
#include <vector>
#include <iostream>
// Define commonly used classes for discrete event system simultation
const double INFINITE = 1e30;
const int END = 0;
class Event
{
	private:
		int _type;
		double _time;
	public:
		Event(int type = END, double time = INFINITE) : _type(type), _time(time) {}
		Event(const Event& event) {std::cout << "Copy Event" << std::endl;}
		const int getType() const {return _type;}
		const double getTime() const {return _time;}
		void setType(const int& type) {_type = type;}
		void setTime(const double& time) {_time = time;}
};

class EventList
{
	private:
		std::vector<Event*> _eventList;
	public:
		EventList();
		~EventList();
		unsigned int nextEvent() const;
		void removeEvent(unsigned int event);
		void clearList();
		void schedule(Event* event);
		Event* operator[](unsigned int i) {return _eventList[i];}
}; 
		
#endif /*DESS_H_*/
