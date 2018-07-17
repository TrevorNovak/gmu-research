#include "DES.h"
#include <iostream>
#include <limits>
using namespace std;

EventList::EventList() 
{
	_eventList.clear();
	// when initializing event list, push an END event so that no dead loop will occur in simulation
	_eventList.push_back(new Event());
}

EventList::~EventList() 
{
	for (unsigned int i = 0; i < _eventList.size(); i++)
	{
		delete _eventList[i];
	}
}
		
unsigned int EventList::nextEvent() const
{
	unsigned int index = 0;
	double minimum = numeric_limits<double>::max();
	for (unsigned int i = 0; i < _eventList.size(); i++)
	{
		if ( _eventList[i]->getTime() < minimum)
		{
			minimum =  _eventList[i]->getTime();
			index = i;
		}
	}
	return index;
}
		
void EventList::removeEvent(unsigned int event)
{
	vector<Event*> tempList;
	if (event >= _eventList.size())
	{
		cout << "Invalid event index " << endl;
	}
	else
	{
		for (unsigned int i = 0; i < _eventList.size(); i++)
		{
			if (i != event)
			{
				tempList.push_back(_eventList[i]);
			}
			else
			{
				// Remove this event
				delete _eventList[i];
			}
		}
		_eventList.clear();
		for (unsigned int i = 0; i < tempList.size(); i++)
		{
			_eventList.push_back(tempList[i]);
		}
	}
}

void EventList::clearList()
{
	_eventList.clear();
}
		
void EventList::schedule(Event* event)
{
	_eventList.push_back(event);
} 




