#ifndef EVENT_NAV_
#define EVENT_NAV_ 1

#include <string>
//#include "TQObject.h"


namespace lcfiplus {

class Event;
class LCIOStorer;

/**
	Event display for LCFIPlus.

	The graphics is based on the TEve framework for ROOT.
	Adapted from the <a href="http://llr.in2p3.fr/~ruan/ILDDisplay/>Druid event display</a>.

	@author T. Tanabe, ICEPP, The University of Tokyo
	@version $Id$
	*/
class EventNavigator {

 public:
  /**
  	Constructor.

  	@param[in] input input file which contains events to inspect.
  	@param[in] start skips given number of events.
  	*/
  EventNavigator(const char* input, int start=0);
  ~EventNavigator();
  /**
  	Display next event.
  	*/
  void Fwd();
  void Bck();
  /**
  	Draws given event.
  	*/
  void drawEvent(Event* event);

 private:
  LCIOStorer* _ls;
  int _start;
};
}


#endif //EVENT_NAV_

