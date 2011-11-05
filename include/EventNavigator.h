#ifndef EVENT_NAV_
#define EVENT_NAV_ 1

#include <string>
//#include "TQObject.h"


namespace lcfiplus {

  class Event;
	class LCIOStorer;

  //class EventNavigator : protected TQObject {
  class EventNavigator {
#ifdef BUILD_EVE
    public:
      EventNavigator(const char* input, int start=0);
			~EventNavigator();
      void Fwd();
      void Bck();
      void drawEvent(Event* event);

    private:
			LCIOStorer* _ls;
			int _start;
#else
    public:
      EventNavigator(const char* input, int start=0) {}
			~EventNavigator() {}

#endif
  };
}


#endif //EVENT_NAV_

