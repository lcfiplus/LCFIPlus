#ifndef EVENT_NAV_
#define EVENT_NAV_ 1

#include <string>
//#include "TQObject.h"


namespace flavtag {

  class Event;
	class LCIOStorer;

  //class EventNavigator : protected TQObject {
  class EventNavigator {
#ifndef NO_EVE
    public:
      EventNavigator(const char* input, int start=0);
			~EventNavigator();
      void Fwd();
      void Bck();
      void drawEvent(Event* event);

    private:
			LCIOStorer* _ls;
			int _start;
#endif
  };
}


#endif //EVENT_NAV_

