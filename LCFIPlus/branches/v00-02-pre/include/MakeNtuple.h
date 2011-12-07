// MakeNtuple.h

#ifndef MakeNtuple_h
#define MakeNtuple_h 1

class TFile;
class TTree;

#include "lcfiplus.h"

namespace lcfiplus{

	class MakeNtuple : public Algorithm {
		public:
			MakeNtuple(){}
			virtual ~MakeNtuple(){}

			void init(Parameters *param);
			void process() {}
			void end();

		private:

			ClassDef(MakeNtuple,1);
	};

}

#endif
