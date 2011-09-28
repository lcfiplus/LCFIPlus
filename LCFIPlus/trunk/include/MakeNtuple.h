// MakeNtuple.h

#ifndef MakeNtuple_h
#define MakeNtuple_h 1

class TFile;
class TTree;

#include "interface.h"
#include "lcfiplus.h"

namespace lcfiplus{

	class MakeNtuple : public LcfiplusAlgorithm {
		public:
			MakeNtuple(){}
			virtual ~MakeNtuple(){}

			void init(LcfiplusParameters *param);
			void process() {}
			void end();

		private:

			ClassDef(MakeNtuple,1);
	};

}

#endif
