// FlavorTag.h

#ifndef FlavorTag_h
#define FlavorTag_h 1

class TFile;
class TTree;

#include "interface.h"
#include "lcfiplus.h"

namespace lcfiplus{

	class FlavorTag : public LcfiplusAlgorithm {
		public:
			FlavorTag(){}
			virtual ~FlavorTag(){}

			void init(LcfiplusParameters *param);
			void process();
			void end();

		private:
			TFile* _file;
			TTree* _t;
			int _nJet;

			ClassDef(FlavorTag,1);
	};

}

#endif
