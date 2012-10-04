// FlavorTag.h

#ifndef FlavorTag_h
#define FlavorTag_h 1

class TFile;
class TTree;

#include "lcfiplus.h"

namespace lcfiplus{

	class FtIPProbHolder;

	/**
		Controls the event data and registers and holds algorithms for
		the computation of flavor tagging variables.

		This algorithm must be included in the
		algorithm tag of the XML steering file
		when using any of the flavor tagging features,
		e.g. ntuple production, TMVA training, and evaluation.

		@author T. Tanabe, ICEPP, The University of Tokyo
		@version $Id$
	 */
	class FlavorTag : public Algorithm {
		public:
			FlavorTag(){}
			virtual ~FlavorTag(){}

			void init(Parameters *param);
			void process();
			void end();

		private:
			TFile* _file;
			TTree* _t;

			int _auxiliaryInfo;
			string _jetcolname;
			string _primvtxcolname;

			FtIPProbHolder *_holder; //!

			int _nhitsJointProbD0;
			int _nhitsJointProbZ0;
			int _nhitsMostSignificantTrack;

			ClassDef(FlavorTag,1);
	};

}

#endif
