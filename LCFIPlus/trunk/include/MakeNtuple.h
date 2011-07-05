// MakeNtuple.h

#ifndef MakeNtuple_h
#define MakeNtuple_h 1

class TFile;
class TTree;

#include "interface.h"
#include "flavtag.h"

namespace flavtag{

	struct FlavtagData {
		double mcnb;
		double mcnc;
		double flavor;
		double px;
		double py;
		double pz;
		double e;
		double nvtx;
		double vtxlen1;
		double vtxsig1;
		double vtxdirdot1;
		double vtxlen2;
		double vtxsig2;
		double vtxdirdot2;
		double vtxlen12;
		double vtxsig12;
		double vtxdirdot12;
		double vtxmom;
		double vtxmass;
		double vtxmasspc;
		double vtxmult;
		double vtxprob;
		double vtxmom1;
		double vtxmass1;
		double vtxmult1;
		double vtxmom2;
		double vtxmass2;
		double vtxmult2;
		double trk1d0sig;
		double trk2d0sig;
		double trk1z0sig;
		double trk2z0sig;
		double trk1pt;
		double trk2pt;
		double jprobr;
		double jprobz;
		double sphericity;
	};

	class MakeNtuple : public FlavtagAlgorithm {
		public:
			MakeNtuple(){}
			virtual ~MakeNtuple(){}

			void init(FlavtagParameters *param);
			void process();
			void end();

		private:
			TFile* _file;
			TTree* _t;

			FlavtagData _data;
			ClassDef(MakeNtuple,1);
	};

}

#endif
