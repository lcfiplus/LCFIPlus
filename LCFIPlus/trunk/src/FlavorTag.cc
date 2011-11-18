#include "FlavorTag.h"

#include <assert.h>
#include "EventStore.h"
#include "LcfiInterface.h"
#include "JetFinder.h"
#include "TreeStorer.h"
#include "VertexFitterLCFI.h"
#include "VertexFinderTearDown.h"
#include "algoSigProb.h"
#include "algoEtc.h"
#include "Driver.h"

#include "TROOT.h"
#include "TInterpreter.h"
#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TCut.h"
#include <string>
#include "TRandom3.h"
#include "flavtag.h"

#include <sstream>

using namespace lcfiplus;
using namespace lcfiplus::algoSigProb;
using namespace lcfiplus::algoEtc;

namespace lcfiplus {

	class FtNvtx : public FTAlgo {
		public:
			FtNvtx() : FTAlgo("nvtx") {}
			void process() {
				_result = _jet->getVertices().size();
			}
	};

	class FtVtxLen1 : public FTAlgo {
		public:
			FtVtxLen1() : FTAlgo("vtxlen1") {}
			void process() {
				_result = 0;
				if (_jet->getVertices().size()>0)
					_result = _jet->getVertices()[0]->length( _event->getPrimaryVertex() );
			}
	};

	class FtVtxLen2 : public FTAlgo {
		public:
			FtVtxLen2() : FTAlgo("vtxlen2") {}
			void process() {
				_result = 0;
				if (_jet->getVertices().size()>1)
					_result = _jet->getVertices()[1]->length( _event->getPrimaryVertex() );
			}
	};

	class FtVtxLen12 : public FTAlgo {
		public:
			FtVtxLen12() : FTAlgo("vtxlen12") {}
			void process() {
				_result = 0;
				if (_jet->getVertices().size()>1)
					_result = _jet->getVertices()[1]->length( _jet->getVertices()[0] );
			}
	};

	class FtVtxSig1 : public FTAlgo {
		public:
			FtVtxSig1() : FTAlgo("vtxsig1") {}
			void process() {
				_result = 0;
				if (_jet->getVertices().size()>0)
					_result = _jet->getVertices()[0]->significance( _event->getPrimaryVertex() );
			}
	};

	class FtVtxSig2 : public FTAlgo {
		public:
			FtVtxSig2() : FTAlgo("vtxsig2") {}
			void process() {
				_result = 0;
				if (_jet->getVertices().size()>1)
					_result = _jet->getVertices()[1]->significance( _event->getPrimaryVertex() );
			}
	};

	class FtVtxSig12 : public FTAlgo {
		public:
			FtVtxSig12() : FTAlgo("vtxsig12") {}
			void process() {
				_result = 0;
				if (_jet->getVertices().size()>1)
					_result = _jet->getVertices()[1]->significance( _jet->getVertices()[0] );
			}
	};

	class FtVtxDirDot1 : public FTAlgo {
		public:
			static const char* name() { return "vtxdirdot1"; }
			FtVtxDirDot1() : FTAlgo("vtxdirdot1") {}
			void process() {
				_result = 0;
				if (_jet->getVertices().size()>0)
					_result = _jet->getVertices()[0]->dirdot( _event->getPrimaryVertex() );
			}
	};

	class FtVtxDirDot2 : public FTAlgo {
		public:
			FtVtxDirDot2() : FTAlgo("vtxdirdot1") {}
			void process() {
				_result = 0;
				if (_jet->getVertices().size()>1)
					_result = _jet->getVertices()[1]->dirdot( _event->getPrimaryVertex() );
			}
	};

	class FtVtxDirDot12 : public FTAlgo {
		public:
			FtVtxDirDot12() : FTAlgo("vtxdirdot12") {}
			void process() {
				_result = 0;
				if (_jet->getVertices().size()>1)
					_result = _jet->getVertices()[1]->dirdot( _jet->getVertices()[0] );
			}
	};


	class FtVtxMom : public FTAlgo {
		public:
			FtVtxMom() : FTAlgo("vtxmom") {}
			void process() {
				_result = 0;
				TLorentzVector vtxp4;
				const vector<const Vertex*>& vtxList = _jet->getVertices();
				for (unsigned int j=0; j<vtxList.size(); ++j) {
					vtxp4 += vtxList[j]->getFourMomentum();
				}
				_result = vtxp4.Vect().Mag();
			}
	};

	class FtVtxMom1 : public FTAlgo {
		public:
			FtVtxMom1() : FTAlgo("vtxmom1") {}
			void process() {
				_result = 0;
				if (_jet->getVertices().size()>0)
					_result = _jet->getVertices()[0]->getFourMomentum().Vect().Mag();
			}
	};

	class FtVtxMom2 : public FTAlgo {
		public:
			FtVtxMom2() : FTAlgo("vtxmom2") {}
			void process() {
				_result = 0;
				if (_jet->getVertices().size()>1)
					_result = _jet->getVertices()[1]->getFourMomentum().Vect().Mag();
			}
	};

	class FtVtxMassPtCorr : public FTAlgo {
		public:
			FtVtxMassPtCorr() : FTAlgo("vtxmasspc") {}
			void process() {
				_result = 0;
				if (_jet->getVertices().size()>0) {
					TLorentzVector vtxp4;
					const vector<const Vertex*>& vtxList = _jet->getVertices();
					for (unsigned int j=0; j<vtxList.size(); ++j) {
						vtxp4 += vtxList[j]->getFourMomentum();
					}

					LcfiInterface interface(_event,_event->getPrimaryVertex());
					double pt = interface.vertexMassPtCorrection(_jet->getVertices()[0],_event->getPrimaryVertex(),vtxp4.Vect(),2);
					double vm = vtxp4.M();
					_result = sqrt( vm*vm+pt*pt ) + pt;
				}
			}
	};

	class FtVtxProb : public FTAlgo {
		public:
			FtVtxProb() : FTAlgo("vtxprob") {}
			void process() {
				_result = 0;
				if (_jet->getVertices().size()>0) {
					double oneMinusProb = 1.;
					VertexVec & vtxList = _jet->getVertices();
					for (unsigned int j=0; j<vtxList.size(); ++j) {
						TrackVec & vtxTracks = vtxList[j]->getTracks();
						int ndf = 2*vtxTracks.size()-3;
						double prob = TMath::Prob(vtxList[j]->getChi2(),ndf);
						oneMinusProb *= (1-prob);
					}
					_result = 1-oneMinusProb;
				}
			}
	};

	class FtVtxMass : public FTAlgo {
		public:
			FtVtxMass() : FTAlgo("vtxmass") {}
			void process() {
				_result = 0;
				if (_jet->getVertices().size()>0) {
					TLorentzVector vtxp4;
					VertexVec & vtxList = _jet->getVertices();
					for (unsigned int j=0; j<vtxList.size(); ++j) {
						vtxp4 += vtxList[j]->getFourMomentum();
					}
					_result = vtxp4.M();
				}
			}
	};

	class FtVtxMass1 : public FTAlgo {
		public:
			FtVtxMass1() : FTAlgo("vtxmass1") {}
			void process() {
				_result = 0;
				if (_jet->getVertices().size()>0)
					_result = _jet->getVertices()[0]->getFourMomentum().M();
			}
	};

	class FtVtxMass2 : public FTAlgo {
		public:
			FtVtxMass2() : FTAlgo("vtxmass2") {}
			void process() {
				_result = 0;
				if (_jet->getVertices().size()>1)
					_result = _jet->getVertices()[1]->getFourMomentum().M();
			}
	};

	class FtVtxMult : public FTAlgo {
		public:
			FtVtxMult() : FTAlgo("vtxmult") {}
			void process() {
				_result = 0;
				if (_jet->getVertices().size()>0) {
					VertexVec & vtxList = _jet->getVertices();
					_result += vtxList.size();
				}
			}
	};

	class FtVtxMult1 : public FTAlgo {
		public:
			FtVtxMult1() : FTAlgo("vtxmult1") {}
			void process() {
				_result = 0;
				if (_jet->getVertices().size()>0)
					_result = _jet->getVertices()[0]->getTracks().size();
			}
	};

	class FtVtxMult2 : public FTAlgo {
		public:
			FtVtxMult2() : FTAlgo("vtxmult2") {}
			void process() {
				_result = 0;
				if (_jet->getVertices().size()>1)
					_result = _jet->getVertices()[1]->getTracks().size();
			}
	};

	class FtTrk1D0Sig : public FTAlgo {
		public:
			FtTrk1D0Sig() : FTAlgo("trk1d0sig") {}
			void process() {
				float sigVec[6];
				findMostSignificantTrack(_jet,_event->getPrimaryVertex(),sigVec);
				_result = sigVec[0];
			}
	};

	class FtTrk2D0Sig : public FTAlgo {
		public:
			FtTrk2D0Sig() : FTAlgo("trk2d0sig") {}
			void process() {
				float sigVec[6];
				findMostSignificantTrack(_jet,_event->getPrimaryVertex(),sigVec);
				_result = sigVec[1];
			}
	};

	class FtTrk1Z0Sig : public FTAlgo {
		public:
			FtTrk1Z0Sig() : FTAlgo("trk1z0sig") {}
			void process() {
				float sigVec[6];
				findMostSignificantTrack(_jet,_event->getPrimaryVertex(),sigVec);
				_result = sigVec[2];
			}
	};

	class FtTrk2Z0Sig : public FTAlgo {
		public:
			FtTrk2Z0Sig() : FTAlgo("trk2z0sig") {}
			void process() {
				float sigVec[6];
				findMostSignificantTrack(_jet,_event->getPrimaryVertex(),sigVec);
				_result = sigVec[3];
			}
	};

	class FtTrk1Pt : public FTAlgo {
		public:
			FtTrk1Pt() : FTAlgo("trk1pt") {}
			void process() {
				float sigVec[6];
				findMostSignificantTrack(_jet,_event->getPrimaryVertex(),sigVec);
				_result = sigVec[4];
			}
	};

	class FtTrk2Pt : public FTAlgo {
		public:
			FtTrk2Pt() : FTAlgo("trk2pt") {}
			void process() {
				float sigVec[6];
				findMostSignificantTrack(_jet,_event->getPrimaryVertex(),sigVec);
				_result = sigVec[5];
			}
	};
	
	class FtJProbR : public FTAlgo {
		public:
			FtJProbR() : FTAlgo("jprobr") {}
			void process() {
				_result = jointProbD0(_jet,_event->getPrimaryVertex());
			}
	};

	class FtJProbZ : public FTAlgo {
		public:
			FtJProbZ() : FTAlgo("jprobz") {}
			void process() {
				_result = jointProbZ0(_jet,_event->getPrimaryVertex());
			}
	};

	class FtSphericity : public FTAlgo {
		public:
			FtSphericity() : FTAlgo("sphericity") {}
			void process() {
				_result = _jet->sphericity();
			}
	};

	void FlavorTag::init(Parameters *param) {
		Algorithm::init(param);

		string outputFilename = param->get("TrainNtupleFile",string("lcfiplus.root"));
		_nJet = (int)param->get("TrainNJet",float(2));

		cout << "FlavorTag: Ntuple file set to " << outputFilename << endl;
		cout << "FlavorTag: Number of jet set to " << _nJet << endl;

		FTManager& mgr = FTManager::getInstance();

		mgr.add( new FtNvtx() );
		mgr.add( new FtVtxLen1() );
		mgr.add( new FtVtxSig1() );
		mgr.add( new FtVtxLen2() );
		mgr.add( new FtVtxSig2() );
		mgr.add( new FtVtxLen12() );
		mgr.add( new FtVtxSig12() );
		mgr.add( new FtVtxDirDot1() );
		mgr.add( new FtVtxDirDot2() );
		mgr.add( new FtVtxDirDot12() );
		mgr.add( new FtVtxMom() );
		mgr.add( new FtVtxMom1() );
		mgr.add( new FtVtxMom2() );
		mgr.add( new FtVtxMass() );
		mgr.add( new FtVtxMass1() );
		mgr.add( new FtVtxMass2() );
		mgr.add( new FtVtxMassPtCorr() );
		mgr.add( new FtVtxMult() );
		mgr.add( new FtVtxMult1() );
		mgr.add( new FtVtxMult2() );
		mgr.add( new FtVtxProb() );
		mgr.add( new FtTrk1D0Sig() );
		mgr.add( new FtTrk2D0Sig() );
		mgr.add( new FtTrk1Z0Sig() );
		mgr.add( new FtTrk2Z0Sig() );
		mgr.add( new FtTrk1Pt() );
		mgr.add( new FtTrk2Pt() );
		mgr.add( new FtJProbR() );
		mgr.add( new FtJProbZ() );
		mgr.add( new FtSphericity() );
	}

	void FlavorTag::process() {

		Event *event = Event::Instance();
		if (event->getTracks().size() == 0) return;

		//TrackVec & tracks = event->getTracks();
		JetVec *jetsPtr(0);
		bool success = event->Get("VertexJets", jetsPtr);
		if (!success) {
			cout << "jets could not be found" << endl;
			return;
		}
		JetVec& jets = *jetsPtr;

		FTManager &mgr = FTManager::getInstance();
		mgr.process(event, jets);
	}

	void FlavorTag::end() {
	}
}
