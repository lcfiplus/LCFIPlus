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
#include "TrackSelector.h"

#include "TROOT.h"
#include "TInterpreter.h"
#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TCut.h"
#include <string>
#include "TRandom3.h"
#include "flavtag.h"

#include "TH1F.h"

#include <sstream>

using namespace lcfiplus;
using namespace lcfiplus::algoSigProb;
using namespace lcfiplus::algoEtc;

namespace lcfiplus {

	class FtAuxiliary : public FTAlgo {
		private:
			int _aux;
		public:
			FtAuxiliary(const char *auxname, int auxval) : FTAlgo(auxname), _aux(auxval) {}
			void process() {
				_result = _aux;
			}
	};

	class FtNtrkWithoutV0 : public FTAlgo {
		public:
			FtNtrkWithoutV0() : FTAlgo("ntrkwithoutv0") {}
			void process() {
				_result = _jet->getAllTracks(true).size();
			}
	};

	class FtNtrk : public FTAlgo {
		public:
			FtNtrk() : FTAlgo("ntrk") {}
			void process() {
				_result = _jet->getAllTracks().size();
			}
	};

	class FtNvtxAll : public FTAlgo {
		public:
			FtNvtxAll() : FTAlgo("nvtxall") {}
			void process() {
				_result = _jet->getVertices().size();
			}
	};

	class FtVtxMassAll : public FTAlgo {
		public:
			FtVtxMassAll() : FTAlgo("vtxmassall") {}
			void process() {
				_result = 0;
				if (_jet->getVertices().size()>0) {
					TLorentzVector vtxp4;
					const VertexVec & vtxList = _jet->getVertices();
					for (unsigned int j=0; j<vtxList.size(); ++j) {
						vtxp4 += vtxList[j]->getFourMomentum();
					}
					_result = vtxp4.M();
				}
			}
	};

	class FtVtxLen12All : public FTAlgo {
		public:
			FtVtxLen12All() : FTAlgo("vtxlen12all") {}
			void process() {
				_result = 0;
				if (_jet->getVertices().size()>1) {
					_result = (_jet->getVertices()[1]->getPos() - _jet->getVertices()[0]->getPos()).Mag();
				}
			}
	};

	class FtVtxLen12AllByJetE : public FTAlgo {
		public:
			FtVtxLen12AllByJetE() : FTAlgo("vtxlen12all_jete") {}
			void process() {
				_result = 0;
				if (_jet->getVertices().size()>1) {
					_result = (_jet->getVertices()[1]->getPos() - _jet->getVertices()[0]->getPos()).Mag() / _jet->Energy();
				}
			}
	};

	class Ft1VtxProb : public FTAlgo {
		public:
			Ft1VtxProb() : FTAlgo("1vtxprob") {}
			void process() {
				_result = 0;
				const VertexVec & vtcs = _jet->getVertices();

				if (_jet->getVertices().size() == 1){
					_result = _jet->getVertices()[0]->getProb();
				}
				else if (_jet->getVertices().size()>=2) {
					if(_jet->params().count("RefinedVertex") > 0){
						_result = _jet->params().find("RefinedVertex")->second.get<double>("SingleVertexProbability");
					}else{

						vector<const Track *> tracks;
						for(unsigned int v=0;v<vtcs.size();v++){
							tracks.insert(tracks.end(), vtcs[v]->getTracks().begin(), vtcs[v]->getTracks().end());
						}
						// run vertex fitter
						Vertex *single = VertexFitterSimple_V()(tracks.begin(), tracks.end());
						_result = single->getProb();
						delete single;
					}
				}
			}
	};

	class FtNvtx : public FTAlgo {
		public:
			FtNvtx() : FTAlgo("nvtx") {}
			void process() {
				_result = _jet->getVerticesForFT().size();
			}
	};

	class FtJetE : public FTAlgo {
		public:
			FtJetE() : FTAlgo("jete") {}
			void process() {
				_result = _jet->Energy();
			}
	};

	class FtVtxLen1 : public FTAlgo {
		public:
			FtVtxLen1() : FTAlgo("vtxlen1") {}
			void process() {
				_result = 0;
				if (_jet->getVerticesForFT().size()>0)
					_result = _jet->getVerticesForFT()[0]->length( _privtx );
			}
	};

	class FtVtxLen2 : public FTAlgo {
		public:
			FtVtxLen2() : FTAlgo("vtxlen2") {}
			void process() {
				_result = 0;
				if (_jet->getVerticesForFT().size()>1) {
					_result = _jet->getVerticesForFT()[1]->length( _privtx );
				}
			}
	};

	class FtVtxLen12 : public FTAlgo {
		public:
			FtVtxLen12() : FTAlgo("vtxlen12") {}
			void process() {
				_result = 0;
				if (_jet->getVerticesForFT().size()>1)
					_result = _jet->getVerticesForFT()[1]->length( _jet->getVerticesForFT()[0] );
			}
	};

	class FtVtxLen1ByJetE : public FTAlgo {
		public:
			FtVtxLen1ByJetE() : FTAlgo("vtxlen1_jete") {}
			void process() {
				_result = 0;
				if (_jet->getVerticesForFT().size()>0)
					_result = _jet->getVerticesForFT()[0]->length( _privtx ) / _jet->Energy();
			}
	};

	class FtVtxLen2ByJetE : public FTAlgo {
		public:
			FtVtxLen2ByJetE() : FTAlgo("vtxlen2_jete") {}
			void process() {
				_result = 0;
				if (_jet->getVerticesForFT().size()>1)
					_result = _jet->getVerticesForFT()[1]->length( _privtx ) / _jet->Energy();
			}
	};

	class FtVtxLen12ByJetE : public FTAlgo {
		public:
			FtVtxLen12ByJetE() : FTAlgo("vtxlen12_jete") {}
			void process() {
				_result = 0;
				if (_jet->getVerticesForFT().size()>1)
					_result = _jet->getVerticesForFT()[1]->length( _jet->getVerticesForFT()[0] ) / _jet->Energy();
			}
	};

	class FtVtxSig1 : public FTAlgo {
		public:
			FtVtxSig1() : FTAlgo("vtxsig1") {}
			void process() {
				_result = 0;
				if (_jet->getVerticesForFT().size()>0)
					_result = _jet->getVerticesForFT()[0]->significance( _privtx );
			}
	};

	class FtVtxSig2 : public FTAlgo {
		public:
			FtVtxSig2() : FTAlgo("vtxsig2") {}
			void process() {
				_result = 0;
				if (_jet->getVerticesForFT().size()>1)
					_result = _jet->getVerticesForFT()[1]->significance( _privtx );
			}
	};

	class FtVtxSig12 : public FTAlgo {
		public:
			FtVtxSig12() : FTAlgo("vtxsig12") {}
			void process() {
				_result = 0;
				if (_jet->getVerticesForFT().size()>1)
					_result = _jet->getVerticesForFT()[1]->significance( _jet->getVerticesForFT()[0] );
			}
	};

	class FtVtxSig1ByJetE : public FTAlgo {
		public:
			FtVtxSig1ByJetE() : FTAlgo("vtxsig1_jete") {}
			void process() {
				_result = 0;
				if (_jet->getVerticesForFT().size()>0)
					_result = _jet->getVerticesForFT()[0]->significance( _privtx ) / _jet->Energy();
			}
	};

	class FtVtxSig2ByJetE : public FTAlgo {
		public:
			FtVtxSig2ByJetE() : FTAlgo("vtxsig2_jete") {}
			void process() {
				_result = 0;
				if (_jet->getVerticesForFT().size()>1)
					_result = _jet->getVerticesForFT()[1]->significance( _privtx ) / _jet->Energy();
			}
	};

	class FtVtxSig12ByJetE : public FTAlgo {
		public:
			FtVtxSig12ByJetE() : FTAlgo("vtxsig12_jete") {}
			void process() {
				_result = 0;
				if (_jet->getVerticesForFT().size()>1)
					_result = _jet->getVerticesForFT()[1]->significance( _jet->getVerticesForFT()[0] ) / _jet->Energy();
			}
	};

	class FtVtxDirAng1 : public FTAlgo {
		public:
			FtVtxDirAng1() : FTAlgo("vtxdirang1") {}
			void process() {
				_result = 0;
				if (_jet->getVerticesForFT().size()>0) {
					_result = _jet->getVerticesForFT()[0]->dirdot( _privtx );
					_result = acos(_result);
				}
			}
	};

	class FtVtxDirAng2 : public FTAlgo {
		public:
			FtVtxDirAng2() : FTAlgo("vtxdirang2") {}
			void process() {
				_result = 0;
				if (_jet->getVerticesForFT().size()>1) {
					_result = _jet->getVerticesForFT()[1]->dirdot( _privtx );
					_result = acos(_result);
				}
			}
	};

	class FtVtxDirAng12 : public FTAlgo {
		public:
			FtVtxDirAng12() : FTAlgo("vtxdirang12") {}
			void process() {
				_result = 0;
				if (_jet->getVerticesForFT().size()>1) {
					_result = _jet->getVerticesForFT()[1]->dirdot( _jet->getVerticesForFT()[0] );
					_result = acos(_result);
				}
			}
	};

	class FtVtxDirAng1TimesJetE : public FTAlgo {
		public:
			FtVtxDirAng1TimesJetE() : FTAlgo("vtxdirang1_jete") {}
			void process() {
				_result = 0;
				if (_jet->getVerticesForFT().size()>0) {
					_result = _jet->getVerticesForFT()[0]->dirdot( _privtx );
					_result = acos(_result)*_jet->Energy();
				}
			}
	};

	class FtVtxDirAng2TimesJetE : public FTAlgo {
		public:
			FtVtxDirAng2TimesJetE() : FTAlgo("vtxdirang2_jete") {}
			void process() {
				_result = 0;
				if (_jet->getVerticesForFT().size()>1) {
					_result = _jet->getVerticesForFT()[1]->dirdot( _privtx );
					_result = acos(_result)*_jet->Energy();
				}
			}
	};

	class FtVtxDirAng12TimesJetE : public FTAlgo {
		public:
			FtVtxDirAng12TimesJetE() : FTAlgo("vtxdirang12_jete") {}
			void process() {
				_result = 0;
				if (_jet->getVerticesForFT().size()>1) {
					_result = _jet->getVerticesForFT()[1]->dirdot( _jet->getVerticesForFT()[0] );
					_result = acos(_result)*_jet->Energy();
				}
			}
	};

	class FtVtxMom : public FTAlgo {
		public:
			FtVtxMom() : FTAlgo("vtxmom") {}
			void process() {
				_result = 0;
				TLorentzVector vtxp4;
				const vector<const Vertex*>& vtxList = _jet->getVerticesForFT();
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
				if (_jet->getVerticesForFT().size()>0)
					_result = _jet->getVerticesForFT()[0]->getFourMomentum().Vect().Mag();
			}
	};

	class FtVtxMom2 : public FTAlgo {
		public:
			FtVtxMom2() : FTAlgo("vtxmom2") {}
			void process() {
				_result = 0;
				if (_jet->getVerticesForFT().size()>1)
					_result = _jet->getVerticesForFT()[1]->getFourMomentum().Vect().Mag();
			}
	};

	class FtVtxMomByJetE : public FTAlgo {
		public:
			FtVtxMomByJetE() : FTAlgo("vtxmom_jete") {}
			void process() {
				_result = 0;
				TLorentzVector vtxp4;
				const vector<const Vertex*>& vtxList = _jet->getVerticesForFT();
				for (unsigned int j=0; j<vtxList.size(); ++j) {
					vtxp4 += vtxList[j]->getFourMomentum();
				}
				_result = vtxp4.Vect().Mag() / _jet->Energy();
			}
	};

	class FtVtxMom1ByJetE : public FTAlgo {
		public:
			FtVtxMom1ByJetE() : FTAlgo("vtxmom1_jete") {}
			void process() {
				_result = 0;
				if (_jet->getVerticesForFT().size()>0)
					_result = _jet->getVerticesForFT()[0]->getFourMomentum().Vect().Mag() / _jet->Energy();
			}
	};

	class FtVtxMom2ByJetE : public FTAlgo {
		public:
			FtVtxMom2ByJetE() : FTAlgo("vtxmom2_jete") {}
			void process() {
				_result = 0;
				if (_jet->getVerticesForFT().size()>1)
					_result = _jet->getVerticesForFT()[1]->getFourMomentum().Vect().Mag() / _jet->Energy();
			}
	};

	class FtVtxMassPtCorr : public FTAlgo {
		public:
			FtVtxMassPtCorr() : FTAlgo("vtxmasspc") {}
			void process() {
				_result = 0;
				if (_jet->getVerticesForFT().size()>0) {
					TLorentzVector vtxp4;
					const vector<const Vertex*>& vtxList = _jet->getVerticesForFT();
					for (unsigned int j=0; j<vtxList.size(); ++j) {
						vtxp4 += vtxList[j]->getFourMomentum();
					}

					LcfiInterface interface(_event,_privtx);
					double pt = interface.vertexMassPtCorrection(_jet->getVerticesForFT()[0],_privtx,vtxp4.Vect(),2);
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
				if (_jet->getVerticesForFT().size()>0) {
					double oneMinusProb = 1.;
					VertexVec & vtxList = _jet->getVerticesForFT();
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
				if (_jet->getVerticesForFT().size()>0) {
					TLorentzVector vtxp4;
					VertexVec & vtxList = _jet->getVerticesForFT();
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
				if (_jet->getVerticesForFT().size()>0)
					_result = _jet->getVerticesForFT()[0]->getFourMomentum().M();
			}
	};

	class FtVtxMass2 : public FTAlgo {
		public:
			FtVtxMass2() : FTAlgo("vtxmass2") {}
			void process() {
				_result = 0;
				if (_jet->getVerticesForFT().size()>1)
					_result = _jet->getVerticesForFT()[1]->getFourMomentum().M();
			}
	};

	class FtVtxMult : public FTAlgo {
		public:
			FtVtxMult() : FTAlgo("vtxmult") {}
			void process() {
				_result = 0;
				if (_jet->getVerticesForFT().size()>0) {
					VertexVec & vtxList = _jet->getVerticesForFT();
					VertexVec::const_iterator iter ;
					for (iter = vtxList.begin(); iter != vtxList.end(); ++iter) {
						_result += (*iter)->getTracks().size();
					}
				}
			}
	};

	class FtVtxMult1 : public FTAlgo {
		public:
			FtVtxMult1() : FTAlgo("vtxmult1") {}
			void process() {
				_result = 0;
				if (_jet->getVerticesForFT().size()>0)
					_result = _jet->getVerticesForFT()[0]->getTracks().size();
			}
	};

	class FtVtxMult2 : public FTAlgo {
		public:
			FtVtxMult2() : FTAlgo("vtxmult2") {}
			void process() {
				_result = 0;
				if (_jet->getVerticesForFT().size()>1)
					_result = _jet->getVerticesForFT()[1]->getTracks().size();
			}
	};

	class FtTrk1D0Sig : public FTAlgo {
		public:
			FtTrk1D0Sig() : FTAlgo("trk1d0sig") {}
			void process() {
				double sigVec[6];
				findMostSignificantTrack(_jet,_privtx,sigVec);
				_result = sigVec[0];
			}
	};

	class FtTrk2D0Sig : public FTAlgo {
		public:
			FtTrk2D0Sig() : FTAlgo("trk2d0sig") {}
			void process() {
				double sigVec[6];
				findMostSignificantTrack(_jet,_privtx,sigVec);
				_result = sigVec[1];
			}
	};

	class FtTrk1Z0Sig : public FTAlgo {
		public:
			FtTrk1Z0Sig() : FTAlgo("trk1z0sig") {}
			void process() {
				double sigVec[6];
				findMostSignificantTrack(_jet,_privtx,sigVec);
				_result = sigVec[2];
			}
	};

	class FtTrk2Z0Sig : public FTAlgo {
		public:
			FtTrk2Z0Sig() : FTAlgo("trk2z0sig") {}
			void process() {
				double sigVec[6];
				findMostSignificantTrack(_jet,_privtx,sigVec);
				_result = sigVec[3];
			}
	};

	class FtTrk1Pt : public FTAlgo {
		public:
			FtTrk1Pt() : FTAlgo("trk1pt") {}
			void process() {
				double sigVec[6];
				findMostSignificantTrack(_jet,_privtx,sigVec);
				_result = sigVec[4];
			}
	};

	class FtTrk2Pt : public FTAlgo {
		public:
			FtTrk2Pt() : FTAlgo("trk2pt") {}
			void process() {
				double sigVec[6];
				findMostSignificantTrack(_jet,_privtx,sigVec);
				_result = sigVec[5];
			}
	};
	
	class FtTrk1PtByJetE : public FTAlgo {
		public:
			FtTrk1PtByJetE() : FTAlgo("trk1pt_jete") {}
			void process() {
				double sigVec[6];
				findMostSignificantTrack(_jet,_privtx,sigVec);
				_result = sigVec[4] / _jet->Energy();
			}
	};

	class FtTrk2PtByJetE : public FTAlgo {
		public:
			FtTrk2PtByJetE() : FTAlgo("trk2pt_jete") {}
			void process() {
				double sigVec[6];
				findMostSignificantTrack(_jet,_privtx,sigVec);
				_result = sigVec[5] / _jet->Energy();
			}
	};
	
	class FtJProbR : public FTAlgo {
		public:
			FtJProbR() : FTAlgo("jprobr") {}
			void process() {
				_result = jointProbD0(_jet,_privtx);
			}
	};

	class FtJProbZ : public FTAlgo {
		public:
			FtJProbZ() : FTAlgo("jprobz") {}
			void process() {
				_result = jointProbZ0(_jet,_privtx);
			}
	};

	class FtJProbR5Sigma : public FTAlgo {
		public:
			FtJProbR5Sigma() : FTAlgo("jprobr5sigma") {}
			void process() {
				_result = jointProbD0(_jet,_privtx,5.);
			}
	};

	class FtJProbZ5Sigma : public FTAlgo {
		public:
			FtJProbZ5Sigma() : FTAlgo("jprobz5sigma") {}
			void process() {
				_result = jointProbZ0(_jet,_privtx,5.);
			}
	};

	class FtSphericity : public FTAlgo {
		public:
			FtSphericity() : FTAlgo("sphericity") {}
			void process() {
				_result = _jet->sphericity();
			}
	};

	class FtTrkMass : public FTAlgo {
		public:
			FtTrkMass() : FTAlgo("trkmass"){}
			void process() {
				vector<const Track *> tracks = _jet->getAllTracks(true);

				TrackSelectorConfig trsel;
				trsel.minD0Sig = 5.;
				trsel.minZ0Sig = 5.;
				trsel.maxD0 = 2.;
				trsel.maxZ0 = 3.;
				vector<const Track *> usedTracks = TrackSelector()(tracks, trsel);

				if(usedTracks.size() < 2)_result = 0.;
				else{
					TLorentzVector v;
					for(unsigned int n=0;n<usedTracks.size(); n++){
						v += *(TLorentzVector *)(usedTracks[n]);
					}
					_result = v.M();
				}
			}
	};

	class FtNMuon : public FTAlgo {
		public:
			FtNMuon() : FTAlgo("nmuon"){}
			void process(){
				TrackVec &tracks = _jet->getTracks(); // do not use vertexed tracks
//				TrackVec tracks = _jet->getAllTracks(true);
				_result = 0;
				for(unsigned int n=0; n<tracks.size();n++){
					if(algoEtc::SimpleSecMuonFinder(tracks[n], 5., 5., 5., -0.1, 0.2, 0.8, 1.5, 4., 0.5)
					|| algoEtc::SimpleSecMuonFinder(tracks[n], 5., 5., 5., 0.05, 0., 10., 0., 10., 10.))
						_result += 1;
				}
			}
	};

	class FtNElectron : public FTAlgo {
		public:
			FtNElectron() : FTAlgo("nelectron"){}
			void process(){
				TrackVec &tracks = _jet->getTracks(); // do not use vertexed tracks
//				TrackVec tracks = _jet->getAllTracks(true);
				_result = 0;
				for(unsigned int n=0; n<tracks.size();n++){
					if(algoEtc::SimpleSecElectronFinder(tracks[n], 5., 5., 5., 5., 0.98, 0.9, 1.15))
						_result += 1;
				}
			}
	};

	// historgram holder for d0/z0 probability
	class FtD0bProb;
	class FtD0cProb;
	class FtD0qProb;
	class FtD0bProbSigned;
	class FtD0cProbSigned;
	class FtD0qProbSigned;
	class FtZ0bProb;
	class FtZ0cProb;
	class FtZ0qProb;

	class FtIPProbHolder {
		friend class FtD0bProb;
		friend class FtD0cProb;
		friend class FtD0qProb;
		friend class FtD0bProbSigned;
		friend class FtD0cProbSigned;
		friend class FtD0qProbSigned;
		friend class FtZ0bProb;
		friend class FtZ0cProb;
		friend class FtZ0qProb;

		public: 
			FtIPProbHolder(const char *d0probfile, const char *z0probfile)
			{
				TDirectory *dir = gDirectory;
				_fd0 = TFile::Open(d0probfile);
				if(!_fd0)throw(Exception("FlavorTag: D0 probability file open error!"));
				_fz0 = TFile::Open(z0probfile);
				if(!_fz0)throw(Exception("FlavorTag: Z0 probability file open error!"));
				dir->cd();

				_hd0[0] = dynamic_cast<TH1F *>(_fd0->Get("hb"));
				_hd0[1] = dynamic_cast<TH1F *>(_fd0->Get("hc"));
				_hd0[2] = dynamic_cast<TH1F *>(_fd0->Get("hq"));

				_hd0p[0] = dynamic_cast<TH1F *>(_fd0->Get("hbp"));
				_hd0p[1] = dynamic_cast<TH1F *>(_fd0->Get("hcp"));
				_hd0p[2] = dynamic_cast<TH1F *>(_fd0->Get("hqp"));

				_hd0n[0] = dynamic_cast<TH1F *>(_fd0->Get("hbn"));
				_hd0n[1] = dynamic_cast<TH1F *>(_fd0->Get("hcn"));
				_hd0n[2] = dynamic_cast<TH1F *>(_fd0->Get("hqn"));

				if(!(_hd0[0] && _hd0[1] && _hd0[2]))throw(Exception("FlavorTag: D0 probability histogram open error!"));

				_hz0[0] = dynamic_cast<TH1F *>(_fz0->Get("hb"));
				_hz0[1] = dynamic_cast<TH1F *>(_fz0->Get("hc"));
				_hz0[2] = dynamic_cast<TH1F *>(_fz0->Get("hq"));

				if(!(_hz0[0] && _hz0[1] && _hz0[2]))throw(Exception("FlavorTag: D0 probability histogram open error!"));
			}

			~FtIPProbHolder(){
				if(_fd0)_fd0->Close();
				if(_fz0)_fz0->Close();
			}

		private:
			TFile *_fd0;
			TFile *_fz0;
			TH1F *_hd0[3];
			TH1F *_hd0p[3]; // signed-positive
			TH1F *_hd0n[3]; // signed-negative
			TH1F *_hz0[3];
	};

	class FtD0bProb : public FTAlgo {
		public:
			FtD0bProb(FtIPProbHolder *holder) : FTAlgo("d0bprob") {_holder = holder;}
			void process() {
				double prob = 1.;
				TrackVec tracks = _jet->getAllTracks(true);
				for(unsigned int n=0;n<tracks.size();n++){
					double d0sig = trackD0Significance(tracks[n], _privtx);
					if( d0sig > 5){
						double ld0 = log10(fabs(tracks[n]->getD0()));

						prob *= _holder->_hd0[0]->GetBinContent(_holder->_hd0[0]->FindBin(ld0)) * 3.; // 3 = b,c,q, averaging effect
					}
				}

				_result = prob;
			}
		private:
			FtIPProbHolder *_holder;
	};

	class FtD0cProb : public FTAlgo {
		public:
			FtD0cProb(FtIPProbHolder *holder) : FTAlgo("d0cprob") {_holder = holder;}
			void process() {
				double prob = 1.;
				TrackVec tracks = _jet->getAllTracks(true);
				for(unsigned int n=0;n<tracks.size();n++){
					double d0sig = trackD0Significance(tracks[n], _privtx);
					if( d0sig > 5){
						double ld0 = log10(fabs(tracks[n]->getD0()));

						prob *= _holder->_hd0[1]->GetBinContent(_holder->_hd0[1]->FindBin(ld0)) * 3.; // 3 = b,c,q, averaging effect
					}
				}

				_result = prob;
			}
		private:
			FtIPProbHolder *_holder;
	};

	class FtD0qProb : public FTAlgo {
		public:
			FtD0qProb(FtIPProbHolder *holder) : FTAlgo("d0qprob") {_holder = holder;}
			void process() {
				double prob = 1.;
				TrackVec tracks = _jet->getAllTracks(true);
				for(unsigned int n=0;n<tracks.size();n++){
					double d0sig = trackD0Significance(tracks[n], _privtx);
					if( d0sig > 5){
						double ld0 = log10(fabs(tracks[n]->getD0()));

						prob *= _holder->_hd0[2]->GetBinContent(_holder->_hd0[2]->FindBin(ld0)) * 3.; // 3 = b,c,q, averaging effect

//						cout << "d0 = " << tracks[n]->getD0() << ", pq = " <<  _holder->_hd0[2]->GetBinContent(_holder->_hd0[2]->FindBin(ld0)) * 3. << endl;
					}
				}

//				cout << "prob = " << prob << endl;
				_result = prob;
			}
		private:
			FtIPProbHolder *_holder;
	};

	class FtD0bProbSigned : public FTAlgo {
		public:
			FtD0bProbSigned(FtIPProbHolder *holder) : FTAlgo("d0bprobsigned") {_holder = holder;}
			void process() {
				double prob = 1.;
				TrackVec tracks = _jet->getAllTracks(true);
				for(unsigned int n=0;n<tracks.size();n++){
					double d0sig = trackD0Significance(tracks[n], _privtx);
					if( d0sig > 5){
						double sd0 = signedD0(tracks[n], _jet, _privtx, true);
						if(sd0>0){
							double ld0 = log10(sd0);
							prob *= _holder->_hd0p[0]->GetBinContent(_holder->_hd0p[0]->FindBin(ld0)) * 3.; // 3 = b,c,q, averaging effect
						}else{
							double ld0 = log10(-sd0);
							prob *= _holder->_hd0n[0]->GetBinContent(_holder->_hd0n[0]->FindBin(ld0)) * 3.; // 3 = b,c,q, averaging effect
						}
					}
				}

				_result = prob;
			}
		private:
			FtIPProbHolder *_holder;
	};

	class FtD0cProbSigned : public FTAlgo {
		public:
			FtD0cProbSigned(FtIPProbHolder *holder) : FTAlgo("d0cprobsigned") {_holder = holder;}
			void process() {
				double prob = 1.;
				TrackVec tracks = _jet->getAllTracks(true);
				for(unsigned int n=0;n<tracks.size();n++){
					double d0sig = trackD0Significance(tracks[n], _privtx);
					if( d0sig > 5){
						double sd0 = signedD0(tracks[n], _jet, _privtx, true);
						if(sd0>0){
							double ld0 = log10(sd0);
							prob *= _holder->_hd0p[1]->GetBinContent(_holder->_hd0p[1]->FindBin(ld0)) * 3.; // 3 = b,c,q, averaging effect
						}else{
							double ld0 = log10(-sd0);
							prob *= _holder->_hd0n[1]->GetBinContent(_holder->_hd0n[1]->FindBin(ld0)) * 3.; // 3 = b,c,q, averaging effect
						}
					}
				}

				_result = prob;
			}
		private:
			FtIPProbHolder *_holder;
	};

	class FtD0qProbSigned : public FTAlgo {
		public:
			FtD0qProbSigned(FtIPProbHolder *holder) : FTAlgo("d0qprobsigned") {_holder = holder;}
			void process() {
				double prob = 1.;
				TrackVec tracks = _jet->getAllTracks(true);
				for(unsigned int n=0;n<tracks.size();n++){
					double d0sig = trackD0Significance(tracks[n], _privtx);
					if( d0sig > 5){
						double sd0 = signedD0(tracks[n], _jet, _privtx, true);
						if(sd0>0){
							double ld0 = log10(sd0);
							prob *= _holder->_hd0p[2]->GetBinContent(_holder->_hd0p[2]->FindBin(ld0)) * 3.; // 3 = b,c,q, averaging effect
						}else{
							double ld0 = log10(-sd0);
							prob *= _holder->_hd0n[2]->GetBinContent(_holder->_hd0n[2]->FindBin(ld0)) * 3.; // 3 = b,c,q, averaging effect
						}
					}
				}

				_result = prob;
			}
		private:
			FtIPProbHolder *_holder;
	};

	class FtZ0bProb : public FTAlgo {
		public:
			FtZ0bProb(FtIPProbHolder *holder) : FTAlgo("z0bprob") {_holder = holder;}
			void process() {
				double prob = 1.;
				TrackVec tracks = _jet->getAllTracks(true);
				for(unsigned int n=0;n<tracks.size();n++){
					double z0sig = trackZ0Significance(tracks[n], _privtx);
					if( z0sig > 5){
						double lz0 = log10(fabs(tracks[n]->getZ0()));

						prob *= _holder->_hz0[0]->GetBinContent(_holder->_hz0[0]->FindBin(lz0)) * 3.; // 3 = b,c,q, averaging effect
					}
				}

				_result = prob;
			}
		private:
			FtIPProbHolder *_holder;
	};

	class FtZ0cProb : public FTAlgo {
		public:
			FtZ0cProb(FtIPProbHolder *holder) : FTAlgo("z0cprob") {_holder = holder;}
			void process() {
				double prob = 1.;
				TrackVec tracks = _jet->getAllTracks(true);
				for(unsigned int n=0;n<tracks.size();n++){
					double z0sig = trackZ0Significance(tracks[n], _privtx);
					if( z0sig > 5){
						double lz0 = log10(fabs(tracks[n]->getZ0()));

						prob *= _holder->_hz0[1]->GetBinContent(_holder->_hz0[1]->FindBin(lz0)) * 3.; // 3 = b,c,q, averaging effect
					}
				}

				_result = prob;
			}
		private:
			FtIPProbHolder *_holder;
	};

	class FtZ0qProb : public FTAlgo {
		public:
			FtZ0qProb(FtIPProbHolder *holder) : FTAlgo("z0qprob") {_holder = holder;}
			void process() {
				double prob = 1.;
				TrackVec tracks = _jet->getAllTracks(true);
				for(unsigned int n=0;n<tracks.size();n++){
					double z0sig = trackZ0Significance(tracks[n], _privtx);
					if( z0sig > 5){
						double lz0 = log10(fabs(tracks[n]->getZ0()));

						prob *= _holder->_hz0[2]->GetBinContent(_holder->_hz0[2]->FindBin(lz0)) * 3.; // 3 = b,c,q, averaging effect
					}
				}

				_result = prob;
			}
		private:
			FtIPProbHolder *_holder;
	};

	void FlavorTag::init(Parameters *param) {
		Algorithm::init(param);

		_primvtxcolname = param->get("PrimaryVertexCollectionName",string("PrimaryVertex"));
		_jetcolname = param->get("FlavorTag.JetCollectionName",string("VertexJets"));
		Event::Instance()->setDefaultPrimaryVertex(_primvtxcolname.c_str()); // backward compatibility

		_auxiliaryInfo = param->get("MakeNtuple.AuxiliaryInfo",int(-1));

		string d0probfilename = param->get("FlavorTag.D0ProbFileName",string("data/vtxprob/d0prob_zpole.root"));
		string z0probfilename = param->get("FlavorTag.Z0ProbFileName",string("data/vtxprob/z0prob_zpole.root"));

		_holder = new FtIPProbHolder(d0probfilename.c_str(), z0probfilename.c_str());

		//string outputFilename = param->get("TrainNtupleFile",string("lcfiplus.root"));
		//_nJet = (int)param->get("TrainNJet",float(2));

		//cout << "FlavorTag: Ntuple file set to " << outputFilename << endl;
		//cout << "FlavorTag: Number of jet set to " << _nJet << endl;

		FTManager& mgr = FTManager::getInstance();

		mgr.add( new FtAuxiliary("aux", _auxiliaryInfo) );

		mgr.add( new FtNtrkWithoutV0() );
		mgr.add( new FtNtrk() );

		mgr.add( new FtNvtxAll() );
		mgr.add( new FtVtxMassAll() );
		mgr.add( new FtVtxLen12All() );
		mgr.add( new FtVtxLen12AllByJetE() );
		mgr.add( new Ft1VtxProb() );

		mgr.add( new FtNvtx() );
		mgr.add( new FtJetE() );
		mgr.add( new FtVtxLen1() );
		mgr.add( new FtVtxLen2() );
		mgr.add( new FtVtxLen12() );
		mgr.add( new FtVtxLen1ByJetE() );
		mgr.add( new FtVtxLen2ByJetE() );
		mgr.add( new FtVtxLen12ByJetE() );
		mgr.add( new FtVtxSig1() );
		mgr.add( new FtVtxSig2() );
		mgr.add( new FtVtxSig12() );
		mgr.add( new FtVtxSig1ByJetE() );
		mgr.add( new FtVtxSig2ByJetE() );
		mgr.add( new FtVtxSig12ByJetE() );
		mgr.add( new FtVtxDirAng1() );
		mgr.add( new FtVtxDirAng2() );
		mgr.add( new FtVtxDirAng12() );
		mgr.add( new FtVtxDirAng1TimesJetE() );
		mgr.add( new FtVtxDirAng2TimesJetE() );
		mgr.add( new FtVtxDirAng12TimesJetE() );
		mgr.add( new FtVtxMom() );
		mgr.add( new FtVtxMom1() );
		mgr.add( new FtVtxMom2() );
		mgr.add( new FtVtxMomByJetE() );
		mgr.add( new FtVtxMom1ByJetE() );
		mgr.add( new FtVtxMom2ByJetE() );
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
		mgr.add( new FtTrk1PtByJetE() );
		mgr.add( new FtTrk2PtByJetE() );
		mgr.add( new FtJProbR() );
		mgr.add( new FtJProbZ() );
		mgr.add( new FtJProbR5Sigma() );
		mgr.add( new FtJProbZ5Sigma() );
		mgr.add( new FtSphericity() );
		mgr.add( new FtTrkMass() );
		mgr.add( new FtNMuon() );
		mgr.add( new FtNElectron() );

		// d0/z0 probs
		mgr.add( new FtD0bProb(_holder));
		mgr.add( new FtD0cProb(_holder));
		mgr.add( new FtD0qProb(_holder));
		mgr.add( new FtZ0bProb(_holder));
		mgr.add( new FtZ0cProb(_holder));
		mgr.add( new FtZ0qProb(_holder));

	}

	void FlavorTag::process() {

		Event *event = Event::Instance();
		//if (event->getTracks().size() == 0) return;

		const Vertex *privtx = event->getPrimaryVertex(_primvtxcolname.c_str());

		//TrackVec & tracks = event->getTracks();
		JetVec *jetsPtr(0);
		bool success = event->Get(_jetcolname.c_str(), jetsPtr);
		if (!success) {
			cout << "jets could not be found" << endl;
			return;
		}
		JetVec& jets = *jetsPtr;

		FTManager &mgr = FTManager::getInstance();
		mgr.process(event, privtx, jets);
	}

	void FlavorTag::end() {
		delete _holder;
	}
}
