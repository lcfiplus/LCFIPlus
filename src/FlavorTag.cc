#include "FlavorTag.h"

#include <assert.h>
#include "EventStore.h"
#include "LcfiInterface.h"
#include "JetFinder.h"
#include "TreeStorer.h"
#include "VertexFitterLCFI.h"
#include "VertexFinderTearDown.h"
#include "VertexFinderPerfect.h"
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
#include "geometry.h"

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
  FtAuxiliary(const char* auxname, int auxval) : FTAlgo(auxname), _aux(auxval) {}
  void process() {
    _result = _aux;
  }
};

class FtAuxiliaryM : public FTAlgo {
 public:
  FtAuxiliaryM(const char* auxname) : FTAlgo(auxname) {}
  void process() {
    _result = FTManager::getInstance().getAuxiliary();
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
      const VertexVec& vtxList = _jet->getVertices();
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
    const VertexVec& vtcs = _jet->getVertices();

    if (_jet->getVertices().size() == 1) {
      _result = _jet->getVertices()[0]->getProb();
    } else if (_jet->getVertices().size()>=2) {
      if (_jet->params().count("RefinedVertex") > 0) {
        _result = _jet->params().find("RefinedVertex")->second.get<double>("SingleVertexProbability");
      } else {

        vector<const Track*> tracks;
        for (unsigned int v=0; v<vtcs.size(); v++) {
          tracks.insert(tracks.end(), vtcs[v]->getTracks().begin(), vtcs[v]->getTracks().end());
        }
        // run vertex fitter
        Vertex* single = VertexFitterSimple_V()(tracks.begin(), tracks.end());
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
      VertexVec& vtxList = _jet->getVerticesForFT();
      for (unsigned int j=0; j<vtxList.size(); ++j) {
        TrackVec& vtxTracks = vtxList[j]->getTracks();
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
      VertexVec& vtxList = _jet->getVerticesForFT();
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
      VertexVec& vtxList = _jet->getVerticesForFT();
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
    findMostSignificantTrack(_jet,_privtx,_nhitsMostSignificantTrack,sigVec);
    _result = sigVec[0];
  }
};

class FtTrk2D0Sig : public FTAlgo {
 public:
  FtTrk2D0Sig() : FTAlgo("trk2d0sig") {}
  void process() {
    double sigVec[6];
    findMostSignificantTrack(_jet,_privtx,_nhitsMostSignificantTrack,sigVec);
    _result = sigVec[1];
  }
};

class FtTrk1Z0Sig : public FTAlgo {
 public:
  FtTrk1Z0Sig() : FTAlgo("trk1z0sig") {}
  void process() {
    double sigVec[6];
    findMostSignificantTrack(_jet,_privtx,_nhitsMostSignificantTrack,sigVec);
    _result = sigVec[2];
  }
};

class FtTrk2Z0Sig : public FTAlgo {
 public:
  FtTrk2Z0Sig() : FTAlgo("trk2z0sig") {}
  void process() {
    double sigVec[6];
    findMostSignificantTrack(_jet,_privtx,_nhitsMostSignificantTrack,sigVec);
    _result = sigVec[3];
  }
};

class FtTrk1Pt : public FTAlgo {
 public:
  FtTrk1Pt() : FTAlgo("trk1pt") {}
  void process() {
    double sigVec[6];
    findMostSignificantTrack(_jet,_privtx,_nhitsMostSignificantTrack,sigVec);
    _result = sigVec[4];
  }
};

class FtTrk2Pt : public FTAlgo {
 public:
  FtTrk2Pt() : FTAlgo("trk2pt") {}
  void process() {
    double sigVec[6];
    findMostSignificantTrack(_jet,_privtx,_nhitsMostSignificantTrack,sigVec);
    _result = sigVec[5];
  }
};

class FtTrk1PtByJetE : public FTAlgo {
 public:
  FtTrk1PtByJetE() : FTAlgo("trk1pt_jete") {}
  void process() {
    double sigVec[6];
    findMostSignificantTrack(_jet,_privtx,_nhitsMostSignificantTrack,sigVec);
    _result = sigVec[4] / _jet->Energy();
  }
};

class FtTrk2PtByJetE : public FTAlgo {
 public:
  FtTrk2PtByJetE() : FTAlgo("trk2pt_jete") {}
  void process() {
    double sigVec[6];
    findMostSignificantTrack(_jet,_privtx,_nhitsMostSignificantTrack,sigVec);
    _result = sigVec[5] / _jet->Energy();
  }
};

class FtJProbR : public FTAlgo {
 public:
  FtJProbR() : FTAlgo("jprobr") {}
  void process() {
    _result = jointProbD0(_jet,_privtx,_nhitsJointProbD0);
  }
};

class FtJProbZ : public FTAlgo {
 public:
  FtJProbZ() : FTAlgo("jprobz") {}
  void process() {
    _result = jointProbZ0(_jet,_privtx,_nhitsJointProbZ0);
  }
};

class FtJProbR5Sigma : public FTAlgo {
 public:
  FtJProbR5Sigma(bool usevt) : FTAlgo(usevt ? "jprobr5sigma" : "jprobr5sigmanv"), _useVertexTracks(usevt) {}
  void process() {
    _result = jointProbD0(_jet,_privtx,_nhitsJointProbD0,5., _useVertexTracks);
  }
 private:
  bool _useVertexTracks;
};

class FtJProbZ5Sigma : public FTAlgo {
 public:
  FtJProbZ5Sigma(bool usevt) : FTAlgo(usevt ? "jprobz5sigma" : "jprobz5sigmanv"), _useVertexTracks(usevt) {}
  void process() {
    _result = jointProbZ0(_jet,_privtx,_nhitsJointProbZ0,5., _useVertexTracks);
  }
  bool _useVertexTracks;
};

class FtJProbR2 : public FTAlgo {
 public:
  FtJProbR2() : FTAlgo("jprobr2") {}
  void process() {
    const FtIPProbHolder* holder = FTManager::getInstance().getIPProbHolder();
    if(holder->_hd0jprob && holder->_hd0jprob2)
      _result = jointProb2D0(_jet,_privtx,_nhitsJointProbD0, 200, true, holder->_hd0jprob, holder->_hd0jprob2);
  }
};

class FtJProbZ2 : public FTAlgo {
 public:
  FtJProbZ2() : FTAlgo("jprobz2") {}
  void process() {
    const FtIPProbHolder* holder = FTManager::getInstance().getIPProbHolder();
    if(holder->_hz0jprob && holder->_hz0jprob2)
      _result = jointProb2Z0(_jet,_privtx,_nhitsJointProbZ0, 200, true, holder->_hz0jprob, holder->_hz0jprob2);
  }
};

class FtJProbR25Sigma : public FTAlgo {
 public:
  FtJProbR25Sigma(bool usevt) : FTAlgo(usevt ? "jprobr25sigma" : "jprobr25sigmanv"), _useVertexTracks(usevt) {}
  void process() {
    const FtIPProbHolder* holder = FTManager::getInstance().getIPProbHolder();
    if(holder->_hd0jprob && holder->_hd0jprob2)
      _result = jointProb2D0(_jet,_privtx,_nhitsJointProbD0,5., _useVertexTracks, holder->_hd0jprob, holder->_hd0jprob2);
  }
 private:
  bool _useVertexTracks;
};

class FtJProbZ25Sigma : public FTAlgo {
 public:
  FtJProbZ25Sigma(bool usevt) : FTAlgo(usevt ? "jprobz25sigma" : "jprobz25sigmanv"), _useVertexTracks(usevt) {}
  void process() {
    const FtIPProbHolder* holder = FTManager::getInstance().getIPProbHolder();
    if(holder->_hz0jprob && holder->_hz0jprob2)
      _result = jointProb2Z0(_jet,_privtx,_nhitsJointProbZ0,5., _useVertexTracks, holder->_hz0jprob, holder->_hz0jprob2);
  }
  bool _useVertexTracks;
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
  FtTrkMass() : FTAlgo("trkmass") {}
  void process() {
    vector<const Track*> tracks = _jet->getAllTracks(true);

    TrackSelectorConfig trsel;
    trsel.minD0Sig = 5.;
    trsel.minZ0Sig = 5.;
    trsel.maxD0 = 2.;
    trsel.maxZ0 = 3.;
    vector<const Track*> usedTracks = TrackSelector()(tracks, trsel);

    if (usedTracks.size() < 2)_result = 0.;
    else {
      TLorentzVector v;
      for (unsigned int n=0; n<usedTracks.size(); n++) {
        v += *(TLorentzVector*)(usedTracks[n]);
      }
      _result = v.M();
      /*					cout << "Track mass = " << v.M() << endl;
      					for(unsigned int n=0;n<usedTracks.size(); n++){
      						const MCParticle *mcp = usedTracks[n]->getMcp();
      						if(!mcp)continue;
      						int ppdg = (mcp->getParent() ? mcp->getParent()->getPDG() : 0);
      						int cpdg = (mcp->getSemiStableCParent() ? mcp->getSemiStableCParent()->getPDG() : 0);
      						int bpdg = (mcp->getSemiStableBParent() ? mcp->getSemiStableBParent()->getPDG() : 0);
      						TLorentzVector v2 = v;
      						v2 -= *(TLorentzVector *)(usedTracks[n]);
      						cout << "  M = " << v2.M() << " without track " << n << " pdg = " << mcp->getPDG() << ", ppdg = " << ppdg;
      						cout << ", cpdg = " << cpdg << ", bpdg = " << bpdg << ", e = " << usedTracks[n]->E() << endl;
      					}
      */
    }
  }
};

class FtTrkMass2 : public FTAlgo {
 public:
  FtTrkMass2() : FTAlgo("trkmass2") {}
  void process() {
    vector<const Track*> tracks = _jet->getAllTracks(true);

    TrackSelectorConfig trsel;
    trsel.minD0Sig = 5.;
    trsel.minZ0Sig = 5.;
    trsel.minD0Z0Sig = 7.;
    trsel.maxD0 = 2.;
    trsel.maxZ0 = 3.;
    vector<const Track*> usedTracks = TrackSelector()(tracks, trsel);

    if (usedTracks.size() < 2)_result = 0.;
    else {
      TLorentzVector v;
      for (unsigned int n=0; n<usedTracks.size(); n++) {
        v += *(TLorentzVector*)(usedTracks[n]);
      }

      // mass selection
      bool modify;
      do {
        modify = false;
        TLorentzVector v2;
        for (unsigned int n=0; n<usedTracks.size(); n++) {
          v2 = v;
          v2 -= *(TLorentzVector*)(usedTracks[n]);

          if (v.M() - v2.M() > usedTracks[n]->E()) {
            //cout << "TrkMass2: track removed by mass selection: mass = " << v.M() << ", mass2 = " << v2.M() << ", e = " << usedTracks[n]->E() << endl;
            usedTracks.erase(usedTracks.begin() + n);
            v = v2;
            modify = true;
            break;
          }
        }
      } while (modify == true);

      _result = v.M();
    }
  }
};

class FtNSecTracks : public FTAlgo {
 public:
  FtNSecTracks(bool usevt) : FTAlgo(usevt ? "nsectracks" : "nsectracksnv"), _useVertexTracks(usevt) {}
  void process() {
    vector<const Track*> tracks = (_useVertexTracks ? _jet->getAllTracks(true) : _jet->getTracks());

    TrackSelectorConfig trsel;
    trsel.minD0Sig = 5.;
    trsel.minZ0Sig = 5.;
    trsel.maxD0 = 2.;
    trsel.maxZ0 = 3.;
    vector<const Track*> usedTracks = TrackSelector()(tracks, trsel);

    _result = usedTracks.size();
  }
 private:
  bool _useVertexTracks;
};

class FtVtxLongitudinalDeviation : public FTAlgo {
 public:
  FtVtxLongitudinalDeviation() : FTAlgo("vtxldev") {}
  void process() {
    _result = 0;
    if (_jet->getVertices().size()==1) {
      const Vertex* vtx = _jet->getVertices()[0];

      //cout << "LongitudinalDeviation: vpos " << vtx->getX() << " " << vtx->getY() << " " << vtx->getZ() << endl;

      double devall = 0;
      double devmax = 0;
      for (unsigned int i=0; i<vtx->getTracks().size(); i++) {
        const Track* tr = vtx->getTracks()[i];

        //const MCParticle* mcp = tr->getMcp();
        //int cpdg = 0, bpdg = 0;
        //if (mcp) {
        //  cpdg = (mcp->getSemiStableCParent() ? mcp->getSemiStableCParent()->getPDG() : 0);
        //  bpdg = (mcp->getSemiStableBParent() ? mcp->getSemiStableBParent()->getPDG() : 0);
        //}

        Helix hel(tr);
        double dev = hel.LongitudinalDeviation(_privtx,vtx);
        //cout << "LongitudinalDeviation: track " << i << ", cpdg " << cpdg << ", bpdg " << bpdg << ", dev " << dev << endl;
        devall += dev;

        if (devmax < dev)devmax = dev;
      }

      _result = devmax;
      //cout << "LongitudinalDeviation: devall " << devall << ", devmax " << devmax << endl;
    }
  }
};

class FtNMuon : public FTAlgo {
 public:
  FtNMuon(bool usevt) : FTAlgo(usevt ? "nmuonall" : "nmuon") , _useVertexTracks(usevt) {}
  void process() {
    TrackVec tracks = (_useVertexTracks ? _jet->getAllTracks(true) : _jet->getTracks());
    _result = 0;
    for (unsigned int n=0; n<tracks.size(); n++) {
      if (algoEtc::SimpleSecMuonFinder(tracks[n], 5., 5., 5., -0.1, 0.2, 0.8, 1.5, 4., 0.5)
          || algoEtc::SimpleSecMuonFinder(tracks[n], 5., 5., 5., 0.05, 0., 10., 0., 10., 10.))
        _result += 1;
    }
  }
 private:
  bool _useVertexTracks;
};

class FtNElectron : public FTAlgo {
 public:
  FtNElectron(bool usevt) : FTAlgo(usevt ? "nelectronall" : "nelectron") , _useVertexTracks(usevt) {}
  void process() {
    TrackVec tracks = (_useVertexTracks ? _jet->getAllTracks(true) : _jet->getTracks());
    _result = 0;
    for (unsigned int n=0; n<tracks.size(); n++) {
      if (algoEtc::SimpleSecElectronFinder(tracks[n], 5., 5., 5., 5., 0.98, 0.9, 1.15))
        _result += 1;
    }
  }
 private:
  bool _useVertexTracks;
};

class FtNElectronPID : public FTAlgo {
 public:
  FtNElectronPID(bool usevt) : FTAlgo(usevt ? "nelectronPIDall" : "nelectronPID") , _useVertexTracks(usevt) {}
  void process(){
    TrackVec tracks = (_useVertexTracks ? _jet->getAllTracks(true) : _jet->getTracks());
    _result = 0;
    for(unsigned int n=0; n<tracks.size();n++){
      if(tracks[n]->getPDG()==11 && tracks[n]->P()>=2.0)
	_result += 1;
    }
  }
 private:
  bool _useVertexTracks;
};
  
class FtNMuonPID : public FTAlgo {
 public:
  FtNMuonPID(bool usevt) : FTAlgo(usevt ? "nmuonPIDall" : "nmuonPID") , _useVertexTracks(usevt) {}
  void process(){
    TrackVec tracks = (_useVertexTracks ? _jet->getAllTracks(true) : _jet->getTracks());
    _result = 0;
    for(unsigned int n=0; n<tracks.size();n++){
      if(tracks[n]->getPDG()==13 && tracks[n]->P()>=5.0)
	_result += 1;
    }
  }
 private:
  bool _useVertexTracks;
};

class FtMCNMuon : public FTAlgo {
 public:
  FtMCNMuon() : FTAlgo("MCnmuon") {}
  void process() {
    TrackVec& tracks = _jet->getTracks(); // do not use vertexed tracks
//				TrackVec tracks = _jet->getAllTracks(true);
    _result = 0;
    for (unsigned int n=0; n<tracks.size(); n++) {
      if (tracks[n]->getMcp() && abs(tracks[n]->getMcp()->getPDG())==13)
        _result += 1;
    }
  }
};

class FtMCNElectron : public FTAlgo {
 public:
  FtMCNElectron() : FTAlgo("MCnelectron") {}
  void process() {
    TrackVec& tracks = _jet->getTracks(); // do not use vertexed tracks
//				TrackVec tracks = _jet->getAllTracks(true);
    _result = 0;
    for (unsigned int n=0; n<tracks.size(); n++) {
      if (tracks[n]->getMcp() && abs(tracks[n]->getMcp()->getPDG())==11)
        _result += 1;
    }
  }
};

class FtMCNB : public FTAlgo {
 public:
  FtMCNB() : FTAlgo("MCnb") {}

 private:
  vector<MCVertex*> _mcvs;

 public:
  void processEvent() {
    // run perfect vertex finder - every event
    if (_event->IsExist(_event->getDefaultMCParticles()))
      VertexFinderPerfect::findPerfectVertices(_event->getTracks(), _event->getMCParticles(), _mcvs, 1, 0.1);
  }

  void process() {
    TrackVec tracks = _jet->getAllTracks(false);

    map<const MCParticle*, vector<const Track*> > mclist;

    for (unsigned int n=0; n<tracks.size(); n++) {
      const MCParticle* mcp = tracks[n]->getMcp();
      if (!mcp)continue;
      const MCParticle* mcpb = mcp->getSemiStableBParent();
      if (!mcpb)continue;

      if (mclist.find(mcpb) == mclist.end()) {
        vector<const Track*> trlist;
        trlist.push_back(tracks[n]);
        mclist[mcpb] = trlist;
      } else {
        mclist[mcpb].push_back(tracks[n]);
      }
    }
    _result = mclist.size();

    if (_result) {
      map<const MCParticle*, vector<const Track*> >::iterator it;
      for (it = mclist.begin(); it != mclist.end(); it++) {
        double ejet = 0.;
//						TVector3 ev = it->first->getEndVertex();

//						cout << "mce = " << it->first->E() << ", pdg = " << it->first->getPDG();
//						cout << ", vpos = (" << ev.x() << " " << ev.y() << " " << ev.z() << ") ";
//						cout << "reco tracks = ";
        for (unsigned int ntr=0; ntr<it->second.size(); ntr++)
          ejet += it->second[ntr]->E();
        /*
        							const MCParticle *mcp = it->second[ntr]->getMcp();
        							const MCParticle *mcp = it->second[ntr]->getMcp();
        							cout << "(" << mcp->getPDG() << " " << mcp->E() << ") ";
        							if(find(_jet->getTracks().begin(), _jet->getTracks().end(), it->second[ntr]) == _jet->getTracks().end())cout << "v ";
        						}
        						cout << endl;*/
        double emc = 0.;
        unsigned int nmc = 0;
        // looking for MCVertex
        for (unsigned int n = 0; n< _mcvs.size(); n++) {
          if (_mcvs[n]->getDaughters().size() == 0) {
            cout << "MCVertex has no daughters!!" << endl;
            continue;
          }
          const MCParticle* mcp = _mcvs[n]->getDaughters()[0];
          if (mcp->isParent(it->first)) {
            nmc += _mcvs[n]->getRecoTracks().size();
            for (unsigned int n2 = 0; n2 < _mcvs[n]->getRecoTracks().size(); n2++) {
              emc += _mcvs[n]->getRecoTracks()[n2]->E();
            }


//							  cout << "MCVertex found at (" << mcvs[n]->getPos().x() << " " << mcvs[n]->getPos().y() << " " << mcvs[n]->getPos().z();
//							  cout << "). # tracks = " << mcvs[n]->getRecoTracks().size() << endl;
          }
        }
//						cout << "ejet = " << ejet << ", njet = " << it->second.size() << ", emc = " << emc << ", nmc = " << nmc << endl;
//						if(it->second.size() * 2 < nmc || ejet * 2 < emc){
        if (ejet * 2 < emc) {
          // not accepted
          _result --;
        }
      }
//					cout << "nvtx = " << _jet->getVerticesForFT().size() << ", nvtxall = " << _jet->getVertices().size() << endl;
    }
  }
};

class FtMCNC : public FTAlgo {
 public:
  FtMCNC() : FTAlgo("MCnc") {}

 private:
  vector<MCVertex*> _mcvs;

 public:
  void processEvent() {
    // run perfect vertex finder - every event
    if (_event->IsExist(_event->getDefaultMCParticles()))
      VertexFinderPerfect::findPerfectVertices(_event->getTracks(), _event->getMCParticles(), _mcvs, 1, 0.1);
  }

  void process() {
    TrackVec tracks = _jet->getAllTracks(false);

    map<const MCParticle*, vector<const Track*> > mclist;

    for (unsigned int n=0; n<tracks.size(); n++) {
      const MCParticle* mcp = tracks[n]->getMcp();
      if (!mcp)continue;
      const MCParticle* mcpc = mcp->getSemiStableCParent();
      if (!mcpc)continue;

      // exclude c from b
      const MCParticle* mcpbq = mcp->getSemiStableBParent();
      if (mcpbq)continue;

      if (mclist.find(mcpc) == mclist.end()) {
        vector<const Track*> trlist;
        trlist.push_back(tracks[n]);
        mclist[mcpc] = trlist;
      } else {
        mclist[mcpc].push_back(tracks[n]);
      }
    }
    _result = mclist.size();

    if (_result) {

      map<const MCParticle*, vector<const Track*> >::iterator it;

      for (it = mclist.begin(); it != mclist.end(); it++) {
        double ejet = 0.;
//						TVector3 ev = it->first->getEndVertex();

//						cout << "mce = " << it->first->E() << ", pdg = " << it->first->getPDG();
//						cout << ", vpos = (" << ev.x() << " " << ev.y() << " " << ev.z() << ") ";
//						cout << "reco tracks = ";
        for (unsigned int ntr=0; ntr<it->second.size(); ntr++)
          ejet += it->second[ntr]->E();
        /*
        							const MCParticle *mcp = it->second[ntr]->getMcp();
        							const MCParticle *mcp = it->second[ntr]->getMcp();
        							cout << "(" << mcp->getPDG() << " " << mcp->E() << ") ";
        							if(find(_jet->getTracks().begin(), _jet->getTracks().end(), it->second[ntr]) == _jet->getTracks().end())cout << "v ";
        						}
        						cout << endl;*/
        double emc = 0.;
        unsigned int nmc = 0;
        // looking for MCVertex
        for (unsigned int n = 0; n< _mcvs.size(); n++) {
          if (_mcvs[n]->getDaughters().size() == 0) {
            cout << "MCVertex has no daughters!!" << endl;
            continue;
          }
          const MCParticle* mcp = _mcvs[n]->getDaughters()[0];
          if (mcp->isParent(it->first)) {
            nmc += _mcvs[n]->getRecoTracks().size();
            for (unsigned int n2 = 0; n2 < _mcvs[n]->getRecoTracks().size(); n2++) {
              emc += _mcvs[n]->getRecoTracks()[n2]->E();
            }


//							  cout << "MCVertex found at (" << _mcvs[n]->getPos().x() << " " << _mcvs[n]->getPos().y() << " " << _mcvs[n]->getPos().z();
//							  cout << "). # tracks = " << _mcvs[n]->getRecoTracks().size() << endl;
          }
        }
//						cout << "ejet = " << ejet << ", njet = " << it->second.size() << ", emc = " << emc << ", nmc = " << nmc << endl;
//						if(it->second.size() * 2 < nmc || ejet * 2 < emc){
        if (ejet * 2 < emc) {
          // not accepted
          _result --;
        }

      }

//					cout << "nvtx = " << _jet->getVerticesForFT().size() << ", nvtxall = " << _jet->getVertices().size() << endl;
    }
  }
};

class FtD0bProb : public FTAlgo {
 public:
  FtD0bProb() : FTAlgo("d0bprob") {}
  void process() {
    double prob = 1.;
    TrackVec tracks = _jet->getAllTracks(true);
    for (unsigned int n=0; n<tracks.size(); n++) {
      double d0sig = trackD0Significance(tracks[n], _privtx);
      if ( d0sig > 5) {
        double ld0 = log10(fabs(tracks[n]->getD0()));

        prob *= FTManager::getInstance().getIPProbHolder()->_hd0[0]->GetBinContent(FTManager::getInstance().getIPProbHolder()->_hd0[0]->FindBin(ld0)) * 3.; // 3 = b,c,q, averaging effect
      }
    }

    _result = prob;
  }
};

class FtD0cProb : public FTAlgo {
 public:
  FtD0cProb() : FTAlgo("d0cprob") {}
  void process() {
    double prob = 1.;
    TrackVec tracks = _jet->getAllTracks(true);
    for (unsigned int n=0; n<tracks.size(); n++) {
      double d0sig = trackD0Significance(tracks[n], _privtx);
      if ( d0sig > 5) {
        double ld0 = log10(fabs(tracks[n]->getD0()));

        prob *= FTManager::getInstance().getIPProbHolder()->_hd0[1]->GetBinContent(FTManager::getInstance().getIPProbHolder()->_hd0[1]->FindBin(ld0)) * 3.; // 3 = b,c,q, averaging effect
      }
    }

    _result = prob;
  }
};

class FtD0qProb : public FTAlgo {
 public:
  FtD0qProb() : FTAlgo("d0qprob") {}
  void process() {
    double prob = 1.;
    TrackVec tracks = _jet->getAllTracks(true);
    for (unsigned int n=0; n<tracks.size(); n++) {
      double d0sig = trackD0Significance(tracks[n], _privtx);
      if ( d0sig > 5) {
        double ld0 = log10(fabs(tracks[n]->getD0()));

        prob *= FTManager::getInstance().getIPProbHolder()->_hd0[2]->GetBinContent(FTManager::getInstance().getIPProbHolder()->_hd0[2]->FindBin(ld0)) * 3.; // 3 = b,c,q, averaging effect

//						cout << "d0 = " << tracks[n]->getD0() << ", pq = " <<  FTManager::getInstance().getIPProbHolder()->_hd0[2]->GetBinContent(FTManager::getInstance().getIPProbHolder()->_hd0[2]->FindBin(ld0)) * 3. << endl;
      }
    }

//				cout << "prob = " << prob << endl;
    _result = prob;
  }
};

class FtD0bProbSigned : public FTAlgo {
 public:
  FtD0bProbSigned(bool usevtxtracks = true) : FTAlgo(usevtxtracks ? "d0bprobsigned" : "d0bprobsignednv") {
    _useVertexTracks = usevtxtracks;
  }
  void process() {
    double prob = 1.;
    TrackVec tracks = (_useVertexTracks ? _jet->getAllTracks(true) : _jet->getTracks());
    for (unsigned int n=0; n<tracks.size(); n++) {
      double d0sig = trackD0Significance(tracks[n], _privtx);
      if ( d0sig > 5) {
        double sd0 = signedD0(tracks[n], _jet, _privtx, true);
        if (sd0>0) {
          double ld0 = log10(sd0);
          prob *= FTManager::getInstance().getIPProbHolder()->_hd0p[0]->GetBinContent(FTManager::getInstance().getIPProbHolder()->_hd0p[0]->FindBin(ld0)) * 3.; // 3 = b,c,q, averaging effect
        } else {
          double ld0 = log10(-sd0);
          prob *= FTManager::getInstance().getIPProbHolder()->_hd0n[0]->GetBinContent(FTManager::getInstance().getIPProbHolder()->_hd0n[0]->FindBin(ld0)) * 3.; // 3 = b,c,q, averaging effect
        }
      }
    }

    _result = prob;
  }
 private:
  bool _useVertexTracks;
};

class FtD0cProbSigned : public FTAlgo {
 public:
  FtD0cProbSigned(bool usevtxtracks = true) : FTAlgo(usevtxtracks ? "d0cprobsigned" : "d0cprobsignednv") {
    _useVertexTracks = usevtxtracks;
  }
  void process() {
    double prob = 1.;
    TrackVec tracks = (_useVertexTracks ? _jet->getAllTracks(true) : _jet->getTracks());
    for (unsigned int n=0; n<tracks.size(); n++) {
      double d0sig = trackD0Significance(tracks[n], _privtx);
      if ( d0sig > 5) {
        double sd0 = signedD0(tracks[n], _jet, _privtx, true);
        if (sd0>0) {
          double ld0 = log10(sd0);
          prob *= FTManager::getInstance().getIPProbHolder()->_hd0p[1]->GetBinContent(FTManager::getInstance().getIPProbHolder()->_hd0p[1]->FindBin(ld0)) * 3.; // 3 = b,c,q, averaging effect
        } else {
          double ld0 = log10(-sd0);
          prob *= FTManager::getInstance().getIPProbHolder()->_hd0n[1]->GetBinContent(FTManager::getInstance().getIPProbHolder()->_hd0n[1]->FindBin(ld0)) * 3.; // 3 = b,c,q, averaging effect
        }
      }
    }

    _result = prob;
  }
 private:
  bool _useVertexTracks;
};

class FtD0qProbSigned : public FTAlgo {
 public:
  FtD0qProbSigned(bool usevtxtracks = true) : FTAlgo(usevtxtracks ? "d0qprobsigned" : "d0qprobsignednv") {
    _useVertexTracks = usevtxtracks;
  }
  void process() {
    double prob = 1.;
    TrackVec tracks = (_useVertexTracks ? _jet->getAllTracks(true) : _jet->getTracks());
    for (unsigned int n=0; n<tracks.size(); n++) {
      double d0sig = trackD0Significance(tracks[n], _privtx);
      if ( d0sig > 5) {
        double sd0 = signedD0(tracks[n], _jet, _privtx, true);
        if (sd0>0) {
          double ld0 = log10(sd0);
          prob *= FTManager::getInstance().getIPProbHolder()->_hd0p[2]->GetBinContent(FTManager::getInstance().getIPProbHolder()->_hd0p[2]->FindBin(ld0)) * 3.; // 3 = b,c,q, averaging effect
        } else {
          double ld0 = log10(-sd0);
          prob *= FTManager::getInstance().getIPProbHolder()->_hd0n[2]->GetBinContent(FTManager::getInstance().getIPProbHolder()->_hd0n[2]->FindBin(ld0)) * 3.; // 3 = b,c,q, averaging effect
        }
      }
    }

    _result = prob;
  }
 private:
  bool _useVertexTracks;
};

class FtD0bProbIP : public FTAlgo {
 public:
  FtD0bProbIP() : FTAlgo("d0nonbprobip") {}
  void process() {
    if(!FTManager::getInstance().getIPProbHolder()->_hd0ip[0])return;
    double prob = 1.;
    TrackVec tracks = _jet->getAllTracks(true);
    for (unsigned int n=0; n<tracks.size(); n++) {
      double sd0sig = signedD0Significance(tracks[n], _jet, _privtx, true);
      if ( fabs(sd0sig) < 5) {
        double sbin = FTManager::getInstance().getIPProbHolder()->_hd0ip[0]->GetBinContent(FTManager::getInstance().getIPProbHolder()->_hd0ip[0]->FindBin(sd0sig));
        double qbin = FTManager::getInstance().getIPProbHolder()->_hd0ip[2]->GetBinContent(FTManager::getInstance().getIPProbHolder()->_hd0ip[2]->FindBin(sd0sig));
//						cout << sd0sig << " " << sbin << " " << qbin << " " << qbin/sbin*FTManager::getInstance().getIPProbHolder()->_normd0ip[0] << endl;
        prob *= qbin / sbin * FTManager::getInstance().getIPProbHolder()->_normd0ip[0];
      }
    }

    _result = prob;
  }
};

class FtD0cProbIP : public FTAlgo {
 public:
  FtD0cProbIP() : FTAlgo("d0noncprobip") {}
  void process() {
    if(!FTManager::getInstance().getIPProbHolder()->_hd0ip[0])return;
    double prob = 1.;
    TrackVec tracks = _jet->getAllTracks(true);
    for (unsigned int n=0; n<tracks.size(); n++) {
      double sd0sig = signedD0Significance(tracks[n], _jet, _privtx, true);
      if ( fabs(sd0sig) < 5) {
        double sbin = FTManager::getInstance().getIPProbHolder()->_hd0ip[1]->GetBinContent(FTManager::getInstance().getIPProbHolder()->_hd0ip[1]->FindBin(sd0sig));
        double qbin = FTManager::getInstance().getIPProbHolder()->_hd0ip[2]->GetBinContent(FTManager::getInstance().getIPProbHolder()->_hd0ip[2]->FindBin(sd0sig));
        prob *= qbin / sbin * FTManager::getInstance().getIPProbHolder()->_normd0ip[1];
      }
    }

    _result = prob;
  }
};

class FtZ0bProb : public FTAlgo {
 public:
  FtZ0bProb(bool usevtxtracks = true) : FTAlgo(usevtxtracks ? "z0bprob" : "z0bprobnv") {
    _useVertexTracks = usevtxtracks;
  }
  void process() {
    double prob = 1.;
    TrackVec tracks = (_useVertexTracks ? _jet->getAllTracks(true) : _jet->getTracks());
    for (unsigned int n=0; n<tracks.size(); n++) {
      double z0sig = trackZ0Significance(tracks[n], _privtx);
      if ( z0sig > 5) {
        double lz0 = log10(fabs(tracks[n]->getZ0()));

        prob *= FTManager::getInstance().getIPProbHolder()->_hz0[0]->GetBinContent(FTManager::getInstance().getIPProbHolder()->_hz0[0]->FindBin(lz0)) * 3.; // 3 = b,c,q, averaging effect
      }
    }

    _result = prob;
  }
 private:
  bool _useVertexTracks;
};

class FtZ0cProb : public FTAlgo {
 public:
  FtZ0cProb(bool usevtxtracks = true) : FTAlgo(usevtxtracks ? "z0cprob" : "z0cprobnv") {
    _useVertexTracks = usevtxtracks;
  }
  void process() {
    double prob = 1.;
    TrackVec tracks = (_useVertexTracks ? _jet->getAllTracks(true) : _jet->getTracks());
    for (unsigned int n=0; n<tracks.size(); n++) {
      double z0sig = trackZ0Significance(tracks[n], _privtx);
      if ( z0sig > 5) {
        double lz0 = log10(fabs(tracks[n]->getZ0()));

        prob *= FTManager::getInstance().getIPProbHolder()->_hz0[1]->GetBinContent(FTManager::getInstance().getIPProbHolder()->_hz0[1]->FindBin(lz0)) * 3.; // 3 = b,c,q, averaging effect
      }
    }

    _result = prob;
  }
 private:
  bool _useVertexTracks;
};

class FtZ0qProb : public FTAlgo {
 public:
  FtZ0qProb(bool usevtxtracks = true) : FTAlgo(usevtxtracks ? "z0qprob" : "z0qprobnv") {
    _useVertexTracks = usevtxtracks;
  }
  void process() {
    double prob = 1.;
    TrackVec tracks = (_useVertexTracks ? _jet->getAllTracks(true) : _jet->getTracks());
    for (unsigned int n=0; n<tracks.size(); n++) {
      double z0sig = trackZ0Significance(tracks[n], _privtx);
      if ( z0sig > 5) {
        double lz0 = log10(fabs(tracks[n]->getZ0()));

        prob *= FTManager::getInstance().getIPProbHolder()->_hz0[2]->GetBinContent(FTManager::getInstance().getIPProbHolder()->_hz0[2]->FindBin(lz0)) * 3.; // 3 = b,c,q, averaging effect
      }
    }

    _result = prob;
  }
 private:
  bool _useVertexTracks;
};

class FtZ0bProbIP : public FTAlgo {
 public:
  FtZ0bProbIP() : FTAlgo("z0nonbprobip") {}
  void process() {
    if(!FTManager::getInstance().getIPProbHolder()->_hd0ip[0])return;
    double prob = 1.;
    TrackVec tracks = _jet->getAllTracks(true);
    for (unsigned int n=0; n<tracks.size(); n++) {
      double sz0sig = signedZ0Significance(tracks[n], _jet, _privtx, true);
      if ( fabs(sz0sig) < 5) {
        double sbin = FTManager::getInstance().getIPProbHolder()->_hz0ip[0]->GetBinContent(FTManager::getInstance().getIPProbHolder()->_hz0ip[0]->FindBin(sz0sig));
        double qbin = FTManager::getInstance().getIPProbHolder()->_hz0ip[2]->GetBinContent(FTManager::getInstance().getIPProbHolder()->_hz0ip[2]->FindBin(sz0sig));
        prob *= qbin / sbin * FTManager::getInstance().getIPProbHolder()->_normz0ip[0];
      }
    }

    _result = prob;
  }
};

class FtZ0cProbIP : public FTAlgo {
 public:
  FtZ0cProbIP() : FTAlgo("z0noncprobip") {}
  void process() {
    if(!FTManager::getInstance().getIPProbHolder()->_hd0ip[0])return;
    double prob = 1.;
    TrackVec tracks = _jet->getAllTracks(true);
    for (unsigned int n=0; n<tracks.size(); n++) {
      double sz0sig = signedZ0Significance(tracks[n], _jet, _privtx, true);
      if ( fabs(sz0sig) < 5) {
        double sbin = FTManager::getInstance().getIPProbHolder()->_hz0ip[1]->GetBinContent(FTManager::getInstance().getIPProbHolder()->_hz0ip[1]->FindBin(sz0sig));
        double qbin = FTManager::getInstance().getIPProbHolder()->_hz0ip[2]->GetBinContent(FTManager::getInstance().getIPProbHolder()->_hz0ip[2]->FindBin(sz0sig));
        prob *= qbin / sbin * FTManager::getInstance().getIPProbHolder()->_normz0ip[1];
      }
    }

    _result = prob;
  }
};

void FTManager::initVars() {
  if (_initialized)return;
  _initialized = true;

  cout << "Initializing FTManager variables." << endl;

  add( new FtAuxiliaryM("aux") );

  add( new FtNtrkWithoutV0() );
  add( new FtNtrk() );

  add( new FtNvtxAll() );
  add( new FtVtxMassAll() );
  add( new FtVtxLen12All() );
  add( new FtVtxLen12AllByJetE() );
  add( new Ft1VtxProb() );

  add( new FtNvtx() );
  add( new FtJetE() );
  add( new FtVtxLen1() );
  add( new FtVtxLen2() );
  add( new FtVtxLen12() );
  add( new FtVtxLen1ByJetE() );
  add( new FtVtxLen2ByJetE() );
  add( new FtVtxLen12ByJetE() );
  add( new FtVtxSig1() );
  add( new FtVtxSig2() );
  add( new FtVtxSig12() );
  add( new FtVtxSig1ByJetE() );
  add( new FtVtxSig2ByJetE() );
  add( new FtVtxSig12ByJetE() );
  add( new FtVtxDirAng1() );
  add( new FtVtxDirAng2() );
  add( new FtVtxDirAng12() );
  add( new FtVtxDirAng1TimesJetE() );
  add( new FtVtxDirAng2TimesJetE() );
  add( new FtVtxDirAng12TimesJetE() );
  add( new FtVtxMom() );
  add( new FtVtxMom1() );
  add( new FtVtxMom2() );
  add( new FtVtxMomByJetE() );
  add( new FtVtxMom1ByJetE() );
  add( new FtVtxMom2ByJetE() );
  add( new FtVtxMass() );
  add( new FtVtxMass1() );
  add( new FtVtxMass2() );
  add( new FtVtxMassPtCorr() );
  add( new FtVtxMult() );
  add( new FtVtxMult1() );
  add( new FtVtxMult2() );
  add( new FtVtxProb() );
  add( new FtTrk1D0Sig() );
  add( new FtTrk2D0Sig() );
  add( new FtTrk1Z0Sig() );
  add( new FtTrk2Z0Sig() );
  add( new FtTrk1Pt() );
  add( new FtTrk2Pt() );
  add( new FtTrk1PtByJetE() );
  add( new FtTrk2PtByJetE() );
  add( new FtJProbR() );
  add( new FtJProbZ() );
  add( new FtJProbR5Sigma(true) );
  add( new FtJProbZ5Sigma(true) );
  add( new FtJProbR5Sigma(false) );
  add( new FtJProbZ5Sigma(false) );
  add( new FtJProbR2() );
  add( new FtJProbZ2() );
  add( new FtJProbR25Sigma(true) );
  add( new FtJProbZ25Sigma(true) );
  add( new FtJProbR25Sigma(false) );
  add( new FtJProbZ25Sigma(false) );
  add( new FtSphericity() );
  add( new FtTrkMass() );
  add( new FtTrkMass2() );
  add( new FtNSecTracks(true) );
  add( new FtNSecTracks(false) );
  add( new FtNMuon(true) );
  add( new FtNElectron(true) );
  add( new FtNMuon(false) );
  add( new FtNElectron(false) );
  add( new FtNMuonPID(true) );
  add( new FtNElectronPID(true) );
  add( new FtNMuonPID(false) );
  add( new FtNElectronPID(false) );
  add( new FtMCNMuon() );
  add( new FtMCNElectron() );
  add( new FtMCNB() );
  add( new FtMCNC() );

  add( new FtVtxLongitudinalDeviation() );

  // d0/z0 probs
  add( new FtD0bProb());
  add( new FtD0cProb());
  add( new FtD0qProb());
  add( new FtD0bProbSigned(true));
  add( new FtD0cProbSigned(true));
  add( new FtD0qProbSigned(true));
  add( new FtD0bProbSigned(false));
  add( new FtD0cProbSigned(false));
  add( new FtD0qProbSigned(false));
  add( new FtD0bProbIP());
  add( new FtD0cProbIP());
  add( new FtZ0bProb(true));
  add( new FtZ0cProb(true));
  add( new FtZ0qProb(true));
  add( new FtZ0bProb(false));
  add( new FtZ0cProb(false));
  add( new FtZ0qProb(false));
  add( new FtZ0bProbIP());
  add( new FtZ0cProbIP());
}

void FlavorTag::init(Parameters* param) {
  Algorithm::init(param);

  _primvtxcolname = param->get("PrimaryVertexCollectionName",string("PrimaryVertex"));
  _jetcolname = param->get("FlavorTag.JetCollectionName",string("VertexJets"));
  Event::Instance()->setDefaultPrimaryVertex(_primvtxcolname.c_str()); // backward compatibility

  _auxiliaryInfo = param->get("MakeNtuple.AuxiliaryInfo",int(-1));

  string d0probfilename = param->get("FlavorTag.D0ProbFileName",string("data/vtxprob/d0prob_zpole.root"));
  string z0probfilename = param->get("FlavorTag.Z0ProbFileName",string("data/vtxprob/z0prob_zpole.root"));

  _holder = new FtIPProbHolder(d0probfilename.c_str(), z0probfilename.c_str());

  _nhitsJointProbD0 = param->get("FlavourTag.NVTXhitsJointProbD0", int(4));
  _nhitsJointProbZ0 = param->get("FlavourTag.NVTXhitsJointProbZ0", int(4));
  _nhitsMostSignificantTrack = param->get("FlavourTag.NhitsMostSignificantTrack", int(4));

  //string outputFilename = param->get("TrainNtupleFile",string("lcfiplus.root"));
  //_nJet = (int)param->get("TrainNJet",float(2));

  //cout << "FlavorTag: Ntuple file set to " << outputFilename << endl;
  //cout << "FlavorTag: Number of jet set to " << _nJet << endl;

  FTManager& mgr = FTManager::getInstance();
  mgr.initVars();

}

void FlavorTag::process() {

  Event* event = Event::Instance();
  //if (event->getTracks().size() == 0) return;

  const Vertex* privtx = event->getPrimaryVertex(_primvtxcolname.c_str());

  //TrackVec & tracks = event->getTracks();
  JetVec* jetsPtr(0);
  bool success = event->Get(_jetcolname.c_str(), jetsPtr);
  if (!success) {
    cout << "jets could not be found" << endl;
    return;
  }
  JetVec& jets = *jetsPtr;

  FTManager& mgr = FTManager::getInstance();

  mgr.setAuxiliary(_auxiliaryInfo);
  mgr.setIPProbHolder(_holder);

  mgr.process(event, privtx, _nhitsJointProbD0, _nhitsJointProbZ0, _nhitsMostSignificantTrack, jets);
}

void FlavorTag::end() {
  delete _holder;
}
}
