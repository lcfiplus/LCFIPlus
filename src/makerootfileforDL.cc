#include <string>

#include "TFile.h"
#include "TROOT.h"
#include "TNtuple.h"
#include "TTree.h"
#include "TNtupleD.h"
#include "TSystem.h"
#include "TPad.h"
#include "TStyle.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TMatrixFSym.h"
#include <TMatrixTSym.h>
#include "TMatrixDSym.h"
#include "HelixClass.h"
#include <HelixClass.h>
#include "TVectorD.h"
#include "TMethodCall.h"

#include "lcfiplus.h"
#include "process.h"
#include "geometry.h"
#include "VertexSelector.h"
#include "algoEtc.h"
#include "algoSigProb.h"
#include "VertexFinderSuehara.h"
#include "TrackSelector.h"
#include "VertexFitterSimple.h"
#include "TrackNtuple.h"
#include "EventStore.h"
#include "LcfiplusProcessor.h"
#include <EVENT/LCEvent.h>

#include "math.h"

#include "makerootfileforDL.h"

using namespace lcfiplus;
using namespace lcfiplus::algoSigProb;
using namespace lcfiplus::algoEtc;

using namespace lcfiplus;

namespace lcfiplus {

// ================================================================================================================ //
// ================================================================================================================ //
void MakeROOTFileCC::init(Parameters* param) {
  Algorithm::init(param);
  
  string filename = param->get("MakeROOTFileCC.RootFileName", string("my_track1.root"));
  
  _file = new TFile(filename.c_str(), "RECREATE");
  
  _ntp = new TTree("track0", "track0");

  nEvt = 0;
  ntr1Trk = 0;
  ntr2Trk = 0;

  TracksData& d = _data;
  _ntp->Branch("nevent", &nEvt, "nevent/I"); // Event number
  _ntp->Branch("ntr1track", &ntr1Trk, "ntr1track/I"); // Event number
  _ntp->Branch("ntr2track", &ntr2Trk, "ntr2track/I"); // Event number
  // Track1 low data
  _ntp->Branch("tr1d0", &d.tr1d0, "tr1d0/D"); // Impact paramter of the track in (r-phi)
  _ntp->Branch("tr1z0", &d.tr1z0, "tr1z0/D"); // mpact paramter of the track in (r-z)
  _ntp->Branch("tr1phi", &d.tr1phi, "tr1phi/D"); // Phi of the track at the reference point
  _ntp->Branch("tr1omega", &d.tr1omega, "tr1omega/D"); // Omega is the signed curvature of the track in [1/mm]
  _ntp->Branch("tr1tanlam", &d.tr1tanlam, "tr1tanlam/D"); // Lambda is the dip angle of the track in r-z at the reference point
  _ntp->Branch("tr1x", &d.tr1x, "tr1x/D");
  _ntp->Branch("tr1y", &d.tr1y, "tr1y/D");
  _ntp->Branch("tr1z", &d.tr1z, "tr1z/D");
  _ntp->Branch("tr1charge", &d.tr1charge, "tr1charge/I");
  _ntp->Branch("tr1energy", &d.tr1energy, "tr1energy/D");
  _ntp->Branch("tr1covmatrixd0d0", &d.tr1covmatrixd0d0, "tr1covmatrixd0d0/D"); // Covariance matrix of the track parameters
  _ntp->Branch("tr1covmatrixd0z0", &d.tr1covmatrixd0z0, "tr1covmatrixd0z0/D"); // Covariance matrix of the track parameters
  _ntp->Branch("tr1covmatrixd0ph", &d.tr1covmatrixd0ph, "tr1covmatrixd0ph/D"); // Covariance matrix of the track parameters
  _ntp->Branch("tr1covmatrixd0om", &d.tr1covmatrixd0om, "tr1covmatrixd0om/D"); // Covariance matrix of the track parameters
  _ntp->Branch("tr1covmatrixd0tl", &d.tr1covmatrixd0tl, "tr1covmatrixd0tl/D"); // Covariance matrix of the track parameters
  _ntp->Branch("tr1covmatrixz0z0", &d.tr1covmatrixz0z0, "tr1covmatrixz0z0/D"); // Covariance matrix of the track parameters
  _ntp->Branch("tr1covmatrixz0ph", &d.tr1covmatrixz0ph, "tr1covmatrixz0ph/D"); // Covariance matrix of the track parameters
  _ntp->Branch("tr1covmatrixz0om", &d.tr1covmatrixz0om, "tr1covmatrixz0om/D"); // Covariance matrix of the track parameters
  _ntp->Branch("tr1covmatrixz0tl", &d.tr1covmatrixz0tl, "tr1covmatrixz0tl/D"); // Covariance matrix of the track parameters
  _ntp->Branch("tr1covmatrixphph", &d.tr1covmatrixphph, "tr1covmatrixphph/D"); // Covariance matrix of the track parameters
  _ntp->Branch("tr1covmatrixphom", &d.tr1covmatrixphom, "tr1covmatrixphom/D"); // Covariance matrix of the track parameters
  _ntp->Branch("tr1covmatrixphtl", &d.tr1covmatrixphtl, "tr1covmatrixphtl/D"); // Covariance matrix of the track parameters
  _ntp->Branch("tr1covmatrixomom", &d.tr1covmatrixomom, "tr1covmatrixomom/D"); // Covariance matrix of the track parameters
  _ntp->Branch("tr1covmatrixomtl", &d.tr1covmatrixomtl, "tr1covmatrixomtl/D"); // Covariance matrix of the track parameters
  _ntp->Branch("tr1covmatrixtltl", &d.tr1covmatrixtltl, "tr1covmatrixtltl/D"); // Covariance matrix of the track parameters
  // ======================================================= //
  // Track2 low data
  _ntp->Branch("tr2d0", &d.tr2d0, "tr2d0/D"); // Impact paramter of the track in (r-phi)
  _ntp->Branch("tr2z0", &d.tr2z0, "tr2z0/D"); // mpact paramter of the track in (r-z)
  _ntp->Branch("tr2phi", &d.tr2phi, "tr2phi/D"); // Phi of the track at the reference point
  _ntp->Branch("tr2omega", &d.tr2omega, "tr2omega/D"); // Omega is the signed curvature of the track in [1/mm]
  _ntp->Branch("tr2tanlam", &d.tr2tanlam, "tr2tanlam/D"); // Lambda is the dip angle of the track in r-z at the reference point
  _ntp->Branch("tr2x", &d.tr2x, "tr2x/D");
  _ntp->Branch("tr2y", &d.tr2y, "tr2y/D");
  _ntp->Branch("tr2z", &d.tr2z, "tr2z/D");
  _ntp->Branch("tr2charge", &d.tr2charge, "tr2charge/I");
  _ntp->Branch("tr2energy", &d.tr2energy, "tr2energy/D");
  _ntp->Branch("tr2covmatrixd0d0", &d.tr2covmatrixd0d0, "tr2covmatrixd0d0/D"); // Covariance matrix of the track parameters
  _ntp->Branch("tr2covmatrixd0z0", &d.tr2covmatrixd0z0, "tr2covmatrixd0z0/D"); // Covariance matrix of the track parameters
  _ntp->Branch("tr2covmatrixd0ph", &d.tr2covmatrixd0ph, "tr2covmatrixd0ph/D"); // Covariance matrix of the track parameters
  _ntp->Branch("tr2covmatrixd0om", &d.tr2covmatrixd0om, "tr2covmatrixd0om/D"); // Covariance matrix of the track parameters
  _ntp->Branch("tr2covmatrixd0tl", &d.tr2covmatrixd0tl, "tr2covmatrixd0tl/D"); // Covariance matrix of the track parameters
  _ntp->Branch("tr2covmatrixz0z0", &d.tr2covmatrixz0z0, "tr2covmatrixz0z0/D"); // Covariance matrix of the track parameters
  _ntp->Branch("tr2covmatrixz0ph", &d.tr2covmatrixz0ph, "tr2covmatrixz0ph/D"); // Covariance matrix of the track parameters
  _ntp->Branch("tr2covmatrixz0om", &d.tr2covmatrixz0om, "tr2covmatrixz0om/D"); // Covariance matrix of the track parameters
  _ntp->Branch("tr2covmatrixz0tl", &d.tr2covmatrixz0tl, "tr2covmatrixz0tl/D"); // Covariance matrix of the track parameters
  _ntp->Branch("tr2covmatrixphph", &d.tr2covmatrixphph, "tr2covmatrixphph/D"); // Covariance matrix of the track parameters
  _ntp->Branch("tr2covmatrixphom", &d.tr2covmatrixphom, "tr2covmatrixphom/D"); // Covariance matrix of the track parameters
  _ntp->Branch("tr2covmatrixphtl", &d.tr2covmatrixphtl, "tr2covmatrixphtl/D"); // Covariance matrix of the track parameters
  _ntp->Branch("tr2covmatrixomom", &d.tr2covmatrixomom, "tr2covmatrixomom/D"); // Covariance matrix of the track parameters
  _ntp->Branch("tr2covmatrixomtl", &d.tr2covmatrixomtl, "tr2covmatrixomtl/D"); // Covariance matrix of the track parameters
  _ntp->Branch("tr2covmatrixtltl", &d.tr2covmatrixtltl, "tr2covmatrixtltl/D"); // Covariance matrix of the track parameters
  // ======================================================= //
  // Track1 MC particle data
  _ntp->Branch("tr1mcx", &d.tr1mcx, "tr1mcx/D"); // the production vertex of the particle in [mm].
  _ntp->Branch("tr1mcy", &d.tr1mcy, "tr1mcy/D"); // the production vertex of the particle in [mm].
  _ntp->Branch("tr1mcz", &d.tr1mcz, "tr1mcz/D"); // the production vertex of the particle in [mm].
  _ntp->Branch("tr1id", &d.tr1id, "tr1id/I"); // the number of parents of this particle - 0 if mother
  _ntp->Branch("tr1pdg", &d.tr1pdg, "tr1pdg/I"); // the number of parents of this particle - 0 if mother
  _ntp->Branch("tr1sscid", &tr1sscid, "tr1sscid/I"); // other : 0, semistable c : xxx, primary : 2
  _ntp->Branch("tr1sscpdg", &tr1sscpdg, "tr1sscpdg/I"); // other : 0, semistable c : 1, primary : 2
  _ntp->Branch("tr1ssc", &tr1ssc, "tr1ssc/I"); // other : 0, semistable c : 1, primary : 2
  // ======================================================= //
  // Track2 MC particle data
  _ntp->Branch("tr2mcx", &d.tr2mcx, "tr2mcx/D"); // the production vertex of the particle in [mm].
  _ntp->Branch("tr2mcy", &d.tr2mcy, "tr2mcy/D"); // the production vertex of the particle in [mm].
  _ntp->Branch("tr2mcz", &d.tr2mcz, "tr2mcz/D"); // the production vertex of the particle in [mm].
  _ntp->Branch("tr2id", &d.tr2id, "tr2id/I"); // the number of parents of this particle - 0 if mother
  _ntp->Branch("tr2pdg", &d.tr2pdg, "tr2pdg/I"); // the number of parents of this particle - 0 if mother
  _ntp->Branch("tr2sscid", &tr2sscid, "tr2sscid/I"); // semistable exist : 1 other : 0
  _ntp->Branch("tr2sscpdg", &tr2sscpdg, "tr2sscpdg/I"); // semistable exist : 1 other : 0
  _ntp->Branch("tr2ssc", &tr2ssc, "tr2ssc/I"); // semistable exist : 1 other : 0
  // ======================================================= //
  _ntp->Branch("coSemiStableC", &d.cosemistablec, "coSemiStableC/I"); // cosemistableC is exist : 1 other : 0
  _ntp->Branch("connect", &d.connect, "connect/I"); // 0 : not connect, 1 : semistablecparent, 2 : primary vertex
  _ntp->Branch("lcfiplustag", &d.lcfiplustag, "lcfiplustag/I"); // 0 : not connect, 1 : semistablecparent, 2 : primary vertex
  // ======================================================= //
  _ntp->Branch("chi2", &d.chi2, "chi2/D"); // chi2
  _ntp->Branch("vchi2", &d.vchi2, "vchi2/D"); // vertex chi2
  _ntp->Branch("vposx", &d.vposx, "vposx/D"); // vertex position x
  _ntp->Branch("vposy", &d.vposy, "vposy/D"); // vertex position y
  _ntp->Branch("vposz", &d.vposz, "vposz/D"); // vertex position z
  _ntp->Branch("mass", &d.mass, "mass/D"); // mass
  _ntp->Branch("mag", &d.mag, "mag/D"); // Mag
  _ntp->Branch("vec", &d.vec, "vec/D"); // Vec direction
  _ntp->Branch("tr1selection", &d.tr1selection, "tr1selection/I"); // track selection
  _ntp->Branch("tr2selection", &d.tr2selection, "tr2selection/I"); // track selection
  _ntp->Branch("v0selection", &d.v0selection, "v0selection/I"); // v0 selection
}

void MakeROOTFileCC::process() {

  Event* event = Event::Instance();
  const TrackVec& track_list = event->getTracks();
	  
  for (unsigned int i=0; i < track_list.size(); ++i) {
    ntr1Trk = i;
    const Track* track1 = track_list[i];
    if (!track1) continue;
    const MCParticle* mcpc1 = event->getMCParticle(track1);
    _hel1 = new lcfiplus::Helix(track1, PointBase::NOTUSED);

    memset(&_data,0,sizeof(_data));

    _data.tr1d0 = track1->getD0();
    _data.tr1z0 = track1->getZ0();
    _data.tr1phi = track1->getPhi();
    _data.tr1omega = track1->getOmega();
    _data.tr1tanlam = track1->getTanLambda();

    _data.tr1charge = (int)track1->getCharge();

    _data.tr1x = _hel1->GetPos(0.).X();
    _data.tr1y = _hel1->GetPos(0.).Y();
    _data.tr1z = _hel1->GetPos(0.).Z();

    const double* cov1 = track1->getCovMatrix();
    
    _data.tr1covmatrixd0d0 = cov1[0]; // d0d0
    _data.tr1covmatrixd0z0 = cov1[1]; // d0z0
    _data.tr1covmatrixd0ph = cov1[2]; // d0ph
    _data.tr1covmatrixd0om = cov1[3]; // d0om
    _data.tr1covmatrixd0tl = cov1[4]; // d0tl
    _data.tr1covmatrixz0z0 = cov1[5]; // z0z0
    _data.tr1covmatrixz0ph = cov1[6]; // z0ph
    _data.tr1covmatrixz0om = cov1[7]; // z0om
    _data.tr1covmatrixz0tl = cov1[8]; // z0tl
    _data.tr1covmatrixphph = cov1[9]; // phph
    _data.tr1covmatrixphom = cov1[10]; // phom
    _data.tr1covmatrixphtl = cov1[11]; // phtl
    _data.tr1covmatrixomom = cov1[12]; // omom
    _data.tr1covmatrixomtl = cov1[13]; // omtl
    _data.tr1covmatrixtltl = cov1[14]; // tltl
    
    if (mcpc1){
      _data.tr1mcx = mcpc1->getVertex().X();
      _data.tr1mcy = mcpc1->getVertex().Y();
      _data.tr1mcz = mcpc1->getVertex().Z();
      _data.tr1id = mcpc1->getId();
      _data.tr1pdg = mcpc1->getPDG();

      const MCParticle* ssmcpc1 = mcpc1->getSemiStableParent();

      if (ssmcpc1==0){
	tr1ssc = 2;
	tr1sscid = 0;
	tr1sscid = 0;
      }
      if (ssmcpc1){
        tr1sscid = ssmcpc1->getId();
        tr1sscpdg = ssmcpc1->getPDG();
        if (ssmcpc1->getParent()->getId()==0){
          tr1ssc = 2;
        }
	else if (ssmcpc1->isSemiStableC()){
	  tr1ssc = 1;
        }
        else {
	  tr1ssc = 0;
        }
      }
    }

    for (unsigned int j=0; j < i; ++j) {
      ntr2Trk = j;
      const Track* track2 = track_list[j];
      if (!track2) continue;

      const MCParticle* mcpc2 = event->getMCParticle(track2);
      _hel2 = new lcfiplus::Helix(track2, PointBase::NOTUSED);

      _data.tr2d0 = track2->getD0();
      _data.tr2z0 = track2->getZ0();
      _data.tr2phi = track2->getPhi();
      _data.tr2omega = track2->getOmega();
      _data.tr2tanlam = track2->getTanLambda();

      _data.tr2charge = (int)track2->getCharge();

      _data.tr2x = _hel2->GetPos(0.).X();
      _data.tr2y = _hel2->GetPos(0.).Y();
      _data.tr2z = _hel2->GetPos(0.).Z();

      const double* cov2 = track2->getCovMatrix();
    
      _data.tr2covmatrixd0d0 = cov2[0]; // d0d0
      _data.tr2covmatrixd0z0 = cov2[1]; // d0z0
      _data.tr2covmatrixd0ph = cov2[2]; // d0ph
      _data.tr2covmatrixd0om = cov2[3]; // d0om
      _data.tr2covmatrixd0tl = cov2[4]; // d0tl
      _data.tr2covmatrixz0z0 = cov2[5]; // z0z0
      _data.tr2covmatrixz0ph = cov2[6]; // z0ph
      _data.tr2covmatrixz0om = cov2[7]; // z0om
      _data.tr2covmatrixz0tl = cov2[8]; // z0tl
      _data.tr2covmatrixphph = cov2[9]; // phph
      _data.tr2covmatrixphom = cov2[10]; // phom
      _data.tr2covmatrixphtl = cov2[11]; // phtl
      _data.tr2covmatrixomom = cov2[12]; // omom
      _data.tr2covmatrixomtl = cov2[13]; // omtl
      _data.tr2covmatrixtltl = cov2[14]; // tltl
      
      if (mcpc2){
        _data.tr2mcx = mcpc2->getVertex().X();
        _data.tr2mcy = mcpc2->getVertex().Y();
        _data.tr2mcz = mcpc2->getVertex().Z();
        _data.tr2id = mcpc2->getId();
        _data.tr2pdg = mcpc2->getPDG();

        const MCParticle* ssmcpc2 = mcpc2->getSemiStableParent();

        if (ssmcpc2==0){
	  tr2ssc = 2;
	  tr2sscid = 0;
	  tr2sscid = 0;
        }
        if (ssmcpc2){
          tr2sscid = ssmcpc2->getId();
          tr2sscpdg = ssmcpc2->getPDG();
          if (ssmcpc2->getParent()->getId()==0){
            tr2ssc = 2;
          }
	  else if (ssmcpc2->isSemiStableC()){
	    tr2ssc = 1;
          }
          else {
	    tr2ssc = 0;
          }
        }
      }

      if (((tr1ssc==1)&&(tr2ssc==1))&&(tr1sscid==tr2sscid)) _data.connect = 1;
      else if ((tr1ssc==2)&&(tr2ssc==2)) _data.connect = 2;
      else _data.connect = 0;

      vector<const Track*> vttmp;
      vttmp.push_back(track1);
      vttmp.push_back(track2);

      Vertex* vtx = VertexFitterSimple_V() (vttmp.begin(), vttmp.end(), 0);

      _data.chi2 = max(vtx->getChi2Track(track1), vtx->getChi2Track(track2));
      _data.vchi2 = vtx->getChi2();

      _data.vposx = vtx->getPos().X();
      _data.vposy = vtx->getPos().Y();
      _data.vposz = vtx->getPos().Z();

      TLorentzVector tr1tlv = *track1;
      TLorentzVector tr2tlv = *track2;

      _data.mass = (tr1tlv+tr2tlv).M();
      _data.mag = vtx->getPos().Mag();
      _data.vec = vtx->getPos().Dot((tr1tlv+tr2tlv).Vect());

      _data.tr1energy = tr1tlv.E();
      _data.tr2energy = tr2tlv.E();

      
      // suehara 200319
      // TrackSelector and VertexSelector config
      TrackSelectorConfig cfgtrack;
      cfgtrack.maxD0 = 10.;
      cfgtrack.maxZ0 = 20.;
      cfgtrack.minPt = 0.1;
      cfgtrack.maxInnermostHitRadius = 1e10;
    
      cfgtrack.maxD0Err = 0.1;
      cfgtrack.maxZ0Err = 0.1;

      cfgtrack.minTpcHits = 20;
      cfgtrack.minFtdHits = 3;
      cfgtrack.minVtxHits = 3;
      cfgtrack.minVtxPlusFtdHits = 1;
    
      VertexFinderSuehara::VertexFinderSueharaConfig cfgvtx;
      // chi2 threshold = 9.
      // mass threshold = 10.
      // minimum distance = 0.3;
      
      TrackSelector sel;
    
      _data.lcfiplustag = 0; // 0: NC, 1: secondary, 2: primary
      const Vertex *ip = event->getPrimaryVertex();
      TrackVec &iptracks = ip->getTracks();
      if(std::find(iptracks.begin(), iptracks.end(), track1) != iptracks.end()
        && std::find(iptracks.begin(), iptracks.end(), track2) != iptracks.end())
          _data.lcfiplustag = 2; // primary tracks
      else{
        // track selection
        if(sel.passesCut(track1, cfgtrack) && sel.passesCut(track2, cfgtrack)){
          // mass, chi2, vpos, v0 selection
          TLorentzVector v1 = *track1;
          TLorentzVector v2 = *track2;
          double mass = (v1+v2).M();
          double chi2 = max(vtx->getChi2Track(track1), vtx->getChi2Track(track2));
          
          if(mass < min(v1.E(), v2.E()) && mass < 10. && chi2 < 9. && vtx->getPos().Mag() > 0.3
            && vtx->getPos().Dot((v1+v2).Vect()) > 0 && !VertexSelector().passesCut(vtx, cfgvtx.v0selVertex, ip)){
            _data.lcfiplustag = 1; // secondary tracks
          }
        }
      }

      if(sel.passesCut(track1, cfgtrack)) _data.tr1selection = 1;
      else _data.tr1selection = 0;

      if(sel.passesCut(track2, cfgtrack)) _data.tr2selection = 1;
      else _data.tr2selection = 0;

      if(VertexSelector().passesCut(vtx, cfgvtx.v0selVertex, ip))  _data.v0selection = 1;
      else _data.v0selection = 0;
      
      _ntp->Fill();
    }

  }
  nEvt++;
}

void MakeROOTFileCC::end() {
  _file->Write();
  _file->Close();
}

// ================================================================================================================ //
// ================================================================================================================ //
// ================================================================================================================ //
void MakeROOTFileBB::init(Parameters* param) {
  Algorithm::init(param);
  
  string filename = param->get("MakeROOTFileBB.RootFileName", string("my_track1.root"));
  
  _file = new TFile(filename.c_str(), "RECREATE");
  
  _ntp = new TTree("track0", "track0");

  nEvt = 0;
  ntr1Trk = 0;
  ntr2Trk = 0;

  TracksData& d = _data;
  _ntp->Branch("nevent", &nEvt, "nevent/I"); // Event number
  _ntp->Branch("ntr1track", &ntr1Trk, "ntr1track/I"); // Event number
  _ntp->Branch("ntr2track", &ntr2Trk, "ntr2track/I"); // Event number
  // Track1 low data
  _ntp->Branch("tr1d0", &d.tr1d0, "tr1d0/D"); // Impact paramter of the track in (r-phi)
  _ntp->Branch("tr1z0", &d.tr1z0, "tr1z0/D"); // mpact paramter of the track in (r-z)
  _ntp->Branch("tr1phi", &d.tr1phi, "tr1phi/D"); // Phi of the track at the reference point
  _ntp->Branch("tr1omega", &d.tr1omega, "tr1omega/D"); // Omega is the signed curvature of the track in [1/mm]
  _ntp->Branch("tr1tanlam", &d.tr1tanlam, "tr1tanlam/D"); // Lambda is the dip angle of the track in r-z at the reference point
  _ntp->Branch("tr1x", &d.tr1x, "tr1x/D");
  _ntp->Branch("tr1y", &d.tr1y, "tr1y/D");
  _ntp->Branch("tr1z", &d.tr1z, "tr1z/D");
  _ntp->Branch("tr1charge", &d.tr1charge, "tr1charge/I");
  _ntp->Branch("tr1energy", &d.tr1energy, "tr1energy/D");
  _ntp->Branch("tr1covmatrixd0d0", &d.tr1covmatrixd0d0, "tr1covmatrixd0d0/D"); // Covariance matrix of the track parameters
  _ntp->Branch("tr1covmatrixd0z0", &d.tr1covmatrixd0z0, "tr1covmatrixd0z0/D"); // Covariance matrix of the track parameters
  _ntp->Branch("tr1covmatrixd0ph", &d.tr1covmatrixd0ph, "tr1covmatrixd0ph/D"); // Covariance matrix of the track parameters
  _ntp->Branch("tr1covmatrixd0om", &d.tr1covmatrixd0om, "tr1covmatrixd0om/D"); // Covariance matrix of the track parameters
  _ntp->Branch("tr1covmatrixd0tl", &d.tr1covmatrixd0tl, "tr1covmatrixd0tl/D"); // Covariance matrix of the track parameters
  _ntp->Branch("tr1covmatrixz0z0", &d.tr1covmatrixz0z0, "tr1covmatrixz0z0/D"); // Covariance matrix of the track parameters
  _ntp->Branch("tr1covmatrixz0ph", &d.tr1covmatrixz0ph, "tr1covmatrixz0ph/D"); // Covariance matrix of the track parameters
  _ntp->Branch("tr1covmatrixz0om", &d.tr1covmatrixz0om, "tr1covmatrixz0om/D"); // Covariance matrix of the track parameters
  _ntp->Branch("tr1covmatrixz0tl", &d.tr1covmatrixz0tl, "tr1covmatrixz0tl/D"); // Covariance matrix of the track parameters
  _ntp->Branch("tr1covmatrixphph", &d.tr1covmatrixphph, "tr1covmatrixphph/D"); // Covariance matrix of the track parameters
  _ntp->Branch("tr1covmatrixphom", &d.tr1covmatrixphom, "tr1covmatrixphom/D"); // Covariance matrix of the track parameters
  _ntp->Branch("tr1covmatrixphtl", &d.tr1covmatrixphtl, "tr1covmatrixphtl/D"); // Covariance matrix of the track parameters
  _ntp->Branch("tr1covmatrixomom", &d.tr1covmatrixomom, "tr1covmatrixomom/D"); // Covariance matrix of the track parameters
  _ntp->Branch("tr1covmatrixomtl", &d.tr1covmatrixomtl, "tr1covmatrixomtl/D"); // Covariance matrix of the track parameters
  _ntp->Branch("tr1covmatrixtltl", &d.tr1covmatrixtltl, "tr1covmatrixtltl/D"); // Covariance matrix of the track parameters
  // ======================================================= //
  // Track2 low data
  _ntp->Branch("tr2d0", &d.tr2d0, "tr2d0/D"); // Impact paramter of the track in (r-phi)
  _ntp->Branch("tr2z0", &d.tr2z0, "tr2z0/D"); // mpact paramter of the track in (r-z)
  _ntp->Branch("tr2phi", &d.tr2phi, "tr2phi/D"); // Phi of the track at the reference point
  _ntp->Branch("tr2omega", &d.tr2omega, "tr2omega/D"); // Omega is the signed curvature of the track in [1/mm]
  _ntp->Branch("tr2tanlam", &d.tr2tanlam, "tr2tanlam/D"); // Lambda is the dip angle of the track in r-z at the reference point
  _ntp->Branch("tr2x", &d.tr2x, "tr2x/D");
  _ntp->Branch("tr2y", &d.tr2y, "tr2y/D");
  _ntp->Branch("tr2z", &d.tr2z, "tr2z/D");
  _ntp->Branch("tr2charge", &d.tr2charge, "tr2charge/I");
  _ntp->Branch("tr2energy", &d.tr2energy, "tr2energy/D");
  _ntp->Branch("tr2covmatrixd0d0", &d.tr2covmatrixd0d0, "tr2covmatrixd0d0/D"); // Covariance matrix of the track parameters
  _ntp->Branch("tr2covmatrixd0z0", &d.tr2covmatrixd0z0, "tr2covmatrixd0z0/D"); // Covariance matrix of the track parameters
  _ntp->Branch("tr2covmatrixd0ph", &d.tr2covmatrixd0ph, "tr2covmatrixd0ph/D"); // Covariance matrix of the track parameters
  _ntp->Branch("tr2covmatrixd0om", &d.tr2covmatrixd0om, "tr2covmatrixd0om/D"); // Covariance matrix of the track parameters
  _ntp->Branch("tr2covmatrixd0tl", &d.tr2covmatrixd0tl, "tr2covmatrixd0tl/D"); // Covariance matrix of the track parameters
  _ntp->Branch("tr2covmatrixz0z0", &d.tr2covmatrixz0z0, "tr2covmatrixz0z0/D"); // Covariance matrix of the track parameters
  _ntp->Branch("tr2covmatrixz0ph", &d.tr2covmatrixz0ph, "tr2covmatrixz0ph/D"); // Covariance matrix of the track parameters
  _ntp->Branch("tr2covmatrixz0om", &d.tr2covmatrixz0om, "tr2covmatrixz0om/D"); // Covariance matrix of the track parameters
  _ntp->Branch("tr2covmatrixz0tl", &d.tr2covmatrixz0tl, "tr2covmatrixz0tl/D"); // Covariance matrix of the track parameters
  _ntp->Branch("tr2covmatrixphph", &d.tr2covmatrixphph, "tr2covmatrixphph/D"); // Covariance matrix of the track parameters
  _ntp->Branch("tr2covmatrixphom", &d.tr2covmatrixphom, "tr2covmatrixphom/D"); // Covariance matrix of the track parameters
  _ntp->Branch("tr2covmatrixphtl", &d.tr2covmatrixphtl, "tr2covmatrixphtl/D"); // Covariance matrix of the track parameters
  _ntp->Branch("tr2covmatrixomom", &d.tr2covmatrixomom, "tr2covmatrixomom/D"); // Covariance matrix of the track parameters
  _ntp->Branch("tr2covmatrixomtl", &d.tr2covmatrixomtl, "tr2covmatrixomtl/D"); // Covariance matrix of the track parameters
  _ntp->Branch("tr2covmatrixtltl", &d.tr2covmatrixtltl, "tr2covmatrixtltl/D"); // Covariance matrix of the track parameters
  // ======================================================= //
  // Track1 MC particle data
  _ntp->Branch("tr1mcx", &d.tr1mcx, "tr1mcx/D"); // the production vertex of the particle in [mm].
  _ntp->Branch("tr1mcy", &d.tr1mcy, "tr1mcy/D"); // the production vertex of the particle in [mm].
  _ntp->Branch("tr1mcz", &d.tr1mcz, "tr1mcz/D"); // the production vertex of the particle in [mm].
  _ntp->Branch("tr1id", &d.tr1id, "tr1id/I"); // the number of MC ID for this particle - 0 if other
  _ntp->Branch("tr1pdg", &d.tr1pdg, "tr1pdg/I"); // the number of MC PDG for this particle - 0 if other
  _ntp->Branch("tr1ssid", &tr1ssid, "tr1ssid/I"); // the number of MC ID for ss parent particle - 0 if other
  _ntp->Branch("tr1sspdg", &tr1sspdg, "tr1sspdg/I"); // the number of MC PDG for ss parent particle - 0 if other
  _ntp->Branch("tr1ssc", &tr1ssc, "tr1ssc/I"); // semistable c : 1, not : 0
  _ntp->Branch("tr1ssb", &tr1ssb, "tr1ssb/I"); // semistable b : 1, not : 0
  _ntp->Branch("tr1oth", &tr1oth, "tr1oth/I"); // other : 1, not : 0
  _ntp->Branch("tr1pri", &tr1pri, "tr1pri/I"); // primary : 1, not : 0
  // ======================================================= //
  // Track2 MC particle data
  _ntp->Branch("tr2mcx", &d.tr2mcx, "tr2mcx/D"); // the production vertex of the particle in [mm].
  _ntp->Branch("tr2mcy", &d.tr2mcy, "tr2mcy/D"); // the production vertex of the particle in [mm].
  _ntp->Branch("tr2mcz", &d.tr2mcz, "tr2mcz/D"); // the production vertex of the particle in [mm].
  _ntp->Branch("tr2id", &d.tr2id, "tr2id/I"); // the number of MC ID for this particle - 0 if other
  _ntp->Branch("tr2pdg", &d.tr2pdg, "tr2pdg/I"); // the number of MC PDG for this particle - 0 if other
  _ntp->Branch("tr2ssid", &tr2ssid, "tr2ssid/I"); // the number of MC ID for ss parent particle - 0 if other
  _ntp->Branch("tr2sspdg", &tr2sspdg, "tr2sspdg/I"); // the number of MC PDG for ss parent particle - 0 if other
  _ntp->Branch("tr2ssc", &tr2ssc, "tr2ssc/I"); // semistable c : 1, not : 0
  _ntp->Branch("tr2ssb", &tr2ssb, "tr2ssb/I"); // semistable b : 1, not : 0
  _ntp->Branch("tr2oth", &tr2oth, "tr2oth/I"); // other : 1, not : 0
  _ntp->Branch("tr2pri", &tr2pri, "tr2pri/I"); // primary : 1, not : 0
  // ======================================================= //
  _ntp->Branch("connect", &d.connect, "connect/I"); // (0 : not connect), (1 : primary), (2 : cc), (3 : bb), (4 : bc)
  _ntp->Branch("lcfiplustag", &d.lcfiplustag, "lcfiplustag/I"); // 0 : not connect, 1 : semistablecparent, 2 : primary vertex
  // ======================================================= //
  // Parameters From Vertex Finder Suehara
  // Fitter Feature Value
  _ntp->Branch("chi2", &d.chi2, "chi2/D"); // chi2
  _ntp->Branch("vchi2", &d.vchi2, "vchi2/D"); // vertex chi2
  _ntp->Branch("vposx", &d.vposx, "vposx/D"); // vertex position x
  _ntp->Branch("vposy", &d.vposy, "vposy/D"); // vertex position y
  _ntp->Branch("vposz", &d.vposz, "vposz/D"); // vertex position z
  _ntp->Branch("mass", &d.mass, "mass/D"); // mass
  _ntp->Branch("mag", &d.mag, "mag/D"); // Mag
  _ntp->Branch("vec", &d.vec, "vec/D"); // Vec direction
  _ntp->Branch("tr1selection", &d.tr1selection, "tr1selection/I"); // track selection
  _ntp->Branch("tr2selection", &d.tr2selection, "tr2selection/I"); // track selection
  _ntp->Branch("v0selection", &d.v0selection, "v0selection/I"); // v0 selection
}

void MakeROOTFileBB::process() {

  Event* event = Event::Instance();
  const TrackVec& track_list = event->getTracks();
	  
  for (unsigned int i=0; i < track_list.size(); ++i) {

    int tr1ssssid = -999;
    int tr1sssspdg = -999;
    int tr2ssssid = -999;
    int tr2sssspdg = -999;
    
    ntr1Trk = i;
    const Track* track1 = track_list[i];
    if (!track1) continue;
    const MCParticle* mcpc1 = event->getMCParticle(track1);
    _hel1 = new lcfiplus::Helix(track1, PointBase::NOTUSED);

    memset(&_data,0,sizeof(_data));

    _data.tr1d0 = track1->getD0();
    _data.tr1z0 = track1->getZ0();
    _data.tr1phi = track1->getPhi();
    _data.tr1omega = track1->getOmega();
    _data.tr1tanlam = track1->getTanLambda();

    _data.tr1charge = (int)track1->getCharge();

    _data.tr1x = _hel1->GetPos(0.).X();
    _data.tr1y = _hel1->GetPos(0.).Y();
    _data.tr1z = _hel1->GetPos(0.).Z();

    const double* cov1 = track1->getCovMatrix();
    
    _data.tr1covmatrixd0d0 = cov1[0]; // d0d0
    _data.tr1covmatrixd0z0 = cov1[1]; // d0z0
    _data.tr1covmatrixd0ph = cov1[2]; // d0ph
    _data.tr1covmatrixd0om = cov1[3]; // d0om
    _data.tr1covmatrixd0tl = cov1[4]; // d0tl
    _data.tr1covmatrixz0z0 = cov1[5]; // z0z0
    _data.tr1covmatrixz0ph = cov1[6]; // z0ph
    _data.tr1covmatrixz0om = cov1[7]; // z0om
    _data.tr1covmatrixz0tl = cov1[8]; // z0tl
    _data.tr1covmatrixphph = cov1[9]; // phph
    _data.tr1covmatrixphom = cov1[10]; // phom
    _data.tr1covmatrixphtl = cov1[11]; // phtl
    _data.tr1covmatrixomom = cov1[12]; // omom
    _data.tr1covmatrixomtl = cov1[13]; // omtl
    _data.tr1covmatrixtltl = cov1[14]; // tltl
    
    if (mcpc1){
      _data.tr1mcx = mcpc1->getVertex().X();
      _data.tr1mcy = mcpc1->getVertex().Y();
      _data.tr1mcz = mcpc1->getVertex().Z();
      _data.tr1id = mcpc1->getId();
      _data.tr1pdg = mcpc1->getPDG();

      const MCParticle* ssmcpc1 = mcpc1->getSemiStableParent();

      if (ssmcpc1==0){
	tr1ssc = 0;
	tr1ssb = 0;
	tr1oth = 0;
	tr1pri = 1;
	tr1ssid = 0;
	tr1sspdg = 0;
      }
      if (ssmcpc1){
        tr1ssid = ssmcpc1->getId();
        tr1sspdg = ssmcpc1->getPDG();
        if (ssmcpc1->getParent()->getId()==0){
	  tr1ssc = 0;
	  tr1ssb = 0;
	  tr1oth = 0;
	  tr1pri = 1;
        }
	else if (ssmcpc1->isSemiStableC()){
	  tr1ssc = 1;
	  tr1ssb = 0;
	  tr1oth = 0;
	  tr1pri = 0;
	  tr1ssssid = ssmcpc1->getSemiStableParent()->getId();
	  tr1sssspdg = ssmcpc1->getSemiStableParent()->getPDG();
        }
	else if (ssmcpc1->isSemiStableB()){
	  tr1ssc = 0;
	  tr1ssb = 1;
	  tr1oth = 0;
	  tr1pri = 0;
        }
        else {
	  tr1ssc = 0;
	  tr1ssb = 0;
	  tr1oth = 1;
	  tr1pri = 0;
        }
	/*
	cout << "SSC : " << tr1ssc << ", ID : " << ssmcpc1->getId() << ", PDG : " << ssmcpc1->getPDG() << 
	", Parent : " << ssmcpc1->getParent()->getId() << endl;
	*/
      }
    }

    for (unsigned int j=0; j < i; ++j) {
      ntr2Trk = j;
      const Track* track2 = track_list[j];
      if (!track2) continue;

      const MCParticle* mcpc2 = event->getMCParticle(track2);
      _hel2 = new lcfiplus::Helix(track2, PointBase::NOTUSED);

      _data.tr2d0 = track2->getD0();
      _data.tr2z0 = track2->getZ0();
      _data.tr2phi = track2->getPhi();
      _data.tr2omega = track2->getOmega();
      _data.tr2tanlam = track2->getTanLambda();

      _data.tr2charge = (int)track2->getCharge();

      _data.tr2x = _hel2->GetPos(0.).X();
      _data.tr2y = _hel2->GetPos(0.).Y();
      _data.tr2z = _hel2->GetPos(0.).Z();

      const double* cov2 = track2->getCovMatrix();
    
      _data.tr2covmatrixd0d0 = cov2[0]; // d0d0
      _data.tr2covmatrixd0z0 = cov2[1]; // d0z0
      _data.tr2covmatrixd0ph = cov2[2]; // d0ph
      _data.tr2covmatrixd0om = cov2[3]; // d0om
      _data.tr2covmatrixd0tl = cov2[4]; // d0tl
      _data.tr2covmatrixz0z0 = cov2[5]; // z0z0
      _data.tr2covmatrixz0ph = cov2[6]; // z0ph
      _data.tr2covmatrixz0om = cov2[7]; // z0om
      _data.tr2covmatrixz0tl = cov2[8]; // z0tl
      _data.tr2covmatrixphph = cov2[9]; // phph
      _data.tr2covmatrixphom = cov2[10]; // phom
      _data.tr2covmatrixphtl = cov2[11]; // phtl
      _data.tr2covmatrixomom = cov2[12]; // omom
      _data.tr2covmatrixomtl = cov2[13]; // omtl
      _data.tr2covmatrixtltl = cov2[14]; // tltl
      
      if (mcpc2){
        _data.tr2mcx = mcpc2->getVertex().X();
        _data.tr2mcy = mcpc2->getVertex().Y();
        _data.tr2mcz = mcpc2->getVertex().Z();
        _data.tr2id = mcpc2->getId();
        _data.tr2pdg = mcpc2->getPDG();

        const MCParticle* ssmcpc2 = mcpc2->getSemiStableParent();

        if (ssmcpc2==0){
	  tr2ssc = 0;
	  tr2ssb = 0;
	  tr2oth = 0;
	  tr2pri = 1;
	  tr2ssid = 0;
	  tr2sspdg = 0;
        }
        if (ssmcpc2){
          tr2ssid = ssmcpc2->getId();
          tr2sspdg = ssmcpc2->getPDG();
          if (ssmcpc2->getParent()->getId()==0){
	    tr2ssc = 0;
	    tr2ssb = 0;
	    tr2oth = 0;
	    tr2pri = 1;
          }
	  else if (ssmcpc2->isSemiStableC()){
	    tr2ssc = 1;
	    tr2ssb = 0;
	    tr2oth = 0;
	    tr2pri = 0;
	    tr2ssssid = ssmcpc2->getSemiStableParent()->getId();
	    tr2sssspdg = ssmcpc2->getSemiStableParent()->getPDG();
          }
	  else if (ssmcpc2->isSemiStableB()){
	    tr2ssc = 0;
	    tr2ssb = 1;
	    tr2oth = 0;
	    tr2pri = 0;
          }
          else {
	    tr2ssc = 0;
	    tr2ssb = 0;
	    tr2oth = 1;
	    tr2pri = 0;
          }
        }
      }

      if (((tr1ssc==1)&&(tr2ssc==1))&&(tr1ssid==tr2ssid)) _data.connect = 2; // cc
      else if (((tr1ssb==1)&&(tr2ssb==1))&&(tr1ssid==tr2ssid)) _data.connect = 3; // bb
      else if ((((tr1ssc==1)&&(tr2ssb==1))&&(tr1ssssid==tr2ssid))
	     ||(((tr1ssb==1)&&(tr2ssc==1))&&(tr1ssid==tr2ssssid))) _data.connect = 4; // bc
      else if ((tr1pri==1)&&(tr2pri==1)) _data.connect = 1; // primary
      else if (((tr1oth==1)&&(tr2oth==1))&&(tr1ssid==tr2ssid)) _data.connect = 5; // v0
      else _data.connect = 0; // others

      vector<const Track*> vttmp;
      vttmp.push_back(track1);
      vttmp.push_back(track2);

      Vertex* vtx = VertexFitterSimple_V() (vttmp.begin(), vttmp.end(), 0);

      _data.chi2 = max(vtx->getChi2Track(track1), vtx->getChi2Track(track2));
      _data.vchi2 = vtx->getChi2();
      /*
      cout << "max(vtx->getChi2Track()) " << _data.chi2 << " : " << max(vtx->getChi2Track(track1), vtx->getChi2Track(track2)) << 
	      " : vtx->getChi2 " << _data.vchi2 << " : " << vtx->getChi2() << endl;
      */      
      _data.vposx = vtx->getPos().X();
      _data.vposy = vtx->getPos().Y();
      _data.vposz = vtx->getPos().Z();

      TLorentzVector tr1tlv = *track1;
      TLorentzVector tr2tlv = *track2;

      _data.mass = (tr1tlv+tr2tlv).M();
      _data.mag = vtx->getPos().Mag();
      _data.vec = vtx->getPos().Dot((tr1tlv+tr2tlv).Vect());

      _data.tr1energy = tr1tlv.E();
      _data.tr2energy = tr2tlv.E();

      
      // suehara 200319
      // TrackSelector and VertexSelector config
      TrackSelectorConfig cfgtrack;
      cfgtrack.maxD0 = 10.;
      cfgtrack.maxZ0 = 20.;
      cfgtrack.minPt = 0.1;
      cfgtrack.maxInnermostHitRadius = 1e10;
    
      cfgtrack.maxD0Err = 0.1;
      cfgtrack.maxZ0Err = 0.1;

      cfgtrack.minTpcHits = 20;
      cfgtrack.minFtdHits = 3;
      cfgtrack.minVtxHits = 3;
      cfgtrack.minVtxPlusFtdHits = 1;
    
      VertexFinderSuehara::VertexFinderSueharaConfig cfgvtx;
      // chi2 threshold = 9.
      // mass threshold = 10.
      // minimum distance = 0.3;
      
      TrackSelector sel;
    
      _data.lcfiplustag = 0; // 0: NC, 1: secondary, 2: primary
      const Vertex *ip = event->getPrimaryVertex();
      TrackVec &iptracks = ip->getTracks();
      if(std::find(iptracks.begin(), iptracks.end(), track1) != iptracks.end()
        && std::find(iptracks.begin(), iptracks.end(), track2) != iptracks.end())
          _data.lcfiplustag = 2; // primary tracks
      else{
        // track selection
        if(sel.passesCut(track1, cfgtrack) && sel.passesCut(track2, cfgtrack)){
          // mass, chi2, vpos, v0 selection
          TLorentzVector v1 = *track1;
          TLorentzVector v2 = *track2;
          double mass = (v1+v2).M();
          double chi2 = max(vtx->getChi2Track(track1), vtx->getChi2Track(track2));
          
          if(mass < min(v1.E(), v2.E()) && mass < 10. && chi2 < 9. && vtx->getPos().Mag() > 0.3
            && vtx->getPos().Dot((v1+v2).Vect()) > 0 && !VertexSelector().passesCut(vtx, cfgvtx.v0selVertex, ip)){
            _data.lcfiplustag = 1; // secondary tracks
          }
        }
      }

      if(sel.passesCut(track1, cfgtrack)) _data.tr1selection = 1;
      else _data.tr1selection = 0;

      if(sel.passesCut(track2, cfgtrack)) _data.tr2selection = 1;
      else _data.tr2selection = 0;

      if(VertexSelector().passesCut(vtx, cfgvtx.v0selVertex, ip))  _data.v0selection = 1;
      else _data.v0selection = 0;
      
      _ntp->Fill();
    }

  }
  nEvt++;
}

void MakeROOTFileBB::end() {
  _file->Write();
  _file->Close();
}

}
