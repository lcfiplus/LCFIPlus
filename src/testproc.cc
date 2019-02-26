#include <string>

#include "TFile.h"
#include "TNtuple.h"
#include "TNtupleD.h"
#include "TSystem.h"
#include "TPad.h"
#include "TStyle.h"

#include "lcfiplus.h"
#include "process.h"
#include "testproc.h"
#include "VertexSelector.h"
#include "algoEtc.h"
#include "VertexFinderSuehara.h"
#include "VertexFitterSimple.h"

using namespace lcfiplus;

namespace lcfiplus {

const Jet* JetMCMatch(JetVec& jets, const MCParticle* mcp, vector<const Track*>& assignedTracks, vector<const Track*>& residualTracks) {
  const vector<const Track*>* pTracks;
  pTracks = &(Event::Instance()->getTracks());

  vector<const Track*> bTracks;

  vector<int> nTrackInJet;
  nTrackInJet.resize(jets.size());
  vector<int> nVertexTrackInJet;
  nVertexTrackInJet.resize(jets.size());

  int nvtx = 0;
  // get tracks
  for (unsigned int i=0; i<pTracks->size(); i++) {
    const MCParticle* mcpc = (*pTracks)[i]->getMcp();

    if (mcpc==0)continue;
    if (mcpc->isParent(mcp)) {
      bTracks.push_back((*pTracks)[i]);
      for (unsigned int j=0; j<jets.size(); j++) {
        for (unsigned int k=0; k<jets[j]->getVertices().size(); k++) {
          const vector<const Track*>& vtr = jets[j]->getVertices()[k]->getTracks();
          if (find(vtr.begin(), vtr.end(), (*pTracks)[i]) != vtr.end()) {
            // the track matched to this jet
            nTrackInJet[j] ++;
// 						if(nvtx == 0 && k > 0){cout << "CAUTION: vertices in the same jet might be from different semistables!" << endl;}
// 						if(nvtx > 0 && k == 0){cout << "CAUTION: vertices in different jets might be from the same samistable!" << endl;}
            nVertexTrackInJet[j] ++;
            nvtx ++;
          }
        }
        if (find(jets[j]->getTracks().begin(), jets[j]->getTracks().end(), (*pTracks)[i]) != jets[j]->getTracks().end()) {
          // the track matched to this jet
          nTrackInJet[j] ++;
        }
      }
    }
  }

  int ntijMax = 0;
  int ntijMaxIndex = -1;

  int ntrsum = 0;
  int ntrvtxsum = 0;
  // determine best-match jet
  for (unsigned int j=0; j<jets.size(); j++) {
    if (ntijMax < nTrackInJet[j]) {
      ntijMax = nTrackInJet[j];
      ntijMaxIndex = j;
    }
    ntrsum += nTrackInJet[j];
    ntrvtxsum += nVertexTrackInJet[j];
  }

  if (ntijMaxIndex == -1)
    return 0;

  vector<const Track*> jetTracks = jets[ntijMaxIndex]->getTracks();
  for (unsigned int i=0; i<jets[ntijMaxIndex]->getVertices().size(); i++) {
    const Vertex* vtx = jets[ntijMaxIndex]->getVertices()[i];
    jetTracks.insert(jetTracks.end(), vtx->getTracks().begin(), vtx->getTracks().end());
  }

  // obtain assignedtracks
  for (unsigned int i=0; i<bTracks.size(); i++) {
    if (find(jetTracks.begin(), jetTracks.end(), bTracks[i]) != jetTracks.end())
      assignedTracks.push_back(bTracks[i]);
    else
      residualTracks.push_back(bTracks[i]);
  }

  cout << "Assigned jet " << ntijMaxIndex << ", Vertex tracks: ";
  for (unsigned int j=0; j<jets.size(); j++)cout << nVertexTrackInJet[j] << ",";
  cout << "/" << ntrvtxsum << ", all tracks: ";
  for (unsigned int j=0; j<jets.size(); j++)cout << nTrackInJet[j] << ",";
  cout << "/" << bTracks.size() << ", PDG: " << mcp->getPDG() << endl;//"," ;
//	cout << " ( " << mcp->getVertex().x() << " " << mcp->getVertex().y() << " " << mcp->getVertex().z() << ")" << endl;

  return jets[ntijMaxIndex];
}

void TestAlgoV0::init(Parameters* param) {
  Algorithm::init(param);
  string filename = param->get("FileName",string("testv0.root"));
  _vtxname = param->get("VertexCollectionName",string("BuildUpVertex"));
  _file = new TFile(filename.c_str(),"RECREATE");
  _ntp = new TTree("v0","v0");
  _vertices = 0;

  VtxData& d = _data;
  _ntp->Branch("x",&d.x,"x/D");
  _ntp->Branch("y",&d.y,"y/D");
  _ntp->Branch("z",&d.z,"z/D");
  _ntp->Branch("r",&d.r,"r/D");
  _ntp->Branch("cs",&d.cs,"cs/D");
  _ntp->Branch("phi",&d.phi,"phi/D");
  _ntp->Branch("chrg",&d.chrg,"chrg/D");
  _ntp->Branch("dirdot",&d.dirdot,"dirdot/D");
  _ntp->Branch("dirdot2",&d.dirdot2,"dirdot2/D");
  _ntp->Branch("ntrk",&d.ntrk,"ntrk/I");
  _ntp->Branch("mks",&d.mks,"mks/D");
  _ntp->Branch("ml0",&d.ml0,"ml0/D");
  _ntp->Branch("mconv",&d.mconv,"mconv/D");
  _ntp->Branch("mks2",&d.mks2,"mks2/D");
  _ntp->Branch("ml02",&d.ml02,"ml02/D");
  _ntp->Branch("v0",&d.v0,"v0/I");
  _ntp->Branch("ks",&d.ks,"ks/I");
  _ntp->Branch("l0",&d.l0,"l0/I");
  _ntp->Branch("conv",&d.conv,"conv/I");
  _ntp->Branch("mcpdg1",&d.mcpdg1,"mcpdg1/I");
  _ntp->Branch("mcpdg2",&d.mcpdg2,"mcpdg2/I");
  _ntp->Branch("mcppdg1",&d.mcppdg1,"mcppdg1/I");
  _ntp->Branch("mcppdg2",&d.mcppdg2,"mcppdg2/I");
  _ntp->Branch("mcpp1",&d.mcpp1,"mcpp1/I");
  _ntp->Branch("mcpp2",&d.mcpp2,"mcpp2/I");
}

void TestAlgoV0::process() {
  if (!_vertices) {
    Event::Instance()->Get(_vtxname.c_str(), _vertices);
  }

  const VertexVec& vtx_list = *_vertices;
  for (unsigned int i=0; i < vtx_list.size(); ++i) {
    const Vertex* vtx = vtx_list[i];
    if (vtx->isPrimary()) continue;

    memset(&_data,0,sizeof(_data));

    _data.x = vtx->getX();
    _data.y = vtx->getY();
    _data.z = vtx->getZ();
    TVector3 pos = vtx->getPos();
    _data.r = pos.Mag();
    _data.cs = pos.CosTheta();
    _data.phi = pos.Phi();

    TVector3 mom;
    TVector3 mom2;
    for (unsigned int j=0; j<vtx->getTracks().size(); ++j) {
      mom += vtx->getTracks()[j]->Vect();
      mom2 += vtx->getTracks()[j]->momentumAtVertex(vtx);
      _data.chrg += vtx->getTracks()[j]->getCharge();
    }

    _data.dirdot = mom.Unit().Dot( pos.Unit() );
    _data.dirdot2 = mom2.Unit().Dot( pos.Unit() );
    _data.ntrk = vtx->getTracks().size();

    // compute ks mass
    if (_data.ntrk == 2 && _data.chrg == 0) {
      const Track* trk1 = vtx->getTracks()[0];
      const Track* trk2 = vtx->getTracks()[1];
      TVector3 mom1 = trk1->Vect();
      TVector3 mom2 = trk2->Vect();
      TVector3 mom1v = trk1->momentumAtVertex(vtx);
      TVector3 mom2v = trk2->momentumAtVertex(vtx);

      TLorentzVector lvec1;
      TLorentzVector lvec2;
      lvec1.SetVectM( mom1, 0.1396 );
      lvec2.SetVectM( mom2, 0.1396 );
      _data.mks = (lvec1+lvec2).M();

      lvec1.SetVectM( mom1v, 0.1396 );
      lvec2.SetVectM( mom2v, 0.1396 );
      _data.mks2 = (lvec1+lvec2).M();

      // compute l0 mass
      TLorentzVector protonForLambda;
      TLorentzVector pionForLambda;
      if (mom1.Mag() > mom2.Mag()) {
        protonForLambda.SetVectM( mom1, 0.9383 );
        pionForLambda.SetVectM( mom2, 0.1396 );
      } else {
        protonForLambda.SetVectM( mom2, 0.9383 );
        pionForLambda.SetVectM( mom1, 0.1396 );
      }
      _data.ml0 = (protonForLambda+pionForLambda).M();

      if (mom1v.Mag() > mom2v.Mag()) {
        protonForLambda.SetVectM( mom1v, 0.9383 );
        pionForLambda.SetVectM( mom2v, 0.1396 );
      } else {
        protonForLambda.SetVectM( mom2v, 0.9383 );
        pionForLambda.SetVectM( mom1v, 0.1396 );
      }
      _data.ml02 = (protonForLambda+pionForLambda).M();

      // compute photon mass
      double ang1 = atan( trk1->getTanLambda() );
      double ang2 = atan( trk2->getTanLambda() );
      _data.mconv = sqrt( mom1.Mag()*mom2.Mag()*(1-cos(ang1-ang2)) );


      const MCParticle* mcp1 = trk1->getMcp();
      const MCParticle* mcp2 = trk2->getMcp();

      if (mcp1 && mcp2) {
        _data.mcpdg1 = mcp1->getPDG();
        _data.mcpdg2 = mcp2->getPDG();
        _data.mcppdg1 = mcp1->getParent()->getPDG();
        _data.mcppdg2 = mcp2->getParent()->getPDG();
      }

      if (mcp1 && mcp2) {
        const MCParticle* parent1 = mcp1->getParent();
        const MCParticle* parent2 = mcp2->getParent();
        const MCParticle* parent = mcp1->getParent();

        _data.mcpp1 = (int)((long long) parent1 );
        _data.mcpp2 = (int)((long long) parent2 );

        if ( abs(mcp1->getPDG())==11 && abs(mcp2->getPDG())==11 && parent->getPDG()==22 ) _data.conv = 1;
        if ( abs(mcp1->getPDG())==211 && abs(mcp2->getPDG())==211 && parent->getPDG()==310 ) _data.ks = 1;
        if ( abs(mcp1->getPDG())==211 && abs(mcp2->getPDG())==2212 && abs(parent->getPDG())==3122 ) _data.l0 = 1;
        if ( abs(mcp2->getPDG())==211 && abs(mcp1->getPDG())==2212 && abs(parent->getPDG())==3122 ) _data.l0 = 1;
      }

      _data.v0 = _data.ks || _data.l0 || _data.conv;
    }
    _ntp->Fill();
  }
}

void TestAlgoV0::end() {
  _file->Write();
  _file->Close();
}

void ZHHAlgo::init(Parameters* param) {
  Algorithm::init(param);

  string filename = param->get("FileName",string("test.root"));
  _jetname8 = param->get("JetCollectionName8",string("RefinedJets_8"));
  _jetname7 = param->get("JetCollectionName7",string("RefinedJets_7"));
  _jetname  = param->get("JetCollectionName6",string("RefinedJets_6"));
  _jetname5 = param->get("JetCollectionName5",string("RefinedJets_5"));
  _jetname4 = param->get("JetCollectionName4",string("RefinedJets_4"));

  _jetnamenv8 = param->get("JetCollectionNameNV8",string("RefinedNJets_8"));
  _jetnamenv7 = param->get("JetCollectionNameNV7",string("RefinedNJets_7"));
  _jetnamenv6 = param->get("JetCollectionNameNV6",string("RefinedNJets_6"));
  _jetnamenv5 = param->get("JetCollectionNameNV5",string("RefinedNJets_5"));
  _jetnamenv4 = param->get("JetCollectionNameNV4",string("RefinedNJets_4"));

  _file = new TFile(filename.c_str(),"RECREATE");
  _tree = new TTree("tree","tree");

  // mc info
  _tree->Branch("mchdecaypdg",&_d.mchdecaypdg,"mchdecaypdg[2]/I");
  _tree->Branch("mchbb",&_d.mchbb,"mchbb/I");
  _tree->Branch("mcnb",&_d.mcnb,"mcnb/I");

  // non-jet variables
  _tree->Branch("ycuts",&_d.ycuts,"ycuts[10]/D");

  _tree->Branch("thrust",&_d.thrust,"thrust/D");
  _tree->Branch("thaxis",&_d.thaxis,"thaxis[3]/D");
  _tree->Branch("ntr", &_d.ntr, "ntr/I");
  _tree->Branch("npfo", &_d.npfo, "ntr/I");

  // combined variables for compatibility
  _tree->Branch("mass",&_d.mass,"mass[15]/D");
  _tree->Branch("ntrjetmin",&_d.ntrjetmin,"ntrjetmin/D");
  _tree->Branch("pmiss",&_d.pmiss,"pmiss[3]/D");
  _tree->Branch("emiss",&_d.emiss,"emiss/D");

  // 6-jet variables
  _tree->Branch("bcat",&_d.bcat,"bcat[6]/D");
  _tree->Branch("btag",&_d.btag,"btag[6]/D");
  _tree->Branch("ctag",&_d.ctag,"ctag[6]/D");
  _tree->Branch("ejet",&_d.ejet,"ejet[6]/D");
  _tree->Branch("pxjet",&_d.pxjet,"pxjet[6]/D");
  _tree->Branch("pyjet",&_d.pyjet,"pyjet[6]/D");
  _tree->Branch("pzjet",&_d.pzjet,"pzjet[6]/D");
  _tree->Branch("ntrjet",&_d.ntrjet,"ntrjet[6]/D");
  _tree->Branch("mcnb6",&_d.mcnb6,"mcnb6[6]/D");
  _tree->Branch("mcnc6",&_d.mcnc6,"mcnc6[6]/D");
  _tree->Branch("twovtxprobjet",&_d.twovtxprobjet,"twovtxprobjet[6]/D");
  _tree->Branch("vtxangle",&_d.vtxangle,"vtxangle[6]/D");

  // 4-jet variables
  _tree->Branch("bcat4",&_d.bcat4,"bcat4[4]/D");
  _tree->Branch("btag4",&_d.btag4,"btag4[4]/D");
  _tree->Branch("ctag4",&_d.ctag4,"ctag4[4]/D");
  _tree->Branch("ejet4",&_d.ejet4,"ejet4[4]/D");
  _tree->Branch("pxjet4",&_d.pxjet4,"pxjet4[4]/D");
  _tree->Branch("pyjet4",&_d.pyjet4,"pyjet4[4]/D");
  _tree->Branch("pzjet4",&_d.pzjet4,"pzjet4[4]/D");
  _tree->Branch("ntrjet4",&_d.ntrjet4,"ntrjet4[4]/D");
  _tree->Branch("twovtxprobjet4",&_d.twovtxprobjet4,"twovtxprobjet4[4]/D");
  _tree->Branch("vtxangle4",&_d.vtxangle4,"vtxangle4[4]/D");

  _tree->Branch("bcat5",&_d.bcat5,"bcat5[5]/D");
  _tree->Branch("btag5",&_d.btag5,"btag5[5]/D");
  _tree->Branch("ctag5",&_d.ctag5,"ctag5[5]/D");
  _tree->Branch("ejet5",&_d.ejet5,"ejet5[5]/D");
  _tree->Branch("pxjet5",&_d.pxjet5,"pxjet5[5]/D");
  _tree->Branch("pyjet5",&_d.pyjet5,"pyjet5[5]/D");
  _tree->Branch("pzjet5",&_d.pzjet5,"pzjet5[5]/D");
  _tree->Branch("ntrjet5",&_d.ntrjet5,"ntrjet5[5]/D");
  _tree->Branch("twovtxprobjet5",&_d.twovtxprobjet5,"twovtxprobjet5[5]/D");
  _tree->Branch("vtxangle5",&_d.vtxangle5,"vtxangle5[5]/D");

  _tree->Branch("bcat7",&_d.bcat7,"bcat7[7]/D");
  _tree->Branch("btag7",&_d.btag7,"btag7[7]/D");
  _tree->Branch("ctag7",&_d.ctag7,"ctag7[7]/D");
  _tree->Branch("ejet7",&_d.ejet7,"ejet7[7]/D");
  _tree->Branch("pxjet7",&_d.pxjet7,"pxjet7[7]/D");
  _tree->Branch("pyjet7",&_d.pyjet7,"pyjet7[7]/D");
  _tree->Branch("pzjet7",&_d.pzjet7,"pzjet7[7]/D");
  _tree->Branch("ntrjet7",&_d.ntrjet7,"ntrjet7[7]/D");
  _tree->Branch("twovtxprobjet7",&_d.twovtxprobjet7,"twovtxprobjet7[7]/D");
  _tree->Branch("vtxangle7",&_d.vtxangle7,"vtxangle7[7]/D");

  _tree->Branch("bcat8",&_d.bcat8,"bcat8[8]/D");
  _tree->Branch("btag8",&_d.btag8,"btag8[8]/D");
  _tree->Branch("ctag8",&_d.ctag8,"ctag8[8]/D");
  _tree->Branch("ejet8",&_d.ejet8,"ejet8[8]/D");
  _tree->Branch("pxjet8",&_d.pxjet8,"pxjet8[8]/D");
  _tree->Branch("pyjet8",&_d.pyjet8,"pyjet8[8]/D");
  _tree->Branch("pzjet8",&_d.pzjet8,"pzjet8[8]/D");
  _tree->Branch("ntrjet8",&_d.ntrjet8,"ntrjet8[8]/D");
  _tree->Branch("twovtxprobjet8",&_d.twovtxprobjet8,"twovtxprobjet8[8]/D");
  _tree->Branch("vtxangle8",&_d.vtxangle8,"vtxangle8[8]/D");

  // jet clustering with no vertex
  _tree->Branch("bcatnv4",&_d.bcatnv4,"bcatnv4[4]/D");
  _tree->Branch("btagnv4",&_d.btagnv4,"btagnv4[4]/D");
  _tree->Branch("ctagnv4",&_d.ctagnv4,"ctagnv4[4]/D");
  _tree->Branch("ejetnv4",&_d.ejetnv4,"ejetnv4[4]/D");
  _tree->Branch("pxjetnv4",&_d.pxjetnv4,"pxjetnv4[4]/D");
  _tree->Branch("pyjetnv4",&_d.pyjetnv4,"pyjetnv4[4]/D");
  _tree->Branch("pzjetnv4",&_d.pzjetnv4,"pzjetnv4[4]/D");
  _tree->Branch("ntrjetnv4",&_d.ntrjetnv4,"ntrjetnv4[4]/D");
  _tree->Branch("twovtxprobjetnv4",&_d.twovtxprobjetnv4,"twovtxprobjetnv4[4]/D");
  _tree->Branch("vtxanglenv4",&_d.vtxanglenv4,"vtxanglenv4[4]/D");

  _tree->Branch("bcatnv5",&_d.bcatnv5,"bcatnv5[5]/D");
  _tree->Branch("btagnv5",&_d.btagnv5,"btagnv5[5]/D");
  _tree->Branch("ctagnv5",&_d.ctagnv5,"ctagnv5[5]/D");
  _tree->Branch("ejetnv5",&_d.ejetnv5,"ejetnv5[5]/D");
  _tree->Branch("pxjetnv5",&_d.pxjetnv5,"pxjetnv5[5]/D");
  _tree->Branch("pyjetnv5",&_d.pyjetnv5,"pyjetnv5[5]/D");
  _tree->Branch("pzjetnv5",&_d.pzjetnv5,"pzjetnv5[5]/D");
  _tree->Branch("ntrjetnv5",&_d.ntrjetnv5,"ntrjetnv5[5]/D");
  _tree->Branch("twovtxprobjetnv5",&_d.twovtxprobjetnv5,"twovtxprobjetnv5[5]/D");
  _tree->Branch("vtxanglenv5",&_d.vtxanglenv5,"vtxanglenv5[5]/D");

  _tree->Branch("bcatnv6",&_d.bcatnv6,"bcatnv6[6]/D");
  _tree->Branch("btagnv6",&_d.btagnv6,"btagnv6[6]/D");
  _tree->Branch("ctagnv6",&_d.ctagnv6,"ctagnv6[6]/D");
  _tree->Branch("ejetnv6",&_d.ejetnv6,"ejetnv6[6]/D");
  _tree->Branch("pxjetnv6",&_d.pxjetnv6,"pxjetnv6[6]/D");
  _tree->Branch("pyjetnv6",&_d.pyjetnv6,"pyjetnv6[6]/D");
  _tree->Branch("pzjetnv6",&_d.pzjetnv6,"pzjetnv6[6]/D");
  _tree->Branch("ntrjetnv6",&_d.ntrjetnv6,"ntrjetnv6[6]/D");
  _tree->Branch("twovtxprobjetnv6",&_d.twovtxprobjetnv6,"twovtxprobjetnv6[6]/D");
  _tree->Branch("vtxanglenv6",&_d.vtxanglenv6,"vtxanglenv6[6]/D");

  _tree->Branch("bcatnv7",&_d.bcatnv7,"bcatnv7[7]/D");
  _tree->Branch("btagnv7",&_d.btagnv7,"btagnv7[7]/D");
  _tree->Branch("ctagnv7",&_d.ctagnv7,"ctagnv7[7]/D");
  _tree->Branch("ejetnv7",&_d.ejetnv7,"ejetnv7[7]/D");
  _tree->Branch("pxjetnv7",&_d.pxjetnv7,"pxjetnv7[7]/D");
  _tree->Branch("pyjetnv7",&_d.pyjetnv7,"pyjetnv7[7]/D");
  _tree->Branch("pzjetnv7",&_d.pzjetnv7,"pzjetnv7[7]/D");
  _tree->Branch("ntrjetnv7",&_d.ntrjetnv7,"ntrjetnv7[7]/D");
  _tree->Branch("twovtxprobjetnv7",&_d.twovtxprobjetnv7,"twovtxprobjetnv7[7]/D");
  _tree->Branch("vtxanglenv7",&_d.vtxanglenv7,"vtxanglenv7[7]/D");

  _tree->Branch("bcatnv8",&_d.bcatnv8,"bcatnv8[8]/D");
  _tree->Branch("btagnv8",&_d.btagnv8,"btagnv8[8]/D");
  _tree->Branch("ctagnv8",&_d.ctagnv8,"ctagnv8[8]/D");
  _tree->Branch("ejetnv8",&_d.ejetnv8,"ejetnv8[8]/D");
  _tree->Branch("pxjetnv8",&_d.pxjetnv8,"pxjetnv8[8]/D");
  _tree->Branch("pyjetnv8",&_d.pyjetnv8,"pyjetnv8[8]/D");
  _tree->Branch("pzjetnv8",&_d.pzjetnv8,"pzjetnv8[8]/D");
  _tree->Branch("ntrjetnv8",&_d.ntrjetnv8,"ntrjetnv8[8]/D");
  _tree->Branch("twovtxprobjetnv8",&_d.twovtxprobjetnv8,"twovtxprobjetnv8[8]/D");
  _tree->Branch("vtxanglenv8",&_d.vtxanglenv8,"vtxanglenv8[8]/D");

  _jets = 0;
  _jets4 = 0;
  _jets5 = 0;
  _jets7 = 0;
  _jets8 = 0;

  _jetsnv4 = 0;
  _jetsnv5 = 0;
  _jetsnv6 = 0;
  _jetsnv7 = 0;
  _jetsnv8 = 0;
}



void ZHHAlgo::process() {
  if (!_jets4) {
    Event::Instance()->Get(_jetname4.c_str(), _jets4);
  }
  if (!_jets5) {
    Event::Instance()->Get(_jetname5.c_str(), _jets5);
  }
  if (!_jets) {
    Event::Instance()->Get(_jetname.c_str(), _jets);
  }
  if (!_jets7) {
    Event::Instance()->Get(_jetname7.c_str(), _jets7);
  }
  if (!_jets8) {
    Event::Instance()->Get(_jetname8.c_str(), _jets8);
  }
  if (!_jetsnv4) {
    Event::Instance()->Get(_jetnamenv4.c_str(), _jetsnv4);
  }
  if (!_jetsnv5) {
    Event::Instance()->Get(_jetnamenv5.c_str(), _jetsnv5);
  }
  if (!_jetsnv6) {
    Event::Instance()->Get(_jetnamenv6.c_str(), _jetsnv6);
  }
  if (!_jetsnv7) {
    Event::Instance()->Get(_jetnamenv7.c_str(), _jetsnv7);
  }
  if (!_jetsnv8) {
    Event::Instance()->Get(_jetnamenv8.c_str(), _jetsnv8);
  }

  const Vertex* privtx = Event::Instance()->getPrimaryVertex();

  // check higgs decay & nbs
  const MCParticleVec& mcps = Event::Instance()->getMCParticles();

  _d.mcnb = 0;
  _d.mchdecaypdg[0] = _d.mchdecaypdg[1] = 0;
  _d.mchbb = 0;

  int hcount = 0;
  for (unsigned int i=0; i<mcps.size(); i++) {
    int abspdg = abs(mcps[i]->getPDG());
    int parpdg = 0;
    if (mcps[i]->getParent())parpdg = abs(mcps[i]->getParent()->getPDG());
    if (((abspdg > 500 && abspdg < 600) || (abspdg > 5000 && abspdg < 6000)) && parpdg < 100)
      _d.mcnb ++;

    if (mcps[i]->getPDG() == 25) {
      // higgs
      if (mcps[i]->getDaughters().size() != 2) {
        cout << "ERR: # of higgs daughters = " << mcps[i]->getDaughters().size() << endl;
        break;
      }
      if (hcount == 2) {
        cout << "Too many higgs found!, ignore decay" << endl;
        break;
      }
      int apdg = abs(mcps[i]->getDaughters()[0]->getPDG());
      _d.mchdecaypdg[hcount++] = apdg;
      if (apdg == 4)_d.mchbb ++;
    }
    if (mcps[i]->getPDG() == 5 && parpdg == 21 && mcps[i]->getParent()->getDaughters().size() == 2) {
      TVector3 v1 = mcps[i]->getParent()->getDaughters()[0]->Vect();
      TVector3 v2 = mcps[i]->getParent()->getDaughters()[1]->Vect();
      cout << "g->bb angle: " << v1.Angle(v2) << endl;
    }
  }

  // thrust
  vector<TVector3> v;
  const TrackVec& tracks = Event::Instance()->getTracks();
  const NeutralVec& neutrals = Event::Instance()->getNeutrals();

  for (unsigned int n=0; n<tracks.size(); n++) {
    v.push_back(tracks[n]->Vect());
  }
  for (unsigned int n=0; n<neutrals.size(); n++) {
    v.push_back(neutrals[n]->Vect());
  }

  TVector3 taxis;
  _d.thrust = algoEtc::calcThrust(v, taxis);
  _d.thaxis[0] = taxis.x();
  _d.thaxis[1] = taxis.y();
  _d.thaxis[2] = taxis.z();

  _d.ntr = tracks.size();
  _d.npfo = tracks.size() + neutrals.size();

  // sorting btag
  vector<const Jet*> jets;
  jets = *_jets;
  sort(jets.begin(), jets.end(), sortBtag);

  int nmass = 0;
  TLorentzVector totp;
  _d.ntrjetmin = 10000.;

  for (unsigned int nj = 0; nj < 6; nj ++) {
    Jet* j = const_cast<Jet*>(jets[nj]);  // TODO: bad boy...
    j->recalcFourMomentum();

    _d.btag[nj] = j->getParam("lcfiplus")->get<double>("BTag");
    _d.bcat[nj] = j->getParam("lcfiplus")->get<double>("Category");
    _d.ctag[nj] = j->getParam("lcfiplus")->get<double>("CTag");
    _d.ejet[nj] = j->E();
    _d.pxjet[nj] = j->Px();
    _d.pyjet[nj] = j->Py();
    _d.pzjet[nj] = j->Pz();

    _d.mcnb6[nj] = j->getParam("lcfiplus")->get<double>("MCnb");
    _d.mcnc6[nj] = j->getParam("lcfiplus")->get<double>("MCnc");

    totp += *j;

    unsigned int ntr = j->getAllTracks().size();
    _d.ntrjet[nj] = ntr;
    if (_d.ntrjetmin > ntr)_d.ntrjetmin = ntr;

    unsigned int nv = j->getVertices().size();
    _d.twovtxprobjet[nj] = 1.;
    _d.vtxangle[nj] = 0.;
    if (nv == 2) {
      // two vtx prob
      _d.twovtxprobjet[nj] *= j->getVertices()[0]->getProb();
      _d.twovtxprobjet[nj] *= j->getVertices()[1]->getProb();

      // vertex angle
      TVector3 pripos = privtx->getPos();
      TVector3 pos1 = j->getVertices()[0]->getPos() - pripos;
      TVector3 pos2 = j->getVertices()[1]->getPos() - pripos;

      _d.vtxangle[nj] = pos1.Angle(pos2);
    }

    // ycut values
    if (nj == 0) {
      for (int i=0; i<10; i++) {
        TString s;
        s.Form("y%d%d",i,i+1);
        _d.ycuts[i] = j->getParam("yth")->get<double>(s);
      }
    }

    // masses
    if (nj == 5)continue;
    for (unsigned int nj2 = nj + 1; nj2 < 6; nj2 ++) {
      const Jet* j2 = jets[nj2];
      TLorentzVector v = *j;
      v += *j2;
      _d.mass[nmass++] = v.M();
    }

  }
  _d.emiss = 500 - totp.E();
  _d.pmiss[0] = -totp.Px();
  _d.pmiss[1] = -totp.Py();
  _d.pmiss[2] = -totp.Pz();

  // sorting btag for
  vector<const Jet*> jets4;
  jets4 = *_jets4;
  sort(jets4.begin(), jets4.end(), sortBtag);

  for (unsigned int nj = 0; nj < 4; nj ++) {
    Jet* j = const_cast<Jet*>(jets4[nj]);  // TODO: bad boy...
    j->recalcFourMomentum();

    _d.btag4[nj] = j->getParam("lcfiplus")->get<double>("BTag");
    _d.bcat4[nj] = j->getParam("lcfiplus")->get<double>("Category");
    _d.ctag4[nj] = j->getParam("lcfiplus")->get<double>("CTag");
    _d.ejet4[nj] = j->E();
    _d.pxjet4[nj] = j->Px();
    _d.pyjet4[nj] = j->Py();
    _d.pzjet4[nj] = j->Pz();

    unsigned int ntr = j->getAllTracks().size();
    _d.ntrjet4[nj] = ntr;

    unsigned int nv = j->getVertices().size();
    _d.twovtxprobjet4[nj] = 1.;
    _d.vtxangle4[nj] = 0.;
    if (nv == 2) {
      // two vtx prob
      _d.twovtxprobjet4[nj] *= j->getVertices()[0]->getProb();
      _d.twovtxprobjet4[nj] *= j->getVertices()[1]->getProb();

      // vertex angle
      TVector3 pripos = privtx->getPos();
      TVector3 pos1 = j->getVertices()[0]->getPos() - pripos;
      TVector3 pos2 = j->getVertices()[1]->getPos() - pripos;

      _d.vtxangle4[nj] = pos1.Angle(pos2);
    }
  }

  // sorting btag for
  vector<const Jet*> jets5;
  jets5 = *_jets5;
  sort(jets5.begin(), jets5.end(), sortBtag);

  for (unsigned int nj = 0; nj < 5; nj ++) {
    Jet* j = const_cast<Jet*>(jets5[nj]);  // TODO: bad boy...
    j->recalcFourMomentum();

    _d.btag5[nj] = j->getParam("lcfiplus")->get<double>("BTag");
    _d.bcat5[nj] = j->getParam("lcfiplus")->get<double>("Category");
    _d.ctag5[nj] = j->getParam("lcfiplus")->get<double>("CTag");
    _d.ejet5[nj] = j->E();
    _d.pxjet5[nj] = j->Px();
    _d.pyjet5[nj] = j->Py();
    _d.pzjet5[nj] = j->Pz();

    unsigned int ntr = j->getAllTracks().size();
    _d.ntrjet5[nj] = ntr;

    unsigned int nv = j->getVertices().size();
    _d.twovtxprobjet5[nj] = 1.;
    _d.vtxangle5[nj] = 0.;
    if (nv == 2) {
      // two vtx prob
      _d.twovtxprobjet5[nj] *= j->getVertices()[0]->getProb();
      _d.twovtxprobjet5[nj] *= j->getVertices()[1]->getProb();

      // vertex angle
      TVector3 pripos = privtx->getPos();
      TVector3 pos1 = j->getVertices()[0]->getPos() - pripos;
      TVector3 pos2 = j->getVertices()[1]->getPos() - pripos;

      _d.vtxangle5[nj] = pos1.Angle(pos2);
    }
  }

  // sorting btag for
  vector<const Jet*> jets7;
  jets7 = *_jets7;
  sort(jets7.begin(), jets7.end(), sortBtag);

  for (unsigned int nj = 0; nj < 7; nj ++) {
    Jet* j = const_cast<Jet*>(jets7[nj]);  // TODO: bad boy...
    j->recalcFourMomentum();

    _d.btag7[nj] = j->getParam("lcfiplus")->get<double>("BTag");
    _d.bcat7[nj] = j->getParam("lcfiplus")->get<double>("Category");
    _d.ctag7[nj] = j->getParam("lcfiplus")->get<double>("CTag");
    _d.ejet7[nj] = j->E();
    _d.pxjet7[nj] = j->Px();
    _d.pyjet7[nj] = j->Py();
    _d.pzjet7[nj] = j->Pz();

    unsigned int ntr = j->getAllTracks().size();
    _d.ntrjet7[nj] = ntr;

    unsigned int nv = j->getVertices().size();
    _d.twovtxprobjet7[nj] = 1.;
    _d.vtxangle7[nj] = 0.;
    if (nv == 2) {
      // two vtx prob
      _d.twovtxprobjet7[nj] *= j->getVertices()[0]->getProb();
      _d.twovtxprobjet7[nj] *= j->getVertices()[1]->getProb();

      // vertex angle
      TVector3 pripos = privtx->getPos();
      TVector3 pos1 = j->getVertices()[0]->getPos() - pripos;
      TVector3 pos2 = j->getVertices()[1]->getPos() - pripos;

      _d.vtxangle7[nj] = pos1.Angle(pos2);
    }
  }

  // sorting btag for
  vector<const Jet*> jets8;
  jets8 = *_jets8;
  sort(jets8.begin(), jets8.end(), sortBtag);

  for (unsigned int nj = 0; nj < 8; nj ++) {
    Jet* j = const_cast<Jet*>(jets8[nj]);  // TODO: bad boy...
    j->recalcFourMomentum();

    _d.btag8[nj] = j->getParam("lcfiplus")->get<double>("BTag");
    _d.bcat8[nj] = j->getParam("lcfiplus")->get<double>("Category");
    _d.ctag8[nj] = j->getParam("lcfiplus")->get<double>("CTag");
    _d.ejet8[nj] = j->E();
    _d.pxjet8[nj] = j->Px();
    _d.pyjet8[nj] = j->Py();
    _d.pzjet8[nj] = j->Pz();

    unsigned int ntr = j->getAllTracks().size();
    _d.ntrjet8[nj] = ntr;

    unsigned int nv = j->getVertices().size();
    _d.twovtxprobjet8[nj] = 1.;
    _d.vtxangle8[nj] = 0.;
    if (nv == 2) {
      // two vtx prob
      _d.twovtxprobjet8[nj] *= j->getVertices()[0]->getProb();
      _d.twovtxprobjet8[nj] *= j->getVertices()[1]->getProb();

      // vertex angle
      TVector3 pripos = privtx->getPos();
      TVector3 pos1 = j->getVertices()[0]->getPos() - pripos;
      TVector3 pos2 = j->getVertices()[1]->getPos() - pripos;

      _d.vtxangle8[nj] = pos1.Angle(pos2);
    }
  }

  vector<const Jet*> jetsnv4;
  jetsnv4 = *_jetsnv4;
  sort(jetsnv4.begin(), jetsnv4.end(), sortBtag);
  for (unsigned int nj = 0; nj < 4; nj ++) {
    Jet* j = const_cast<Jet*>(jetsnv4[nj]);  // TODO: bad boy...
    j->recalcFourMomentum();

    _d.btagnv4[nj] = j->getParam("lcfiplus")->get<double>("BTag");
    _d.bcatnv4[nj] = j->getParam("lcfiplus")->get<double>("Category");
    _d.ctagnv4[nj] = j->getParam("lcfiplus")->get<double>("CTag");
    _d.ejetnv4[nj] = j->E();
    _d.pxjetnv4[nj] = j->Px();
    _d.pyjetnv4[nj] = j->Py();
    _d.pzjetnv4[nj] = j->Pz();

    unsigned int ntr = j->getAllTracks().size();
    _d.ntrjetnv4[nj] = ntr;

    unsigned int nv = j->getVertices().size();
    _d.twovtxprobjetnv4[nj] = 1.;
    _d.vtxanglenv4[nj] = 0.;
    if (nv == 2) {
      // two vtx prob
      _d.twovtxprobjetnv4[nj] *= j->getVertices()[0]->getProb();
      _d.twovtxprobjetnv4[nj] *= j->getVertices()[1]->getProb();

      // vertex angle
      TVector3 pripos = privtx->getPos();
      TVector3 pos1 = j->getVertices()[0]->getPos() - pripos;
      TVector3 pos2 = j->getVertices()[1]->getPos() - pripos;

      _d.vtxanglenv4[nj] = pos1.Angle(pos2);
    }
  }

  vector<const Jet*> jetsnv5;
  jetsnv5 = *_jetsnv5;
  sort(jetsnv5.begin(), jetsnv5.end(), sortBtag);
  for (unsigned int nj = 0; nj < 5; nj ++) {
    Jet* j = const_cast<Jet*>(jetsnv5[nj]);  // TODO: bad boy...
    j->recalcFourMomentum();

    _d.btagnv5[nj] = j->getParam("lcfiplus")->get<double>("BTag");
    _d.bcatnv5[nj] = j->getParam("lcfiplus")->get<double>("Category");
    _d.ctagnv5[nj] = j->getParam("lcfiplus")->get<double>("CTag");
    _d.ejetnv5[nj] = j->E();
    _d.pxjetnv5[nj] = j->Px();
    _d.pyjetnv5[nj] = j->Py();
    _d.pzjetnv5[nj] = j->Pz();

    unsigned int ntr = j->getAllTracks().size();
    _d.ntrjetnv5[nj] = ntr;

    unsigned int nv = j->getVertices().size();
    _d.twovtxprobjetnv5[nj] = 1.;
    _d.vtxanglenv5[nj] = 0.;
    if (nv == 2) {
      // two vtx prob
      _d.twovtxprobjetnv5[nj] *= j->getVertices()[0]->getProb();
      _d.twovtxprobjetnv5[nj] *= j->getVertices()[1]->getProb();

      // vertex angle
      TVector3 pripos = privtx->getPos();
      TVector3 pos1 = j->getVertices()[0]->getPos() - pripos;
      TVector3 pos2 = j->getVertices()[1]->getPos() - pripos;

      _d.vtxanglenv5[nj] = pos1.Angle(pos2);
    }
  }

  vector<const Jet*> jetsnv6;
  jetsnv6 = *_jetsnv6;
  sort(jetsnv6.begin(), jetsnv6.end(), sortBtag);
  for (unsigned int nj = 0; nj < 6; nj ++) {
    Jet* j = const_cast<Jet*>(jetsnv6[nj]);  // TODO: bad boy...
    j->recalcFourMomentum();

    _d.btagnv6[nj] = j->getParam("lcfiplus")->get<double>("BTag");
    _d.bcatnv6[nj] = j->getParam("lcfiplus")->get<double>("Category");
    _d.ctagnv6[nj] = j->getParam("lcfiplus")->get<double>("CTag");
    _d.ejetnv6[nj] = j->E();
    _d.pxjetnv6[nj] = j->Px();
    _d.pyjetnv6[nj] = j->Py();
    _d.pzjetnv6[nj] = j->Pz();

    unsigned int ntr = j->getAllTracks().size();
    _d.ntrjetnv6[nj] = ntr;

    unsigned int nv = j->getVertices().size();
    _d.twovtxprobjetnv6[nj] = 1.;
    _d.vtxanglenv6[nj] = 0.;
    if (nv == 2) {
      // two vtx prob
      _d.twovtxprobjetnv6[nj] *= j->getVertices()[0]->getProb();
      _d.twovtxprobjetnv6[nj] *= j->getVertices()[1]->getProb();

      // vertex angle
      TVector3 pripos = privtx->getPos();
      TVector3 pos1 = j->getVertices()[0]->getPos() - pripos;
      TVector3 pos2 = j->getVertices()[1]->getPos() - pripos;

      _d.vtxanglenv6[nj] = pos1.Angle(pos2);
    }
  }

  vector<const Jet*> jetsnv7;
  jetsnv7 = *_jetsnv7;
  sort(jetsnv7.begin(), jetsnv7.end(), sortBtag);
  for (unsigned int nj = 0; nj < 7; nj ++) {
    Jet* j = const_cast<Jet*>(jetsnv7[nj]);  // TODO: bad boy...
    j->recalcFourMomentum();

    _d.btagnv7[nj] = j->getParam("lcfiplus")->get<double>("BTag");
    _d.bcatnv7[nj] = j->getParam("lcfiplus")->get<double>("Category");
    _d.ctagnv7[nj] = j->getParam("lcfiplus")->get<double>("CTag");
    _d.ejetnv7[nj] = j->E();
    _d.pxjetnv7[nj] = j->Px();
    _d.pyjetnv7[nj] = j->Py();
    _d.pzjetnv7[nj] = j->Pz();

    unsigned int ntr = j->getAllTracks().size();
    _d.ntrjetnv7[nj] = ntr;

    unsigned int nv = j->getVertices().size();
    _d.twovtxprobjetnv7[nj] = 1.;
    _d.vtxanglenv7[nj] = 0.;
    if (nv == 2) {
      // two vtx prob
      _d.twovtxprobjetnv7[nj] *= j->getVertices()[0]->getProb();
      _d.twovtxprobjetnv7[nj] *= j->getVertices()[1]->getProb();

      // vertex angle
      TVector3 pripos = privtx->getPos();
      TVector3 pos1 = j->getVertices()[0]->getPos() - pripos;
      TVector3 pos2 = j->getVertices()[1]->getPos() - pripos;

      _d.vtxanglenv7[nj] = pos1.Angle(pos2);
    }
  }

  vector<const Jet*> jetsnv8;
  jetsnv8 = *_jetsnv8;
  sort(jetsnv8.begin(), jetsnv8.end(), sortBtag);
  for (unsigned int nj = 0; nj < 8; nj ++) {
    Jet* j = const_cast<Jet*>(jetsnv8[nj]);  // TODO: bad boy...
    j->recalcFourMomentum();

    _d.btagnv8[nj] = j->getParam("lcfiplus")->get<double>("BTag");
    _d.bcatnv8[nj] = j->getParam("lcfiplus")->get<double>("Category");
    _d.ctagnv8[nj] = j->getParam("lcfiplus")->get<double>("CTag");
    _d.ejetnv8[nj] = j->E();
    _d.pxjetnv8[nj] = j->Px();
    _d.pyjetnv8[nj] = j->Py();
    _d.pzjetnv8[nj] = j->Pz();

    unsigned int ntr = j->getAllTracks().size();
    _d.ntrjetnv8[nj] = ntr;

    unsigned int nv = j->getVertices().size();
    _d.twovtxprobjetnv8[nj] = 1.;
    _d.vtxanglenv8[nj] = 0.;
    if (nv == 2) {
      // two vtx prob
      _d.twovtxprobjetnv8[nj] *= j->getVertices()[0]->getProb();
      _d.twovtxprobjetnv8[nj] *= j->getVertices()[1]->getProb();

      // vertex angle
      TVector3 pripos = privtx->getPos();
      TVector3 pos1 = j->getVertices()[0]->getPos() - pripos;
      TVector3 pos2 = j->getVertices()[1]->getPos() - pripos;

      _d.vtxanglenv8[nj] = pos1.Angle(pos2);
    }
  }

  _tree->Fill();
}

void ZHHAlgo::end() {
  _file->Write();
  _file->Close();
}

void TestAlgo::init(Parameters* param) {
  Algorithm::init(param);

  string filename = param->get("FileName",string("test.root"));
  _jetname = param->get("JetCollectionName",string("Durham_2Jets"));
  string primvtxcolname = param->get("PrimaryVertexCollectionName",string("PrimaryVertex"));
  string secvtxcolname = param->get("SecondaryVertexCollectionName",string("BuildUpVertex"));
  Event::Instance()->setDefaultPrimaryVertex(primvtxcolname.c_str());
  Event::Instance()->setDefaultSecondaryVertices(secvtxcolname.c_str());

  _file = new TFile(filename.c_str(),"RECREATE");
//		_nt = new TNtupleD("nt","nt","nev:pdg:d0:d0sig:z0:z0sig:e:pt:pz:chi2:sd0:sz0:ecaldep:hcaldep");
  _nt = new TNtupleD("nt","nt","nev:pdg:parpdg:d0:d0sig:z0:z0sig:e:pt:pz:chi2:invtx:minchi2:nvtx");

  _jets = 0;
  _nev = 0;
}

void TestAlgo::process() {

  if (!_jets) {
    Event::Instance()->Get(_jetname.c_str(), _jets);
  }
//		TrackVec &tracks = Event::Instance()->getTracks();
  const Vertex* privtx = Event::Instance()->getPrimaryVertex();
//		VertexVec &vtcs = Event::Instance()->getSecondaryVertices();

  for (unsigned int nj = 0; nj < _jets->size(); nj++) {
    const Jet* j = (*_jets)[nj];
    TrackVec tracks = j->getAllTracks(true);
    VertexVec& vtcs = j->getVertices();

    int nvtx = 0;
    for (unsigned int nv=0; nv<vtcs.size(); nv++)
      if (vtcs[nv]->getTracks().size() >=2) nvtx ++;

    for (unsigned int n=0; n<tracks.size(); n++) {
      const Track* tr = tracks[n];

      //double sd0 = signedD0(tr, j, privtx, true);
      //double sd0sig = signedD0Significance(tr, j, privtx, true);
      //double sz0 = signedZ0(tr, j, privtx, true);
      //double sz0sig = signedZ0Significance(tr, j, privtx, true);

      // vertex-track association
      int invtx = 0;
      double minchi2 = 1e+300;

      if (find(privtx->getTracks().begin(), privtx->getTracks().end(), tr) != privtx->getTracks().end())invtx = 1;
      else {
        for (unsigned int nv = 0; nv < vtcs.size(); nv ++) {
          double chi2 = 1e+300;
          if (find(vtcs[nv]->getTracks().begin(), vtcs[nv]->getTracks().end(), tr) != vtcs[nv]->getTracks().end()) {
            invtx = 2;
            chi2 = vtcs[nv]->getChi2Track(tr);
          } else { // if(vtcs[nv]->getTracks().size() >= 2){
            chi2 = vtcs[nv]->getChi2TrackFit(tr,3);
          }

          if (minchi2 > chi2)minchi2 = chi2;
        }
      }

      const MCParticle* mcp = tracks[n]->getMcp();
      const MCParticle* pmcp = (mcp ? mcp->getSemiStableBParent() : 0);
      _nt->Fill(_nev, mcp ? mcp->getPDG() : 0, pmcp ? pmcp->getPDG() : 0, fabs(tr->getD0()), fabs(tr->getD0() / sqrt(tr->getCovMatrix()[tpar::d0d0])),
                fabs(tr->getZ0()), fabs(tr->getZ0() / sqrt(tr->getCovMatrix()[tpar::z0z0])),
                tr->E(), tr->Pt(), tr->Pz(), tr->getChi2(), invtx, minchi2, nvtx * 10 + vtcs.size());
      //sd0, sz0, tr->getCaloEdep()[tpar::ecal], tr->getCaloEdep()[tpar::hcal]);

    }
  }

  _nev ++;
}

void TestAlgo::end() {
  _file->Write();
  _file->Close();
}

void FlavtagReader::init(Parameters* param) {
  Algorithm::init(param);

  string filename = param->get("FileName",string("test.root"));
  _jetname = param->get("JetCollectionName",string("Durham_2Jets"));
  string primvtxcolname = param->get("PrimaryVertexCollectionName",string("PrimaryVertex"));
  Event::Instance()->setDefaultPrimaryVertex(primvtxcolname.c_str());

  _file = new TFile(filename.c_str(),"RECREATE");
  _nt = new TNtupleD("nt","nt","nev:nj:e:px:py:pz:btag:ctag:otag:bbtag:bctag:cctag");
  _ntev = new TNtupleD("ntev","ntev","nev:btag1:btag2:btag3:btag4:btag5:btag6:ctag1:ctag2:ctag3:ctag4:ctag5:ctag6");

  _jets = 0;
  _nev = 0;
}

void FlavtagReader::process() {
  if (!_jets) {
    Event::Instance()->Get(_jetname.c_str(), _jets);
  }
  //const Vertex * privtx = Event::Instance()->getPrimaryVertex();

  vector<double> btags, ctags;

  for (unsigned int nj = 0; nj < _jets->size(); nj++) {
    const Jet* j = (*_jets)[nj];

    const Parameters* para = j->getParam("lcfiplus");
    _nt->Fill(_nev, nj, j->E(), j->Px(), j->Py(), j->Pz(),
              para->get<double>("BTag"), para->get<double>("CTag"), para->get<double>("OTag"),  para->get<double>("BBTag"),  para->get<double>("CCTag"),  para->get<double>("BCTag"));
    btags.push_back(para->get<double>("BTag"));
    ctags.push_back(para->get<double>("CTag"));

    cout << "nvtx = " << para->get<double>("nvtx") << ", nvtxall = " << para->get<double>("nvtxall") << endl;
  }

  std::sort(btags.begin(), btags.end());
  std::sort(ctags.begin(), ctags.end());
  if (_jets->size() >= 6)
    _ntev->Fill(_nev,  btags[0],  btags[1],  btags[2],  btags[3],  btags[4],  btags[5],  ctags[0],  ctags[1],  ctags[2],  ctags[3],  ctags[4],  ctags[5]);

  _nev ++;
}

void FlavtagReader::end() {
  _file->Write();
  _file->Close();
}

#if 0

void TestAlgo::init(Parameters* param) {
  Algorithm::init(param);

  gStyle->SetPalette(1);
  _h = new TH2D("h","h",200,-2,2,200,-2,2);
  _he = new TH2D("he","he",200,-2,2,200,-2,2);
}

void TestAlgo::process() {
  // check bbbbbb (reject H->WW etc.)
  const MCParticleVec& mcps = Event::Instance()->getMCParticles();
  const MCColorSingletVec& mccss = Event::Instance()->getMCColorSinglets();
  cout << "# mccs = " << mccss.size() << endl;
  if (mccss.size() != 3) return;

  int nq[3];
  for (unsigned int i=0; i<3; i++) {
    nq[i] = mccss[i]->_initials.size();
  }
  cout << "# qs = " << nq[0] << " " << nq[1] << " " << nq[2] << endl;
  if (nq[0]!=2 ||nq[1]!=2 ||nq[2]!=2)return;

  for (unsigned int i=0; i<mcps.size(); i++) {
    const MCColorSinglet* mccs = mcps[i]->getColorSinglet(&mccss);
    const MCParticle* p1 = mccs->_initials[0];
    const MCParticle* p2 = mccs->_initials[1];

    TVector3 normal = p1->Vect().Cross(p2->Vect()).Unit();
    double ndp = (normal.Dot(mcps[i]->Vect().Unit()));
    TVector3 pplane = mcps[i]->Vect().Unit() - (normal * ndp);
    double nxp1 = p1->Vect().Unit().Dot(pplane);
    double nxp2 = p2->Vect().Unit().Dot(pplane);
    double nxp12 = fabs(p1->Vect().Unit().Dot(p2->Vect().Unit()));
    double nxp = (nxp1 - nxp2) / (1 - nxp12);

    _h->Fill(nxp, ndp);
    cout << nxp << " " << ndp << endl;
  }
}

void TestAlgo::end() {
  _h->Draw("colz");
  gPad->Update();
  gSystem->ProcessEvents();
}


void TestAlgo::init(Parameters* param) {
  Algorithm::init(param);

  string filename = param->get("FileName",string("test.root"));
  _privtxname = param->get("PrimaryVertexCollectionName",string("PrimaryVertex"));
  _vtxname = param->get("BuildUpVertexCollectionName",string("BuildUpVertex"));
  _v0vtxname = param->get("V0VertexCollectionName",string("BuildUpVertex_V0"));
  _jetname = param->get("JetCollectionName",string("Durham_2Jets"));
  _vtxsel = param->get("VertexSelection",(int)0);
  _refine = param->get("PerformRefining",(int)0);
  _bbhh = param->get("IsBBHH",int(0));

  _file = new TFile(filename.c_str(),"RECREATE");
  _ntJet = new TNtupleD("ntJet","ntJet","nvtx:1vtxprob:2vtxprob:cflt:ecvtx:ejet:vangle:vmass:esingle");

  _vertices = 0;
  _jets = 0;
}



void TestAlgo::process() {
  if (!_vertices) {
    Event::Instance()->Get(_vtxname.c_str(), _vertices);
    //cout << "vtx name: " << _vtxname << ", pointer = " << (unsigned int)_vertices << endl;
  }
  if (!_v0vertices) {
    Event::Instance()->Get(_v0vtxname.c_str(), _v0vertices);
    //cout << "vtx name: " << _vtxname << ", pointer = " << (unsigned int)_vertices << endl;
  }

  if (!_jets) {
    Event::Instance()->Get(_jetname.c_str(), _jets);
    //cout << "jet name: " << _jetname << ", pointer = " << (unsigned int)_jets << endl;
  }

  // check bbbbbb (reject H->WW etc.)
  const MCParticleVec& mcps = Event::Instance()->getMCParticles();
  if (_bbhh) {
    int hcount = 0;
    for (unsigned int i=0; i<mcps.size(); i++) {
      if (mcps[i]->getPDG() != 25)continue;

      // higgs
      if (mcps[i]->getDaughters().size() != 2) {
        cout << "ERR: # of higgs daughters = " << mcps[i]->getDaughters().size() << endl;
        break;
      }
      if (abs(mcps[i]->getDaughters()[0]->getPDG()) != 5)break;
      if (abs(mcps[i]->getDaughters()[1]->getPDG()) != 5)break;
      hcount ++;
    }
    if (hcount < 2)return;
  }

  // select vertices
  vector<const Track*> residualTracks;
  vector<const Vertex*> selectedVertices;
  const vector<const Vertex*>* pVertices;
  if (_vertices && _vtxsel) {
    VertexSelectorConfig vscfg;
    vscfg.rejectdist = true;
    vscfg.minpos = .3;
    vscfg.maxpos = 30.;
    vscfg.rejectk0 = true;
    vscfg.k0width = .01;

    selectedVertices = VertexSelector()(*_vertices, vscfg, residualTracks,false);
    pVertices = &selectedVertices;
  } else {
    pVertices = _vertices;
  }

  cout << "# jet = " << _jets->size() << endl;

  // copy vertices
  vector<Vertex*> vtcs2;
  const Vertex* ip = Event::Instance()->getPrimaryVertex(_privtxname.c_str());
  if (ip == 0)throw(Exception("IP not found!"));

  for (unsigned int v=0; v<pVertices->size(); v++) {
    if (!(*pVertices)[v]->isPrimary())
      vtcs2.push_back(new Vertex(*(*pVertices)[v]));
  }

  cout << "# sec vtx = " << vtcs2.size() << endl;

  vector<vector<Vertex*> > jetVertices;
  vector<vector<const Track*> > jetResidualTracks;
  // jet-vtx association
  lcfiplus::algoEtc::connectVerticesToJets(*_jets, vtcs2, jetVertices, jetResidualTracks,ip);

  VertexFinderSuehara::VertexFinderSueharaConfig cfg;

  for (unsigned int j=0; j<_jets->size(); j++) {

    // single track probability
    double singleprob = 0;
    double twoprob = 0;
    double cflt = 0;
    double ecvtx = 0;
    double vangle = 0;
    double vmass = 0;
    double esingle = 0;

    if (_refine) {
      vector<Vertex*> singleVtcs = VertexFinderSuehara::makeSingleTrackVertices(constVector(jetVertices[j]), jetResidualTracks[j], *_v0vertices, ip, cfg);

      if (jetVertices[j].size() + singleVtcs.size() >= 2) {
        cout << "Before recombination:" << endl;
        for (unsigned int k=0; k<jetVertices[j].size(); k++)
          jetVertices[j][k]->Print();
        for (unsigned int k=0; k<singleVtcs.size(); k++)
          singleVtcs[k]->Print();
      }

      VertexFinderSuehara::recombineVertices(jetVertices[j], singleVtcs);

      // v0 selection again
      VertexSelector()(jetVertices[j], cfg.v0selVertex);

      vector<const Track*> singletracklist;
      if (jetVertices[j].size() > 1) {
        twoprob = jetVertices[j][0]->getProb() * jetVertices[j][1]->getProb();
        singletracklist.resize(jetVertices[j][0]->getTracks().size() + jetVertices[j][1]->getTracks().size());
        std::copy(jetVertices[j][0]->getTracks().begin(), jetVertices[j][0]->getTracks().end(), singletracklist.begin());
        std::copy(jetVertices[j][1]->getTracks().begin(), jetVertices[j][1]->getTracks().end(), singletracklist.begin() + jetVertices[j][0]->getTracks().size());

        Vertex* single = VertexFitterSimple_V()(singletracklist.begin(), singletracklist.end());

        singleprob = single->getProb();

        cout << "twoprob: " << jetVertices[j][0]->getProb() << " " << jetVertices[j][1]->getProb() << " oneprob: " << singleprob << endl;

        delete single;

        cflt = (jetVertices[j][1]->getPos() - jetVertices[j][0]->getPos()).Mag();

        // looking for near vertex
        int nnear = (jetVertices[j][0]->getPos().Mag() > jetVertices[j][1]->getPos().Mag() ? 0 : 1);

        for (unsigned int ntr=0; ntr<jetVertices[j][nnear]->getTracks().size(); ntr++) {
          const Track* tr = jetVertices[j][nnear]->getTracks()[ntr];
          ecvtx += tr->E();
        }
        cout << "cflt = " << cflt << ", ecvtx = " << ecvtx << ", ejet = " << (*_jets)[j]->E() << endl;

        // single track investigation
        int idx = -1;
        if (jetVertices[j][0]->getTracks().size() == 1) idx = 0;
        if (jetVertices[j][1]->getTracks().size() == 1) idx = 1;

        if (idx >= 0) {
          const Track* tr = jetVertices[j][idx]->getTracks()[0];
          TVector3 vpos = jetVertices[j][idx]->getPos();
          vangle = vpos.Angle(tr->Vect());
          vmass = 2 * tr->E() * tr->E() * (1 - cos(vpos.Angle(tr->Vect())));
          esingle = tr->E();
        }

      } else if (jetVertices[j].size() == 1) {
        singleprob = jetVertices[j][0]->getProb();
      }
    }

    if (jetVertices[j].size() >= 2) {
      for (unsigned int k=0; k<jetVertices[j].size(); k++)
        jetVertices[j][k]->Print();
    }

    _file->cd();
    _ntJet->Fill((int)jetVertices[j].size(),singleprob, twoprob, cflt, ecvtx, (*_jets)[j]->E(), vangle, vmass, esingle);
  }
}

void TestAlgo::init(Parameters* param) {
  Algorithm::init(param);


  string filename = param->get("FileName",string("test.root"));
  _jetname = param->get("JetCollectionName",string("Durham_6Jets2"));
  _bbhh = param->get("IsBBHH",int(0));

  _file = new TFile(filename.c_str(),"RECREATE");

  _ntJet2 = new TNtuple("ntJet2","ntJet2","nev:njet:nbjetmc:nvtx:nvtxjet:ngoodvtx:fracgoodvtxtrack:ycut:nbjet:fracgoodtrack");
  _nbJet = new TNtuple("nbJet", "number of b tracks in eachjet", "nev:nb1:nb2:nb3:nb4:nb5:nb6:nb11:nb12:nb13:nb14:nb15:nb16");

  _jets = 0;
}

void TestAlgo::process() {
  if (!_jets) {
    Event::Instance()->Get(_jetname.c_str(), _jets);
    //cout << "jet name: " << _jetname << ", pointer = " << (unsigned int)_jets << endl;
  }
  JetVec& jets = *_jets;
  unsigned int nj = jets.size();

  Event* event = Event::Instance();
  MCParticleVec& mcps = event->getMCParticles();

  // check bbbbbb (reject H->WW etc.)
  if (_bbhh) {
    int hcount = 0;
    for (unsigned int i=0; i<mcps.size(); i++) {
      if (mcps[i]->getPDG() != 25)continue;

      // higgs
      if (mcps[i]->getDaughters().size() != 2) {
        cout << "ERR: # of higgs daughters = " << mcps[i]->getDaughters().size() << endl;
        break;
      }
      if (abs(mcps[i]->getDaughters()[0]->getPDG()) != 5)break;
      if (abs(mcps[i]->getDaughters()[1]->getPDG()) != 5)break;
      hcount ++;
    }
    if (hcount < 2)return;
  }

  // semistable B
  vector<const MCParticle*> blist;
  for (unsigned int i=0; i<mcps.size(); i++) {
    if (mcps[i]->isSemiStableB()) {
      blist.push_back(mcps[i]);
    }
  }
  cout << "Number of semistable B: " << blist.size() << endl;

//	TNtuple *ntResidual = new TNtuple("ntResidual","ResidualTracks","nev:bid:btracks:mcvx:mcvy:mcvz:d0:d0err:z0:z0err:tre");
  // calculate btracks
  /*		int *btracks = new int[blist.size()];
  		memset(btracks,0,sizeof(int)*blist.size());
  		for(unsigned int i=0;i<tracks.size();i++){
  			for(unsigned int k=0;k<blist.size();k++){
  				if(tracks[i]->getMcp()->isParent(blist[k]))btracks[k] ++;
  			}
  		}
  */
  vector<const Track*> assignedTracks;
  vector<const Track*> residualTracks;

  map<const Jet*, int > nbmap;
  int nvtxjet = 0;
  for (unsigned int nj=0; nj<jets.size(); nj++) {
    nbmap[jets[nj]] = 0;
    if (jets[nj]->getVertices().size()>0)nvtxjet ++;
  }

  // nbjet fill
  vector<int> nbjet0(max(int(jets.size()),6));
  vector<int> nbjet1(max(int(jets.size()),6));
  for (unsigned int nj=0; nj<jets.size(); nj++) {
    nbjet0[nj] = nbjet1[nj] = 0;

    for (unsigned int i=0; i<jets[nj]->getTracks().size(); i++) {
      const Track* tr = jets[nj]->getTracks()[i];
      if (tr->getMcp() == 0)continue;
      if (tr->getMcp()->getSemiStableParent() == 0)continue;
      int pdg = tr->getMcp()->getSemiStableParent()->getPDG();
      if ((abs(pdg)>400&&abs(pdg)<600) || (abs(pdg)>4000&&abs(pdg)<6000)) {
        for (unsigned int k=0; k<blist.size(); k++) {
          if (tr->getMcp()->isParent(blist[k])) {
            nbjet0[nj] ++;
            if (tr->E()>1.)nbjet1[nj] ++;
            break;
          }
        }
      }
    }
    for (unsigned int nv=0; nv<jets[nj]->getVertices().size(); nv++) {
      for (unsigned int i=0; i<jets[nj]->getVertices()[nv]->getTracks().size(); i++) {
        const Track* tr = jets[nj]->getVertices()[nv]->getTracks()[i];
        if (tr->getMcp() == 0)continue;
        if (tr->getMcp()->getSemiStableParent() == 0)continue;
        int pdg = tr->getMcp()->getSemiStableParent()->getPDG();
        if ((abs(pdg)>400&&abs(pdg)<600) || (abs(pdg)>4000&&abs(pdg)<6000)) {
          for (unsigned int k=0; k<blist.size(); k++) {
            if (tr->getMcp()->isParent(blist[k])) {
              nbjet0[nj] ++;
              if (tr->E()>1.)nbjet1[nj] ++;
              break;
            }
          }
        }
      }
    }
  }

  // sort by decending order
  sort(nbjet0.begin(), nbjet0.end(), greater<int>());
  sort(nbjet1.begin(), nbjet1.end(), greater<int>());
  _nbJet->Fill(0,nbjet0[0], nbjet0[1], nbjet0[2], nbjet0[3], nbjet0[4], nbjet0[5],
               nbjet1[0], nbjet1[1], nbjet1[2], nbjet1[3], nbjet0[4], nbjet0[5]);

  for (unsigned int ib=0; ib<blist.size(); ib++) {
    vector<const Track*> aTracks, rTracks;
//Jet * JetMCMatch(vector<Jet *> &jets, MCParticle *mcp, vector<Track *> &assignedTracks, vector<Track *> &residualTracks)
    const Jet* jet = JetMCMatch(jets, blist[ib], aTracks, rTracks);
    if (jet) {
      nbmap[jet] ++;
      assignedTracks.insert(assignedTracks.end(), aTracks.begin(), aTracks.end());
      residualTracks.insert(residualTracks.end(), rTracks.begin(), rTracks.end());
    }
  }

  int nbjet = 0;
  for (unsigned int nj=0; nj<jets.size(); nj++) {
    if (nbmap[jets[nj]] > 0)nbjet ++;
  }
//	TNtuple *ntJet2 = new TNtuple("ntJet2","ntJet2","nev:njet:nbjetmc:nvtx:nvtxjet:ngoodvtx:ngoodvtxtrack:ycut:nbjet:fracgoodtrack");
//	TNtuple *ntJet2 = new TNtuple("ntJet2","ntJet2","nev:njet:ycut:nbjet:fracgoodtrack");
  double fracgoodtrack = double(residualTracks.size()) / double(assignedTracks.size()+residualTracks.size());
  _ntJet2->Fill(0, jets.size(), blist.size(), 0, nvtxjet, 0, 0, 0., nbjet, fracgoodtrack);
}


void TestAlgo::end() {
  _file->Write();
  _file->Close();
}
#endif

void VertexAnalysis::init(Parameters* param) {
  Algorithm::init(param);
  _privtxname = param->get("PrimaryVertexCollectionName",string("PrimaryVertex"));

  string filename = param->get("VertexAnalysis.FileName",string("VertexAnalysis.root"));
  _secvtxname = param->get("VertexAnalysis.SecondaryVertexCollectionName",string("RefinedVertex"));
  _jetname = param->get("VertexAnalysis.JetCollectionName",string("RefinedJets"));
  _file = new TFile(filename.c_str(),"RECREATE");
  _nt = new TNtupleD("vtxtree","vtxtree","track:invtx:d0:d0err:z0:z0err:e:pt:chi2:ndf:vtxftdhits:parent:bhadron:bvtxdist:blepton:chadron:cvtxdist:clepton:recvtxdist:recvtxdist2");
}

void VertexAnalysis::process() {
  const Vertex* privtx = Event::Instance()->getPrimaryVertex(_privtxname.c_str());
  TrackVec& tracks = Event::Instance()->getTracks();
  MCParticleVec& mcps = Event::Instance()->getMCParticles();
  VertexVec* psecvtx = 0;
  if (_secvtxname != "")
    Event::Instance()->Get(_secvtxname.c_str(), psecvtx);

  JetVec *pjet = 0;
  if (_jetname != "")
    Event::Instance()->Get(_jetname.c_str(), pjet);

  if (mcps.size() == 0) {
    cout << "MCParticle collection not specified. We need the MCParticle collection for vertex analysis." << endl;
    return;
  }


  if(pjet){
    for(unsigned int j=0;j<pjet->size();j++){
      const Jet *jet = (*pjet)[j];
      VertexVec &v = jet->getVertices();
      cout << "Jet " << j << " has " << v.size() << " vertices. Direction (" << jet->Px() << ", " << jet->Py() << ", " << jet->Pz() << ")" << endl;
      
      for(unsigned int k=0;k<v.size();k++){
	TrackVec &vtrs = v[k]->getTracks();
	cout << "Vertex " << k << " has " << vtrs.size() << " tracks. position (" << v[k]->getX() << ", " << v[k]->getY() << ", " << v[k]->getZ() << ")" << endl;
	
	for(unsigned int l=0;l<vtrs.size();l++){
	  const Track *vtr = vtrs[l];
	  const MCParticle *vmcp = vtr->getMcp();
	  if(!vmcp){
	    cout << "Track " << l << " cannot be associated with MCP. Energy = " << vtr->E() << endl;
	  }else{
	    const MCParticle *parb, *parc, *par;
	    parb = vmcp->getSemiStableBParent();
	    parc = vmcp->getSemiStableCParent();
	    par = vmcp->getSemiStableParent();
	    cout << "Track " << l << ": PDG " << vmcp->getPDG() << ", B parent " << (parb ? parb->getPDG() : 0)
		 << ", C parent " << (parc ? parc->getPDG() : 0) << ", parent " << (par ? par->getPDG() : 0) << endl;
	    cout << "Start point (" << vmcp->getVertex().x() << ", " << vmcp->getVertex().y() << ", " << vmcp->getVertex().z() << ")" << endl;
	  } 
	}
      }
    }
  }
  

  for (unsigned int i=0; i < tracks.size(); ++i) {
    const Track* tr = tracks[i];
    const MCParticle* mcp = tr->getMcp();

    struct {
      double trackseed;
      double invtx;
      double d0;
      double d0err;
      double z0;
      double z0err;
      double e;
      double pt;
      double chi2;
      double ndf;
      double vtxftdhits;
      double parent;
      double bhadron;
      double bvtxdist;
      double blepton;
      double chadron;
      double cvtxdist;
      double clepton;
      double recvtxdist;
      double recvtxdist2;
    } d;

    memset(&d,0,sizeof(d));

    if(pjet){
      bool found = false;
      for(unsigned int j=0;j<pjet->size();j++){
        TrackVec& vtr = (*pjet)[j]->getAllTracks();
        if (find(vtr.begin(), vtr.end(), tr) != vtr.end()){
	  const Jet *jet = (*pjet)[j];
	  TVector3 pospri = privtx->getPos();
	  if(jet->getVertices().size() > 0){
	    TVector3 posvtx = jet->getVertices()[0]->getPos();
	    d.recvtxdist = (posvtx - pospri).Mag();
	    //cout << "privtx " << pospri.Mag() << " sec " << posvtx.Mag() << " dist " << d.recvtxdist << endl;
	  }
	  if(jet->getVertices().size() > 1){
	    TVector3 posvtx = jet->getVertices()[1]->getPos();
	    d.recvtxdist2 = (posvtx - pospri).Mag();
	  }
	  found = true;
	  break;
	}
      }
      if(!found)continue;
    }


    const MCParticle *mcd;
    //    if (mcp->getSemiStableParent() == 0) d.trackseed = 0.;
    if (mcp->getSemiStableParent() == 0 || (mcp->getVertex().x() == 0 && mcp->getVertex().y() == 0)) d.trackseed = 0.;
    else{
      d.parent = mcp->getSemiStableParent()->getPDG();
      d.trackseed = 3.;
      if ((mcd = mcp->getSemiStableBParent()) != 0){
	d.trackseed = 1.;
	d.bhadron = mcd->getPDG();
	d.bvtxdist = mcd->decayDistance();
	const MCParticle *mcdl = mcd->semileptonicDecay();
	d.blepton = (mcdl ? mcdl->getPDG() : 0);
      }
      if ((mcd = mcp->getSemiStableCParent()) != 0){
	d.trackseed = 2.;
	d.chadron = mcd->getPDG();
	d.cvtxdist = mcd->decayDistance();
	const MCParticle *mcdl = mcd->semileptonicDecay();
	d.clepton = (mcdl ? mcdl->getPDG() : 0);
      }
    }

    d.invtx = (find(privtx->getTracks().begin(), privtx->getTracks().end(), tr) != privtx->getTracks().end());
    if (d.invtx == 0. && psecvtx) {
      // looking for secondary vertices
      for (unsigned int j=0; j<psecvtx->size(); j++) {
        TrackVec& vtr = (*psecvtx)[j]->getTracks();
        if (find(vtr.begin(), vtr.end(), tr) != vtr.end())
          d.invtx = 2.;
      }
    }

    d.d0 = tr->getD0();
    d.d0err = sqrt(tr->getCovMatrix()[tpar::d0d0]);
    d.z0 = tr->getZ0();
    d.z0err = sqrt(tr->getCovMatrix()[tpar::z0z0]);

    d.e = tr->E();
    d.pt = tr->Pt();

    d.chi2 = tr->getChi2();
    d.ndf = tr->getNdf();
    d.vtxftdhits = tr->getVtxHits() + tr->getFtdHits();

    _nt->Fill((double *)&d);

  }
}

void VertexAnalysis::end() {
  _file->Write();
  _file->Close();
}


}

