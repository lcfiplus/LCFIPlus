#include <string>

#include "TFile.h"
#include "TNtuple.h"
#include "TNtupleD.h"
#include "TSystem.h"
#include "TPad.h"
#include "TStyle.h"

#include "lcfiplus.h"
#include "process.h"
#include "DNNProvider2.h"
#include "VertexSelector.h"
#include "algoEtc.h"
#include "VertexFinderSuehara.h"
#include "VertexFitterSimple.h"

#include <utility>
#include <algorithm>

using namespace lcfiplus;

namespace lcfiplus {

void DNNProvider2::init(Parameters* param) {
  Algorithm::init(param);
  string filename = param->get("DNNProvider2.FileName",string("DNNProvider2.root"));
  _jetname = param->get("DNNProvider2.JetCollectionName",string("RefinedJets"));
  string privtx = param->get("PrimaryVertexCollectionName",string("PrimaryVertex"));
  Event::Instance()->setDefaultPrimaryVertex(privtx.c_str());

  _mcIsB = param->get("DNNProvider2.MC_IsB",int(0));
  _mcIsC = param->get("DNNProvider2.MC_IsC",int(0));
  _mcIsQ = param->get("DNNProvider2.MC_IsQ",int(0));

  _jets = 0;

  if(_mcIsB + _mcIsC + _mcIsQ != 1){
    cout << "Label parameter is wrong! Put 1 on one of MC_Is(BCQ)." << endl;
  }

  cout << filename << endl;

  _file = new TFile(filename.c_str(),"RECREATE");
  _ntp = new TTree("tree","tree");

  DNNData & d = _data;
  // jet variables
  _ntp->Branch("jet_px",&d.jet_px,"jet_px/F");
  _ntp->Branch("jet_py",&d.jet_py,"jet_py/F");
  _ntp->Branch("jet_pz",&d.jet_pz,"jet_pz/F");
  _ntp->Branch("jet_e",&d.jet_e,"jet_e/F");
  _ntp->Branch("jet_mass",&d.jet_mass,"jet_mass/F");
  _ntp->Branch("jet_ntracks",&d.jet_ntracks,"jet_ntracks/I");
  _ntp->Branch("jet_nneutrals",&d.jet_nneutrals,"jet_nneutrals/I");

  // particle kinematics
  _ntp->Branch("px",&d.px);
  _ntp->Branch("py",&d.py);
  _ntp->Branch("pz",&d.pz);
  _ntp->Branch("e",&d.e);
  _ntp->Branch("efrac",&d.efrac);
  _ntp->Branch("dtheta",&d.dtheta);
  _ntp->Branch("dphi",&d.dphi);

  // particle displacements
  _ntp->Branch("d0",&d.d0);
  _ntp->Branch("d0sig",&d.d0sig);
  _ntp->Branch("z0",&d.z0);
  _ntp->Branch("z0sig",&d.z0sig);
  _ntp->Branch("ip3d",&d.ip3d);
  _ntp->Branch("ip3dsig",&d.ip3dsig);

  // particle ID
  _ntp->Branch("charge",&d.charge);
  _ntp->Branch("isMuon",&d.ismuon);
  _ntp->Branch("isElectron",&d.iselectron);
  _ntp->Branch("isPhoton",&d.isphoton);
  _ntp->Branch("isChargedHadron",&d.ischargedhadron);
  _ntp->Branch("isNeutralHadron",&d.isneutralhadron);

  // label
  _ntp->Branch("mc_b",&d.mc_b,"mc_b/I");
  _ntp->Branch("mc_c",&d.mc_c,"mc_c/I");
  _ntp->Branch("mc_q",&d.mc_q,"mc_q/I");
}

void DNNProvider2::process() {
  if (!_jets) {
    Event::Instance()->Get(_jetname.c_str(), _jets);
  }

  const JetVec& jets = *_jets;
  const Vertex *privtx = Event::Instance()->getPrimaryVertex();

  DNNData &d = _data;

  for (unsigned int j=0; j < jets.size(); ++j) {
    const Jet* jet = jets[j];

    memset(&_data,0,sizeof(_data));

    d.jet_px = jet->Px();
    d.jet_py = jet->Py();
    d.jet_pz = jet->Pz();
    d.jet_e = jet->E();
    d.jet_mass = jet->M();
    TrackVec &tracks = jet->getAllTracks();
    NeutralVec &neutrals = jet->getNeutrals();

    d.jet_ntracks = tracks.size();
    d.jet_nneutrals = neutrals.size();

    float jet_theta = jet->Theta();
    float jet_phi = jet->Phi();

    d.mc_b = _mcIsB;
    d.mc_c = _mcIsC;
    d.mc_q = _mcIsQ;

    // probably order of tracks/netural does not matter...
    int nall = d.jet_ntracks + d.jet_nneutrals;
    d.px.resize(nall);
    d.py.resize(nall);
    d.pz.resize(nall);
    d.e.resize(nall);
    d.efrac.resize(nall);
    d.dtheta.resize(nall);
    d.dphi.resize(nall);

    d.d0.resize(nall);
    d.d0sig.resize(nall);
    d.z0.resize(nall);
    d.z0sig.resize(nall);
    d.ip3d.resize(nall);
    d.ip3dsig.resize(nall);

    d.charge.resize(nall);
    d.ismuon.resize(nall);
    d.iselectron.resize(nall);
    d.isphoton.resize(nall);
    d.ischargedhadron.resize(nall);
    d.isneutralhadron.resize(nall);

    vector<std::pair<float, int> > order;
    order.resize(nall);

    int i;

    for(i=0;i<d.jet_ntracks;i++){
      const Track *tr = tracks[i];
      order[i] = std::pair<float, int>(tr->E(), i);
    }
    for(;i<nall;i++){
      const Neutral *neu = neutrals[i-d.jet_ntracks];
      order[i] = std::pair<float, int>(neu->E(), i);
    }

    std::sort(order.begin(), order.end(), [](std::pair<float, int>a, std::pair<float, int> b){
	return a.first > b.first;
      });

    for(i=0;i<nall;i++){
      if(order[i].second >= d.jet_ntracks) continue;
      //cout << i << " " << order[i].second << " " << d.jet_ntracks << " " << nall << endl;
      const Track *tr = tracks[order[i].second];
      d.px[i] = tr->Px();
      d.py[i] = tr->Py();
      d.pz[i] = tr->Pz();
      d.e[i] = tr->E();
      d.efrac[i] = tr->E() / jet->E();
      d.dtheta[i] = tr->Theta() - jet_theta;
      d.dphi[i] = tr->Phi() - jet_phi;

      d.d0[i] = tr->getD0();
      d.d0sig[i] = tr->getD0() / sqrt(tr->getCovMatrix()[tpar::cov::d0d0]);
      d.z0[i] = tr->getZ0();
      d.z0sig[i] = tr->getZ0() / sqrt(tr->getCovMatrix()[tpar::cov::z0z0]);

      d.ip3d[i] = sqrt(tr->getD0() * tr->getD0() + tr->getZ0() * tr->getZ0());
      d.ip3dsig[i] = d.ip3d[i] / sqrt(tr->getCovMatrix()[tpar::cov::d0d0] + tr->getCovMatrix()[tpar::cov::z0z0] + 2 * tr->getCovMatrix()[tpar::cov::d0z0]);

      d.charge[i] = tr->getCharge();
      // tracing LCFIPlus default
      d.ismuon[i] = algoEtc::SimpleSecMuonFinder(tr, 5,5,5, -0.1, 0.2, 0.8, 1.5, 4, 0.5, privtx);
      d.iselectron[i] = algoEtc::SimpleSecElectronFinder(tr, 5,5,5,5,0.98,0.9, 1.15, privtx);
      d.isphoton[i] = 0;
      d.ischargedhadron[i] = !(d.ismuon[i] || d.iselectron[i]);
      d.isneutralhadron[i] = 0;

    }

    for(i=0;i<nall;i++){
      if(order[i].second < d.jet_ntracks) continue;
      const Neutral *neu = neutrals[order[i].second-d.jet_ntracks];
      d.px[i] = neu->Px();
      d.py[i] = neu->Py();
      d.pz[i] = neu->Pz();
      d.e[i] = neu->E();
      d.efrac[i] = neu->E() / jet->E();
      d.dtheta[i] = neu->Theta() - jet_theta;
      d.dphi[i] = neu->Phi() - jet_phi;

      d.d0[i] = 0.;
      d.d0sig[i] = 0.;
      d.z0[i] = 0.;
      d.z0sig[i] = 0.;
      d.ip3d[i] = 0.;
      d.ip3dsig[i] = 0.;

      d.charge[i] = 0;
      d.ismuon[i] = 0;
      d.iselectron[i] = 0;
      // simple photon finder
      double ecaldep = neu->getCaloEdep()[tpar::ecal];
      double hcaldep = neu->getCaloEdep()[tpar::hcal];
      d.isphoton[i] = (ecaldep / (ecaldep + hcaldep) > 0.98);
      d.ischargedhadron[i] = 0;
      d.isneutralhadron[i] = !d.isphoton[i];
    }

    _ntp->Fill();
  }
}

void DNNProvider2::end() {
  _file->Write();
  _file->Close();
}

}

