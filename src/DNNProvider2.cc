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
  _mcIsU = param->get("DNNProvider2.MC_IsU",int(0));
  _mcIsD = param->get("DNNProvider2.MC_IsD",int(0));
  _mcIsS = param->get("DNNProvider2.MC_IsS",int(0));
  _mcIsG = param->get("DNNProvider2.MC_IsG",int(0)); 

  _jets = 0;

  if(_mcIsB + _mcIsC + _mcIsU + _mcIsD + _mcIsS + _mcIsG != 1){
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
  _ntp->Branch("jet_mass",&d.jet_mass,"jet_mass/F");
  _ntp->Branch("jet_ntracks",&d.jet_ntracks,"jet_ntracks/I");
  _ntp->Branch("jet_nneutrals",&d.jet_nneutrals,"jet_nneutrals/I");

  _ntp->Branch("jet_phi",&d.jet_phi,"jet_phi/F");
  _ntp->Branch("jet_theta",&d.jet_theta,"jet_theta/F");

  // particle kinematics
  _ntp->Branch("pfcand_px",&d.px);
  _ntp->Branch("pfcand_py",&d.py);
  _ntp->Branch("pfcand_pz",&d.pz);
  _ntp->Branch("pfcand_e",&d.e);
  _ntp->Branch("pfcand_efrac",&d.efrac);
  _ntp->Branch("pfcand_erel_log",&d.erel_log);
  _ntp->Branch("pfcand_thetarel",&d.dtheta);
  _ntp->Branch("pfcand_phirel",&d.dphi);
  _ntp->Branch("pfcand_thetarel_ilc",&d.dtheta_ilc);
  _ntp->Branch("pfcand_phirel_ilc",&d.dphi_ilc);

  // track errors
  _ntp->Branch("pfcand_dptdpt",&d.cov_omega);
  _ntp->Branch("pfcand_detadeta",&d.cov_tanlambda);
  _ntp->Branch("pfcand_dphidphi",&d.cov_phi);
  _ntp->Branch("pfcand_dxydxy",&d.cov_d0);
  _ntp->Branch("pfcand_dzdz",&d.cov_z0);
  _ntp->Branch("pfcand_dxydz",&d.cov_d0_z0);
  _ntp->Branch("pfcand_dphidxy",&d.cov_d0_phi);
  _ntp->Branch("pfcand_dlambdadz",&d.cov_z0_tanlambda);
  _ntp->Branch("pfcand_dxyc",&d.cov_d0_omega);
  _ntp->Branch("pfcand_dxyctgtheta",&d.cov_d0_tanlambda);
  _ntp->Branch("pfcand_phic",&d.cov_phi_omega);
  _ntp->Branch("pfcand_phidz",&d.cov_z0_phi);
  _ntp->Branch("pfcand_phictgtheta",&d.cov_phi_tanlambda);
  _ntp->Branch("pfcand_cdz",&d.cov_z0_omega);
  _ntp->Branch("pfcand_cctgtheta",&d.cov_omega_tanlambda);

  // particle displacements
  _ntp->Branch("d0",&d.d0);
  _ntp->Branch("d0sig",&d.d0sig);
  _ntp->Branch("z0",&d.z0);
  _ntp->Branch("z0sig",&d.z0sig);
  _ntp->Branch("ip3d",&d.ip3d);
  _ntp->Branch("ip3dsig",&d.ip3dsig);

  _ntp->Branch("pfcand_dxy",&d.dxy);
  _ntp->Branch("pfcand_dz",&d.dz);
  _ntp->Branch("pfcand_btagSip2dVal",&d.ip2d_fcc);
  _ntp->Branch("pfcand_btagSip2dSig",&d.ip2dsig_fcc);
  _ntp->Branch("pfcand_btagSip3dVal",&d.ip3d_fcc);
  _ntp->Branch("pfcand_btagSip3dSig",&d.ip3dsig_fcc);
  _ntp->Branch("pfcand_btagJetDistVal",&d.jetdist_fcc);
  _ntp->Branch("pfcand_btagJetDistSig",&d.jetdistsig_fcc);

  // particle ID
  _ntp->Branch("pfcand_charge",&d.charge);
  _ntp->Branch("pfcand_isMu",&d.ismuon);
  _ntp->Branch("pfcand_isEl",&d.iselectron);
  _ntp->Branch("pfcand_isGamma",&d.isphoton);
  _ntp->Branch("pfcand_isChargedHad",&d.ischargedhadron);
  _ntp->Branch("pfcand_isNeutralHad",&d.isneutralhadron);
  _ntp->Branch("pfcand_type",&d.pdg_pfa);

  // label
  _ntp->Branch("mc_b",&d.mc_b,"mc_b/I");
  _ntp->Branch("mc_c",&d.mc_c,"mc_c/I");
  _ntp->Branch("mc_u",&d.mc_u,"mc_u/I");
  _ntp->Branch("mc_d",&d.mc_d,"mc_d/I");
  _ntp->Branch("mc_s",&d.mc_s,"mc_s/I");
  _ntp->Branch("mc_g",&d.mc_g,"mc_g/I");
  _ntp->Branch("mc_q",&d.mc_q,"mc_q/I"); // usdg
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

    d.jet_theta = jet->Theta();
    d.jet_phi = jet->Phi();

    d.jet_ntracks = tracks.size();
    d.jet_nneutrals = neutrals.size();

    //float jet_theta = jet->Theta();
    //float jet_phi = jet->Phi();

    d.mc_b = _mcIsB;
    d.mc_c = _mcIsC;
    d.mc_u = _mcIsU;
    d.mc_d = _mcIsD;
    d.mc_s = _mcIsS;
    d.mc_g = _mcIsG;
    d.mc_q = _mcIsU || _mcIsD || _mcIsS || _mcIsG;

    // probably order of tracks/netural does not matter...
    int nall = d.jet_ntracks + d.jet_nneutrals;
    d.px.resize(nall);
    d.py.resize(nall);
    d.pz.resize(nall);
    d.e.resize(nall);
    d.efrac.resize(nall);
    d.erel_log.resize(nall);
    d.dtheta.resize(nall);
    d.dphi.resize(nall);
    d.dtheta_ilc.resize(nall);
    d.dphi_ilc.resize(nall);

    d.cov_d0.resize(nall);
    d.cov_z0.resize(nall);
    d.cov_phi.resize(nall);
    d.cov_omega.resize(nall);
    d.cov_tanlambda.resize(nall);

    d.cov_d0_z0.resize(nall);
    d.cov_d0_phi.resize(nall);
    d.cov_d0_omega.resize(nall);
    d.cov_d0_tanlambda.resize(nall);

    d.cov_z0_phi.resize(nall);
    d.cov_z0_omega.resize(nall);
    d.cov_z0_tanlambda.resize(nall);

    d.cov_phi_omega.resize(nall);
    d.cov_phi_tanlambda.resize(nall);
    d.cov_omega_tanlambda.resize(nall);

    d.d0.resize(nall);
    d.d0sig.resize(nall);
    d.z0.resize(nall);
    d.z0sig.resize(nall);
    d.ip3d.resize(nall);
    d.ip3dsig.resize(nall);

    d.dxy.resize(nall);
    d.dz.resize(nall);
    d.ip2d_fcc.resize(nall);
    d.ip2dsig_fcc.resize(nall);
    d.ip3d_fcc.resize(nall);
    d.ip3dsig_fcc.resize(nall);
    d.jetdist_fcc.resize(nall);
    d.jetdistsig_fcc.resize(nall);

    d.charge.resize(nall);
    d.ismuon.resize(nall);
    d.iselectron.resize(nall);
    d.isphoton.resize(nall);
    d.ischargedhadron.resize(nall);
    d.isneutralhadron.resize(nall);
    d.pdg_pfa.resize(nall);

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
      d.erel_log[i] = log10(d.efrac[i]);
      d.dtheta_ilc[i] = tr->Theta() - jet->Theta();
      d.dphi_ilc[i] = tr->Phi() - jet->Phi();
      if(d.dphi_ilc[i] < -TMath::Pi())d.dphi_ilc[i] += TMath::Pi() * 2;
      if(d.dphi_ilc[i] > TMath::Pi())d.dphi_ilc[i] -= TMath::Pi() * 2;
      calc_thetaphi(jet->Vect(), tr->Vect(), d.dtheta[i], d.dphi[i]);

      // track covmatrix
      d.cov_d0[i] = tr->getCovMatrix()[tpar::d0d0];
      d.cov_z0[i] = tr->getCovMatrix()[tpar::z0z0];
      d.cov_phi[i] = tr->getCovMatrix()[tpar::phph];
      d.cov_omega[i] = tr->getCovMatrix()[tpar::omom];
      d.cov_tanlambda[i] = tr->getCovMatrix()[tpar::tdtd];

      d.cov_d0_z0[i] = tr->getCovMatrix()[tpar::d0z0];
      d.cov_d0_phi[i] = tr->getCovMatrix()[tpar::d0ph];
      d.cov_d0_omega[i] = tr->getCovMatrix()[tpar::d0om];
      d.cov_d0_tanlambda[i] = tr->getCovMatrix()[tpar::d0td];

      d.cov_z0_phi[i] = tr->getCovMatrix()[tpar::z0ph];
      d.cov_z0_omega[i] = tr->getCovMatrix()[tpar::z0om];
      d.cov_z0_tanlambda[i] = tr->getCovMatrix()[tpar::z0td];

      d.cov_phi_omega[i] = tr->getCovMatrix()[tpar::phom];
      d.cov_phi_tanlambda[i] = tr->getCovMatrix()[tpar::phtd];
      d.cov_omega_tanlambda[i] = tr->getCovMatrix()[tpar::omtd];

      d.d0[i] = tr->getD0();
      d.d0sig[i] = tr->getD0() / sqrt(tr->getCovMatrix()[tpar::cov::d0d0]);
      d.z0[i] = tr->getZ0();
      d.z0sig[i] = tr->getZ0() / sqrt(tr->getCovMatrix()[tpar::cov::z0z0]);

      d.ip3d[i] = sqrt(tr->getD0() * tr->getD0() + tr->getZ0() * tr->getZ0());
      d.ip3dsig[i] = d.ip3d[i] / sqrt(tr->getCovMatrix()[tpar::cov::d0d0] + tr->getCovMatrix()[tpar::cov::z0z0] + 2 * tr->getCovMatrix()[tpar::cov::d0z0]);
      
      d.dxy[i] = calc_dxy(tr->getD0(), tr->getZ0(), tr->getPhi(), tr->Vect(), privtx->getPos(), tr->getCharge());
      d.dz[i] = calc_dz(tr->getD0(), tr->getZ0(), tr->getPhi(), tr->Vect(), privtx->getPos(), tr->getCharge());
      d.ip2d_fcc[i] = calc_sip2d(tr->getD0(), tr->getPhi(), jet->Px(), jet->Py());
      d.ip2dsig_fcc[i] = d.ip2d_fcc[i] / sqrt(d.cov_d0[i]);
      d.ip3d_fcc[i] = calc_sip3d(tr->getD0(), tr->getZ0(), tr->getPhi(), jet->Vect());
      d.ip3dsig_fcc[i] = d.ip3d_fcc[i] / sqrt(d.cov_d0[i] + d.cov_z0[i]);
      d.jetdist_fcc[i] = calc_jetDist(tr->getD0(), tr->getZ0(), tr->getPhi(), tr->Vect(), jet->Vect());
      d.jetdistsig_fcc[i] = d.jetdist_fcc[i] / sqrt(d.cov_d0[i] + d.cov_z0[i]);

      d.charge[i] = tr->getCharge();
      // tracing LCFIPlus default
      d.ismuon[i] = algoEtc::SimpleSecMuonFinder(tr, 5,5,5, -0.1, 0.2, 0.8, 1.5, 4, 0.5, privtx);
      d.iselectron[i] = algoEtc::SimpleSecElectronFinder(tr, 5,5,5,5,0.98,0.9, 1.15, privtx);
      d.isphoton[i] = 0;
      d.ischargedhadron[i] = !(d.ismuon[i] || d.iselectron[i]);
      d.isneutralhadron[i] = 0;
      d.pdg_pfa[i] = tr->getPDG();

    }

    for(i=0;i<nall;i++){
      if(order[i].second < d.jet_ntracks) continue;
      const Neutral *neu = neutrals[order[i].second-d.jet_ntracks];
      d.px[i] = neu->Px();
      d.py[i] = neu->Py();
      d.pz[i] = neu->Pz();
      d.e[i] = neu->E();
      d.efrac[i] = neu->E() / jet->E();
      d.erel_log[i] = log10(d.efrac[i]);
      d.dtheta_ilc[i] = neu->Theta() - jet->Theta();
      d.dphi_ilc[i] = neu->Phi() - jet->Phi();
      if(d.dphi_ilc[i] < -TMath::Pi())d.dphi_ilc[i] += TMath::Pi() * 2;
      if(d.dphi_ilc[i] > TMath::Pi())d.dphi_ilc[i] -= TMath::Pi() * 2;
      calc_thetaphi(jet->Vect(), neu->Vect(), d.dtheta[i], d.dphi[i]);

      d.cov_d0[i] = -9;
      d.cov_z0[i] = -9;
      d.cov_phi[i] = -9;
      d.cov_omega[i] = -9;
      d.cov_tanlambda[i] = -9;

      d.cov_d0_z0[i] = -9;
      d.cov_d0_phi[i] = -9;
      d.cov_d0_omega[i] = -9;
      d.cov_d0_tanlambda[i] = -9;

      d.cov_z0_phi[i] = -9;
      d.cov_z0_omega[i] = -9;
      d.cov_z0_tanlambda[i] = -9;

      d.cov_phi_omega[i] = -9;
      d.cov_phi_tanlambda[i] = -9;
      d.cov_omega_tanlambda[i] = -9;

      d.dxy[i] = -9;
      d.dz[i] = -9;
      d.ip2d_fcc[i] = -9;
      d.ip2dsig_fcc[i] = -9;
      d.ip3d_fcc[i] = -9;
      d.ip3dsig_fcc[i] = -9;
      d.jetdist_fcc[i] = -9;
      d.jetdistsig_fcc[i] = -9;

      d.charge[i] = 0;
      d.ismuon[i] = 0;
      d.iselectron[i] = 0;
      // simple photon finder
      double ecaldep = neu->getCaloEdep()[tpar::ecal];
      double hcaldep = neu->getCaloEdep()[tpar::hcal];
      d.isphoton[i] = (ecaldep / (ecaldep + hcaldep) > 0.98);
      d.ischargedhadron[i] = 0;
      d.isneutralhadron[i] = !d.isphoton[i];
      d.pdg_pfa[i] = neu->getPDG();
    }

    _ntp->Fill();
  }
}

void DNNProvider2::end() {
  _file->Write();
  _file->Close();
}

// copied from FCCANalyses/analyzers/dataframe/src/ReconstructedParticle2Track.cc

float DNNProvider2::calc_dxy(float D0_wrt0, float Z0_wrt0, float phi0_wrt0, TVector3 p, TVector3 privtx, int charge){
  double Bz = 3.5;
  double cSpeed = 2.99792458e8 * 1e-9;

  TVector3 X( - D0_wrt0 * TMath::Sin(phi0_wrt0) , D0_wrt0 * TMath::Cos(phi0_wrt0) , Z0_wrt0);
  TVector3 x = X - privtx;
  //std::cout<<"vertex: "<<V.Vect().X()<<", "<<V.Vect().Y()<<", "<<V.Vect().Z()<<", "<<std::endl;

  double a = - charge * Bz * cSpeed;
  double pt = p.Pt();
  double r2 = x(0) * x(0) + x(1) * x(1);
  double cross = x(0) * p(1) - x(1) * p(0);
  double D=-9;
  if (pt * pt - 2 * a * cross + a * a * r2 > 0) {
    double T = TMath::Sqrt(pt * pt - 2 * a * cross + a * a * r2);
    if (pt < 10.0) D = (T - pt) / a;
    else D = (-2 * cross + a * r2) / (T + pt);
  }
  return D;
}

float DNNProvider2::calc_dz(float D0_wrt0, float Z0_wrt0, float phi0_wrt0, TVector3 p, TVector3 privtx, int charge){
  double Bz = 3.5;
  double cSpeed = 2.99792458e8 * 1e-9;

  TVector3 X( - D0_wrt0 * TMath::Sin(phi0_wrt0) , D0_wrt0 * TMath::Cos(phi0_wrt0) , Z0_wrt0);
  TVector3 x = X - privtx;

  double a = - charge * Bz * cSpeed;
  double pt = p.Pt();
  double C = a/(2 * pt);
  double r2 = x(0) * x(0) + x(1) * x(1);
  double cross = x(0) * p(1) - x(1) * p(0);
  double T = TMath::Sqrt(pt * pt - 2 * a * cross + a * a * r2);
  double D;
  if (pt < 10.0) D = (T - pt) / a;
  else D = (-2 * cross + a * r2) / (T + pt);
  double B = C * TMath::Sqrt(TMath::Max(r2 - D * D, 0.0) / (1 + 2 * C * D));
  if ( TMath::Abs(B) > 1.) B = TMath::Sign(1, B);
  double st = TMath::ASin(B) / C;
  double ct = p(2) / pt;
  double z0;
  double dot = x(0) * p(0) + x(1) * p(1);
  if (dot > 0.0) z0 = x(2) - ct * st;
  else z0 = x(2) + ct * st;

  return z0;
}

float DNNProvider2::calc_sip2d(float D0, float phi0, float jetpx, float jetpy)
{
  TVector2 p(jetpx, jetpy);

  TVector2 d0(-D0 * TMath::Sin(phi0), D0 * TMath::Cos(phi0));
  return TMath::Sign(1, d0 * p) * fabs(D0);
}

float DNNProvider2::calc_sip3d(float D0, float Z0, float phi0, TVector3 p_jet)
{
  TVector3 d(-D0 * TMath::Sin(phi0), D0 * TMath::Cos(phi0), Z0);
  return TMath::Sign(1, d * p_jet) * fabs(sqrt(D0 * D0 + Z0 * Z0));
}

float DNNProvider2::calc_jetDist(float D0, float Z0, float phi0, TVector3 p_ct, TVector3 p_jet)
{
  TVector3 d(-D0 * TMath::Sin(phi0), D0 * TMath::Cos(phi0), Z0);
  TVector3 r_jet(0.0, 0.0, 0.0);
  TVector3 n = p_ct.Cross(p_jet).Unit(); // What if they are parallel?
  return n.Dot(d - r_jet);
}

void DNNProvider2::calc_thetaphi(TVector3 jet, TVector3 part, float &theta, float &phi){
  part.RotateZ(-jet.Phi());
  part.RotateY(-jet.Theta());

  theta = part.Theta();
  phi = part.Phi();  
}

}

