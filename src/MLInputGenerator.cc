#include <string>

#include "TFile.h"
#include "TNtuple.h"
#include "TNtupleD.h"
#include "TSystem.h"
#include "TPad.h"
#include "TStyle.h"

#include "lcfiplus.h"
#include "process.h"
#include "MLInputGenerator.h"
#include "VertexSelector.h"
#include "algoEtc.h"
#include "VertexFinderSuehara.h"
#include "VertexFitterSimple.h"

#include <utility>
#include <algorithm>

#include <map>
#include <functional>

using namespace lcfiplus;
using namespace std;

namespace lcfiplus {

// use namespace instead of class for static data
namespace MLInputGenerator {
  void init(){

    // jet inputs
    varJet["jet_px"] = [](const Jet* jet){ return jet->Px(); };
    varJet["jet_py"] = [](const Jet* jet){ return jet->Py(); };
    varJet["jet_pz"] = [](const Jet* jet){ return jet->Pz(); };
    varJet["jet_mass"] = [](const Jet* jet){ return jet->M(); };
    varJet["jet_ntracks"] = [](const Jet* jet){ return jet->getAllTracks().size(); };
    varJet["jet_nneutrals"] = [](const Jet* jet){ return jet->getNeutrals().size(); };
    varJet["jet_phi"] = [](const Jet* jet){ return jet->Phi(); };
    varJet["jet_theta"] = [](const Jet* jet){ return jet->Theta(); };

    // track inputs
    // particle kinematics
    varTrack["tr_px"] = [](const Track* tr){ return tr->Px(); };
    varTrack["tr_py"] = [](const Track* tr){ return tr->Py(); };
    varTrack["tr_pz"] = [](const Track* tr){ return tr->Pz(); };
    varTrack["tr_e"] = [](const Track* tr){ return tr->E(); };
    varTrackJet["tr_efrac"] = [](const Track* tr, const Jet* jet){ return tr->E() / jet->E(); };
    varTrackJet["tr_erel_log"] = [](const Track* tr, const Jet* jet){ return log10( tr->E() / jet->E() ); };
    varTrackJet["tr_thetarel"] = [](const Track* tr, const Jet* jet){
      float theta, phi;
      calc_thetaphi( tr->Vect(), jet->Vect(), theta, phi);
      return theta;
    };
    varTrackJet["tr_phirel"] = [](const Track* tr, const Jet* jet){
      float theta, phi;
      calc_thetaphi( tr->Vect(), jet->Vect(), theta, phi);
      return phi;
    };
    varTrackJet["tr_thetarel_ilc"] = [](const Track* tr, const Jet* jet){ return tr->Theta() - jet->Theta(); };
    varTrackJet["tr_phirel_ilc"] = [](const Track* tr, const Jet* jet){
      auto ret = tr->Phi() - jet->Phi();
      if(ret < -TMath::Pi()) ret += TMath::Pi() * 2;
      if(ret > TMath::Pi()) ret -= TMath::Pi() * 2;
      return ret;
    };

    // track errors
    varTrack["tr_dptdpt"] = [](const Track* tr){ return tr->getCovMatrix()[tpar::omom]; };
    varTrack["tr_detadeta"] = [](const Track* tr){ return tr->getCovMatrix()[tpar::tdtd]; };
    varTrack["tr_dphidphi"] = [](const Track* tr){ return tr->getCovMatrix()[tpar::phph]; };
    varTrack["tr_dxydxy"] = [](const Track* tr){ return tr->getCovMatrix()[tpar::d0d0]; };
    varTrack["tr_dzdz"] = [](const Track* tr){ return tr->getCovMatrix()[tpar::z0z0]; };
    varTrack["tr_dxydz"] = [](const Track* tr){ return tr->getCovMatrix()[tpar::d0z0]; };
    varTrack["tr_dphidxy"] = [](const Track* tr){ return tr->getCovMatrix()[tpar::d0ph]; };
    varTrack["tr_dlambdadz"] = [](const Track* tr){ return tr->getCovMatrix()[tpar::z0td]; };
    varTrack["tr_dxyc"] = [](const Track* tr){ return tr->getCovMatrix()[tpar::d0om]; };
    varTrack["tr_dxyctgtheta"] = [](const Track* tr){ return tr->getCovMatrix()[tpar::d0td]; };
    varTrack["tr_phic"] = [](const Track* tr){ return tr->getCovMatrix()[tpar::phom]; };
    varTrack["tr_phidz"] = [](const Track* tr){ return tr->getCovMatrix()[tpar::z0ph]; };
    varTrack["tr_phictgtheta"] = [](const Track* tr){ return tr->getCovMatrix()[tpar::phtd]; };
    varTrack["tr_cdz"] = [](const Track* tr){ return tr->getCovMatrix()[tpar::z0om]; };
    varTrack["tr_cctgtheta"] = [](const Track* tr){ return tr->getCovMatrix()[tpar::omtd]; };

    // particle displacements
    varTrack["tr_d0"] = [](const Track* tr){ return tr->getD0(); };
    varTrack["tr_d0sig"] = [](const Track* tr){ return tr->getD0() / sqrt(tr->getCovMatrix()[tpar::cov::d0d0]); };
    varTrack["tr_z0"] = [](const Track* tr){ return tr->getZ0(); };
    varTrack["tr_z0sig"] = [](const Track* tr){ return tr->getZ0() / sqrt(tr->getCovMatrix()[tpar::cov::z0z0]); };

    varTrack["tr_ip3d"] = [](const Track* tr){
      return sqrt(tr->getD0() * tr->getD0() + tr->getZ0() * tr->getZ0());
    };

    varTrack["tr_ip3dsig"] = [](const Track* tr){
      auto ip3d = sqrt(tr->getD0() * tr->getD0() + tr->getZ0() * tr->getZ0());
      return ip3d / sqrt(tr->getCovMatrix()[tpar::cov::d0d0] + tr->getCovMatrix()[tpar::cov::z0z0] + 2 * tr->getCovMatrix()[tpar::cov::d0z0]);
    };

    varTrackVertex["tr_dxy"] = [](const Track* tr, const Vertex* privtx){
      return calc_dxy(tr->getD0(), tr->getZ0(), tr->getPhi(), tr->Vect(), privtx->getPos(), tr->getCharge());
    };

    varTrackVertex["tr_dz"] = [](const Track* tr, const Vertex* privtx){
      return calc_dz(tr->getD0(), tr->getZ0(), tr->getPhi(), tr->Vect(), privtx->getPos(), tr->getCharge());
    };

    varTrackJet["tr_btagSip2dVal"] = [](const Track* tr, const Jet* jet){
      return calc_sip2d(tr->getD0(), tr->getPhi(), jet->Px(), jet->Py());
    };

    varTrackJet["tr_btagSip2dSig"] = [](const Track* tr, const Jet* jet){

      return calc_sip2d(tr->getD0(), tr->getPhi(), jet->Px(), jet->Py()) / sqrt( tr->getCovMatrix()[tpar::d0d0] );
    };

    varTrackJet["tr_btagSip3dVal"] = [](const Track* tr, const Jet* jet){
      return calc_sip3d(tr->getD0(), tr->getZ0(), tr->getPhi(), jet->Vect());
    };

    varTrackJet["tr_btagSip3dSig"] = [](const Track* tr, const Jet* jet){
      return calc_sip3d(tr->getD0(), tr->getZ0(), tr->getPhi(), jet->Vect()) / sqrt( tr->getCovMatrix()[tpar::d0d0] + tr->getCovMatrix()[tpar::z0z0] );
    };

    varTrackJet["tr_btagJetDistVal"] = [](const Track* tr, const Jet* jet){
      return calc_jetDist(tr->getD0(), tr->getZ0(), tr->getPhi(), tr->Vect(), jet->Vect());
    };

    varTrackJet["tr_btagJetDistSig"] = [](const Track* tr, const Jet* jet){
      return calc_jetDist(tr->getD0(), tr->getZ0(), tr->getPhi(), tr->Vect(), jet->Vect()) / sqrt( tr->getCovMatrix()[tpar::d0d0] + tr->getCovMatrix()[tpar::z0z0] );
    };

    // particle ID
    varTrack["tr_dEdx"] = [](const Track* tr){ return tr->getdEdx(); };
    varTrack["tr_charge"] = [](const Track* tr){ return tr->getCharge(); };
    varTrack["tr_isMu"] = [](const Track* tr){ return tr->getParticleIDProbability("muonProbability"); };
    varTrack["tr_isEl"] = [](const Track* tr){ return tr->getParticleIDProbability("electronProbability"); };
    varTrack["tr_isGamma"] = [](const Track*){ return 0.; };
    varTrack["tr_isChargedHad"] = [](const Track* tr){
      auto isMu = tr->getParticleIDProbability("muonProbability");
      auto isEl = tr->getParticleIDProbability("electronProbability");
      return !(isMu || isEl);
    };
    varTrack["tr_isNeutralHad"] = [](const Track*){ return 0.; };
    varTrack["tr_type"] = [](const Track* tr){ return tr->getPDG(); };
    varTrack["tr_mcpid"] = [](const Track* tr){ return (tr->getMcp()) ? tr->getMcp()->getId() : 0.; };
    varTrack["tr_mcp_pdg"] = [](const Track* tr){ return (tr->getMcp()) ? tr->getMcp()->getPDG() : 0.; };
    varTrack["tr_Ktype"] = [](const Track* tr){
      auto pdg = abs(tr->getPDG());
      if (pdg==321 || pdg==310) return 1;
      return 0;
    };
    varTrack["tr_isPion"] = [](const Track* tr){ return tr->getParticleIDProbability("pionProbability"); };
    varTrack["tr_isKaon"] = [](const Track* tr){ return tr->getParticleIDProbability("kaonProbability"); };
    varTrack["tr_isProton"] = [](const Track* tr){ return tr->getParticleIDProbability("protonProbability"); };

    // for PI3 (20240203)
    varTrack["tr_proton_K"] = [](const Track* tr){ 
      auto isproton = tr->getParticleIDProbability("protonProbability");
      auto iskaon = tr->getParticleIDProbability("kaonProbability");
      return isproton-iskaon;
    };
    varTrack["tr_pion_K"] = [](const Track* tr){ 
      auto ispion = tr->getParticleIDProbability("pionProbability");
      auto iskaon = tr->getParticleIDProbability("kaonProbability");
      return ispion-iskaon;
    };
    varTrack["tr_proton_Klike"] = [](const Track* tr){ 
      auto isproton = tr->getParticleIDProbability("protonLikelihood");
      auto iskaon = tr->getParticleIDProbability("kaonLikelihood");
      return isproton-iskaon;
    };
    varTrack["tr_pion_Klike"] = [](const Track* tr){ 
      auto ispion = tr->getParticleIDProbability("pionLikelihood");
      auto iskaon = tr->getParticleIDProbability("kaonLikelihood");
      return ispion-iskaon;
    };
    varTrack["tr_dEdxEl"] = [](const Track* tr){ return tr->getParticleIDProbability("electron_dEdxdistance"); };
    varTrack["tr_dEdxMu"] = [](const Track* tr){ return tr->getParticleIDProbability("muon_dEdxdistance"); };
    varTrack["tr_dEdxPion"] = [](const Track* tr){ return tr->getParticleIDProbability("pion_dEdxdistance"); };
    varTrack["tr_dEdxKaon"] = [](const Track* tr){ return tr->getParticleIDProbability("kaon_dEdxdistance"); };
    varTrack["tr_dEdxProton"] = [](const Track* tr){ return tr->getParticleIDProbability("proton_dEdxdistance"); };
    varTrack["tr_iselectronlike"] = [](const Track* tr){ return tr->getParticleIDProbability("electronLikelihood"); };
    varTrack["tr_ismuonlike"] = [](const Track* tr){ return tr->getParticleIDProbability("muonLikelihood"); };
    varTrack["tr_ispionlike"] = [](const Track* tr){ return tr->getParticleIDProbability("pionLikelihood"); };
    varTrack["tr_iskaonlike"] = [](const Track* tr){ return tr->getParticleIDProbability("kaonLikelihood"); };
    varTrack["tr_isprotonlike"] = [](const Track* tr){ return tr->getParticleIDProbability("protonLikelihood"); };

    // for neutrals
    // particle kinematics
    varNeutral["neu_px"] = [](const Neutral* neu){ return neu->Px(); };
    varNeutral["neu_py"] = [](const Neutral* neu){ return neu->Py(); };
    varNeutral["neu_pz"] = [](const Neutral* neu){ return neu->Pz(); };
    varNeutral["neu_e"] = [](const Neutral* neu){ return neu->E(); };
    varNeutralJet["neu_efrac"] = [](const Neutral* neu, const Jet* jet){ return neu->E() / jet->E(); };
    varNeutralJet["neu_erel_log"] = [](const Neutral* neu, const Jet* jet){ return log10(neu->E() / jet->E()); };
    varNeutralJet["neu_thetarel"] = [](const Neutral* neu, const Jet* jet){
      float theta, phi;
      calc_thetaphi(jet->Vect(), neu->Vect(), theta, phi);
      return theta;
    };
    varNeutralJet["neu_phirel"] = [](const Neutral* neu, const Jet* jet){
      float theta, phi;
      calc_thetaphi(jet->Vect(), neu->Vect(), theta, phi);
      return phi;
    };
    varNeutralJet["neu_thetarel_ilc"] = [](const Neutral* neu, const Jet* jet){ return neu->Theta() - jet->Theta(); };
    varNeutralJet["neu_phirel_ilc"] = [](const Neutral* neu, const Jet* jet){
      auto ret = neu->Phi() - jet->Phi();
      if(ret < -TMath::Pi()) ret += TMath::Pi() * 2;
      if(ret > TMath::Pi()) ret -= TMath::Pi() * 2;
      return ret;
    };

    // track errors -- why are these implemented for neutrals?
    varNeutral["neu_dptdpt"] = [](const Neutral*){ return -9; };
    varNeutral["neu_detadeta"] = [](const Neutral*){ return -9; };
    varNeutral["neu_dphidphi"] = [](const Neutral*){ return -9; };
    varNeutral["neu_dxydxy"] = [](const Neutral*){ return -9; };
    varNeutral["neu_dzdz"] = [](const Neutral*){ return -9; };
    varNeutral["neu_dxydz"] = [](const Neutral*){ return -9; };
    varNeutral["neu_dphidxy"] = [](const Neutral*){ return -9; };
    varNeutral["neu_dlambdadz"] = [](const Neutral*){ return -9; };
    varNeutral["neu_dxyc"] = [](const Neutral*){ return -9; };
    varNeutral["neu_dxyctgtheta"] = [](const Neutral*){ return -9; };
    varNeutral["neu_phic"] = [](const Neutral*){ return -9; };
    varNeutral["neu_phidz"] = [](const Neutral*){ return -9; };
    varNeutral["neu_phictgtheta"] = [](const Neutral*){ return -9; };
    varNeutral["neu_cdz"] = [](const Neutral*){ return -9; };
    varNeutral["neu_cctgtheta"] = [](const Neutral*){ return -9; };

    // particle displacements -- why are these implemented for neutrals?
    varNeutral["neu_d0"] = [](const Neutral*){ return -9; };
    varNeutral["neu_d0sig"] = [](const Neutral*){ return -9; };
    varNeutral["neu_z0"] = [](const Neutral*){ return -9; };
    varNeutral["neu_z0sig"] = [](const Neutral*){ return -9; };
    varNeutral["neu_ip3d"] = [](const Neutral*){ return -9; };
    varNeutral["neu_ip3dsig"] = [](const Neutral*){ return -9; };

    varNeutral["neu_dxy"] = [](const Neutral*){ return -9; };
    varNeutral["neu_dz"] = [](const Neutral*){ return -9; };
    varNeutral["neu_btagSip2dVal"] = [](const Neutral*){ return -9; };
    varNeutral["neu_btagSip2dSig"] = [](const Neutral*){ return -9; };
    varNeutral["neu_btagSip3dVal"] = [](const Neutral*){ return -9; };
    varNeutral["neu_btagSip3dSig"] = [](const Neutral*){ return -9; };
    varNeutral["neu_btagJetDistVal"] = [](const Neutral*){ return -9; };
    varNeutral["neu_btagJetDistSig"] = [](const Neutral*){ return -9; };

    // particle ID
    varNeutral["neu_charge"] = [](const Neutral*){ return 0; };
    varNeutral["neu_isMu"] = [](const Neutral*){ return 0; };
    varNeutral["neu_isEl"] = [](const Neutral*){ return 0; };
    varNeutral["neu_isGamma"] = [](const Neutral* neu){
      // simple photon finder
      double ecaldep = neu->getCaloEdep()[tpar::ecal];
      double hcaldep = neu->getCaloEdep()[tpar::hcal];
      return (ecaldep / (ecaldep + hcaldep) > 0.98);
    };
    varNeutral["neu_isChargedHad"] = [](const Neutral*){ return 0; };
    varNeutral["neu_isNeutralHad"] = [](const Neutral* neu){ 
      // simple photon finder
      double ecaldep = neu->getCaloEdep()[tpar::ecal];
      double hcaldep = neu->getCaloEdep()[tpar::hcal];
      return !(ecaldep / (ecaldep + hcaldep) > 0.98);
    };
    varNeutral["neu_type"] = [](const Neutral* neu){ return neu->getPDG(); };
    varNeutral["neu_mcpid"] = [](const Neutral* neu){ return (neu->getMcp()) ? neu->getMcp()->getId() : 0.; };
    varNeutral["neu_mcp_pdg"] = [](const Neutral* neu){ return (neu->getMcp()) ? neu->getMcp()->getPDG() : 0.; };

    varNeutral["neu_Ktype"] = [](const Neutral* neu){
      auto pdg = abs(neu->getPDG());
      return ( (pdg==321 || pdg==310) );
    };
    varNeutral["neu_isPion"] = [](const Neutral*){ return 0; };
    varNeutral["neu_isKaon"] = [](const Neutral*){ return 0; };
    varNeutral["neu_isProton"] = [](const Neutral*){ return 0; };

  }

#if 0 // copied from DNNProvider2.cc
  // label
  _ntp->Branch("mc_b",&d.mc_b,"mc_b/I");
  _ntp->Branch("mc_c",&d.mc_c,"mc_c/I");
  _ntp->Branch("mc_u",&d.mc_u,"mc_u/I");
  _ntp->Branch("mc_d",&d.mc_d,"mc_d/I");
  _ntp->Branch("mc_s",&d.mc_s,"mc_s/I");
  _ntp->Branch("mc_g",&d.mc_g,"mc_g/I");
  _ntp->Branch("mc_q",&d.mc_q,"mc_q/I"); // usdg
#endif

  // copied from FCCANalyses/analyzers/dataframe/src/ReconstructedParticle2Track.cc
  float calc_dxy(float D0_wrt0, float Z0_wrt0, float phi0_wrt0, TVector3 p, TVector3 privtx, int charge){
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

  float calc_dz(float D0_wrt0, float Z0_wrt0, float phi0_wrt0, TVector3 p, TVector3 privtx, int charge){
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

  float calc_sip2d(float D0, float phi0, float jetpx, float jetpy){
    TVector2 p(jetpx, jetpy);
    TVector2 d0(-D0 * TMath::Sin(phi0), D0 * TMath::Cos(phi0));
    return TMath::Sign(1, d0 * p) * fabs(D0);
  }

  float calc_sip3d(float D0, float Z0, float phi0, TVector3 p_jet){
    TVector3 d(-D0 * TMath::Sin(phi0), D0 * TMath::Cos(phi0), Z0);
    return TMath::Sign(1, d * p_jet) * fabs(sqrt(D0 * D0 + Z0 * Z0));
  }

  float calc_jetDist(float D0, float Z0, float phi0, TVector3 p_ct, TVector3 p_jet){
    TVector3 d(-D0 * TMath::Sin(phi0), D0 * TMath::Cos(phi0), Z0);
    TVector3 r_jet(0.0, 0.0, 0.0);
    TVector3 n = p_ct.Cross(p_jet).Unit(); // What if they are parallel?
    return n.Dot(d - r_jet);
  }

  void calc_thetaphi(TVector3 jet, TVector3 part, float &theta, float &phi){
    part.RotateZ(-jet.Phi());
    part.RotateY(-jet.Theta());
    theta = part.Theta();
    phi = part.Phi();  
  }

}
}