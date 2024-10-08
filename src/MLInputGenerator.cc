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
#include <variant>

using namespace lcfiplus;
using namespace std;

namespace lcfiplus {

// use namespace instead of class for static data
namespace MLInputGenerator {

  map<string, variant<
    function<double(const Jet*)>,
    function<double(const Track*)>,
    function<double(const Track*, const Jet*)>,
    function<double(const Track*, const Vertex*)>,
    function<double(const Neutral*)>,
    function<double(const Neutral*, const Jet*)>,
    function<double(const Neutral*, const Vertex*)>
  > > calcInput;

  bool _initialized(false);

  void init(){

    if (_initialized) {
      return;
    }
    _initialized = true;
    //cout << "MLInputGenerator: preparing functions for input variables" << endl;

    //const string _jet_prefix = "jet_";
    //const string _trk_prefix = "tr_";
    //const string _neu_prefix = "neu_";

    const string _jet_prefix = "jet_";
    const string _trk_prefix = "pfcand_";
    const string _neu_prefix = "neu_pfcand_";

    calcInput[_jet_prefix+"px"] = [](const Jet* jet){ return jet->Px(); };
    calcInput[_trk_prefix+"efrac"] = [](const Track* tr, const Jet* jet){ return tr->E() / jet->E(); };
    calcInput[_neu_prefix+"px"] = [](const Neutral* neu){ return neu->Px(); };

    // jet inputs
    calcInput[_jet_prefix+"px"] = [](const Jet* jet){ return jet->Px(); };
    calcInput[_jet_prefix+"py"] = [](const Jet* jet){ return jet->Py(); };
    calcInput[_jet_prefix+"pz"] = [](const Jet* jet){ return jet->Pz(); };
    calcInput[_jet_prefix+"mass"] = [](const Jet* jet){ return jet->M(); };
    calcInput[_jet_prefix+"ntracks"] = [](const Jet* jet){ return jet->getAllTracks().size(); };
    calcInput[_jet_prefix+"nneutrals"] = [](const Jet* jet){ return jet->getNeutrals().size(); };
    calcInput[_jet_prefix+"phi"] = [](const Jet* jet){ return jet->Phi(); };
    calcInput[_jet_prefix+"theta"] = [](const Jet* jet){ return jet->Theta(); };

    // track inputs
    // particle kinematics
    calcInput[_trk_prefix+"px"] = [](const Track* tr){ return tr->Px(); };
    calcInput[_trk_prefix+"py"] = [](const Track* tr){ return tr->Py(); };
    calcInput[_trk_prefix+"pz"] = [](const Track* tr){ return tr->Pz(); };
    calcInput[_trk_prefix+"e"] = [](const Track* tr){ return tr->E(); };
    calcInput[_trk_prefix+"efrac"] = [](const Track* tr, const Jet* jet){ return tr->E() / jet->E(); };
    calcInput[_trk_prefix+"erel_log"] = [](const Track* tr, const Jet* jet){ return log10( tr->E() / jet->E() ); };
    calcInput[_trk_prefix+"thetarel"] = [](const Track* tr, const Jet* jet){
      float theta, phi;
      calc_thetaphi( tr->Vect(), jet->Vect(), theta, phi);
      return theta;
    };
    calcInput[_trk_prefix+"phirel"] = [](const Track* tr, const Jet* jet){
      float theta, phi;
      calc_thetaphi( tr->Vect(), jet->Vect(), theta, phi);
      return phi;
    };
    calcInput[_trk_prefix+"thetarel_ilc"] = [](const Track* tr, const Jet* jet){ return tr->Theta() - jet->Theta(); };
    calcInput[_trk_prefix+"phirel_ilc"] = [](const Track* tr, const Jet* jet){
      auto ret = tr->Phi() - jet->Phi();
      if(ret < -TMath::Pi()) ret += TMath::Pi() * 2;
      if(ret > TMath::Pi()) ret -= TMath::Pi() * 2;
      return ret;
    };

    // track errors
    calcInput[_trk_prefix+"dptdpt"] = [](const Track* tr){ return tr->getCovMatrix()[tpar::omom]; };
    calcInput[_trk_prefix+"detadeta"] = [](const Track* tr){ return tr->getCovMatrix()[tpar::tdtd]; };
    calcInput[_trk_prefix+"dphidphi"] = [](const Track* tr){ return tr->getCovMatrix()[tpar::phph]; };
    calcInput[_trk_prefix+"dxydxy"] = [](const Track* tr){ return tr->getCovMatrix()[tpar::d0d0]; };
    calcInput[_trk_prefix+"dzdz"] = [](const Track* tr){ return tr->getCovMatrix()[tpar::z0z0]; };
    calcInput[_trk_prefix+"dxydz"] = [](const Track* tr){ return tr->getCovMatrix()[tpar::d0z0]; };
    calcInput[_trk_prefix+"dphidxy"] = [](const Track* tr){ return tr->getCovMatrix()[tpar::d0ph]; };
    calcInput[_trk_prefix+"dlambdadz"] = [](const Track* tr){ return tr->getCovMatrix()[tpar::z0td]; };
    calcInput[_trk_prefix+"dxyc"] = [](const Track* tr){ return tr->getCovMatrix()[tpar::d0om]; };
    calcInput[_trk_prefix+"dxyctgtheta"] = [](const Track* tr){ return tr->getCovMatrix()[tpar::d0td]; };
    calcInput[_trk_prefix+"phic"] = [](const Track* tr){ return tr->getCovMatrix()[tpar::phom]; };
    calcInput[_trk_prefix+"phidz"] = [](const Track* tr){ return tr->getCovMatrix()[tpar::z0ph]; };
    calcInput[_trk_prefix+"phictgtheta"] = [](const Track* tr){ return tr->getCovMatrix()[tpar::phtd]; };
    calcInput[_trk_prefix+"cdz"] = [](const Track* tr){ return tr->getCovMatrix()[tpar::z0om]; };
    calcInput[_trk_prefix+"cctgtheta"] = [](const Track* tr){ return tr->getCovMatrix()[tpar::omtd]; };

    // particle displacements
    calcInput[_trk_prefix+"d0"] = [](const Track* tr){ return tr->getD0(); };
    calcInput[_trk_prefix+"d0sig"] = [](const Track* tr){ return tr->getD0() / sqrt(tr->getCovMatrix()[tpar::cov::d0d0]); };
    calcInput[_trk_prefix+"z0"] = [](const Track* tr){ return tr->getZ0(); };
    calcInput[_trk_prefix+"z0sig"] = [](const Track* tr){ return tr->getZ0() / sqrt(tr->getCovMatrix()[tpar::cov::z0z0]); };

    calcInput[_trk_prefix+"ip3d"] = [](const Track* tr){
      return sqrt(tr->getD0() * tr->getD0() + tr->getZ0() * tr->getZ0());
    };

    calcInput[_trk_prefix+"ip3dsig"] = [](const Track* tr){
      auto ip3d = sqrt(tr->getD0() * tr->getD0() + tr->getZ0() * tr->getZ0());
      return ip3d / sqrt(tr->getCovMatrix()[tpar::cov::d0d0] + tr->getCovMatrix()[tpar::cov::z0z0] + 2 * tr->getCovMatrix()[tpar::cov::d0z0]);
    };

    calcInput[_trk_prefix+"dxy"] = [](const Track* tr, const Vertex* privtx){
      return calc_dxy(tr->getD0(), tr->getZ0(), tr->getPhi(), tr->Vect(), privtx->getPos(), tr->getCharge());
    };

    calcInput[_trk_prefix+"dz"] = [](const Track* tr, const Vertex* privtx){
      return calc_dz(tr->getD0(), tr->getZ0(), tr->getPhi(), tr->Vect(), privtx->getPos(), tr->getCharge());
    };

    calcInput[_trk_prefix+"btagSip2dVal"] = [](const Track* tr, const Jet* jet){
      return calc_sip2d(tr->getD0(), tr->getPhi(), jet->Px(), jet->Py());
    };

    calcInput[_trk_prefix+"btagSip2dSig"] = [](const Track* tr, const Jet* jet){

      return calc_sip2d(tr->getD0(), tr->getPhi(), jet->Px(), jet->Py()) / sqrt( tr->getCovMatrix()[tpar::d0d0] );
    };

    calcInput[_trk_prefix+"btagSip3dVal"] = [](const Track* tr, const Jet* jet){
      return calc_sip3d(tr->getD0(), tr->getZ0(), tr->getPhi(), jet->Vect());
    };

    calcInput[_trk_prefix+"btagSip3dSig"] = [](const Track* tr, const Jet* jet){
      return calc_sip3d(tr->getD0(), tr->getZ0(), tr->getPhi(), jet->Vect()) / sqrt( tr->getCovMatrix()[tpar::d0d0] + tr->getCovMatrix()[tpar::z0z0] );
    };

    calcInput[_trk_prefix+"btagJetDistVal"] = [](const Track* tr, const Jet* jet){
      return calc_jetDist(tr->getD0(), tr->getZ0(), tr->getPhi(), tr->Vect(), jet->Vect());
    };

    calcInput[_trk_prefix+"btagJetDistSig"] = [](const Track* tr, const Jet* jet){
      return calc_jetDist(tr->getD0(), tr->getZ0(), tr->getPhi(), tr->Vect(), jet->Vect()) / sqrt( tr->getCovMatrix()[tpar::d0d0] + tr->getCovMatrix()[tpar::z0z0] );
    };

    // particle ID

    // no prefix in existing weight files
    //calcInput[_trk_prefix+"dEdx"] = [](const Track* tr){ return tr->getdEdx(); };
    calcInput["dEdx"] = [](const Track* tr){ return tr->getdEdx(); };

    calcInput[_trk_prefix+"charge"] = [](const Track* tr){ return tr->getCharge(); };
    calcInput[_trk_prefix+"isMu"] = [](const Track* tr){ return tr->getParticleIDProbability("muonProbability"); };
    calcInput[_trk_prefix+"isEl"] = [](const Track* tr){ return tr->getParticleIDProbability("electronProbability"); };
    calcInput[_trk_prefix+"isGamma"] = [](const Track*){ return 0.; };
    calcInput[_trk_prefix+"isChargedHad"] = [](const Track* tr){
      auto isMu = tr->getParticleIDProbability("muonProbability");
      auto isEl = tr->getParticleIDProbability("electronProbability");
      return !(isMu || isEl);
    };
    calcInput[_trk_prefix+"isNeutralHad"] = [](const Track*){ return 0.; };
    calcInput[_trk_prefix+"type"] = [](const Track* tr){ return tr->getPDG(); };
    calcInput[_trk_prefix+"mcpid"] = [](const Track* tr){ return (tr->getMcp()) ? tr->getMcp()->getId() : 0.; };
    calcInput[_trk_prefix+"mcp_pdg"] = [](const Track* tr){ return (tr->getMcp()) ? tr->getMcp()->getPDG() : 0.; };
    calcInput[_trk_prefix+"Ktype"] = [](const Track* tr){
      auto pdg = abs(tr->getPDG());
      if (pdg==321 || pdg==310) return 1;
      return 0;
    };
    calcInput[_trk_prefix+"isPion"] = [](const Track* tr){ return tr->getParticleIDProbability("pionProbability"); };
    calcInput[_trk_prefix+"isKaon"] = [](const Track* tr){ return tr->getParticleIDProbability("kaonProbability"); };
    calcInput[_trk_prefix+"isProton"] = [](const Track* tr){ return tr->getParticleIDProbability("protonProbability"); };

    // for PI3 (20240203)
    calcInput[_trk_prefix+"proton_K"] = [](const Track* tr){ 
      auto isproton = tr->getParticleIDProbability("protonProbability");
      auto iskaon = tr->getParticleIDProbability("kaonProbability");
      return isproton-iskaon;
    };
    calcInput[_trk_prefix+"pion_K"] = [](const Track* tr){ 
      auto ispion = tr->getParticleIDProbability("pionProbability");
      auto iskaon = tr->getParticleIDProbability("kaonProbability");
      return ispion-iskaon;
    };
    calcInput[_trk_prefix+"proton_Klike"] = [](const Track* tr){ 
      auto isproton = tr->getParticleIDProbability("protonLikelihood");
      auto iskaon = tr->getParticleIDProbability("kaonLikelihood");
      return isproton-iskaon;
    };
    calcInput[_trk_prefix+"pion_Klike"] = [](const Track* tr){ 
      auto ispion = tr->getParticleIDProbability("pionLikelihood");
      auto iskaon = tr->getParticleIDProbability("kaonLikelihood");
      return ispion-iskaon;
    };
    calcInput[_trk_prefix+"dEdxEl"] = [](const Track* tr){ return tr->getParticleIDProbability("electron_dEdxdistance"); };
    calcInput[_trk_prefix+"dEdxMu"] = [](const Track* tr){ return tr->getParticleIDProbability("muon_dEdxdistance"); };
    calcInput[_trk_prefix+"dEdxPion"] = [](const Track* tr){ return tr->getParticleIDProbability("pion_dEdxdistance"); };
    calcInput[_trk_prefix+"dEdxKaon"] = [](const Track* tr){ return tr->getParticleIDProbability("kaon_dEdxdistance"); };
    calcInput[_trk_prefix+"dEdxProton"] = [](const Track* tr){ return tr->getParticleIDProbability("proton_dEdxdistance"); };
    calcInput[_trk_prefix+"iselectronlike"] = [](const Track* tr){ return tr->getParticleIDProbability("electronLikelihood"); };
    calcInput[_trk_prefix+"ismuonlike"] = [](const Track* tr){ return tr->getParticleIDProbability("muonLikelihood"); };
    calcInput[_trk_prefix+"ispionlike"] = [](const Track* tr){ return tr->getParticleIDProbability("pionLikelihood"); };
    calcInput[_trk_prefix+"iskaonlike"] = [](const Track* tr){ return tr->getParticleIDProbability("kaonLikelihood"); };
    calcInput[_trk_prefix+"isprotonlike"] = [](const Track* tr){ return tr->getParticleIDProbability("protonLikelihood"); };

    // for neutrals
    // particle kinematics
    calcInput[_neu_prefix+"px"] = [](const Neutral* neu){ return neu->Px(); };
    calcInput[_neu_prefix+"py"] = [](const Neutral* neu){ return neu->Py(); };
    calcInput[_neu_prefix+"pz"] = [](const Neutral* neu){ return neu->Pz(); };
    calcInput[_neu_prefix+"e"] = [](const Neutral* neu){ return neu->E(); };
    calcInput[_neu_prefix+"efrac"] = [](const Neutral* neu, const Jet* jet){ return neu->E() / jet->E(); };
    calcInput[_neu_prefix+"erel_log"] = [](const Neutral* neu, const Jet* jet){ return log10(neu->E() / jet->E()); };
    calcInput[_neu_prefix+"thetarel"] = [](const Neutral* neu, const Jet* jet){
      float theta, phi;
      calc_thetaphi(jet->Vect(), neu->Vect(), theta, phi);
      return theta;
    };
    calcInput[_neu_prefix+"phirel"] = [](const Neutral* neu, const Jet* jet){
      float theta, phi;
      calc_thetaphi(jet->Vect(), neu->Vect(), theta, phi);
      return phi;
    };
    calcInput[_neu_prefix+"thetarel_ilc"] = [](const Neutral* neu, const Jet* jet){ return neu->Theta() - jet->Theta(); };
    calcInput[_neu_prefix+"phirel_ilc"] = [](const Neutral* neu, const Jet* jet){
      auto ret = neu->Phi() - jet->Phi();
      if(ret < -TMath::Pi()) ret += TMath::Pi() * 2;
      if(ret > TMath::Pi()) ret -= TMath::Pi() * 2;
      return ret;
    };

    // track errors -- why are these implemented for neutrals?
    calcInput[_neu_prefix+"dptdpt"] = [](const Neutral*){ return -9; };
    calcInput[_neu_prefix+"detadeta"] = [](const Neutral*){ return -9; };
    calcInput[_neu_prefix+"dphidphi"] = [](const Neutral*){ return -9; };
    calcInput[_neu_prefix+"dxydxy"] = [](const Neutral*){ return -9; };
    calcInput[_neu_prefix+"dzdz"] = [](const Neutral*){ return -9; };
    calcInput[_neu_prefix+"dxydz"] = [](const Neutral*){ return -9; };
    calcInput[_neu_prefix+"dphidxy"] = [](const Neutral*){ return -9; };
    calcInput[_neu_prefix+"dlambdadz"] = [](const Neutral*){ return -9; };
    calcInput[_neu_prefix+"dxyc"] = [](const Neutral*){ return -9; };
    calcInput[_neu_prefix+"dxyctgtheta"] = [](const Neutral*){ return -9; };
    calcInput[_neu_prefix+"phic"] = [](const Neutral*){ return -9; };
    calcInput[_neu_prefix+"phidz"] = [](const Neutral*){ return -9; };
    calcInput[_neu_prefix+"phictgtheta"] = [](const Neutral*){ return -9; };
    calcInput[_neu_prefix+"cdz"] = [](const Neutral*){ return -9; };
    calcInput[_neu_prefix+"cctgtheta"] = [](const Neutral*){ return -9; };

    // particle displacements -- why are these implemented for neutrals?
    calcInput[_neu_prefix+"d0"] = [](const Neutral*){ return -9; };
    calcInput[_neu_prefix+"d0sig"] = [](const Neutral*){ return -9; };
    calcInput[_neu_prefix+"z0"] = [](const Neutral*){ return -9; };
    calcInput[_neu_prefix+"z0sig"] = [](const Neutral*){ return -9; };
    calcInput[_neu_prefix+"ip3d"] = [](const Neutral*){ return -9; };
    calcInput[_neu_prefix+"ip3dsig"] = [](const Neutral*){ return -9; };

    calcInput[_neu_prefix+"dxy"] = [](const Neutral*){ return -9; };
    calcInput[_neu_prefix+"dz"] = [](const Neutral*){ return -9; };
    calcInput[_neu_prefix+"btagSip2dVal"] = [](const Neutral*){ return -9; };
    calcInput[_neu_prefix+"btagSip2dSig"] = [](const Neutral*){ return -9; };
    calcInput[_neu_prefix+"btagSip3dVal"] = [](const Neutral*){ return -9; };
    calcInput[_neu_prefix+"btagSip3dSig"] = [](const Neutral*){ return -9; };
    calcInput[_neu_prefix+"btagJetDistVal"] = [](const Neutral*){ return -9; };
    calcInput[_neu_prefix+"btagJetDistSig"] = [](const Neutral*){ return -9; };

    // particle ID
    calcInput[_neu_prefix+"charge"] = [](const Neutral*){ return 0; };
    calcInput[_neu_prefix+"isMu"] = [](const Neutral*){ return 0; };
    calcInput[_neu_prefix+"isEl"] = [](const Neutral*){ return 0; };
    calcInput[_neu_prefix+"isGamma"] = [](const Neutral* neu){
      // simple photon finder
      double ecaldep = neu->getCaloEdep()[tpar::ecal];
      double hcaldep = neu->getCaloEdep()[tpar::hcal];
      return (ecaldep / (ecaldep + hcaldep) > 0.98);
    };
    calcInput[_neu_prefix+"isChargedHad"] = [](const Neutral*){ return 0; };
    calcInput[_neu_prefix+"isNeutralHad"] = [](const Neutral* neu){ 
      // simple photon finder
      double ecaldep = neu->getCaloEdep()[tpar::ecal];
      double hcaldep = neu->getCaloEdep()[tpar::hcal];
      return !(ecaldep / (ecaldep + hcaldep) > 0.98);
    };
    calcInput[_neu_prefix+"type"] = [](const Neutral* neu){ return neu->getPDG(); };
    calcInput[_neu_prefix+"mcpid"] = [](const Neutral* neu){ return (neu->getMcp()) ? neu->getMcp()->getId() : 0.; };
    calcInput[_neu_prefix+"mcp_pdg"] = [](const Neutral* neu){ return (neu->getMcp()) ? neu->getMcp()->getPDG() : 0.; };

    calcInput[_neu_prefix+"Ktype"] = [](const Neutral* neu){
      auto pdg = abs(neu->getPDG());
      return ( (pdg==321 || pdg==310) );
    };
    calcInput[_neu_prefix+"isPion"] = [](const Neutral*){ return 0; };
    calcInput[_neu_prefix+"isKaon"] = [](const Neutral*){ return 0; };
    calcInput[_neu_prefix+"isProton"] = [](const Neutral*){ return 0; };

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