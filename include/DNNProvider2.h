// DNNProvider2.h

#ifndef dnnprovider2_h
#define dnnprovider2_h 1

#include "lcfiplus.h"
#include "TTree.h"
#include "TFile.h"
#include "TVector3.h"

#include <vector>

namespace lcfiplus {

class DNNProvider2 : public Algorithm {
 public:
  DNNProvider2() {}
  virtual ~DNNProvider2() {}

  void init(Parameters* param);
  void process();
  void end();
  ClassDef(DNNProvider2,1);
 private:
  TTree* _ntp;
  TFile* _file;
  JetVec* _jets;  //!
  string _jetname;
  int _mcIsB;
  int _mcIsC;
  int _mcIsU;
  int _mcIsD;
  int _mcIsS;
  int _mcIsG;

  // FCC functions
  float calc_dxy(float, float, float, TVector3, TVector3, int);
  float calc_dz(float, float, float, TVector3, TVector3, int);
  float calc_sip2d(float, float, float, float);
  float calc_sip3d(float, float, float, TVector3);
  float calc_jetDist(float, float, float, TVector3, TVector3);
  void calc_thetaphi(TVector3, TVector3, float&, float&);

  struct DNNData {
    float jet_px;
    float jet_py;
    float jet_pz;
    float jet_e;
    float jet_mass;
    int jet_ntracks;
    int jet_nneutrals;

    float jet_theta;
    float jet_phi;

    // for charged hadrons
    std::vector<float> px;
    std::vector<float> py;
    std::vector<float> pz;
    std::vector<float> e;
    std::vector<float> efrac;
    std::vector<float> erel_log;
    std::vector<float> dtheta;
    std::vector<float> dphi;
    std::vector<float> dtheta_ilc;
    std::vector<float> dphi_ilc;

    std::vector<float> dEdx;

    // covariant matrix of tracks
    std::vector<float> cov_d0;
    std::vector<float> cov_z0;
    std::vector<float> cov_phi;
    std::vector<float> cov_omega;
    std::vector<float> cov_tanlambda;

    std::vector<float> cov_d0_z0;
    std::vector<float> cov_d0_phi;
    std::vector<float> cov_d0_omega;
    std::vector<float> cov_d0_tanlambda;

    std::vector<float> cov_z0_phi;
    std::vector<float> cov_z0_omega;
    std::vector<float> cov_z0_tanlambda;

    std::vector<float> cov_phi_omega;
    std::vector<float> cov_phi_tanlambda;
    std::vector<float> cov_omega_tanlambda;

    std::vector<float> d0;
    std::vector<float> d0sig;
    std::vector<float> z0;
    std::vector<float> z0sig;
    std::vector<float> ip3d;
    std::vector<float> ip3dsig;

    std::vector<float> dxy;
    std::vector<float> dz;
    std::vector<float> ip2d_fcc;
    std::vector<float> ip2dsig_fcc;
    std::vector<float> ip3d_fcc;
    std::vector<float> ip3dsig_fcc;
    std::vector<float> jetdist_fcc;
    std::vector<float> jetdistsig_fcc;

    std::vector<int> charge;
    std::vector<float> ismuon;
    std::vector<float> iselectron;
    std::vector<float> isphoton;
    std::vector<float> ispion;
    std::vector<float> iskaon;
    std::vector<float> isproton;
    std::vector<float> ischargedhadron;
    std::vector<float> isneutralhadron;
    std::vector<int> iskaon0;
    std::vector<int> pdg_pfa;
    std::vector<int> mcpid;
    std::vector<int> mcp_pdg;
    std::vector<int> K_pdg_pfa; //20240203

    // for PI3 (20240203)
    std::vector<float> proton_K;
    std::vector<float> pion_K;
    std::vector<float> proton_Klike;
    std::vector<float> pion_Klike;
    std::vector<float> electron_dEdxdistance;
    std::vector<float> muon_dEdxdistance;
    std::vector<float> kaon_dEdxdistance;
    std::vector<float> pion_dEdxdistance;
    std::vector<float> proton_dEdxdistance;

    std::vector<float> ismuonlike;
    std::vector<float> iselectronlike;
    std::vector<float> ispionlike;
    std::vector<float> iskaonlike;
    std::vector<float> isprotonlike;


    int mc_b;
    int mc_c;
    int mc_u;
    int mc_d;
    int mc_s;
    int mc_g;
    int mc_q;

    // for neutral hadrons
    std::vector<float> neu_px;
    std::vector<float> neu_py;
    std::vector<float> neu_pz;
    std::vector<float> neu_e;
    std::vector<float> neu_efrac;
    std::vector<float> neu_erel_log;
    std::vector<float> neu_dtheta;
    std::vector<float> neu_dphi;
    std::vector<float> neu_dtheta_ilc;
    std::vector<float> neu_dphi_ilc;

    // covariant matrix of tracks
    std::vector<float> neu_cov_d0;
    std::vector<float> neu_cov_z0;
    std::vector<float> neu_cov_phi;
    std::vector<float> neu_cov_omega;
    std::vector<float> neu_cov_tanlambda;

    std::vector<float> neu_cov_d0_z0;
    std::vector<float> neu_cov_d0_phi;
    std::vector<float> neu_cov_d0_omega;
    std::vector<float> neu_cov_d0_tanlambda;

    std::vector<float> neu_cov_z0_phi;
    std::vector<float> neu_cov_z0_omega;
    std::vector<float> neu_cov_z0_tanlambda;

    std::vector<float> neu_cov_phi_omega;
    std::vector<float> neu_cov_phi_tanlambda;
    std::vector<float> neu_cov_omega_tanlambda;

    std::vector<float> neu_d0;
    std::vector<float> neu_d0sig;
    std::vector<float> neu_z0;
    std::vector<float> neu_z0sig;
    std::vector<float> neu_ip3d;
    std::vector<float> neu_ip3dsig;

    std::vector<float> neu_dxy;
    std::vector<float> neu_dz;
    std::vector<float> neu_ip2d_fcc;
    std::vector<float> neu_ip2dsig_fcc;
    std::vector<float> neu_ip3d_fcc;
    std::vector<float> neu_ip3dsig_fcc;
    std::vector<float> neu_jetdist_fcc;
    std::vector<float> neu_jetdistsig_fcc;

    std::vector<int> neu_charge;
    std::vector<float> neu_ismuon;
    std::vector<float> neu_iselectron;
    std::vector<float> neu_isphoton;
    std::vector<float> neu_ispion;
    std::vector<float> neu_iskaon;
    std::vector<float> neu_isproton;
    std::vector<float> neu_ischargedhadron;
    std::vector<float> neu_isneutralhadron;
    std::vector<int> neu_iskaon0;
    std::vector<int> neu_pdg_pfa;
    std::vector<int> neu_mcpid;
    std::vector<int> neu_mcp_pdg;
    std::vector<int> neu_K_pdg_pfa; //20240201




  };
  DNNData _data;
};

}

#endif
