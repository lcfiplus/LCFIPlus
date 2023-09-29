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
    std::vector<int> ismuon;
    std::vector<int> iselectron;
    std::vector<int> isphoton;
    std::vector<int> ischargedhadron;
    std::vector<int> isneutralhadron;
    std::vector<int> pdg_pfa;

    int mc_b;
    int mc_c;
    int mc_u;
    int mc_d;
    int mc_s;
    int mc_g;
    int mc_q;
  };
  DNNData _data;
};

}

#endif
