// DNNProvider2.h

#ifndef dnnprovider2_h
#define dnnprovider2_h 1

#include "lcfiplus.h"
#include "TTree.h"
#include "TFile.h"

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
  int _mcIsQ;

  struct DNNData {
    float jet_px;
    float jet_py;
    float jet_pz;
    float jet_e;
    float jet_mass;
    int jet_ntracks;
    int jet_nneutrals;

    std::vector<float> px;
    std::vector<float> py;
    std::vector<float> pz;
    std::vector<float> e;
    std::vector<float> efrac;
    std::vector<float> dtheta;
    std::vector<float> dphi;

    std::vector<float> d0;
    std::vector<float> d0sig;
    std::vector<float> z0;
    std::vector<float> z0sig;
    std::vector<float> ip3d;
    std::vector<float> ip3dsig;

    std::vector<int> charge;
    std::vector<int> ismuon;
    std::vector<int> iselectron;
    std::vector<int> isphoton;
    std::vector<int> ischargedhadron;
    std::vector<int> isneutralhadron;

    int mc_b;
    int mc_c;
    int mc_q;
  };
  DNNData _data;
};

}

#endif
