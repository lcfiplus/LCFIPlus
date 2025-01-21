// MLInputGenerator.h
#ifndef mlinputgenerator_h
#define mlinputgenerator_h 1

#include "lcfiplus.h"
#include "TTree.h"
#include "TFile.h"
#include "TVector3.h"

#include <vector>
#include <variant>

namespace lcfiplus {
namespace MLInputGenerator {

  extern map<string, variant<
    function<double(const Jet*)>,
    function<double(const Track*)>,
    function<double(const Track*, const Jet*)>,
    function<double(const Track*, const Vertex*)>,
    function<double(const Neutral*)>,
    function<double(const Neutral*, const Jet*)>,
    function<double(const Neutral*, const Vertex*)>
  > > calcInput;

  extern bool _initialized;

  void init();

  // FCC functions
  float calc_dxy(float, float, float, TVector3, TVector3, int);
  float calc_dz(float, float, float, TVector3, TVector3, int);
  float calc_sip2d(float, float, float, float);
  float calc_sip3d(float, float, float, TVector3);
  float calc_jetDist(float, float, float, TVector3, TVector3);
  void calc_thetaphi(TVector3, TVector3, float&, float&);  
}
}

#endif
