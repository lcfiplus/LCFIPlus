// MLInputGenerator.h

#ifndef mlinputgenerator_h
#define mlinputgenerator_h 1

#include "lcfiplus.h"
#include "TTree.h"
#include "TFile.h"
#include "TVector3.h"

#include <vector>

namespace lcfiplus {
namespace MLInputGenerator {
  static map<string, function<double(const Jet*)> > varJet;
  static map<string, function<double(const Track*)> > varTrack;
  static map<string, function<double(const Track*, const Jet*)> > varTrackJet;
  static map<string, function<double(const Track*, const Vertex*)> > varTrackVertex;
  static map<string, function<double(const Neutral*)> > varNeutral;
  static map<string, function<double(const Neutral*, const Jet*)> > varNeutralJet;
  static map<string, function<double(const Neutral*, const Vertex*)> > varNeutralVertex;
  
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
