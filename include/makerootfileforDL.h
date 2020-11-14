// testproc.h

#ifndef makerootfileforDL_h
#define makerootfileforDL_h 1

#include "lcfiplus.h"
#include "LcfiplusProcessor.h"
#include "geometry.h"
#include "TNtuple.h"
#include "TNtupleD.h"
#include "TH2D.h"
#include <EVENT/Track.h>

namespace lcfiplus {

class MakeROOTFileCC : public Algorithm {
 public:
  MakeROOTFileCC() {}
  virtual ~MakeROOTFileCC() {}
  void init(Parameters* param);
  void process();
  void end();
  ClassDef(MakeROOTFileCC,1);
 private:
  TTree* _ntp;
  TFile* _file;
  Helix* _hel1;
  Helix* _hel2;
  int nEvt;
  int ntr1Trk;
  int ntr2Trk;
  int tr1sscid = -999;
  int tr1sscpdg = -999;
  int tr1ssc = -999;
  int tr2sscid = -999;
  int tr2sscpdg = -999;
  int tr2ssc = -999;

  TVector3* xyz;
  //enum par : int;
  //enum tpar : int;
  struct TracksData {
    double tr1d0 = -999;
    double tr1z0 = -999;
    double tr1phi = -999;
    double tr1omega = -999;
    double tr1tanlam = -999;
    double tr1x = -999;
    double tr1y = -999;
    double tr1z = -999;
    int tr1charge = -999;
    double tr1energy = -999;
    double tr1covmatrixd0d0 = -999;
    double tr1covmatrixd0z0 = -999;
    double tr1covmatrixd0ph = -999;
    double tr1covmatrixd0om = -999;
    double tr1covmatrixd0tl = -999;
    double tr1covmatrixz0z0 = -999;
    double tr1covmatrixz0ph = -999;
    double tr1covmatrixz0om = -999;
    double tr1covmatrixz0tl = -999;
    double tr1covmatrixphph = -999;
    double tr1covmatrixphom = -999;
    double tr1covmatrixphtl = -999;
    double tr1covmatrixomom = -999;
    double tr1covmatrixomtl = -999;
    double tr1covmatrixtltl = -999;
    double tr2d0 = -999;
    double tr2z0 = -999;
    double tr2phi = -999;
    double tr2omega = -999;
    double tr2tanlam = -999;
    double tr2x = -999;
    double tr2y = -999;
    double tr2z = -999;
    int tr2charge = -999;
    double tr2energy = -999;
    double tr2covmatrixd0d0 = -999;
    double tr2covmatrixd0z0 = -999;
    double tr2covmatrixd0ph = -999;
    double tr2covmatrixd0om = -999;
    double tr2covmatrixd0tl = -999;
    double tr2covmatrixz0z0 = -999;
    double tr2covmatrixz0ph = -999;
    double tr2covmatrixz0om = -999;
    double tr2covmatrixz0tl = -999;
    double tr2covmatrixphph = -999;
    double tr2covmatrixphom = -999;
    double tr2covmatrixphtl = -999;
    double tr2covmatrixomom = -999;
    double tr2covmatrixomtl = -999;
    double tr2covmatrixtltl = -999;
    double tr1mcx = -999;
    double tr1mcy = -999;
    double tr1mcz = -999;
    int tr1id = -999;
    int tr1pdg = -999;
    double tr2mcx = -999;
    double tr2mcy = -999;
    double tr2mcz = -999;
    int tr2id = -999;
    int tr2pdg = -999;
    int cosemistable = -999;
    int cosemistablec = -999;
    int connect = -999;
    double tr1tlvx = -999;
    double tr1tlvy = -999;
    double tr1tlvz = -999;
    double tr2tlvx = -999;
    double tr2tlvy = -999;
    double tr2tlvz = -999;
    double mass = -999;
    double l0mass = -999;
    double minE = -999;
    double chi2 = -999;
    double vchi2 = -999;
    double vposx = -999;
    double vposy = -999;
    double vposz = -999;
    double mag = -999;
    double vec = -999;
    int tr1selection = -999;
    int tr2selection = -999;
    int v0selection = -999;
    int lcfiplustag = -999;
  };
  TracksData _data;
};

class MakeROOTFileBB : public Algorithm {
 public:
  MakeROOTFileBB() {}
  virtual ~MakeROOTFileBB() {}
  void init(Parameters* param);
  void process();
  void end();
  ClassDef(MakeROOTFileBB,1);
 private:
  TTree* _ntp;
  TFile* _file;
  Helix* _hel1;
  Helix* _hel2;
  int nEvt;
  int ntr1Trk;
  int ntr2Trk;
  int tr1ssid = -999;
  int tr1sspdg = -999;
  int tr1ssc = -999;
  int tr1ssb = -999;
  int tr1oth = -999;
  int tr1pri = -999;
  int tr2ssid = -999;
  int tr2sspdg = -999;
  int tr2ssc = -999;
  int tr2ssb = -999;
  int tr2oth = -999;
  int tr2pri = -999;

  TVector3* xyz;
  //enum par : int;
  //enum tpar : int;
  struct TracksData {
    double tr1d0 = -999;
    double tr1z0 = -999;
    double tr1phi = -999;
    double tr1omega = -999;
    double tr1tanlam = -999;
    double tr1x = -999;
    double tr1y = -999;
    double tr1z = -999;
    int tr1charge = -999;
    double tr1energy = -999;
    double tr1covmatrixd0d0 = -999;
    double tr1covmatrixd0z0 = -999;
    double tr1covmatrixd0ph = -999;
    double tr1covmatrixd0om = -999;
    double tr1covmatrixd0tl = -999;
    double tr1covmatrixz0z0 = -999;
    double tr1covmatrixz0ph = -999;
    double tr1covmatrixz0om = -999;
    double tr1covmatrixz0tl = -999;
    double tr1covmatrixphph = -999;
    double tr1covmatrixphom = -999;
    double tr1covmatrixphtl = -999;
    double tr1covmatrixomom = -999;
    double tr1covmatrixomtl = -999;
    double tr1covmatrixtltl = -999;
    double tr2d0 = -999;
    double tr2z0 = -999;
    double tr2phi = -999;
    double tr2omega = -999;
    double tr2tanlam = -999;
    double tr2x = -999;
    double tr2y = -999;
    double tr2z = -999;
    int tr2charge = -999;
    double tr2energy = -999;
    double tr2covmatrixd0d0 = -999;
    double tr2covmatrixd0z0 = -999;
    double tr2covmatrixd0ph = -999;
    double tr2covmatrixd0om = -999;
    double tr2covmatrixd0tl = -999;
    double tr2covmatrixz0z0 = -999;
    double tr2covmatrixz0ph = -999;
    double tr2covmatrixz0om = -999;
    double tr2covmatrixz0tl = -999;
    double tr2covmatrixphph = -999;
    double tr2covmatrixphom = -999;
    double tr2covmatrixphtl = -999;
    double tr2covmatrixomom = -999;
    double tr2covmatrixomtl = -999;
    double tr2covmatrixtltl = -999;
    double tr1mcx = -999;
    double tr1mcy = -999;
    double tr1mcz = -999;
    int tr1id = -999;
    int tr1pdg = -999;
    double tr2mcx = -999;
    double tr2mcy = -999;
    double tr2mcz = -999;
    int tr2id = -999;
    int tr2pdg = -999;
    int cosemistable = -999;
    int cosemistablec = -999;
    int connect = -999;
    double tr1tlvx = -999;
    double tr1tlvy = -999;
    double tr1tlvz = -999;
    double tr2tlvx = -999;
    double tr2tlvy = -999;
    double tr2tlvz = -999;
    double mass = -999;
    double l0mass = -999;
    double minE = -999;
    double chi2 = -999;
    double vchi2 = -999;
    double vposx = -999;
    double vposy = -999;
    double vposz = -999;
    double mag = -999;
    double vec = -999;
    int tr1selection = -999;
    int tr2selection = -999;
    int v0selection = -999;
    int lcfiplustag = -999;
  };
  TracksData _data;
};

}

#endif
