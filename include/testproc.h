// testproc.h

#ifndef testproc_h
#define testproc_h 1

#include "lcfiplus.h"
#include "TNtuple.h"
#include "TNtupleD.h"
#include "TH2D.h"

namespace lcfiplus {

class ZHHAlgo : public Algorithm {
 public:
  struct data {
    int mchdecaypdg[2];
    int mchbb;
    int mcnb;

    double thrust;
    double thaxis[3];
    double ycuts[10];
    int ntr;
    int npfo;

    double mass[15];
    double ntrjetmin;
    double pmiss[3];
    double emiss;

    double bcat[6];
    double btag[6];
    double ctag[6];
    double ejet[6];
    double pxjet[6];
    double pyjet[6];
    double pzjet[6];
    double ntrjet[6];
    double twovtxprobjet[6];
    double vtxangle[6];
    double mcnb6[6];
    double mcnc6[6];

    double bcat4[4];
    double btag4[4];
    double ctag4[4];
    double ejet4[4];
    double pxjet4[4];
    double pyjet4[4];
    double pzjet4[4];
    double ntrjet4[4];
    double twovtxprobjet4[4];
    double vtxangle4[4];

    double bcat5[5];
    double btag5[5];
    double ctag5[5];
    double ejet5[5];
    double pxjet5[5];
    double pyjet5[5];
    double pzjet5[5];
    double ntrjet5[5];
    double twovtxprobjet5[5];
    double vtxangle5[5];

    double bcat7[7];
    double btag7[7];
    double ctag7[7];
    double ejet7[7];
    double pxjet7[7];
    double pyjet7[7];
    double pzjet7[7];
    double ntrjet7[7];
    double twovtxprobjet7[7];
    double vtxangle7[7];

    double bcat8[8];
    double btag8[8];
    double ctag8[8];
    double ejet8[8];
    double pxjet8[8];
    double pyjet8[8];
    double pzjet8[8];
    double ntrjet8[8];
    double twovtxprobjet8[8];
    double vtxangle8[8];

    double bcatnv4[4];
    double btagnv4[4];
    double ctagnv4[4];
    double ejetnv4[4];
    double pxjetnv4[4];
    double pyjetnv4[4];
    double pzjetnv4[4];
    double ntrjetnv4[4];
    double twovtxprobjetnv4[4];
    double vtxanglenv4[4];

    double bcatnv5[5];
    double btagnv5[5];
    double ctagnv5[5];
    double ejetnv5[5];
    double pxjetnv5[5];
    double pyjetnv5[5];
    double pzjetnv5[5];
    double ntrjetnv5[5];
    double twovtxprobjetnv5[5];
    double vtxanglenv5[5];

    double bcatnv6[6];
    double btagnv6[6];
    double ctagnv6[6];
    double ejetnv6[6];
    double pxjetnv6[6];
    double pyjetnv6[6];
    double pzjetnv6[6];
    double ntrjetnv6[6];
    double twovtxprobjetnv6[6];
    double vtxanglenv6[6];

    double bcatnv7[7];
    double btagnv7[7];
    double ctagnv7[7];
    double ejetnv7[7];
    double pxjetnv7[7];
    double pyjetnv7[7];
    double pzjetnv7[7];
    double ntrjetnv7[7];
    double twovtxprobjetnv7[7];
    double vtxanglenv7[7];

    double bcatnv8[8];
    double btagnv8[8];
    double ctagnv8[8];
    double ejetnv8[8];
    double pxjetnv8[8];
    double pyjetnv8[8];
    double pzjetnv8[8];
    double ntrjetnv8[8];
    double twovtxprobjetnv8[8];
    double vtxanglenv8[8];

  };

  ZHHAlgo() {}
  virtual ~ZHHAlgo() {}

  static bool sortBtag(const Jet* j1, const Jet* j2) {
    double btag1 = j1->getParam("lcfiplus")->get<double>("BTag");
    double btag2 = j2->getParam("lcfiplus")->get<double>("BTag");

    return btag1 > btag2;
  }

  void init(Parameters* param);
  void process();
  void end();

  ClassDef(ZHHAlgo,1);

 private:
  TFile* _file;

  string _jetname;
  string _jetname4;
  string _jetname5;
  string _jetname7;
  string _jetname8;

  string _jetnamenv4;
  string _jetnamenv5;
  string _jetnamenv6;
  string _jetnamenv7;
  string _jetnamenv8;

  JetVec* _jets;  //!
  JetVec* _jets4;  //!
  JetVec* _jets5;  //!
  JetVec* _jets7;  //!
  JetVec* _jets8;  //!

  JetVec* _jetsnv4;  //!
  JetVec* _jetsnv5;  //!
  JetVec* _jetsnv6;  //!
  JetVec* _jetsnv7;  //!
  JetVec* _jetsnv8;  //!

  TTree* _tree;

  data _d;

};

class TestAlgo : public Algorithm {
 public:
  TestAlgo() {}
  virtual ~TestAlgo() {}

  void init(Parameters* param);
  void process();
  void end();

  ClassDef(TestAlgo,1);

 private:
  TNtupleD* _nt;
  int _nev;

  TNtuple* _ntJet2;
  TNtuple* _nbJet;
  TFile* _file;

  string _v0vtxname;
  string _privtxname;
  string _jetname;
  bool _bbhh;

  VertexVec* _vertices;  //!
  VertexVec* _v0vertices;  //!
  JetVec* _jets;  //!

  // for old version

  string _vtxname;
  int _vtxsel;
  int _refine;
  TNtupleD* _ntJet;

  TH2D* _h;
  TH2D* _he;

};

class VertexAnalysis : public Algorithm {
 public:
  VertexAnalysis() {}
  virtual ~VertexAnalysis() {}

  void init(Parameters* param);
  void process();
  void end();

  ClassDef(VertexAnalysis,1);

 private:
  TNtupleD* _nt;
  int _nev;

  TFile* _file;

  string _privtxname;
  string _secvtxname;
};

class FlavtagReader : public Algorithm {
 public:
  FlavtagReader() {}
  virtual ~FlavtagReader() {}

  void init(Parameters* param);
  void process();
  void end();

  ClassDef(FlavtagReader,1);

 private:
  TNtupleD* _nt;
  TNtupleD* _ntev;
  int _nev;

  TNtuple* _ntJet2;
  TNtuple* _nbJet;
  TFile* _file;

  string _v0vtxname;
  string _privtxname;
  string _jetname;
  bool _bbhh;

  VertexVec* _vertices;  //!
  VertexVec* _v0vertices;  //!
  JetVec* _jets;  //!

  // for old version

  string _vtxname;
  int _vtxsel;
  int _refine;
  TNtupleD* _ntJet;

  TH2D* _h;
  TH2D* _he;

};

class TestAlgoV0 : public Algorithm {
 public:
  TestAlgoV0() {}
  virtual ~TestAlgoV0() {}
  void init(Parameters* param);
  void process();
  void end();
  ClassDef(TestAlgoV0,1);
 private:
  TTree* _ntp;
  TFile* _file;
  VertexVec* _vertices;  //!
  string _vtxname;

  struct VtxData {
    double x;
    double y;
    double z;
    double r;
    double cs;
    double phi;
    double chrg;
    double dirdot;
    double dirdot2;
    int ntrk;
    double mks;
    double ml0;
    double mconv;
    double mks2;
    double ml02;
    int v0;
    int ks;
    int l0;
    int conv;
    int mcpdg1;
    int mcpdg2;
    int mcppdg1;
    int mcppdg2;
    int mcpp1;
    int mcpp2;
  };
  VtxData _data;
};

}

#endif
