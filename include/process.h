// process.h

#ifndef process_h
#define process_h 1

#include "lcfiplus.h"
#include "VertexFinderSuehara.h"

namespace lcfiplus {
struct TrackSelectorConfig;

class PrimaryVertexFinder : public Algorithm {
 public:
  PrimaryVertexFinder() {}
  virtual ~PrimaryVertexFinder() {}

  void init(Parameters* param);
  void process();
  void end();

  ClassDef(PrimaryVertexFinder,1);

 private:
  vector<Vertex*>* _vertex;	//!

  // parameters
  double _chi2th;
  bool _beamspotConstraint;
  bool _beamspotSmearing;

  // track cut parameters
  TrackSelectorConfig* _secVtxCfg; //!
};

class BuildUpVertex : public Algorithm {
 public:
  BuildUpVertex() : _vertices(0) {}
  virtual ~BuildUpVertex() {}

  void init(Parameters* param);
  void process();
  void end();

  ClassDef(BuildUpVertex,1);

 private:
  std::vector<Vertex*>* _vertices;	//!
  std::vector<Vertex*>* _v0vertices;	//!

  // parameters
  std::string _primvtxcolname;

  // vertex formation limits
  double _chi2thpri;
  double _chi2thsec;
  double _massth;
  double _posth;
  double _chi2orderinglimit;
  int _v0sel;

  // association parameters
  bool _doassoc;
  double _minimumdist;
  double _chi2ratio;

  // track cut parameters
  TrackSelectorConfig* _secVtxCfg; //!
};

class JetClustering : public Algorithm {
 public:
  JetClustering() {}
  virtual ~JetClustering() {}

  void init(Parameters* param);
  void process();
  void end();

  ClassDef(JetClustering,1);

 private:
  VertexVec* _vertices;  //!
  map<double, vector<Jet*> * > _jetsmap;  //!
  map<double, vector<Vertex*> * > _jetvtxmap;  //!
  vector<int> _njets;
  vector<double> _ycut;
  string _algo;
  int _useBeamJets;
  double _rParameter;
  double _alphaParameter;
  double _betaParameter;
  bool _useMuonID;
  bool _muonIDExternal;
  double _muonIDMinD0Sig;
  double _muonIDMinZ0Sig;
  double _muonIDMaxDist;
  double _muonIDMinProb;

  double _vsMinDist;
  double _vsMaxDist;
  double _vsK0MassWidth;
  bool _outputVertexStoresVertex;
  string _vcolname;
  int _maxYth;

  double _yaddVV;
  double _yaddVL;
  double _yaddLL;
};

class JetVertexRefiner : public Algorithm {
 public:
  JetVertexRefiner() {}
  virtual ~JetVertexRefiner() {}

  void init(Parameters* param);
  void process();
  void end();

  ClassDef(JetVertexRefiner,1);

 private:
  VertexVec* _invertices;  //!
  VertexVec* _v0vertices;  //!
  vector<Vertex*>* _outvertices;   //!

  JetVec* _inputJets;	//!
  vector<Jet*>* _outputJets;   //!

  VertexFinderSuehara::VertexFinderSueharaConfig _cfg;
  double _oneVtxProbTh;
  double _maxCharmFlightLengthPerJetEnergy;

  string _jincolname;
  string _vincolname;
  string _vv0colname;
  string _vprimcolname;
};

}

#endif
