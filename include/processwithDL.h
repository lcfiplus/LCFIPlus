// process.h

#ifndef processwithDL_h
#define processwithDL_h 1

#include "lcfiplus.h"
#include "VertexFinderSuehara.h"
#include "VertexFinderwithDL.h"

#include "TMVA/Reader.h"

#include <tensorflow/cc/saved_model/loader.h>
#include <tensorflow/cc/saved_model/tag_constants.h>


namespace lcfiplus {
struct TrackSelectorConfig;

class VertexFinderLSTM : public Algorithm {
 public:
  VertexFinderLSTM() {}
  virtual ~VertexFinderLSTM() {}

  void init(Parameters* param);
  void process();
  void end();

  ClassDef(VertexFinderLSTM,1);

 private:
  std::vector<Vertex*>* _primary_vertex;	//!
  std::vector<Vertex*>* _secondary_vertices;	//!

  // parameters
  tensorflow::ConfigProto config;

  int NTrackVariable;
  int MaxSample, MaxEvent, MaxTrack, MaxNpyVariable, NPairVariable, NTrackVariable;
  int MaxPrimaryVertexLoop;
  double ThresholdPairSecondaryScoreBBCC, ThresholdPairSecondaryScore, ThresholdPairPosScore;
  double ThresholdPrimaryScore, ThresholdSecondaryScore;
  bool debug;

  tensorflow::string pair_path, pair_pos_path, lstm_path, slstm_path;
  tensorflow::SessionOptions session_options;
  tensorflow::RunOptions run_options;
  tensorflow::SavedModelBundleLite pair_model_bundle, pair_pos_model_bundle, lstm_model_bundle, slstm_model_bundle;

  tensorflow::Status pair_status, pair_pos_status lstm_status, slstm_status;


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

  //AVF parameters
  bool _avf; //flag AVF/chi2
  double _temperature;  //AVF parameter

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
  double _gammaParameter;
  bool _useMuonID;
  bool _muonIDExternal;
  double _muonIDMinEnergy;
  double _muonIDMinD0Sig;
  double _muonIDMinZ0Sig;
  double _muonIDMaxDist;
  double _muonIDMinProb;

  double _vsMinDist;
  double _vsMaxDist;
  double _vsK0MassWidth;
  bool _outputVertexStoresVertex;
  string _vpricolname;
  string _vseccolname;
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

  //BNess tagger
  TMVA::Reader *_bness=nullptr;
  float _var[8];
  string _bnessbookname;
  string _bnessbookname1;
  
};

}

#endif
