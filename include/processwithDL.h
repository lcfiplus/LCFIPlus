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

class VertexFindingwithDL : public Algorithm {
 public:
  VertexFindingwithDL() {}
  virtual ~VertexFindingwithDL() {}

  void init(Parameters* param);
  void process();
  void end();

  ClassDef(VertexFindingwithDL,1);

 private:
  std::vector<Vertex*>* _primary_vertex;	//!
  std::vector<Vertex*>* _secondary_vertices;	//!

  // parameters
  tensorflow::ConfigProto config;

  int NEventNumber;
  int NTrackVariable;
  int MaxTrack;
  int MaxPrimaryVertexLoop;
  double ThresholdPairSecondaryScoreBBCC, ThresholdPairSecondaryScore, ThresholdPairPosScore;
  double ThresholdPrimaryScore, ThresholdSecondaryScore;
  bool debug;

  tensorflow::string pair_path, pair_pos_path, lstm_path, slstm_path;
  tensorflow::SessionOptions session_options;
  tensorflow::RunOptions run_options;
  tensorflow::SavedModelBundleLite pair_model_bundle, pair_pos_model_bundle, lstm_model_bundle, slstm_model_bundle;

  //tensorflow::Status pair_status, pair_pos_status, lstm_status, slstm_status;


};

}

#endif
