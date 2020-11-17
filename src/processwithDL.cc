#include <assert.h>
#include <string>
#include <memory>

#include "EventStore.h"
#include "TrackSelector.h"
#include "VertexSelector.h"
#include "JetFinder.h"
#include "VertexFitterLCFI.h"
#include "VertexFinderTearDown.h"
#include "VertexFinderSuehara.h"
#include "VertexFinderwithDL.h"

#include "TROOT.h"
#include "TInterpreter.h"
#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TCut.h"
#include "TRandom.h"
#include "TNtuple.h"

#include "process.h"
#include "processwithDL.h"
#include "geometry.h"

#include <tensorflow/cc/saved_model/loader.h>
#include <tensorflow/cc/saved_model/tag_constants.h>


using namespace lcfiplus;
using namespace lcfiplus::algoEtc;

namespace lcfiplus {

void VertexFindingwithDL::init(Parameters* param) {
  Algorithm::init(param);

  config.mutable_gpu_options()->set_visible_device_list("1");
  config.mutable_gpu_options()->set_allow_growth(true);

  _primary_vertex = 0;

  std::string vpricolname = param->get("VertexFindingwithDL.PrimaryVertexCollectionName",string("PrimaryVertex"));
  std::string vseccolname = param->get("VertexFindingwithDL.SecondaryVerticesCollectionName",string("BuildUpVertex"));
  Event::Instance()->Register(vpricolname.c_str(), _primary_vertex, EventStore::PERSIST);
  Event::Instance()->Register(vseccolname.c_str(), _secondary_vertices, EventStore::PERSIST);

  // default setting
  Event::Instance()->setDefaultPrimaryVertex(vpricolname.c_str());

  NEventNumber = 0;
  MaxTrack = param->get("VertexFindingwithDL.MaxTrack", int(53));
  NTrackVariable = param->get("VertexFindingwithDL.NTrackVariable", int(22));
  MaxPrimaryVertexLoop = param->get("VertexFindingwithDL.MaxPrimaryVertexLoop", int(3));

  ThresholdPairSecondaryScoreBBCC = param->get("VertexFindingwithDL.ThresholdPairSecondaryScoreBBCC", double(0.6));
  ThresholdPairSecondaryScore = param->get("VertexFindingwithDL.ThresholdPairSecondaryScore", double(0.8));
  ThresholdPairPosScore = param->get("VertexFindingwithDL.ThresholdPairPosScore", double(5));

  ThresholdPrimaryScore = param->get("VertexFindingwithDL.ThresholdPrimaryScore", double(0.5));
  ThresholdSecondaryScore = param->get("VertexFindingwithDL.ThresholdSecondaryScore", double(0.8));

  debug = param->get("VertexFindingwithDL.Debug", bool(0));

  pair_path = param->get("VertexFindingwithDL.PairModelPath", string("/home/goto/ILC/Deep_Learning/model/Pair_Model_vfdnn04_1Msamples_2500epochs"));
  pair_pos_path = param->get("VertexFindingwithDL.PairPosModelPath", string("/home/goto/ILC/Deep_Learning/model/Pair_Pos_Model_vfdnn04_1Msamples_2500epochs"));
  lstm_path = param->get("VertexFindingwithDL.LSTMModelPath", string("/home/goto/ILC/Deep_Learning/model/Attention_Bidirectional_VLSTM_Model_vfdnn06_50000samples_100epochs"));
  slstm_path = param->get("VertexFindingwithDL.SLSTMModelPath", string("/home/goto/ILC/Deep_Learning/model/Attention_Bidirectional_VLSTM_Model_vfdnn06_50000samples_100epochs_ps_100epochs_s"));

  session_options = tensorflow::SessionOptions();
  run_options = tensorflow::RunOptions();

  tensorflow::Status pair_status = LoadSavedModel(session_options, run_options, pair_path, {tensorflow::kSavedModelTagServe}, &pair_model_bundle);
  tensorflow::Status pair_pos_status = LoadSavedModel(session_options, run_options, pair_pos_path, {tensorflow::kSavedModelTagServe}, &pair_pos_model_bundle);
  tensorflow::Status lstm_status = LoadSavedModel(session_options, run_options, lstm_path, {tensorflow::kSavedModelTagServe}, &lstm_model_bundle);
  tensorflow::Status slstm_status = LoadSavedModel(session_options, run_options, slstm_path, {tensorflow::kSavedModelTagServe}, &slstm_model_bundle);

}

void VertexFindingwithDL::process() {
  bool verbose = true;
  Event* event = Event::Instance();
  if(verbose==true){
    std::cout << "########################################################################################################" << std::endl;
    std::cout << "Event Number : " << NEventNumber << std::endl;
    std::cout << "########################################################################################################" << std::endl;
  }

  // clearing old vertices
  if(_primary_vertex->size()>0) {
    delete (*_primary_vertex)[0];
    _primary_vertex->clear();
  }
  for(std::size_t n=0; n<_secondary_vertices->size(); n++)
    delete (*_secondary_vertices)[n];
  _secondary_vertices->clear();

  TrackVec& tracks = event->getTracks();
  std::vector<std::vector<double> > index, pairs, encoder_tracks, decoder_tracks;
  VertexFinderwithDL::GetPairsEncoderDecoderTracks(tracks, NTrackVariable, MaxTrack, index, pairs, encoder_tracks, decoder_tracks);
  if(debug==true) VertexFinderwithDL::DebugPrintVariableSize(pairs, encoder_tracks, decoder_tracks);

  if(verbose==true) std::cout << "Get Event Data ..." << std::endl;
  std::vector<std::vector<double> > event_data = VertexFinderwithDL::GetEventData(debug, index, pairs, pair_model_bundle, pair_pos_model_bundle);
  if(debug==true) VertexFinderwithDL::DebugPrintGetTracks(encoder_tracks);

  // Secondary Seed Selection
  if(verbose==true) std::cout << "Secondary Seed Selection ..." << std::endl;
  std::vector<std::vector<double> > secondary_event_data = VertexFinderwithDL::SecondarySeedSelection(event_data, ThresholdPairSecondaryScore, ThresholdPairPosScore);
  if(debug==true) VertexFinderwithDL::DebugPrintSecondarySort(secondary_event_data);

  // Primary Vertex Finder
  if(verbose==true) std::cout << "Primary Vertex Prediction ..." << std::endl;
  std::vector<int> primary_track_list;
  std::vector<std::vector<std::vector<double> > > primary_scores;
  std::vector<double> bigger_primary_scores;
  VertexFinderwithDL::PrimaryVertexFinder(debug, MaxPrimaryVertexLoop, ThresholdPrimaryScore, event_data, encoder_tracks, decoder_tracks, 
		                          lstm_model_bundle, primary_track_list, primary_scores, bigger_primary_scores);
  if(debug==true) VertexFinderwithDL::DebugPrintVLSTMPrediction(primary_scores);
  if(verbose==true){
    std::cout << "Primary Track List" << std::endl;
    for(std::size_t i=0; i<primary_track_list.size(); i++){
      std::cout << primary_track_list.at(i) << " ";
    }
    std::cout << std::endl;
  }

  // Secondary Vertex Finder
  if(verbose==true) std::cout << "Secondary Vertex Prediction ..." << std::endl;
  std::vector<std::vector<int> > secondary_track_lists;
  VertexFinderwithDL::SecondaryVertexFinder(debug, ThresholdSecondaryScore, bigger_primary_scores, secondary_event_data, encoder_tracks, decoder_tracks,
		                            slstm_model_bundle, primary_track_list, secondary_track_lists);
  if(verbose==true) VertexFinderwithDL::PrintResults(primary_track_list, secondary_track_lists);

  if(verbose==true) std::cout << "Merge Single Track ..." << std::endl;
  std::vector<std::vector<int> > merge_secondary_track_lists = VertexFinderwithDL::MergeSingleTrack(tracks, secondary_track_lists);
  if(verbose==true) VertexFinderwithDL::PrintResults(primary_track_list, merge_secondary_track_lists);

  if(verbose==true) std::cout << "Vertex Making ..." << std::endl;
  VertexFinderwithDL::PrimarySecondaryVertices(verbose, tracks, primary_track_list, merge_secondary_track_lists, _primary_vertex, _secondary_vertices);
  //if(!_primary_vertex->at(0)) std::cout << "PrimaryVertexFinder: No primary vertex found." << std::endl; 
  //if(verbose==true) std::cout << "PrimaryVertexFinder: " << _primary_vertex->at(0)->getTracks().size() << " tracks associated to the primary vertex." << std::endl;
  NEventNumber++;
  if(verbose==true) std::cout << std::endl;
}

void VertexFindingwithDL::end() {
}

}
