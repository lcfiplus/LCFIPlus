// VertexFinderwithDL.h

#ifndef VertexFinderwithDL_h
#define VertexFinderwithDL_h 1

#include "lcfiplus.h"
#include <list>
#include <vector>
#include <functional>
#include <algorithm>

#include <tensorflow/cc/saved_model/loader.h>
#include <tensorflow/cc/saved_model/tag_constants.h>


using namespace std;

namespace lcfiplus {

namespace VertexFinderwithDL {

class VertexFinderwithDLConfig {
 public:
  // main parameters
  int MaxSample, MaxEvent, MaxTrack, MaxNpyVariable, NPairVariable, NTrackVariable;
  int MaxPrimaryVertexLoop;
  double ThresholdPairSecondaryScoreBBCC, ThresholdPairSecondaryScore, ThresholdPairPosScore;
  double ThresholdPrimaryScore, ThresholdSecondaryScore;
  bool debug;

  tensorflow::string pair_path, pair_pos_path, lstm_path, slstm_path;


  // default values
  VertexFinderwithDL(){
    MaxSample = 100000;
    MaxEvent = 100;
    MaxTrack = 53;
    MaxNpyVariable = 59;
    NPairVariable = 44;
    NTrackVariable=22;
    MaxPrimaryVertexLoop = 3;

    ThresholdPairSecondaryScoreBBCC = 0.6;
    ThresholdPairSecondaryScore = 0.8;
    ThresholdPairPosScore = 5;
    ThresholdPrimaryScore = 0.5;
    ThresholdSecondaryScore = 0.8;
  
    debug = false;

    pair_path = "/home/goto/ILC/Deep_Learning/model/Pair_Model_vfdnn04_1Msamples_2500epochs";
    pair_pos_path = "/home/goto/ILC/Deep_Learning/model/Pair_Pos_Model_vfdnn04_1Msamples_2500epochs";
    lstm_path = "/home/goto/ILC/Deep_Learning/model/Attention_Bidirectional_VLSTM_Model_vfdnn06_50000samples_100epochs";
    slstm_path = "/home/goto/ILC/Deep_Learning/model/Attention_Bidirectional_VLSTM_Model_vfdnn06_50000samples_100epochs_ps_100epochs_s";

  }
};

std::vector<std::vector<double> > GetEventTrackPair(std::vector<lcfiplus::Vertex>& vtx, VertexFinderwithDLConfig& cfg);

std::vector<std::vector<double> > SliceN2DVector(std::vector<std::vector<double> > vec, int rowstart, int rowend, int colstart, int colend);
std::vector<std::vector<double> > ConcatN2DVector(std::vector<std::vector<double> > vec1, std::vector<std::vector<double> > vec2);

tensorflow::Tensor N2DVector2Tensor(std::vector<std::vector<double> > vec);
tensorflow::Tensor N3DVector2Tensor(std::vector<std::vector<std::vector<double> > > vec);

std::vector<std::vector<double> > Tensor2N2DVector(tensorflow::Tensor tensor);
std::vector<std::vector<std::vector<double> > > Tensor2N3DVector(tensorflow::Tensor tensor);

void DebugPrintPairPrediction(std::vector<std::vector<double> > data, std::vector<tensorflow::Tensor> tmppair_outputs, std::vector<std::vector<double> > labels);
void DebugPrintGetEventData(std::vector<std::vector<double> > event_data, int ncomb, int NCombination);
void DebugPrintSecondarySort(std::vector<std::vector<double> > secondary_seeds);
void DebugPrintPrimarySort(std::vector<std::vector<double> > primary_seeds);
void DebugPrintGetTracks(std::vector<std::vector<double> > tracks);
void DebugPrintVLSTMPrediction(std::vector<std::vector<std::vector<double> > > scores);

void PrintResults(std::vector<int> primary_track_list, std::vector<std::vector<int> > secondary_track_lists);

void GetPairsEncoderDecoderTracks(TrackVec& tracks, int NTrackVariable, int MaxTrack, std::vector<std::vector<double> >& pairs, 
				  std::vector<std::vector<double> >& encoder_tracks, std::vector<std::vector<double> >& decoder_tracks);

std::vector<std::vector<double> > GetEventData(bool debug, std::vector<std::vector<double> > variables,
                                               tensorflow::SavedModelBundleLite& pair_model_bundle, tensorflow::SavedModelBundleLite& pair_pos_model_bundle);

std::vector<std::vector<double> > GetRemainDecoderTracks(std::vector<std::vector<double> > decoder_tracks, std::vector<int> track_list);

std::vector<std::vector<double> > SecondarySeedSelection(std::vector<std::vector<double> > event_data, 
		                                         int ThresholdPairSecondaryScoreBBCC, int ThresholdPairSecondaryScore, int ThresholdPairPosScore);

void PrimaryVertexFinder(bool debug, int MaxPrimaryVertexLoop, double ThresholdPrimaryScore, std::vector<std::vector<double> > event_data, 
		         std::vector<std::vector<double> > encoder_tracks, std::vector<std::vector<double> > decoder_tracks, 
                         tensorflow::SavedModelBundleLite& lstm_model_bundle,
			 std::vector<int>& primary_track_list, std::vector<std::vector<std::vector<double> > >& primary_scores, std::vector<double>& primary_scores);

void SecondaryVertexFinder(bool debug, double ThresholdSecondaryScore, std::vector<double> primary_scores,
		           std::vector<std::vector<double> > secondary_event_data, 
		           std::vector<std::vector<double> > encoder_tracks, std::vector<std::vector<double> > decoder_tracks, 
                           tensorflow::SavedModelBundleLite& slstm_model_bundle, 
			   std::vector<int>& primary_track_list, std::vector<std::vector<int> >& secondary_track_lists);

void PrimarySecondaryVertices(TrackVec& tracks, std::vector<int> primary_track_list, std::vector<std::vector<int> > secondary_track_lists,
		              Vertex& vtx, std::vector<Vertex*>& vtces)


}
}

#endif

