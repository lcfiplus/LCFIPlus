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

std::vector<std::vector<double> > SliceN2DVector(std::vector<std::vector<double> > vec, int rowstart, int rowend, int colstart, int colend);
std::vector<std::vector<double> > ConcatN2DVector(std::vector<std::vector<double> > vec1, std::vector<std::vector<double> > vec2);

tensorflow::Tensor N2DVector2Tensor(std::vector<std::vector<double> > vec);
tensorflow::Tensor N3DVector2Tensor(std::vector<std::vector<std::vector<double> > > vec);

std::vector<std::vector<double> > Tensor2N2DVector(tensorflow::Tensor tensor);
std::vector<std::vector<std::vector<double> > > Tensor2N3DVector(tensorflow::Tensor tensor);

void DebugPrintGetEventData(std::vector<std::vector<double> > event_data, int ncomb, int NCombination);
void DebugPrintSecondarySort(std::vector<std::vector<double> > secondary_seeds);
void DebugPrintPrimarySort(std::vector<std::vector<double> > primary_seeds);
void DebugPrintGetTracks(std::vector<std::vector<double> > tracks);
void DebugPrintVLSTMPrediction(std::vector<std::vector<std::vector<double> > > scores);

std::vector<std::vector<double> > GetEventData(std::vector<std::vector<double> > variables,
                                               tensorflow::SavedModelBundleLite& pair_model_bundle, tensorflow::SavedModelBundleLite& pair_pos_model_bundle);

std::vector<std::vector<double> > GetRemainDecoderTracks(std::vector<std::vector<double> > decoder_tracks, std::vector<int> track_list);

std::vector<std::vector<double> > SecondarySeedSelection(std::vector<std::vector<double> > event_data, 
		                                         int ThresholdPairSecondaryScore, int ThresholdPairPosScore);

void PrintResults(std::vector<int> primary_track_list, std::vector<std::vector<int> > secondary_track_lists);

void GetPairsEncoderDecoderTracks(TrackVec& tracks, int NTrackVariable, int MaxTrack, std::vector<std::vector<double> >& pairs, 
				  std::vector<std::vector<double> >& encoder_tracks, std::vector<std::vector<double> >& decoder_tracks);

void PrimaryVertexFinder(bool debug, int MaxPrimaryVertexLoop, double ThresholdPrimaryScore, std::vector<std::vector<double> > event_data, 
		         std::vector<std::vector<double> > encoder_tracks, std::vector<std::vector<double> > decoder_tracks, 
                         tensorflow::SavedModelBundleLite& lstm_model_bundle,
			 std::vector<int>& primary_track_list, std::vector<std::vector<std::vector<double> > >& primary_scores, std::vector<double>& bigger_primary_scores);

void SecondaryVertexFinder(bool debug, double ThresholdSecondaryScore, std::vector<double> primary_scores,
		           std::vector<std::vector<double> > secondary_event_data, 
		           std::vector<std::vector<double> > encoder_tracks, std::vector<std::vector<double> > decoder_tracks, 
                           tensorflow::SavedModelBundleLite& slstm_model_bundle, 
			   std::vector<int>& primary_track_list, std::vector<std::vector<int> >& secondary_track_lists);

std::vector<std::vector<int> > MergeSingleTrack(TrackVec& tracks, std::vector<std::vector<int> > secondary_track_lists);

void PrimarySecondaryVertices(TrackVec& tracks, std::vector<int> primary_track_list, std::vector<std::vector<int> > secondary_track_lists,
		              Vertex& vtx, std::vector<Vertex*>& vtces);

}

}


#endif

