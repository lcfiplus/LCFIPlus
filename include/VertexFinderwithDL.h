// VertexFinderwithDL.h

#ifndef VertexFinderwithDL_h
#define VertexFinderwithDL_h 1

#include "lcfiplus.h"
#include <list>
#include <vector>
#include <functional>

#include <tensorflow/cc/saved_model/loader.h>
#include <tensorflow/cc/saved_model/tag_constants.h>


using namespace std;

namespace lcfiplus {

namespace VertexFinderwithDL {
}
}
/*
class VertexFinderwithDLConfig {
 public:
  // main parameters
  VertexSelectorConfig v0selTrack; // selector for tight v0 selection: used for track rejection
  VertexSelectorConfig v0selVertex;// selector for loose v0 selection: used for vertex rejection
  double chi2th;	// track chi2 threshold to accept vertices
  double chi2thV0SelTrack; // track chi2 threshold to reject tracks in v0 rejection
  double massth;	// maximum mass to accept vertices
  double chi2orderinglimit; // chi2 threshold to order distance rather than chi2 value

  // associateIPtracks parameters
  double minimumdistIP; // minimum distance to associate IP tracks
  double chi2ratioIP; 	// bias factor to associate IP rather than secondary

  // singletrackvertex parameters
  double minPosSingle;
  double maxPosSingle;
  double minEnergySingle;
  double maxAngleSingle;
  double maxSeparationPerPosSingle;
  double mind0SigSingle;
  double minz0SigSingle;

  //flg for AVF/chi2 algorithm
  bool avf;
  double temperature;   //parameter for avf

  //flg for BNess tagger fake rejection
  bool useBNess;
  double cutBNess;  //parameter for BNess
  double cutBNessE1;  //parameter for BNess

  // default values
  VertexFinderSueharaConfig() {
    v0selTrack.setV0Tight();
    v0selTrack.rejectdist = true;
    v0selTrack.rejectdistnegative = true;
    v0selTrack.rejectdistor = true;
    v0selTrack.minpos = 0.5;

    chi2thV0SelTrack = 2.;

    v0selVertex.setV0Loose();
    v0selVertex.rejectdist = true;
    v0selVertex.minpos = 0.3;

    chi2th = 9.;
    massth = 10.;
    chi2orderinglimit = 1.;

    minimumdistIP = 0.; // not used in default
    chi2ratioIP = 2.;		// biased twice for primary

    minPosSingle = 0.3;
    maxPosSingle = 30.;
    minEnergySingle = 1.;
    maxAngleSingle = 0.5;
    maxSeparationPerPosSingle = 0.1;
    mind0SigSingle = 5.;
    minz0SigSingle = 5.;

    avf = false;   //default: use chi2 algorithm
    temperature = 5.0;

    useBNess = false; //default: do not use BNess
    cutBNess = -0.80;
    cutBNessE1 = -0.15;
  }
};

std::vector<std::vector<double> > GetEventTrackPair(std::vector<lcfiplus::Vertex>& vtx, VertexFinderwithDLConfig& cfg);

int CountCombination(std::vector<std::vector<double> > data, int ievent);
  int ncomb=0;
  for(int i=0; i<data.size(); i++){
    if(data.at(i).at(0)==ievent) ncomb++;
  }
  return ncomb;
}


std::vector<std::vector<double> > SliceN2DVector(std::vector<std::vector<double> > vec, int rowstart, int rowend, int colstart, int colend){
  std::vector<std::vector<double> > new_vec(rowend-rowstart, std::vector<double>(colend-colstart));
  for(int i=rowstart; i<rowend; i++){
    for(int j=colstart; j<colend; j++){
      new_vec.at(i-rowstart).at(j-colstart) = vec.at(i).at(j);
    }
  }
  return new_vec;
}

std::vector<std::vector<double> > ConcatN2DVector(std::vector<std::vector<double> > vec1, std::vector<std::vector<double> > vec2){
  if(vec1.size()!=vec2.size()){
    std::cout << "error: vector shapes are not correct." << std::endl;
    std::cout << "vec1: " << vec1.size() << "vec2: " << vec2.size() << std::endl;
    exit(1);
  }
  std::vector<std::vector<double> > new_vec = vec1;
  for(int i=0; i<new_vec.size(); i++){
    new_vec.at(i).insert(new_vec.at(i).end(), vec2.at(i).begin(), vec2.at(i).end());
  }
  return new_vec;
}

tensorflow::Tensor N2DVector2Tensor(std::vector<std::vector<double> > vec){
  const int row = vec.size(), col = vec.at(0).size();
  tensorflow::Tensor tensor = tensorflow::Tensor(tensorflow::DT_FLOAT, tensorflow::TensorShape({row, col}));
  tensorflow::TTypes<float, 2>::Tensor tensor_element = tensor.tensor<float, 2>();
  tensor_element.setZero();
  for(int i=0; i<row; i++){
    for(int j=0; j<col; j++){
      tensor_element(i, j) = float(vec.at(i).at(j));
    }
  }
  return tensor;
}

tensorflow::Tensor N3DVector2Tensor(std::vector<std::vector<std::vector<double> > > vec){
  const int height = vec.size(), width = vec.at(0).size(), depth = vec.at(0).at(0).size();
  tensorflow::Tensor tensor = tensorflow::Tensor(tensorflow::DT_FLOAT, tensorflow::TensorShape({height, width, depth}));
  tensorflow::TTypes<float, 3>::Tensor tensor_element = tensor.tensor<float, 3>();
  tensor_element.setZero();
  for(int i=0; i<height; i++){
    for(int j=0; j<width; j++){
      for(int k=0; k<depth; k++){
        tensor_element(i, j, k) = float(vec.at(i).at(j).at(k));
      }
    }
  }
  return tensor;
}

std::vector<std::vector<double> > Tensor2N2DVector(tensorflow::Tensor tensor){
  const int row = tensor.dim_size(0), col = tensor.dim_size(1);
  std::vector<std::vector<double> > vec(row, std::vector<double>(col));
  tensorflow::TTypes<float, 2>::Tensor tensor_element = tensor.tensor<float, 2>();
  for(int i=0; i<row; i++){
    for(int j=0; j<col; j++){
      vec.at(i).at(j) = double(tensor_element(i, j));
    }
  }
  return vec;
}

std::vector<std::vector<std::vector<double> > > Tensor2N3DVector(tensorflow::Tensor tensor){
  const int height = tensor.dim_size(0), width = tensor.dim_size(1), depth = tensor.dim_size(2);
  std::vector<std::vector<std::vector<double> > > vec(height, std::vector<std::vector<double> > (width, std::vector<double>(depth)));
  tensorflow::TTypes<float, 3>::Tensor tensor_element = tensor.tensor<float, 3>();
  for(int i=0; i<height; i++){
    for(int j=0; j<width; j++){
      for(int k=0; k<depth; k++){
        vec.at(i).at(j).at(k) = double(tensor_element(i, j, k));
      }
    }
  }
  return vec;
}

void DebugPrintPairPrediction(std::vector<std::vector<double> > data, std::vector<tensorflow::Tensor> tmppair_outputs, std::vector<std::vector<double> > labels){
  std::cout << "Predict Vertex Finder Shape: " << tmppair_outputs.size() << std::endl;
  std::cout << "Num Elements: " << std::to_string(tmppair_outputs[0].NumElements()) 
	    << " Num Elements Dims 0: " << std::to_string(tmppair_outputs[0].dim_size(0)) 
	    << " Num Elements Dims 0: " << std::to_string(tmppair_outputs[0].dim_size(1)) << std::endl;

  std::cout << "Pair Prediction Checking ..." << std::endl;
  for(int i=0; i<10; i++){
    std::cout << "Predicted Scores: " << labels.at(i).at(0) << " " << labels.at(i).at(1) << " " 
	                              << labels.at(i).at(2) << " " << labels.at(i).at(3) << " " << labels.at(i).at(4)
	      << " ;; True Label: " << data.at(i).at(57) << std::endl;
  }
}

void DebugPrintGetEventData(std::vector<std::vector<double> > event_data, int ncomb, int NCombination){
  std::cout << "Event Data Size Checking ... : " << event_data.size() << " = " << ncomb << " = " << NCombination << std::endl;
}

void DebugPrintSecondarySort(std::vector<std::vector<double> > secondary_seeds){
  int iMax = secondary_seeds.size();
  if(secondary_seeds.size() > 10) iMax = 10;
  for(int i=0; i<iMax; i++){
    std::cout << "Track 1 " << secondary_seeds.at(i).at(1) << " Track 2 " << secondary_seeds.at(i).at(2) 
	      << " SV Score " << secondary_seeds.at(i).at(61)+secondary_seeds.at(i).at(62)+secondary_seeds.at(i).at(63)
	      << " NC " << secondary_seeds.at(i).at(59) << " PV " << secondary_seeds.at(i).at(60) << std::endl;
  }
}

void DebugPrintPrimarySort(std::vector<std::vector<double> > primary_seeds){
  int iMax = primary_seeds.size();
  if(primary_seeds.size() > 10) iMax = 10;
  for(int i=0; i<iMax; i++){
    std::cout << "Track 1 " << primary_seeds.at(i).at(0) << " Track 2 " << primary_seeds.at(i).at(1) << " Score " << primary_seeds.at(i).at(60) << std::endl;
  }
}

void DebugPrintGetTracks(std::vector<std::vector<double> > tracks){
  std::cout << "Track Shape Check ;; NTracks: " << tracks.size() << " NVariables: " << tracks.at(0).size() << std::endl;
  std::cout << "Track Variables Check" << std::endl;
  for(int i=0; i<tracks.size(); i++){
    std::cout << "Track " << i << std::endl;
    for(int j=0; j<tracks.at(0).size(); j++){
      std::cout << j << ": " <<  tracks.at(i).at(j) << "  ";
      if((j+1)%10==0) std::cout << std::endl;
    }
    std::cout << std::endl;
  }
}

void DebugPrintVLSTMPrediction(std::vector<std::vector<std::vector<double> > > scores){
  std::cout << "VLSTM Prediction Checking ..." << std::endl;
  int MaxLoop = scores.size();
  if(scores.size()>5) MaxLoop = 5;
  for(int i=0; i<MaxLoop; i++){
    std::cout << "Predicted Scores" << std::endl;
    for(int j=0; j<scores.at(0).size(); j++){
      std::cout << "Track " << j << ": " << std::setprecision(4) << scores.at(i).at(j).at(0) << " ";
      if((j+1)%5==0) std::cout << std::endl;
    }
    std::cout << std::endl;
  }
}

void PrintResults(std::vector<int> primary_track_list, std::vector<std::vector<int> > secondary_track_lists){
  std::cout << "Finish !!" << std::endl;

  std::cout << "Primary Track List" << std::endl;
  for(int i=0; i<primary_track_list.size(); i++){
    std::cout << primary_track_list.at(i) << " ";
  }
  std::cout << std::endl;
  std::cout << std::endl;

  std::cout << "Secondary Track Lists" << std::endl;
  for(int i=0; i<secondary_track_lists.size(); i++){
    std::cout << "List " << i << std::endl;
    for(int j=0; j<secondary_track_lists.at(i).size(); j++){
      std::cout << secondary_track_lists.at(i).at(j) << " ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
}

void GetEventData(std::vector<std::vector<double> > data, int ievent, int MaxNpyVariable, std::vector<std::vector<double> >& event_data, int& ncomb){
  for(int i=0; i<data.size(); i++){
    if(data.at(i).at(0)==ievent){
      event_data.at(ncomb).resize(MaxNpyVariable);
      event_data.at(ncomb) = data.at(i);
      ncomb++;
    }
  }
}

void GetEncoderDecoderTracks(std::vector<std::vector<double> > event_data, int NTrack, int NTrackVariable,
                             std::vector<std::vector<double> >& encoder_tracks, std::vector<std::vector<double> >& decoder_tracks){
  int offset1, offset2;
  for(int i=0; i<NTrack; i++){
    if(i<NTrack-1) offset1 = 1, offset2 = 24;
    else if(i==NTrack-1) offset1 = 0, offset2 = 2;
    for(int j=1; j<NTrackVariable+1; j++){
      if(j==1){
        encoder_tracks.at(i).at(j-1) = double(1);
        decoder_tracks.at(i).at(j-1) = double(1);
      }
      encoder_tracks.at(i).at(j) = event_data.at(int(event_data.size())-NTrack+offset1+i).at(j+offset2);
      decoder_tracks.at(i).at(j) = event_data.at(int(event_data.size())-NTrack+offset1+i).at(j+offset2);
    }
  }
}

std::vector<std::vector<double> > GetRemainDecoderTracks(std::vector<std::vector<double> > decoder_tracks, std::vector<int> track_list){
  std::vector<std::vector<double> > remain_decoder_tracks;
  for(int i=0; i<decoder_tracks.size(); i++){
    if(std::find(track_list.begin(), track_list.end(), i) != track_list.end()){
      //std::cout << "add track num " << i << std::endl;
      remain_decoder_tracks.push_back(decoder_tracks.at(i));
    }
  }
  return remain_decoder_tracks;
}

std::vector<std::vector<double> > SecondarySeedSelection(std::vector<std::vector<double> > event_data, 
		                                         int ThresholdPairSecondaryScoreBBCC, int ThresholdPairSecondaryScore, int ThresholdPairPosScore){ 
  std::vector<std::vector<double> > secondary_event_data;
  for(int i; i<event_data.size(); i++){
    double tmp_secondary_score = event_data.at(i).at(61) + event_data.at(i).at(62) + event_data.at(i).at(63);
    if((event_data.at(i).at(59) > event_data.at(i).at(60) and event_data.at(i).at(59) > event_data.at(i).at(61) and 
	event_data.at(i).at(59) > event_data.at(i).at(62) and event_data.at(i).at(59) > event_data.at(i).at(63)) or 
       (event_data.at(i).at(60) > event_data.at(i).at(59) and event_data.at(i).at(60) > event_data.at(i).at(61) and
        event_data.at(i).at(60) > event_data.at(i).at(62) and event_data.at(i).at(59) > event_data.at(i).at(63))) continue;
    //if((event_data.at(i).at(62) > ThresholdPairSecondaryScoreBBCC) or (event_data.at(i).at(63) > ThresholdPairSecondaryScoreBBCC)) continue;
    if(tmp_secondary_score < ThresholdPairSecondaryScore) continue;
    if(event_data.at(i).at(64) > ThresholdPairPosScore) continue;
    secondary_event_data.push_back(event_data.at(i));
  }
  sort(secondary_event_data.begin(), secondary_event_data.end(), [](const std::vector<double> &alpha, const std::vector<double> &beta){
    return alpha.at(61)+alpha.at(62)+alpha.at(63) > beta.at(61)+beta.at(62)+beta.at(63);
  });
  return secondary_event_data;
}

void SecondaryFirstVertexFinder(std::vector<std::vector<double> > secondary_event_data, 
		                std::vector<std::vector<double> > encoder_tracks, std::vector<std::vector<double> > decoder_tracks, 
                                tensorflow::SavedModelBundleLite& slstm_model_bundle, std::vector<std::vector<std::vector<double> > >& secondary_first_scores){
  std::vector<std::vector<double> > secondary_pairs = SliceN2DVector(secondary_event_data, 0, secondary_event_data.size(), 3, 47);
  std::vector<std::vector<std::vector<double> > > secondary_encoder_tracks(secondary_pairs.size(), encoder_tracks);
  std::vector<std::vector<std::vector<double> > > secondary_decoder_tracks(secondary_pairs.size(), decoder_tracks);

  tensorflow::Tensor tsecondary_pairs = N2DVector2Tensor(secondary_pairs);
  tensorflow::Tensor tsecondary_encoder_tracks = N3DVector2Tensor(secondary_encoder_tracks);
  tensorflow::Tensor tsecondary_decoder_tracks = N3DVector2Tensor(secondary_decoder_tracks);
    
  std::vector<tensorflow::Tensor> tmpsecondary_first_outputs;
  tensorflow::Status runStatus = slstm_model_bundle.GetSession()->Run({{"serving_default_Decoder_Input:0", tsecondary_decoder_tracks},
		                                                       {"serving_default_Encoder_Input:0", tsecondary_encoder_tracks},
						                       {"serving_default_Pair_Input:0", tsecondary_pairs}}, 
				          	  	              {"StatefulPartitionedCall:0"},
						                      {}, &tmpsecondary_first_outputs);

  secondary_first_scores = Tensor2N3DVector(tmpsecondary_first_outputs[0]);
}


void PrimaryVertexFinder(int MaxPrimaryVertexLoop, double ThresholdPrimaryScore, std::vector<std::vector<double> > event_data, 
		         std::vector<std::vector<double> > encoder_tracks, std::vector<std::vector<double> > decoder_tracks, 
			 bool debug,
                         tensorflow::SavedModelBundleLite& lstm_model_bundle,
			 std::vector<int>& primary_track_list, std::vector<std::vector<std::vector<double> > >& primary_scores){
  std::vector<std::vector<double> > primary_event_data = event_data;
  sort(primary_event_data.begin(), primary_event_data.end(), [](const std::vector<double> &alpha, const std::vector<double> &beta){
    return alpha.at(60) > beta.at(60);
  });

  if(debug==true){
    for(int i=0; i<MaxPrimaryVertexLoop; i++){
      std::cout << "Track 1 " << primary_event_data.at(i).at(1) << " Track 2 " << primary_event_data.at(i).at(2) << " Score " << primary_event_data.at(i).at(60) << std::endl;
    }
  }

  std::vector<std::vector<double> > primary_pairs = SliceN2DVector(primary_event_data, 0, MaxPrimaryVertexLoop, 3, 47);
  std::vector<std::vector<std::vector<double> > > primary_encoder_tracks(primary_pairs.size(), encoder_tracks);
  std::vector<std::vector<std::vector<double> > > primary_decoder_tracks(primary_pairs.size(), decoder_tracks);

  tensorflow::Tensor tprimary_pairs = N2DVector2Tensor(primary_pairs);
  tensorflow::Tensor tprimary_encoder_tracks = N3DVector2Tensor(primary_encoder_tracks);
  tensorflow::Tensor tprimary_decoder_tracks = N3DVector2Tensor(primary_decoder_tracks);
    
  std::vector<tensorflow::Tensor> tmpprimary_outputs;
  tensorflow::Status runStatus = lstm_model_bundle.GetSession()->Run({{"serving_default_Decoder_Input:0", tprimary_decoder_tracks},
                                                                      {"serving_default_Encoder_Input:0", tprimary_encoder_tracks},
                                                                      {"serving_default_Pair_Input:0", tprimary_pairs}}, 
                                                                     {"StatefulPartitionedCall:0"},
							     	     {}, &tmpprimary_outputs);

  primary_scores = Tensor2N3DVector(tmpprimary_outputs[0]);
  for(int i=0; i<primary_scores.size(); i++){
    for(int j=0; j<primary_scores.at(0).size(); j++){
      if(primary_scores.at(i).at(j).at(0) > ThresholdPrimaryScore){
	if(debug==true) std::cout << "Track " << j << " Primary Score: " << primary_scores.at(i).at(j).at(0) << std::endl;
	primary_track_list.push_back(j);
      }
    }
  }
  std::sort(primary_track_list.begin(), primary_track_list.end());
  primary_track_list.erase(std::unique(primary_track_list.begin(), primary_track_list.end()), primary_track_list.end());
}


void SecondaryVertexFinder(double ThresholdSecondaryScore, std::vector<double> primary_scores,
		           std::vector<std::vector<double> > secondary_event_data, 
		           std::vector<std::vector<double> > encoder_tracks, std::vector<std::vector<double> > decoder_tracks, 
			   bool debug,
                           tensorflow::SavedModelBundleLite& slstm_model_bundle, 
			   std::vector<int>& primary_track_list, std::vector<std::vector<int> >& secondary_track_lists){
  std::vector<int> track_list(decoder_tracks.size());
  std::iota(track_list.begin(), track_list.end(), 0);
  for(int i=0; i<secondary_event_data.size(); i++){
    int track1 = secondary_event_data.at(i).at(1), track2 = secondary_event_data.at(i).at(2);
    if(std::find(primary_track_list.begin(), primary_track_list.end(), track1) != primary_track_list.end() or
       std::find(primary_track_list.begin(), primary_track_list.end(), track2) != primary_track_list.end()) continue;

    if(debug==true){
      std::cout << "Track List" << std::endl;
      for(int i=0; i<track_list.size(); i++){
        std::cout << track_list.at(i) << " ";
      }
      std::cout << std::endl;
    
      std::cout << "seed: " << track1 << " " << track2 << std::endl;
    }

    std::vector<std::vector<double> > remain_decoder_tracks = GetRemainDecoderTracks(decoder_tracks, track_list);

    std::vector<std::vector<double> > secondary_pair = SliceN2DVector(secondary_event_data, i, i+1, 3, 47);
    std::vector<std::vector<std::vector<double> > > secondary_encoder_tracks(1, encoder_tracks);
    std::vector<std::vector<std::vector<double> > > secondary_remain_decoder_tracks(1, remain_decoder_tracks);

    tensorflow::Tensor tsecondary_pair = N2DVector2Tensor(secondary_pair);
    tensorflow::Tensor tsecondary_encoder_tracks = N3DVector2Tensor(secondary_encoder_tracks);
    tensorflow::Tensor tsecondary_remain_decoder_tracks = N3DVector2Tensor(secondary_remain_decoder_tracks);
    
    std::vector<tensorflow::Tensor> tmpsecondary_outputs;
    tensorflow::Status runStatus = slstm_model_bundle.GetSession()->Run({{"serving_default_Decoder_Input:0", tsecondary_remain_decoder_tracks},
		                                                         {"serving_default_Encoder_Input:0", tsecondary_encoder_tracks},
						                         {"serving_default_Pair_Input:0", tsecondary_pair}}, 
				          	  	                {"StatefulPartitionedCall:0"},
									{}, &tmpsecondary_outputs);
    std::vector<std::vector<std::vector<double> > > secondary_scores = Tensor2N3DVector(tmpsecondary_outputs[0]);
    if(debug==true) DebugPrintVLSTMPrediction(secondary_scores);

    std::vector<int> tmpsecondary_track_list, tmptrack_list = track_list;
    for(int j=0; j<track_list.size(); j++){
      if(secondary_scores.at(0).at(j).at(0) > ThresholdSecondaryScore){
        if(std::find(primary_track_list.begin(), primary_track_list.end(), track_list.at(j)) == primary_track_list.end()){
          tmpsecondary_track_list.push_back(track_list.at(j));
	  tmptrack_list.erase(std::remove(tmptrack_list.begin(), tmptrack_list.end(), track_list.at(j)), tmptrack_list.end());
	}
	else if(std::find(primary_track_list.begin(), primary_track_list.end(), track_list.at(j)) != primary_track_list.end() and
	        secondary_scores.at(0).at(j).at(0) > primary_scores.at(track_list.at(j))){
          tmpsecondary_track_list.push_back(track_list.at(j));
	  if(debug==true){
	    std::cout << "scramble track num " << track_list.at(j) << " sv score " << secondary_scores.at(0).at(j).at(0) 
		      << " pv score " << primary_scores.at(track_list.at(j)) << std::endl;
	  }
	  primary_track_list.erase(std::remove(primary_track_list.begin(), primary_track_list.end(), track_list.at(j)), primary_track_list.end());
	  tmptrack_list.erase(std::remove(tmptrack_list.begin(), tmptrack_list.end(), track_list.at(j)), tmptrack_list.end());
	}
      }
    }
    if(tmpsecondary_track_list.size()!=0) secondary_track_lists.push_back(tmpsecondary_track_list);
    track_list = tmptrack_list;
  }
}




//vector<lcfiplus::Vertex*> * findSueharaVertices(const Event& evt, const Jet& jet);
}

#endif
*/

