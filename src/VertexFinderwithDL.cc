// VertexFinderSuehara.cc
#include <cmath>
#include <algorithm>

#include "lcfiplus.h"
#include "VertexFinderwithDL.h"
#include "LcfiInterface.h"
#include "VertexFitterSimple.h"

#include <tensorflow/cc/saved_model/loader.h>
#include <tensorflow/cc/saved_model/tag_constants.h>


using namespace lcfiplus;
using namespace lcfiplus::VertexFinderwithDL;


std::vector<std::vector<double> > VertexFinderwithDL::SliceN2DVector(std::vector<std::vector<double> > vec, int rowstart, int rowend, int colstart, int colend){
  std::vector<std::vector<double> > new_vec(rowend-rowstart, std::vector<double>(colend-colstart));
  for(int i=rowstart; i<rowend; i++){
    for(int j=colstart; j<colend; j++){
      new_vec.at(i-rowstart).at(j-colstart) = vec.at(i).at(j);
    }
  }
  return new_vec;
}

std::vector<std::vector<double> > VertexFinderwithDL::ConcatN2DVector(std::vector<std::vector<double> > vec1, std::vector<std::vector<double> > vec2){
  if(vec1.size()!=vec2.size()){
    std::cout << "error: vector shapes are not correct." << std::endl;
    std::cout << "vec1: " << vec1.size() << "vec2: " << vec2.size() << std::endl;
    exit(1);
  }
  std::vector<std::vector<double> > new_vec = vec1;
  for(std::size_t i=0; i<new_vec.size(); i++){
    new_vec.at(i).insert(new_vec.at(i).end(), vec2.at(i).begin(), vec2.at(i).end());
  }
  return new_vec;
}

tensorflow::Tensor VertexFinderwithDL::N2DVector2Tensor(std::vector<std::vector<double> > vec){
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

tensorflow::Tensor VertexFinderwithDL::N3DVector2Tensor(std::vector<std::vector<std::vector<double> > > vec){
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

std::vector<std::vector<double> > VertexFinderwithDL::Tensor2N2DVector(tensorflow::Tensor tensor){
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

std::vector<std::vector<std::vector<double> > > VertexFinderwithDL::Tensor2N3DVector(tensorflow::Tensor tensor){
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

void VertexFinderwithDL::DebugPrintGetEventData(std::vector<std::vector<double> > event_data, int ncomb, int NCombination){
  std::cout << "Event Data Size Checking ... : " << event_data.size() << " = " << ncomb << " = " << NCombination << std::endl;
}

void VertexFinderwithDL::DebugPrintSecondarySort(std::vector<std::vector<double> > secondary_seeds){
  int iMax = secondary_seeds.size();
  if(secondary_seeds.size() > 10) iMax = 10;
  for(int i=0; i<iMax; i++){
    std::cout << "Track 1 " << secondary_seeds.at(i).at(1) << " Track 2 " << secondary_seeds.at(i).at(2)
	      << " SV Score " << secondary_seeds.at(i).at(61)+secondary_seeds.at(i).at(62)+secondary_seeds.at(i).at(63)
	      << " NC " << secondary_seeds.at(i).at(59) << " PV " << secondary_seeds.at(i).at(60) << std::endl;
  }
}

void VertexFinderwithDL::DebugPrintPrimarySort(std::vector<std::vector<double> > primary_seeds){
  int iMax = primary_seeds.size();
  if(primary_seeds.size() > 10) iMax = 10;
  for(int i=0; i<iMax; i++){
    std::cout << "Track 1 " << primary_seeds.at(i).at(0) << " Track 2 " << primary_seeds.at(i).at(1) << " Score " << primary_seeds.at(i).at(60) << std::endl;
  }
}

void VertexFinderwithDL::DebugPrintGetTracks(std::vector<std::vector<double> > tracks){
  std::cout << "Track Shape Check ;; NTracks: " << tracks.size() << " NVariables: " << tracks.at(0).size() << std::endl;
  std::cout << "Track Variables Check" << std::endl;
  for(std::size_t i=0; i<tracks.size(); i++){
    std::cout << "Track " << i << std::endl;
    for(std::size_t j=0; j<tracks.at(0).size(); j++){
      std::cout << j << ": " <<  tracks.at(i).at(j) << "  ";
      if((j+1)%10==0) std::cout << std::endl;
    }
    std::cout << std::endl;
  }
}

void VertexFinderwithDL::DebugPrintVLSTMPrediction(std::vector<std::vector<std::vector<double> > > scores){
  std::cout << "VLSTM Prediction Checking ..." << std::endl;
  int MaxLoop = scores.size();
  if(scores.size()>5) MaxLoop = 5;
  for(int i=0; i<MaxLoop; i++){
    std::cout << "Predicted Scores" << std::endl;
    for(std::size_t j=0; j<scores.at(0).size(); j++){
      std::cout << "Track " << j << ": " << std::setprecision(4) << scores.at(i).at(j).at(0) << " ";
      if((j+1)%5==0) std::cout << std::endl;
    }
    std::cout << std::endl;
  }
}

void VertexFinderwithDL::PrintResults(std::vector<int> primary_track_list, std::vector<std::vector<int> > secondary_track_lists){
  std::cout << "Finish !!" << std::endl;

  std::cout << "Primary Track List" << std::endl;
  for(std::size_t i=0; i<primary_track_list.size(); i++){
    std::cout << primary_track_list.at(i) << " ";
  }
  std::cout << std::endl;
  std::cout << std::endl;

  std::cout << "Secondary Track Lists" << std::endl;
  for(std::size_t i=0; i<secondary_track_lists.size(); i++){
    std::cout << "List " << i << std::endl;
    for(std::size_t j=0; j<secondary_track_lists.at(i).size(); j++){
      std::cout << secondary_track_lists.at(i).at(j) << " ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
}


void VertexFinderwithDL::GetPairsEncoderDecoderTracks(TrackVec& tracks, int NTrackVariable, int MaxTrack,
                                                 std::vector<std::vector<double> >& pairs,
                                                 std::vector<std::vector<double> >& encoder_tracks, 
						 std::vector<std::vector<double> >& decoder_tracks){
  
  std::vector<double> tmpzero_track(NTrackVariable+1);

  // One track
  for(std::size_t i=0; i<tracks.size(); i++){
    const Track* track = tracks.at(i);
    std::vector<double> tmptrack;
    std::vector<double> tmppair;
    double tr1d0 = tanh(track->getD0());
    double tr1z0 = tanh(track->getZ0());
    double tr1phi = (1.0/M_PI) * track->getPhi();
    double tr1omega = tanh(200 * track->getOmega());
    double tr1tanlam = tanh(0.3 * track->getTanLambda());
    double tr1charge = track->getCharge();
    
    tmptrack.push_back(1);
    tmptrack.push_back(tr1d0);
    tmptrack.push_back(tr1z0);
    tmptrack.push_back(tr1phi);
    tmptrack.push_back(tr1omega);
    tmptrack.push_back(tr1tanlam);
    tmptrack.push_back(tr1charge);
    tmppair.push_back(tr1d0);
    tmppair.push_back(tr1z0);
    tmppair.push_back(tr1phi);
    tmppair.push_back(tr1omega);
    tmppair.push_back(tr1tanlam);
    tmppair.push_back(tr1charge);

    TLorentzVector tlv = *track;
    double tr1energy = tanh(0.5 * (tlv.E() - 5.0));
    tmptrack.push_back(tr1energy);
    tmppair.push_back(tr1energy);

    const double* cov = track->getCovMatrix();
    for(int c=0; c<15; c++){
      tmptrack.push_back(tanh(8000 * (cov[c] + 0.0005)));
      tmppair.push_back(tanh(8000 * (cov[c] + 0.0005)));
    }

    // The other track
    for(std::size_t j=0; j<i; j++){
      const Track* other_track = tracks.at(j);
      double tr2d0 = tanh(other_track->getD0());
      double tr2z0 = tanh(other_track->getZ0());
      double tr2phi = (1.0/M_PI) * other_track->getPhi();
      double tr2omega = tanh(200 * other_track->getOmega());
      double tr2tanlam = tanh(0.3 * other_track->getTanLambda());
      double tr2charge = other_track->getCharge();
      tmppair.push_back(tr2d0);
      tmppair.push_back(tr2z0);
      tmppair.push_back(tr2phi);
      tmppair.push_back(tr2omega);
      tmppair.push_back(tr2tanlam);
      tmppair.push_back(tr2charge);

      TLorentzVector other_tlv = *other_track;
      double tr2energy = tanh(0.5 * (other_tlv.E() - 5.0));
      tmptrack.push_back(tr2energy);

      const double* other_cov = other_track->getCovMatrix();
      for(int c=0; c<15; c++){
        tmptrack.push_back(tanh(8000 * (other_cov[c] + 0.0005)));
      }
    }
    // One track
    encoder_tracks.push_back(tmptrack);
    decoder_tracks.push_back(tmptrack);
    pairs.push_back(tmppair);
  }
  for(int i=tracks.size(); i<MaxTrack; i++){
    encoder_tracks.push_back(tmpzero_track);
  }
}

std::vector<std::vector<double> > VertexFinderwithDL::GetEventData(std::vector<std::vector<double> > variables,
                                                                   tensorflow::SavedModelBundleLite& pair_model_bundle,
								   tensorflow::SavedModelBundleLite& pair_pos_model_bundle){
  tensorflow::Tensor tvariables = VertexFinderwithDL::N2DVector2Tensor(variables);
  std::vector<tensorflow::Tensor> tmppair_outputs;
  tensorflow::Status pair_runStatus = pair_model_bundle.GetSession()->Run({{"serving_default_input_1:0", tvariables}},
				          		                  {"StatefulPartitionedCall:0"},
						                          {}, &tmppair_outputs);
  std::vector<std::vector<double> > _labels = VertexFinderwithDL::Tensor2N2DVector(tmppair_outputs[0]);

  std::vector<tensorflow::Tensor> tmppair_pos_outputs;
  tensorflow::Status pair_pos_runStatus = pair_pos_model_bundle.GetSession()->Run({{"serving_default_input_1:0", tvariables}},
				          		                          {"StatefulPartitionedCall:0"},
						                                  {}, &tmppair_pos_outputs);
  std::vector<std::vector<double> > _positions = VertexFinderwithDL::Tensor2N2DVector(tmppair_pos_outputs[0]);

  std::vector<std::vector<double> > _variables_labels = VertexFinderwithDL::ConcatN2DVector(variables, _labels); // 44:NC 45:PV 46CC 47BB 48BC
  std::vector<std::vector<double> > event_data = VertexFinderwithDL::ConcatN2DVector(_variables_labels, _positions); // 49Pos
  return event_data;
}

std::vector<std::vector<double> > VertexFinderwithDL::GetRemainDecoderTracks(std::vector<std::vector<double> > decoder_tracks, std::vector<int> track_list){
  std::vector<std::vector<double> > remain_decoder_tracks;
  for(std::size_t i=0; i<decoder_tracks.size(); i++){
    if(std::find(track_list.begin(), track_list.end(), i) != track_list.end()){
      remain_decoder_tracks.push_back(decoder_tracks.at(i));
    }
  }
  return remain_decoder_tracks;
}

std::vector<std::vector<double> > VertexFinderwithDL::SecondarySeedSelection(std::vector<std::vector<double> > event_data, int ThresholdPairSecondaryScore, int ThresholdPairPosScore){
  std::vector<std::vector<double> > secondary_event_data;
  for(std::size_t i=0; i<event_data.size(); i++){
    double tmp_secondary_score = event_data.at(i).at(46) + event_data.at(i).at(47) + event_data.at(i).at(48);
    if((event_data.at(i).at(44) > event_data.at(i).at(45) and event_data.at(i).at(44) > event_data.at(i).at(46) and
	event_data.at(i).at(44) > event_data.at(i).at(47) and event_data.at(i).at(44) > event_data.at(i).at(48)) or
       (event_data.at(i).at(45) > event_data.at(i).at(44) and event_data.at(i).at(45) > event_data.at(i).at(46) and
        event_data.at(i).at(45) > event_data.at(i).at(47) and event_data.at(i).at(45) > event_data.at(i).at(48))) continue;
    //if((event_data.at(i).at(46) > ThresholdPairSecondaryScoreBBCC) or (event_data.at(i).at(47) > ThresholdPairSecondaryScoreBBCC)) continue;
    if(tmp_secondary_score < ThresholdPairSecondaryScore) continue;
    if(event_data.at(i).at(49) > ThresholdPairPosScore) continue;
    secondary_event_data.push_back(event_data.at(i));
  }
  sort(secondary_event_data.begin(), secondary_event_data.end(), [](const std::vector<double> &alpha, const std::vector<double> &beta){
    return alpha.at(46)+alpha.at(47)+alpha.at(48) > beta.at(46)+beta.at(47)+beta.at(48);
  });
  return secondary_event_data;
}

void VertexFinderwithDL::PrimaryVertexFinder(bool debug, int MaxPrimaryVertexLoop, double ThresholdPrimaryScore, std::vector<std::vector<double> > event_data,
		         std::vector<std::vector<double> > encoder_tracks, std::vector<std::vector<double> > decoder_tracks,
                         tensorflow::SavedModelBundleLite& lstm_model_bundle,
			 std::vector<int>& primary_track_list, std::vector<std::vector<std::vector<double> > >& primary_scores, 
			 std::vector<double>& bigger_primary_scores){
  std::vector<std::vector<double> > primary_event_data = event_data;
  sort(primary_event_data.begin(), primary_event_data.end(), [](const std::vector<double> &alpha, const std::vector<double> &beta){
    return alpha.at(60) > beta.at(60);
  });

  if(debug==true){
    for(int i=0; i<MaxPrimaryVertexLoop; i++){
      std::cout << "Track 1 " << primary_event_data.at(i).at(1) << " Track 2 " << primary_event_data.at(i).at(2) << " Score " << primary_event_data.at(i).at(60) << std::endl;
    }
  }

  std::vector<std::vector<double> > primary_pairs = VertexFinderwithDL::SliceN2DVector(primary_event_data, 0, MaxPrimaryVertexLoop, 3, 47);
  std::vector<std::vector<std::vector<double> > > primary_encoder_tracks(primary_pairs.size(), encoder_tracks);
  std::vector<std::vector<std::vector<double> > > primary_decoder_tracks(primary_pairs.size(), decoder_tracks);

  tensorflow::Tensor tprimary_pairs = VertexFinderwithDL::N2DVector2Tensor(primary_pairs);
  tensorflow::Tensor tprimary_encoder_tracks = VertexFinderwithDL::N3DVector2Tensor(primary_encoder_tracks);
  tensorflow::Tensor tprimary_decoder_tracks = VertexFinderwithDL::N3DVector2Tensor(primary_decoder_tracks);

  std::vector<tensorflow::Tensor> tmpprimary_outputs;
  tensorflow::Status runStatus = lstm_model_bundle.GetSession()->Run({{"serving_default_Decoder_Input:0", tprimary_decoder_tracks},
                                                                      {"serving_default_Encoder_Input:0", tprimary_encoder_tracks},
                                                                      {"serving_default_Pair_Input:0", tprimary_pairs}},
                                                                     {"StatefulPartitionedCall:0"},
							     	     {}, &tmpprimary_outputs);

  primary_scores = VertexFinderwithDL::Tensor2N3DVector(tmpprimary_outputs[0]);
  for(std::size_t i=0; i<primary_scores.size(); i++){
    for(std::size_t j=0; j<primary_scores.at(0).size(); j++){
      if(primary_scores.at(i).at(j).at(0) > ThresholdPrimaryScore){
	if(debug==true) std::cout << "Track " << j << " Primary Score: " << primary_scores.at(i).at(j).at(0) << std::endl;
	primary_track_list.push_back(j);
      }
    }
  }
  std::sort(primary_track_list.begin(), primary_track_list.end());
  primary_track_list.erase(std::unique(primary_track_list.begin(), primary_track_list.end()), primary_track_list.end());
  
  for(std::size_t i=0; i<encoder_tracks.size(); i++){
    double tmpbigger_primary_scores = 0;
    for(int j=0; j<MaxPrimaryVertexLoop; j++){
      if(tmpbigger_primary_scores < primary_scores.at(j).at(i).at(0)) tmpbigger_primary_scores = primary_scores.at(j).at(i).at(0);
    }
    bigger_primary_scores.push_back(tmpbigger_primary_scores);
  }
}


void VertexFinderwithDL::SecondaryVertexFinder(bool debug, double ThresholdSecondaryScore, std::vector<double> primary_scores,
		           std::vector<std::vector<double> > secondary_event_data,
		           std::vector<std::vector<double> > encoder_tracks, std::vector<std::vector<double> > decoder_tracks,
                           tensorflow::SavedModelBundleLite& slstm_model_bundle,
			   std::vector<int>& primary_track_list, std::vector<std::vector<int> >& secondary_track_lists){
  std::vector<int> track_list(decoder_tracks.size());
  std::iota(track_list.begin(), track_list.end(), 0);
  for(std::size_t i=0; i<secondary_event_data.size(); i++){
    int track1 = secondary_event_data.at(i).at(1), track2 = secondary_event_data.at(i).at(2);
    if(std::find(primary_track_list.begin(), primary_track_list.end(), track1) != primary_track_list.end() or
       std::find(primary_track_list.begin(), primary_track_list.end(), track2) != primary_track_list.end()) continue;
    if(debug==true){
      std::cout << "Track List" << std::endl;
      for(std::size_t i=0; i<track_list.size(); i++){
        std::cout << track_list.at(i) << " ";
      }
      std::cout << std::endl;

      std::cout << "seed: " << track1 << " " << track2 << std::endl;
    }

    std::vector<std::vector<double> > remain_decoder_tracks = VertexFinderwithDL::GetRemainDecoderTracks(decoder_tracks, track_list);

    std::vector<std::vector<double> > secondary_pair = VertexFinderwithDL::SliceN2DVector(secondary_event_data, i, i+1, 3, 47);
    std::vector<std::vector<std::vector<double> > > secondary_encoder_tracks(1, encoder_tracks);
    std::vector<std::vector<std::vector<double> > > secondary_remain_decoder_tracks(1, remain_decoder_tracks);

    tensorflow::Tensor tsecondary_pair = VertexFinderwithDL::N2DVector2Tensor(secondary_pair);
    tensorflow::Tensor tsecondary_encoder_tracks = VertexFinderwithDL::N3DVector2Tensor(secondary_encoder_tracks);
    tensorflow::Tensor tsecondary_remain_decoder_tracks = VertexFinderwithDL::N3DVector2Tensor(secondary_remain_decoder_tracks);

    std::vector<tensorflow::Tensor> tmpsecondary_outputs;
    tensorflow::Status runStatus = slstm_model_bundle.GetSession()->Run({{"serving_default_Decoder_Input:0", tsecondary_remain_decoder_tracks},
		                                                         {"serving_default_Encoder_Input:0", tsecondary_encoder_tracks},
						                         {"serving_default_Pair_Input:0", tsecondary_pair}},
				          	  	                {"StatefulPartitionedCall:0"},
									{}, &tmpsecondary_outputs);

    std::vector<std::vector<std::vector<double> > > secondary_scores = VertexFinderwithDL::Tensor2N3DVector(tmpsecondary_outputs[0]);
    if(debug==true) VertexFinderwithDL::DebugPrintVLSTMPrediction(secondary_scores);

    std::vector<int> tmpsecondary_track_list, tmptrack_list = track_list;
    for(std::size_t j=0; j<track_list.size(); j++){
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

std::vector<std::vector<int> > VertexFinderwithDL::MergeSingleTrack(TrackVec& tracks, std::vector<std::vector<int> > secondary_track_lists){
  std::vector<std::vector<int> > merge_secondary_track_lists, tmp_secondary_track_lists = secondary_track_lists;
  std::vector<int> tmp_single_track_list, tmp_not_single_track_list;
  std::vector<const Track*> tmp_tracks(2);
  double chi;
  int merge_track_num;

  for(std::size_t i=0; i<secondary_track_lists.size(); i++){
    if(secondary_track_lists.at(i).size()==1){
      tmp_single_track_list.push_back(secondary_track_lists.at(i).at(0));
      tmp_secondary_track_lists.erase(std::remove(tmp_secondary_track_lists.begin(), tmp_secondary_track_lists.end(), secondary_track_lists.at(i)), tmp_secondary_track_lists.end());
    }
    else{
      for(std::size_t j=0; j<secondary_track_lists.at(i).size(); j++){
        tmp_not_single_track_list.push_back(secondary_track_lists.at(i).at(j));
      }
    }
  }
  for(std::size_t i=0; i<tmp_single_track_list.size(); i++){
    tmp_tracks.at(0) = tracks.at(tmp_single_track_list.at(i));
    for(std::size_t j=0; j<tmp_not_single_track_list.size(); j++){
      tmp_tracks.at(1) = tracks.at(tmp_not_single_track_list.at(j));
      Vertex* vtx = VertexFitterSimple_V()(tmp_tracks.begin(), tmp_tracks.end());
      if(j==0){
        chi = vtx->getChi2();
	merge_track_num = tmp_not_single_track_list.at(j);
      } 
      else if(chi>vtx->getChi2()){
        chi = vtx->getChi2();
	merge_track_num = tmp_not_single_track_list.at(j);
      }
    }
    for(std::size_t j=0; j<tmp_secondary_track_lists.size(); j++){
      if(std::find(tmp_secondary_track_lists.at(j).begin(), tmp_secondary_track_lists.at(j).end(), merge_track_num) != tmp_secondary_track_lists.at(j).end()){
        tmp_secondary_track_lists.at(j).push_back(tmp_single_track_list.at(i));
      }
    }
  }
  merge_secondary_track_lists = tmp_secondary_track_lists;
  return merge_secondary_track_lists;
}

void VertexFinderwithDL::PrimarySecondaryVertices(TrackVec& tracks, std::vector<int> primary_track_list, std::vector<std::vector<int> > secondary_track_lists,
		                                  Vertex& vtx, std::vector<Vertex*>& vtces){
  //TrackVec& primary_tracks;
  std::vector<const Track*> primary_tracks;
  for(std::size_t i=0; i<primary_track_list.size(); i++){
    primary_tracks.push_back(tracks.at(primary_track_list.at(i)));
  }
  Vertex* tmpvtx = VertexFitterSimple_V()(primary_tracks.begin(), primary_tracks.end());
  vtx = *tmpvtx;
  for(std::size_t i=0; i<secondary_track_lists.size(); i++){
    //TrackVec& tmpsecondary_tracks;
    std::vector<const Track*> tmpsecondary_tracks;
    for(std::size_t j=0; j<secondary_track_lists.at(i).size(); j++){
      tmpsecondary_tracks.push_back(tracks.at(secondary_track_lists.at(i).at(j)));
    }
    Vertex* tmpsecondary_vtx = VertexFitterSimple_V()(tmpsecondary_tracks.begin(), tmpsecondary_tracks.end());
    vtces.push_back(tmpsecondary_vtx);
  }
}

