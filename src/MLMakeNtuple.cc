#include "MLInputGenerator.h"
#include "MLMakeNtuple.h"

#include <assert.h>
#include "EventStore.h"
#include "lcfiplus.h"
#include "JetFinder.h"
#include "algoEtc.h"

#include "TROOT.h"
#include <variant>
#include <string>
#include <typeinfo>

using namespace lcfiplus;
using namespace MLInputGenerator;

namespace lcfiplus {

float& MLMakeNtuple::MLData::newData(const string& key) {
  auto const& inserted = _mapData.insert({key, float()});
  auto const& iter = inserted.first;
  return iter->second;
}

vector<float>& MLMakeNtuple::MLData::newDataVec(const string& key) {
  auto const& inserted = _mapDataVec.insert({key, vector<float>()});
  auto const& iter = inserted.first;
  return iter->second;
}

void MLMakeNtuple::MLData::setData(const string& key, const float& val) {
  _mapData[key] = val;
}

void MLMakeNtuple::MLData::addDataVec(const string& key, const float& val) {
  _mapDataVec[key].push_back(val);
}

void MLMakeNtuple::MLData::resetData() {
  for (auto& d: _mapData) d.second = 0;
  for (auto& d: _mapDataVec) d.second.clear();
}

void MLMakeNtuple::init(Parameters* param) {
  Algorithm::init(param);
  _jetname = param->get("MLMakeNtuple.JetCollectionName",string("RefinedJets"));
  string privtx = param->get("PrimaryVertexCollectionName",string("PrimaryVertex"));
  Event::Instance()->setDefaultPrimaryVertex(privtx.c_str());

  string outputFilename = param->get("MLMakeNtuple.OutputRootFileName",string("ml_flavtag.root"));
  string treeName = param->get("MLMakeNtuple.TTreeName",string("ntp"));
  _labelKeep = param->get("MLMakeNtuple.Label",int(0));
  _outEvent = param->get("MLMakeNtuple.EventClassification",int(0));
  _outEventNoJets = param->get("MLMakeNtuple.EventClassificationNoJets",int(0));
  
  //cout << "MLMakeNtuple: Ntuple file set to " << outputFilename << endl;
  _file = new TFile(outputFilename.c_str(),"RECREATE");
  
  //cout << "MLMakeNtuple: setting TTree name '" << treeName.c_str() << "'" << endl;
  _tree = new TTree(treeName.c_str(),"ML input data for flavor tagging");

  MLInputGenerator::init();

  _tree->Branch("label",&_label,"Label/I");

  // TODO: implement labels
#if 0
  // label
  _tree->Branch("mc_b",&d.mc_b,"mc_b/I");
  _tree->Branch("mc_c",&d.mc_c,"mc_c/I");
  _tree->Branch("mc_u",&d.mc_u,"mc_u/I");
  _tree->Branch("mc_d",&d.mc_d,"mc_d/I");
  _tree->Branch("mc_s",&d.mc_s,"mc_s/I");
  _tree->Branch("mc_g",&d.mc_g,"mc_g/I");
  _tree->Branch("mc_q",&d.mc_q,"mc_q/I"); // usdg
#endif

  if (_outEvent){ 
    string key = "jetIndex";
    _tree->Branch( key.c_str(), &_data.newDataVec(key) );

    key = "neu_jetIndex";
    _tree->Branch( key.c_str(), &_data.newDataVec(key) );
  }

  for (const auto& v: calcInput) {
    const auto& key = v.first;

    if(_outEventNoJets){
      if (std::holds_alternative<function<double(const Track*)> >(v.second)
       || std::holds_alternative<function<double(const Track*, const Vertex*)> >(v.second)
       || std::holds_alternative<function<double(const Neutral*)> >(v.second)
       || std::holds_alternative<function<double(const Neutral*, const Vertex*)> >(v.second)){
	_tree->Branch( key.c_str(), &_data.newDataVec(key) );
      }
    }

    if (_outEvent){ 
      if (std::holds_alternative<function<double(const Track*)> >(v.second)
       || std::holds_alternative<function<double(const Track*, const Vertex*)> >(v.second)
       || std::holds_alternative<function<double(const Neutral*)> >(v.second)
       || std::holds_alternative<function<double(const Neutral*, const Vertex*)> >(v.second)
       || std::holds_alternative<function<double(const Track*, const Jet*)> >(v.second)
       || std::holds_alternative<function<double(const Neutral*, const Jet*)> >(v.second)
       )
      {
	_tree->Branch( key.c_str(), &_data.newDataVec(key) );
      }

    } else {

      if (std::holds_alternative<function<double(const Jet*)> >(v.second)) {
        string buf = key + "/F";
        _tree->Branch( key.c_str(), &_data.newData(key), buf.c_str() );
      } else {
        _tree->Branch( key.c_str(), &_data.newDataVec(key) );
      }
    }

  }
}

//template<typename ...Ts> struct overloaded : Ts... { using Ts::operator()...; };
//template<typename ...Ts> overloaded(Ts...) -> overloaded<Ts...>;

void MLMakeNtuple::process() {

  // only one of EventClassification or EventClassificationNoJets is allowed to be 1
  if (_outEvent && _outEventNoJets) {
    cout << "Skipping due to ambiguous setting: MLMakeNtuple.EventClassification and MLMakeNtuple.EventClassificationNoJets are both turned on" << endl;
    return;
  }

  if(_outEvent)
    processEvent();
  else if (_outEventNoJets)
    processEventNoJets();
  else
    processJets();
}

void MLMakeNtuple::processEventNoJets(){
  Event* event = Event::Instance();
  const Vertex* privtx = Event::Instance()->getPrimaryVertex();
  _label = _labelKeep;

  // get tracks and neutrals
  TrackVec &tracks_orig = event->getTracks();
  NeutralVec &neutrals_orig = event->getNeutrals();

  vector<const Track *> tracks(tracks_orig.size());
  vector<const Neutral *> neutrals(neutrals_orig.size());

  std::partial_sort_copy(tracks_orig.begin(),tracks_orig.end(),tracks.begin(), tracks.end(), [](const Track *a, const Track *b){
    return a->E() > b->E();
  });
  std::partial_sort_copy(neutrals_orig.begin(),neutrals_orig.end(),neutrals.begin(), neutrals.end(),[](const Neutral *a, const Neutral *b){
    return a->E() > b->E();
  });

  _data.resetData();
  for (const auto& v: calcInput) {
    const auto& key = v.first;

    if (std::holds_alternative<function<double(const Track*)> >(v.second)) {
      auto f = std::get<function<double(const Track*)> >(v.second);
      for (auto tr: tracks) {
	double ret = f(tr);
	_data.addDataVec(key,ret);
      }
    }

    if (std::holds_alternative<function<double(const Track*, const Vertex*)> >(v.second)) {
      auto f = std::get<function<double(const Track*, const Vertex*)> >(v.second);
      for (auto tr: tracks) {
	double ret = f(tr,privtx);
	_data.addDataVec(key,ret);
      }
    }
    
    if (std::holds_alternative<function<double(const Neutral*)> >(v.second)) {
      auto f = std::get<function<double(const Neutral*)> >(v.second);
      for (auto neu: neutrals) {
	double ret = f(neu);
	_data.addDataVec(key,ret);
      }
    }
    
    if (std::holds_alternative<function<double(const Neutral*, const Vertex*)> >(v.second)) {
      auto f = std::get<function<double(const Neutral*, const Vertex*)> >(v.second);
      for (auto neu: neutrals) {
	double ret = f(neu,privtx);
	_data.addDataVec(key,ret);
      }
    }
  }

  _tree->Fill();

}

void MLMakeNtuple::processEvent(){

  Event* event = Event::Instance();
  const Vertex* privtx = Event::Instance()->getPrimaryVertex();
  JetVec* jetsPtr(0);
  bool success = event->Get(_jetname.c_str(), jetsPtr);
  if (!success) {
    cout << "jets could not be found" << endl;
    return;
  }
  JetVec& jets = *jetsPtr;

  _data.resetData();
  _label = _labelKeep;

  for (unsigned int njet = 0; njet < jets.size(); ++njet) {
    
    const Jet* jet = jets[njet];
    TrackVec tracks = jet->getAllTracksSorted(true);
    NeutralVec neutrals = jet->getNeutralsSorted();

    for (auto tr: tracks) {
      (void) tr;
      _data.addDataVec("jetIndex",(float)njet);
    }
    for (auto neu: neutrals) {
      (void) neu;
      _data.addDataVec("neu_jetIndex",(float)njet);
    }

    for (const auto& v: calcInput) {
      const auto& key = v.first;

      /*
      if (std::holds_alternative<function<double(const Jet*)> >(v.second)) {
        auto f = std::get<function<double(const Jet*)> >(v.second);
        double ret = f(jet);
        _data.setData(key,ret);
      }
      */

      if (std::holds_alternative<function<double(const Track*)> >(v.second)) {
        auto f = std::get<function<double(const Track*)> >(v.second);
        for (auto tr: tracks) {
          double ret = f(tr);
          _data.addDataVec(key,ret);
        }
      }

      if (std::holds_alternative<function<double(const Track*, const Jet*)> >(v.second)) {
        auto f = std::get<function<double(const Track*, const Jet*)> >(v.second);
        for (auto tr: tracks) {
          double ret = f(tr,jet);
          _data.addDataVec(key,ret);
        }
      }

      if (std::holds_alternative<function<double(const Track*, const Vertex*)> >(v.second)) {
        auto f = std::get<function<double(const Track*, const Vertex*)> >(v.second);
        for (auto tr: tracks) {
          double ret = f(tr,privtx);
          _data.addDataVec(key,ret);
        }
      }

      if (std::holds_alternative<function<double(const Neutral*)> >(v.second)) {
        auto f = std::get<function<double(const Neutral*)> >(v.second);
        for (auto neu: neutrals) {
          double ret = f(neu);
          _data.addDataVec(key,ret);
        }
      }

      if (std::holds_alternative<function<double(const Neutral*, const Jet*)> >(v.second)) {
        auto f = std::get<function<double(const Neutral*, const Jet*)> >(v.second);
        for (auto neu: neutrals) {
          double ret = f(neu,jet);
          _data.addDataVec(key,ret);
        }
      }

      if (std::holds_alternative<function<double(const Neutral*, const Vertex*)> >(v.second)) {
        auto f = std::get<function<double(const Neutral*, const Vertex*)> >(v.second);
        for (auto neu: neutrals) {
          double ret = f(neu,privtx);
          _data.addDataVec(key,ret);
        }
      }
    }
  }

  _tree->Fill();

}
  
void MLMakeNtuple::processJets(){
  Event* event = Event::Instance();
  const Vertex* privtx = Event::Instance()->getPrimaryVertex();
  JetVec* jetsPtr(0);
  bool success = event->Get(_jetname.c_str(), jetsPtr);
  if (!success) {
    cout << "jets could not be found" << endl;
    return;
  }
  JetVec& jets = *jetsPtr;

  vector<const MCParticle *> mcpJets;
  if(_labelKeep == 0 && Event::Instance()->isMCExist()){
    algoEtc::AssignJetsToMC(jets, mcpJets);
  }
  for (unsigned int njet = 0; njet < jets.size(); ++njet) {
    _label = _labelKeep;
    if(_label == 0 && Event::Instance()->isMCExist()){
      _label = (mcpJets[njet] ? mcpJets[njet]->getPDG() : 0);
    }
    
    _data.resetData();
    const Jet* jet = jets[njet];
    TrackVec tracks = jet->getAllTracksSorted(true);
    NeutralVec neutrals = jet->getNeutralsSorted();

    for (const auto& v: calcInput) {
      const auto& key = v.first;

      if (std::holds_alternative<function<double(const Jet*)> >(v.second)) {
        auto f = std::get<function<double(const Jet*)> >(v.second);
        double ret = f(jet);
        _data.setData(key,ret);
      }

      if (std::holds_alternative<function<double(const Track*)> >(v.second)) {
        auto f = std::get<function<double(const Track*)> >(v.second);
        for (auto tr: tracks) {
          double ret = f(tr);
          _data.addDataVec(key,ret);
        }
      }

      if (std::holds_alternative<function<double(const Track*, const Jet*)> >(v.second)) {
        auto f = std::get<function<double(const Track*, const Jet*)> >(v.second);
        for (auto tr: tracks) {
          double ret = f(tr,jet);
          _data.addDataVec(key,ret);
        }
      }

      if (std::holds_alternative<function<double(const Track*, const Vertex*)> >(v.second)) {
        auto f = std::get<function<double(const Track*, const Vertex*)> >(v.second);
        for (auto tr: tracks) {
          double ret = f(tr,privtx);
          _data.addDataVec(key,ret);
        }
      }

      if (std::holds_alternative<function<double(const Neutral*)> >(v.second)) {
        auto f = std::get<function<double(const Neutral*)> >(v.second);
        for (auto neu: neutrals) {
          double ret = f(neu);
          _data.addDataVec(key,ret);
        }
      }

      if (std::holds_alternative<function<double(const Neutral*, const Jet*)> >(v.second)) {
        auto f = std::get<function<double(const Neutral*, const Jet*)> >(v.second);
        for (auto neu: neutrals) {
          double ret = f(neu,jet);
          _data.addDataVec(key,ret);
        }
      }

      if (std::holds_alternative<function<double(const Neutral*, const Vertex*)> >(v.second)) {
        auto f = std::get<function<double(const Neutral*, const Vertex*)> >(v.second);
        for (auto neu: neutrals) {
          double ret = f(neu,privtx);
          _data.addDataVec(key,ret);
        }
      }
    }
    _tree->Fill();
  }
}

void MLMakeNtuple::end() {
  _file->Write();
  _file->Close();
}

}
