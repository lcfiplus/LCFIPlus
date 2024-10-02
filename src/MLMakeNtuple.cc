#include "MLInputGenerator.h"
#include "MLMakeNtuple.h"

#include <assert.h>
#include "EventStore.h"
#include "lcfiplus.h"
#include "JetFinder.h"

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
  
  //cout << "MLMakeNtuple: Ntuple file set to " << outputFilename << endl;
  _file = new TFile(outputFilename.c_str(),"RECREATE");
  
  //cout << "MLMakeNtuple: setting TTree name '" << treeName.c_str() << "'" << endl;
  _tree = new TTree(treeName.c_str(),"ML input data for flavor tagging");

  MLInputGenerator::init();

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

  for (const auto& v: calcInput) {
    const auto& key = v.first;
    if (std::holds_alternative<function<double(const Jet*)> >(v.second)) {
      string buf = key + "/F";
      _tree->Branch( key.c_str(), &_data.newData(key), buf.c_str() );
    } else {
      _tree->Branch( key.c_str(), &_data.newDataVec(key) );
    }
  }
}

//template<typename ...Ts> struct overloaded : Ts... { using Ts::operator()...; };
//template<typename ...Ts> overloaded(Ts...) -> overloaded<Ts...>;

void MLMakeNtuple::process() {
  Event* event = Event::Instance();
  const Vertex* privtx = Event::Instance()->getPrimaryVertex();
  JetVec* jetsPtr(0);
  bool success = event->Get(_jetname.c_str(), jetsPtr);
  if (!success) {
    cout << "jets could not be found" << endl;
    return;
  }
  JetVec& jets = *jetsPtr;

  for (unsigned int njet = 0; njet < jets.size(); ++njet) {
    _data.resetData();
    const Jet* jet = jets[njet];
    TrackVec tracks = jet->getAllTracks(true);
    NeutralVec neutrals = jet->getNeutrals();

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
