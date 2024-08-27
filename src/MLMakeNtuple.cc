#include "MLInputGenerator.h"
#include "MLMakeNtuple.h"

#include <assert.h>
#include "EventStore.h"
#include "lcfiplus.h"
#include "JetFinder.h"

#include "TROOT.h"
#include <string>
#include <typeinfo>

using namespace lcfiplus;
namespace lcfiplus {

using namespace MLInputGenerator;

struct MLData {
  vector<float> data;
  vector< vector<float> > dataVec;

  MLData() {
    data.reserve(20);
    dataVec.reserve(100);
  }

  float& newData() {
    float d;
    data.push_back(d);
    return data.back();
  }

  vector<float>& newDataVec() {
    vector<float> vec;
    dataVec.push_back(vec);
    return dataVec.back();
  }
};

// Relation to DNNProvider2:
//   It does the same thing but cleaner implementation
//
// Relation to MakeNtuple + FTManager:
//   It is a modernized implementation for input variable output

void MLMakeNtuple::init(Parameters* param) {
  Algorithm::init(param);
  string outputFilename = param->get("MLMakeNtuple.OutputRootFileName",string("ml_flavtag.root"));
  string treeName = param->get("MLMakeNtuple.TTreeName",string("ntp"));
  
  //cout << "MLMakeNtuple: Ntuple file set to " << outputFilename << endl;
  _file = new TFile(outputFilename.c_str(),"RECREATE");
  
  //cout << "MLMakeNtuple: setting TTree name '" << treeName.c_str() << "'" << endl;
  _tree = new TTree(treeName.c_str(),"ML input data for flavor tagging");

  MLInputGenerator::init();

  MLData d;

  for (auto v : calcInput) {
    auto key = v.first;
    _tree->Branch( key.c_str(), &d.newDataVec() );
  }

#if 0
  char buf[1024];
  for ( vector<FTAlgo*>::iterator iter = _algoList.begin(); iter != _algoList.end(); ++iter ) {
    FTAlgo* algo = *iter;
    snprintf( buf, 1000, "%s/F", algo->getName().c_str() );
    _tree->Branch( algo->getName().c_str(), algo->getValueAddress(), buf );
#endif
}

void MLMakeNtuple::process() {

}

void MLMakeNtuple::end() {
  _file->Write();
  _file->Close();
}

}
