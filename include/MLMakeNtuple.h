// MakeNtuple.h

// Relation to DNNProvider2:
//   It does the same thing but cleaner implementation
//
// Relation to MakeNtuple + FTManager:
//   It is a modernized implementation for input variable output

#ifndef MLMakeNtuple_h
#define MLMakeNtuple_h 1

class TFile;
class TTree;

#include "lcfiplus.h"

namespace lcfiplus {

class MLMakeNtuple : public Algorithm {
 public:
  MLMakeNtuple() {}
  virtual ~MLMakeNtuple() {}

  void init(Parameters* param);
  void process();
  void end();

 private:
  string _jetname;
  TFile* _file;
  TTree* _tree;
  int _label;
  int _labelKeep;

  struct MLData {
    // methods
    MLData() {};
    float& newData(const string& key);
    vector<float>& newDataVec(const string& key);

    void setData(const string& key, const float& val);
    void addDataVec(const string& key, const float& val);

    void resetData();

    // private members
    map<string,float> _mapData;
    map<string,vector<float> > _mapDataVec;
  } _data;

  ClassDef(MLMakeNtuple,1);
};

}

#endif
