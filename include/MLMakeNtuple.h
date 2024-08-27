// MakeNtuple.h

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
  TFile* _file;
  TTree* _tree;

  ClassDef(MLMakeNtuple,1);
};

}

#endif
