// MakeNtuple.h

#ifndef MakeNtuple_h
#define MakeNtuple_h 1

class TFile;
class TTree;

#include "lcfiplus.h"

namespace lcfiplus {

/**
	Lcfiplus algorithm for computing variables, to be used in multivariate analysis.
	@author T. Tanabe, ICEPP, The University of Tokyo
	@version $Id$
 */
class MakeNtuple : public Algorithm {
 public:
  MakeNtuple() {}
  virtual ~MakeNtuple() {}

  void init(Parameters* param);
  void process() {}
  void end();

 private:

  ClassDef(MakeNtuple,1);
};

}

#endif
