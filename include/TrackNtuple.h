// TrackNtuple.h

#ifndef TrackNtuple_h
#define TrackNtuple_h 1

class TFile;
class TNtuple;

#include "lcfiplus.h"

namespace lcfiplus {

/**
	Making track d0/z0 ntuple needed for flavor tagging.
	Should be rerun when vertex configuration is changed.
	Run TrackProb.C with hadded output files of this algorithm.

	@author T. Suehara, Dept. of Physics, Kyushu University
	@version $Id$
 */

class TrackNtuple : public Algorithm {
 public:
  TrackNtuple() {}
  virtual ~TrackNtuple() {}

  void init(Parameters* param);
  void process();
  void end();

 private:
  TFile* _file;
  TNtuple* _tree;

  string _jetcolname;
  string _primvtxcolname;
  int _hitcutJprob;

  ClassDef(TrackNtuple,1);
};

}

#endif
