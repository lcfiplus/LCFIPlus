// VertexNtuple.h

#ifndef VertexNtuple_h
#define VertexNtuple_h 1

class TFile;
class TNtuple;

#include "lcfiplus.h"

namespace lcfiplus {

/**
	Storing the number of tracks associated to vertices..

	@author R.Yonamine, Dept. of Physics, Tohok University
	@version $Id$
 */

class VertexNtuple : public Algorithm {
 public:
  VertexNtuple() {}
  virtual ~VertexNtuple() {}

  void init(Parameters* param);
  void process();
  void end();

 private:
  TFile* _file;
  TNtuple* _tree;

  string _jetcolname;
  string _primvtxcolname;

  ClassDef(VertexNtuple,1)
};

}

#endif
