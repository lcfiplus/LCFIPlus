#ifndef VertexMassRecovery_h
#define VertexMassRecovery_h 1

//pi0 vertex finder
#include "Pi0VertexFinder.h"

namespace lcfiplus {

  //class Pi0VertexFinder;

class VertexMassRecovery : public Algorithm {
 public:	  
  VertexMassRecovery(){}
  virtual ~VertexMassRecovery(){}
	  
  void init(Parameters *param);
  void process();
  void end();

 private:
  Pi0VertexFinder *pi0vtxfinder;
  
  JetVec *_inputJets;
  VertexVec *_invertices;
  vector<Vertex *> * _outvertices; //!
  vector<Jet *> * _outputJets; //!
  
  string _jincolname;
  string _vincolname;
  //string _vv0colname;
  string _vprimcolname;
  
  ClassDef(VertexMassRecovery, 1);
};

}

#endif
