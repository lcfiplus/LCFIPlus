#ifndef VertexMassRecovery_h
#define VertexMassRecovery_h 1

//pi0 vertex finder
//#include "Pi0VertexFinder.h"

namespace lcfiplus {

class Pi0VertexFinder;

class VertexMassRecovery : public Algorithm {
 public:	  
  VertexMassRecovery(){}
  virtual ~VertexMassRecovery(){}
	  
  void init(Parameters *param);
  void process();
  void end();

  ClassDef(VertexMassRecovery, 1);

 private:
  Pi0VertexFinder *pi0vtxfinder; //!
  
  JetVec *_inputJets; //!
  VertexVec *_invertices; //!
  std::vector<Vertex *> * _outvertices; //!
  std::vector<Jet *> * _outputJets; //!
  
  std::string _jincolname;
  std::string _vincolname;
  //string _vv0colname;
  std::string _vprimcolname;
  
};

}

#endif
