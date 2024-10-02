// MLInferenceWeaver.h

#ifndef MLInferenceWeaver_h
#define MLInferenceWeaver_h 1

#include "ROOT/RVec.hxx"
#include "lcfiplus.h"

namespace rv = ROOT::VecOps;
namespace lcfiplus {

class WeaverInterface;

class MLInferenceWeaver : public Algorithm {
 public:
  MLInferenceWeaver() {}
  virtual ~MLInferenceWeaver() {}

  void init(Parameters* param);
  void process();
  void end();
  
 private:
  std::string _jetCollectionName;
  std::string _jsonFileName;
  std::string _onnxFileName;
  rv::RVec<std::string> _variables; //!
  WeaverInterface* _weaver; //!

  ClassDef(MLInferenceWeaver,1);
};

}

#endif
