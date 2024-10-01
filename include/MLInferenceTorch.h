// MLInferenceTorch.h

#ifndef MLInferenceTorch_h
#define MLInferenceTorch_h 1

#include "lcfiplus.h"

#undef ClassDef
#include <torch/script.h>
//namespace torch::jit::script { struct Module; }
//namespace torch::jit { struct Module; }

namespace lcfiplus {

class MLInferenceTorch : public Algorithm {
 public:
  MLInferenceTorch() {}
  virtual ~MLInferenceTorch() {}

  void init(Parameters* param);
  void process();
  void end();
  
 private:
  torch::jit::script::Module _model;
  //torch::jit::script::Module* _model; //!
  //ClassDef(MLInferenceTorch,1);

};

}

#endif
