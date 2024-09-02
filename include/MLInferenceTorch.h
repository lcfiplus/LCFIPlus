// MLInferenceTorch.h

#ifndef MLInferenceTorch_h
#define MLInferenceTorch_h 1

#include "lcfiplus.h"

namespace lcfiplus {

class MLInferenceTorch : public Algorithm {
 public:
  MLInferenceTorch() {}
  virtual ~MLInferenceTorch() {}

  void init(Parameters* param);
  void process();
  void end();
 private:
  ClassDef(MLInferenceTorch,1);
};

}

#endif
