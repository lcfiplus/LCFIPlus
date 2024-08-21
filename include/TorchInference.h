// DNNProvider2.h

#ifndef torchinference_h
#define torchinference_h 1

#include "lcfiplus.h"

namespace lcfiplus {

class TorchInference : public Algorithm {
 public:
  TorchInference() {}
  virtual ~TorchInference() {}

  void init(Parameters* param);
  void process();
  void end();
 private:
  ClassDef(TorchInference,1);
};

}

#endif
