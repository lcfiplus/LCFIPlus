#include <functional>
#include <map>
#include <string>
#include <sstream>

#include "lcfiplus.h"
#include "process.h"
#include "MLInputGenerator.h"
#include "MLInferenceWeaver.h"
#include "WeaverInterface.h"

using namespace lcfiplus;
using namespace MLInputGenerator;

void MLInferenceWeaver::init(Parameters* param){
  Algorithm::init(param);
  _jetCollectionName = param->get("MLInferenceWeaver.JetCollectionName",string("RefinedJets"));
  _jsonFileName = param->get("MLInferenceWeaver.JsonFileName",string("preprocess.json"));
  _onnxFileName = param->get("MLInferenceWeaver.OnnxFileName",string("test.onnx"));
  std::vector<string> varNames;
  param->fetchArray("MLInferenceWeaver.Variables", varNames);
  _variables = varNames; // convert to RVec<string>
  _weaver = new WeaverInterface(_onnxFileName, _jsonFileName, _variables);
}

void MLInferenceWeaver::process() {
  rv::RVec< rv::RVec<float> > out;
  rv::RVec<float> vars;

  for (size_t i=0; i<_variables.size(); ++i) {
    vars.push_back(0.1);
  }
  out.emplace_back(_weaver->run(vars));

  std::cout << out << std::endl;
}

void MLInferenceWeaver::end() {

}