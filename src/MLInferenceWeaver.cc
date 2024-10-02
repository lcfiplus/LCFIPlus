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
  string privtx = param->get("PrimaryVertexCollectionName",string("PrimaryVertex"));
  Event::Instance()->setDefaultPrimaryVertex(privtx.c_str());
  _jsonFileName = param->get("MLInferenceWeaver.JsonFileName",string("preprocess.json"));
  _onnxFileName = param->get("MLInferenceWeaver.OnnxFileName",string("test.onnx"));
  std::vector<string> varNames;
  param->fetchArray("MLInferenceWeaver.Variables", varNames);
  _variables = varNames; // needed to convert to RVec<string>
  _weaver = new WeaverInterface(_onnxFileName, _jsonFileName, _variables);
  MLInputGenerator::init();
}

void MLInferenceWeaver::process() {

  Event* event = Event::Instance();
  const Vertex* privtx = Event::Instance()->getPrimaryVertex();
  JetVec* jetsPtr(0);
  bool success = event->Get(_jetCollectionName.c_str(), jetsPtr);
  if (!success) {
    cout << "jets could not be found" << endl;
    return;
  }
  JetVec& jets = *jetsPtr;

  // output vector to hold the weights for every jet
  rv::RVec< rv::RVec<float> > out;

  // loop over jets
  for (unsigned int njet = 0; njet < jets.size(); ++njet) {
    const Jet* jet = jets[njet];
    TrackVec tracks = jet->getAllTracks(true);
    NeutralVec neutrals = jet->getNeutrals();

    // prepare input
    rv::RVec< rv::RVec<float> > vars;
    for (size_t i=0; i<_variables.size(); ++i) {
      
      // hold computed variable for all the candidates in a jet
      rv::RVec<float> cand_vars;

      // find function provided by MLInputGenerator; first check if it exists
      const auto& name = _variables[i];
      if (calcInput.find(name) == calcInput.end()) {
        cerr << "MLInferenceWeaver: no function to compute input variable " << name << endl;
        exit(1);
      }
      const auto& func = calcInput[name];
      
      for (auto tr: tracks) {
      
        float val(0.);

        if (std::holds_alternative<function<double(const Jet*)> >(func)) {
          auto f = std::get<function<double(const Jet*)> >(func);
          val = f(jet);
        }

        if (std::holds_alternative<function<double(const Track*)> >(func)) {
          auto f = std::get<function<double(const Track*)> >(func);
          val = f(tr);
        }

        if (std::holds_alternative<function<double(const Track*, const Jet*)> >(func)) {
          auto f = std::get<function<double(const Track*, const Jet*)> >(func);
          val = f(tr,jet);
        }

        if (std::holds_alternative<function<double(const Track*, const Vertex*)> >(func)) {
          auto f = std::get<function<double(const Track*, const Vertex*)> >(func);
          val = f(tr,privtx);
        }

        /*
        if (std::holds_alternative<function<double(const Neutral*)> >(func)) {
          auto f = std::get<function<double(const Neutral*)> >(func);
          val = f(neu);
        }

        if (std::holds_alternative<function<double(const Neutral*, const Jet*)> >(func)) {
          auto f = std::get<function<double(const Neutral*, const Jet*)> >(func);
          val = f(neu,jet);
        }

        if (std::holds_alternative<function<double(const Neutral*, const Vertex*)> >(func)) {
          auto f = std::get<function<double(const Neutral*, const Vertex*)> >(func);
          val = f(neu,privtx);
        }
        */
       
        cand_vars.push_back(val);
      }

      vars.push_back(cand_vars);
    }
    out.emplace_back(_weaver->run(vars));
  }
  std::cout << out << std::endl;
}

void MLInferenceWeaver::end() {

}