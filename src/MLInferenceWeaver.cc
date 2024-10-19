#include <functional>
#include <map>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>

#include "lcfiplus.h"
#include "process.h"
#include "MLInputGenerator.h"
#include "MLInferenceWeaver.h"
#include "WeaverInterface.h"
#include "nlohmann/json.hpp"

using namespace lcfiplus;
using namespace MLInputGenerator;

void MLInferenceWeaver::parseJSON(const string& json_filename) {
  std::ifstream json_file(json_filename);
  std::vector<std::string> input_names;
  std::vector<std::string> output_names;
  try {
    const auto json = nlohmann::json::parse(json_file);

    // process output names
    json.at("output_names").get_to(output_names);
    for (const auto& output : output_names) {
      _outputVariables.emplace_back(output);
    }

    // process input names
    json.at("input_names").get_to(input_names);
    for (const auto& input : input_names) {

      // skip if string ends with _mask
      const string ending = "_mask";
      if (input.compare(input.size()-ending.size(), ending.size(), ending) == 0) {
        continue;
      }

      const auto& group_params = json.at(input);
      auto& info = prep_info_map_[input];
      info.name = input;
      // retrieve the variables names
      group_params.at("var_names").get_to(info.var_names);
      for (const auto& v : info.var_names) {
        _variables.emplace_back(v);
      }
    }
  } catch (const nlohmann::json::exception& exc) {
    throw std::runtime_error("Failed to parse input JSON file '" + json_filename + "'.\n" + exc.what());
  }
}

void MLInferenceWeaver::init(Parameters* param) {
  Algorithm::init(param);
  _jetCollectionName = param->get("MLInferenceWeaver.JetCollectionName",string("RefinedJets"));
  string privtx = param->get("PrimaryVertexCollectionName",string("PrimaryVertex"));
  Event::Instance()->setDefaultPrimaryVertex(privtx.c_str());
  _jsonFileName = param->get("MLInferenceWeaver.JsonFileName",string("preprocess.json"));
  _onnxFileName = param->get("MLInferenceWeaver.OnnxFileName",string("test.onnx"));
  _outputParameterName = param->get("MLInferenceWeaver.OutputParameterName",string("weaver"));

  // event classification (true) vs. jet classification (false; default)
  _eventClassification = param->get("MLInferenceWeaver.EventClassification",false);
  
  parseJSON(_jsonFileName);
  _weaver = new WeaverInterface(_onnxFileName, _jsonFileName, _variables);
  MLInputGenerator::init();
}

void MLInferenceWeaver::process() {
  if (_eventClassification) {
    processEvent();
  } else {
    processJet();
  }
}

void MLInferenceWeaver::processJet() {

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
    TrackVec tracks = jet->getAllTracksSorted(true);
    NeutralVec neutrals = jet->getNeutralsSorted();

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

      if (std::holds_alternative<function<double(const Jet*)> >(func)) {
        auto f = std::get<function<double(const Jet*)> >(func);
        cand_vars.push_back( f(jet) );
      }
      else if (std::holds_alternative<function<double(const Track*)> >(func)) {
        auto f = std::get<function<double(const Track*)> >(func);
        for (const auto& tr: tracks) {
          cand_vars.push_back( f(tr) );
        }
      }
      else if (std::holds_alternative<function<double(const Track*, const Jet*)> >(func)) {
        auto f = std::get<function<double(const Track*, const Jet*)> >(func);
        for (const auto& tr: tracks) {
          cand_vars.push_back( f(tr,jet) );
        }
      }
      else if (std::holds_alternative<function<double(const Track*, const Vertex*)> >(func)) {
        auto f = std::get<function<double(const Track*, const Vertex*)> >(func);
        for (const auto& tr: tracks) {
          cand_vars.push_back( f(tr,privtx) );
        }
      }
      else if (std::holds_alternative<function<double(const Neutral*)> >(func)) {
        auto f = std::get<function<double(const Neutral*)> >(func);
        for (const auto& neu: neutrals) {
          cand_vars.push_back( f(neu) );
        }
      }
      else if (std::holds_alternative<function<double(const Neutral*, const Jet*)> >(func)) {
        auto f = std::get<function<double(const Neutral*, const Jet*)> >(func);
        for (const auto& neu: neutrals) {
          cand_vars.push_back( f(neu,jet) );
        }
      }
      else if (std::holds_alternative<function<double(const Neutral*, const Vertex*)> >(func)) {
        auto f = std::get<function<double(const Neutral*, const Vertex*)> >(func);
        for (const auto& neu: neutrals) {
          cand_vars.push_back( f(neu,privtx) );
        }
      }

      vars.push_back(cand_vars);
    }

    rv::RVec<float> res = _weaver->run(vars);
    out.emplace_back(res);

    assert(res.size() == _outputVariables.size());

    Parameters param;
    for (size_t i=0; i<res.size(); ++i) {
      //param.add( _outputVariables[i].data(), res[i] );
      param.add( _outputVariables[i].data(), (double)res[i] );
    }
    //param.add( "Category", (double)category );
    jet->addParam( _outputParameterName.data(), param );
  }

}

void MLInferenceWeaver::processEvent() {

  Event* event = Event::Instance();
  const Vertex* privtx = Event::Instance()->getPrimaryVertex();

  // output vector to hold the weights for every jet
  rv::RVec< rv::RVec<float> > out;

  // loop over jets
  TrackVec tracks = event->getTracks();
  NeutralVec neutrals = event->getNeutrals();

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

    if (std::holds_alternative<function<double(const Track*)> >(func)) {
      auto f = std::get<function<double(const Track*)> >(func);
      for (const auto& tr: tracks) {
        cand_vars.push_back( f(tr) );
      }
    }
    else if (std::holds_alternative<function<double(const Track*, const Vertex*)> >(func)) {
      auto f = std::get<function<double(const Track*, const Vertex*)> >(func);
      for (const auto& tr: tracks) {
        cand_vars.push_back( f(tr,privtx) );
      }
    }
    else if (std::holds_alternative<function<double(const Neutral*)> >(func)) {
      auto f = std::get<function<double(const Neutral*)> >(func);
      for (const auto& neu: neutrals) {
        cand_vars.push_back( f(neu) );
      }
    }
    else if (std::holds_alternative<function<double(const Neutral*, const Vertex*)> >(func)) {
      auto f = std::get<function<double(const Neutral*, const Vertex*)> >(func);
      for (const auto& neu: neutrals) {
        cand_vars.push_back( f(neu,privtx) );
      }
    }

    vars.push_back(cand_vars);
  }

  rv::RVec<float> res = _weaver->run(vars);
  out.emplace_back(res);

  assert(res.size() == _outputVariables.size());

  Parameters param;
  for (size_t i=0; i<res.size(); ++i) {
    //param.add( _outputVariables[i].data(), res[i] );
    param.add( _outputVariables[i].data(), (double)res[i] );
  }
  //param.add( "Category", (double)category );

  // TODO: implement add parameter output for Event class  
  //event->addParam( _outputParameterName.data(), param );

}

void MLInferenceWeaver::end() {
  // do nothing
}
