#include <functional>
#include <map>
#include <string>
#include <sstream>

#include "lcfiplus.h"
#include "process.h"
#include "MLInputGenerator.h"
#include "MLInferenceTorch.h"

using namespace lcfiplus;

void MLInferenceTorch::init(Parameters* param){
    Algorithm::init(param);

    string _torchScriptFileName = "/data/suehara/part/training/unprocessed/ilc_nnqq_default/net_best_epoch_state.pt";

    try {
        //torch::jit::script::Module _model;
        _model = torch::jit::load(_torchScriptFileName);
        _model.to(torch::kCPU);
        _model.eval();
    } catch (const c10::Error& e) {
        std::stringstream message;
        message << "Error loading the model \'" << _torchScriptFileName << "\': " << e.what();
        //streamlog_out(ERROR) << message.str() << std::endl;
        //throw EVENT::Exception(message.str());
        std::cerr << message.str() << std::endl;
        auto x = message.str();
        throw(Exception(x.c_str()));
    }

}

void MLInferenceTorch::process() {

}

void MLInferenceTorch::end() {

}