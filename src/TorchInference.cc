#include <string>
#include <sstream>

#include "lcfiplus.h"
#include "process.h"
#include "MLInputGenerator.h"
#include "TorchInference.h"

#undef ClassDef
#include <torch/script.h>

#include <map>
#include <functional>

using namespace lcfiplus;

void TorchInference::init(Parameters* param){
    Algorithm::init(param);

    auto _torchScriptFileName = "/data/suehara/part/training/unprocessed/ilc_nnqq_default/net_best_epoch_state.pt";

    try {
        auto _model = torch::jit::load(_torchScriptFileName);
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

void TorchInference::process() {

}

void TorchInference::end() {

}