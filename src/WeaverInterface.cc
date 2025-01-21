// Code from https://github.com/HEP-FCC/FCCAnalyses

#include "WeaverInterface.h"

#include "nlohmann/json.hpp"
#include <fstream>
#include <iostream>

//#define DUMP_WEAVER_INPUT 1

using namespace lcfiplus;

WeaverInterface::WeaverInterface(const std::string& onnx_filename,
                                 const std::string& json_filename,
                                 const rv::RVec<std::string>& vars)
    : variables_names_(vars.begin(), vars.end()) {
  if (onnx_filename.empty())
    throw std::runtime_error("ONNX modeld input file not specified!");
  if (json_filename.empty())
    throw std::runtime_error("JSON preprocessed input file not specified!");

  // the preprocessing JSON was found ; extract the variables listing and all useful information
  std::ifstream json_file(json_filename);
  std::vector<std::string> input_names;
  try {
    const auto json = nlohmann::json::parse(json_file);
    json.at("input_names").get_to(input_names);
    for (const auto& input : input_names) {
      const auto& group_params = json.at(input);
      auto& info = prep_info_map_[input];
      info.name = input;
      // retrieve the variables names
      group_params.at("var_names").get_to(info.var_names);
      // retrieve the shapes for all variables
      if (group_params.contains("var_length")) {
        info.min_length = info.max_length = group_params.at("var_length");
        input_shapes_.push_back({1, (int64_t)info.var_names.size(), (int64_t)info.min_length});
      } else {
        info.min_length = group_params.at("min_length");
        info.max_length = group_params.at("max_length");
        input_shapes_.push_back({1, (int64_t)info.var_names.size(), -1});
      }
      // for all variables, retrieve the allowed range
      const auto& var_info_params = group_params.at("var_infos");
      for (const auto& name : info.var_names) {
        const auto& var_params = var_info_params.at(name);
        info.var_info_map[name] =
            PreprocessParams::VarInfo(var_params.at("median"),
                                      var_params.at("norm_factor"),
                                      var_params.at("replace_inf_value"),
                                      var_params.at("lower_bound"),
                                      var_params.at("upper_bound"),
                                      var_params.contains("pad") ? (double)var_params.at("pad") : 0.);
      }
      // create data storage with a fixed size vector initialised with 0's
      const auto& len = input_sizes_.emplace_back(info.max_length * info.var_names.size());
      data_.emplace_back(len, 0);
    }
  } catch (const nlohmann::json::exception& exc) {
    throw std::runtime_error("Failed to parse input JSON file '" + json_filename + "'.\n" + exc.what());
  }
  onnx_ = std::make_unique<ONNXRuntime>(onnx_filename, input_names);
}

std::vector<float> WeaverInterface::center_norm_pad(const rv::RVec<float>& input,
                                                    float center,
                                                    float scale,
                                                    size_t min_length,
                                                    size_t max_length,
                                                    float pad_value,
                                                    float replace_inf_value,
                                                    float min,
                                                    float max) {
  if (min > pad_value || pad_value > max)
    throw std::runtime_error("Pad value not within (min, max) range");
  if (min_length > max_length)
    throw std::runtime_error("Variable length mismatch (min_length >= max_length)");

  auto ensure_finitude = [](const float in, const float replace_val) {
    if (!std::isfinite(in))
      return replace_val;
    return in;
  };

  size_t target_length = std::clamp(input.size(), min_length, max_length);
  std::vector<float> out(target_length, pad_value);
  for (size_t i = 0; i < input.size() && i < target_length; ++i)
    out[i] = std::clamp((ensure_finitude(input[i], replace_inf_value) - center) * scale, min, max);
  return out;
}

size_t WeaverInterface::variablePos(const std::string& var_name) const {
  auto var_it = std::find(variables_names_.begin(), variables_names_.end(), var_name);
  if (var_it == variables_names_.end())
    throw std::runtime_error("Unable to find variable with name '" + var_name +
                             "' in the list of registered variables");
  return var_it - variables_names_.begin();
}

#ifdef DUMP_WEAVER_INPUT
void WeaverInterface::writeToFile(const rv::RVec< rv::RVec<float> >& vars) {
  static size_t cnt(0);
  ++cnt;
  char buf[1024];
  sprintf(buf,"pre%05ld.txt",cnt);
  std::ofstream out(buf);
  for (size_t i=0; i<vars.size(); ++i) {
    const rv::RVec<float>& row = vars[i];
    //const auto& var_name = params.var_names.at(i);
    for (size_t j=0; j<row.size(); ++j) {
      out << row[j];
      if (j!=row.size()-1)
        out << ",";
    }
    out << std::endl;
  }
  out.close();
}

void WeaverInterface::writeToFile(const std::vector< std::vector<float> >& vars) {
  static size_t cnt(0);
  ++cnt;

  char buf[1024];
  sprintf(buf,"norm%05ld.txt",cnt);
  std::ofstream out(buf);

  //out << "input_sizes_.size()=" << input_sizes_.size() << std::endl;
  //out << "vars.size()=" << vars.size() << std::endl;
  for (size_t i=0; i<vars.size(); ++i) {
    const auto& name = onnx_->inputNames()[i];
    const auto& params = prep_info_map_.at(name);
    out << "BLOCK[" << i << "]: " << name;
    const std::vector<float>& row = vars[i];
    //out << "row.size()=" << row.size() << std::endl;
    for (size_t j=0; j<row.size(); ++j) {
      if ( (j%75) == 0) {
        const auto& var_name = params.var_names.at(j/75);
        out << std::endl << "  VARNAME[" << j/75 << "]: " << var_name << std::endl;
      }

      out << "    " << row[j];
      if (j!=row.size()-1) {
        out << ",";
      }
    }
    out << std::endl;
  }
  out.close();
}
#endif


rv::RVec<float> WeaverInterface::run(const rv::RVec<ConstituentVars>& constituents) {
  size_t i = 0;
  for (const auto& name : onnx_->inputNames()) {
    const auto& params = prep_info_map_.at(name);
    auto& values = data_[i];
    values.resize(input_sizes_.at(i));
    std::fill(values.begin(), values.end(), 0);
    size_t it_pos = 0;
    ConstituentVars jc;
    for (size_t j = 0; j < params.var_names.size(); ++j) {  // transform and add the proper amount of padding
      const auto& var_name = params.var_names.at(j);
      //if (std::find(variables_names_.begin(), variables_names_.end(), "pfcand_mask") == variables_names_.end())
      //  jc = ConstituentVars(constituents.at(0).size(), 1.f);

      //if (var_name.find("pfcand_mask") != std::string::npos)
        //jc = ConstituentVars(constituents.at(0).size(), 1.f);

      if (var_name == "pfcand_mask") {
        jc = ConstituentVars(constituents.at(0).size(), 1.f);
      } else if (var_name == "neu_pfcand_mask") {
        jc = ConstituentVars(constituents.at(constituents.size()-1).size(), 1.f);
      } else {
        jc = constituents.at(variablePos(var_name));
      }
      const auto& var_info = params.info(var_name);
      auto val = center_norm_pad(jc,
                                 var_info.center,
                                 var_info.norm_factor,
                                 params.min_length,
                                 params.max_length,
                                 var_info.pad,
                                 var_info.replace_inf_value,
                                 var_info.lower_bound,
                                 var_info.upper_bound);
      std::copy(val.begin(), val.end(), values.begin() + it_pos);
      it_pos += val.size();
    }
    values.resize(it_pos);
    ++i;
  }
#ifdef DUMP_WEAVER_INPUT
  writeToFile(constituents);
  writeToFile(data_);
#endif
  return onnx_->run<float>(data_, input_shapes_)[0];
}

void WeaverInterface::PreprocessParams::dumpVars() const {
  std::cout << "List of variables for preprocessing parameter '" << name
            << "': " << rv::RVec<std::string>(var_names.begin(), var_names.end()) << "." << std::endl;
}