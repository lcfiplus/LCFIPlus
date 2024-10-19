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
    void processJet(); // for jet classification mode
    void processEvent(); // for event classification mode
    std::string _jetCollectionName;
    std::string _jsonFileName;
    std::string _onnxFileName;
    std::string _outputParameterName;
    bool _eventClassification;
    std::vector<std::string> _outputVariables;
    rv::RVec<std::string> _variables; //!
    WeaverInterface* _weaver; //!

    struct PreprocessParams {
      struct VarInfo {
        VarInfo() {}
        VarInfo(float imedian,
                float inorm_factor,
                float ireplace_inf_value,
                float ilower_bound,
                float iupper_bound,
                float ipad)
            : center(imedian),
              norm_factor(inorm_factor),
              replace_inf_value(ireplace_inf_value),
              lower_bound(ilower_bound),
              upper_bound(iupper_bound),
              pad(ipad) {}

        float center{0.};
        float norm_factor{1.};
        float replace_inf_value{0.};
        float lower_bound{-5.};
        float upper_bound{5.};
        float pad{0.};
      };
      std::string name;
      size_t min_length{0}, max_length{0};
      std::vector<std::string> var_names;
      std::unordered_map<std::string, VarInfo> var_info_map;
      VarInfo info(const std::string& name) const { return var_info_map.at(name); }
      void dumpVars() const;
    };
    std::unordered_map<std::string, PreprocessParams> prep_info_map_; //!
    void parseJSON(const string& json_filename);

    ClassDef(MLInferenceWeaver,1);
};

}

#endif
