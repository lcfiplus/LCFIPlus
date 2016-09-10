#ifndef flavtag_h
#define flavtag_h 1

#include "flavtag.h"
#include "JetFinder.h"

#include "TMVA/Reader.h"

class TFile;
class TTree;
class TH1F;

namespace lcfiplus {

// definition for flavor tag category
struct FlavtagCategory {
  TString definition;
  TString preselection;
  std::vector<std::string> vars;
  std::vector<std::string> spec;
  void AddVariable(std::string s, char /*c*/) {
    vars.push_back(s);
  }
  void AddSpectator(std::string s) {
    spec.push_back(s);
  }
};

struct FlavtagType {
  TString name;
  TString cut;
};


// base class for algorithm to compute flavor tagging variables
class FTAlgo {
 public:
  FTAlgo(string name) : _name(name) {}
  virtual ~FTAlgo() {}
  void setEvent(const Event* event, const Vertex* privtx);
  void setJet(const Jet* jet);
  void setNHitsJointProbD0(int value);
  void setNHitsJointProbZ0(int value);
  void setNHitsMostSignificantTrack(int value);
  float getValue();
  const string& getName() const {
    return _name;
  }
  float* getValueAddress() {
    return &_result;
  }

 protected:
  const Event* _event;
  const Vertex* _privtx;
  const Jet* _jet;
  int _nhitsJointProbD0;
  int _nhitsJointProbZ0;
  int _nhitsMostSignificantTrack;
  float _result;
  string _name;

  //////////////////////////////////////////
  // methods to be overloaded by subclass //
  //////////////////////////////////////////
 public:
  virtual void processEvent() {} // called once per event
  virtual void process() {} // called for each jet
};

// forward declaration for singleton
class FTManager;
class FtIPProbHolder;

class FTManager {
 private:
  static FTManager _theInstance;

 public:
  static FTManager& getInstance() {
    return _theInstance;
  }

  void initVars();

  void fillTree();
  void openTree();
  void openFile(const char* filename);
  void closeFile();
  void process(const Event* event, const Vertex* privtx, int nhitsJointProbD0, int nhitsJointProbZ0, int nhitsMostSignificantTrack, JetVec& jets);

  float* getVarAddress(const string& varname);
  void setEval(bool seteval, bool exportAllVars) {
    _evaluate = seteval;
    _exportAllVars = exportAllVars;
  }

  void addReader(TMVA::Reader* reader, const FlavtagCategory& c);
  void setParamName(TString s) {
    _paramName = s;
  }


  // accessor of variables for flavor tagging
  const FtIPProbHolder* getIPProbHolder()const {
    return _holder;
  }
  void setIPProbHolder(FtIPProbHolder* holder) {
    _holder = holder;
  }

  double getAuxiliary()const {
    return _aux;
  }
  void setAuxiliary(double aux) {
    _aux = aux;
  }

 private:
  void add(FTAlgo* v); // 121214 moved to private
  FTManager();

  vector<FTAlgo*> _algoList;

  TFile* _file;
  TTree* _tree;
  string _ntpName;

  bool _initialized;
  bool _evaluate;
  bool _exportAllVars;

  vector<TMVA::Reader*> _readers;
  vector<FlavtagCategory> _categories;
  TString _paramName;

  // variables for flavor tagging
  FtIPProbHolder* _holder;
  double _aux;
};

// historgram holder for d0/z0 probability
class FtD0bProb;
class FtD0cProb;
class FtD0qProb;
class FtD0bProbSigned;
class FtD0cProbSigned;
class FtD0qProbSigned;
class FtD0bProbIP;
class FtD0cProbIP;
class FtZ0bProb;
class FtZ0cProb;
class FtZ0qProb;
class FtZ0bProbIP;
class FtZ0cProbIP;
class FtJProbR2;
class FtJProbZ2;
class FtJProbR25Sigma;
class FtJProbZ25Sigma;

class FtIPProbHolder {
  friend class FtD0bProb;
  friend class FtD0cProb;
  friend class FtD0qProb;
  friend class FtD0bProbSigned;
  friend class FtD0cProbSigned;
  friend class FtD0qProbSigned;
  friend class FtD0bProbIP;
  friend class FtD0cProbIP;
  friend class FtZ0bProb;
  friend class FtZ0cProb;
  friend class FtZ0qProb;
  friend class FtZ0bProbIP;
  friend class FtZ0cProbIP;
  friend class FtJProbR2;
  friend class FtJProbZ2;
  friend class FtJProbR25Sigma;
  friend class FtJProbZ25Sigma;

 public:
  FtIPProbHolder(const char* d0probfile, const char* z0probfile);
  ~FtIPProbHolder();

 private:
  TFile* _fd0;
  TFile* _fz0;
  TH1F* _hd0[3];
  TH1F* _hd0p[3]; // signed-positive
  TH1F* _hd0n[3]; // signed-negative
  TH1F* _hd0ip[3];
  TH1F* _hz0[3];
  TH1F* _hz0ip[3];

  TH1F* _hd0jprob;
  TH1F* _hd0jprob2;
  TH1F* _hz0jprob;
  TH1F* _hz0jprob2;

  // normalization factors
  double _normd0ip[2];
  double _normz0ip[2];
};



}

#endif

