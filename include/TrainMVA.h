// TrainMVA.h

#ifndef TrainMVA_h
#define TrainMVA_h 1

class TFile;
class TTree;

#include "lcfiplus.h"
#include "flavtag.h"

#include "TMVA/Types.h"

namespace lcfiplus {

/**
	Lcfiplus algorithm for training classifications using TMVA.
	@author T. Tanabe, ICEPP, The University of Tokyo
	@version $Id$
 */
class TrainMVA : public Algorithm {
 public:
  TrainMVA() {}
  virtual ~TrainMVA() {}

  void init(Parameters* param);
  void process();
  void end();


  struct InputFileInfo {
    TString name;
    TString fileName;
    TString treeName;
    TString presel;
    TFile* file;
    TTree* tree;

    InputFileInfo() {}
    InputFileInfo( TString n, TString f, TString t, TString p ) :
      name(n), fileName(f), treeName(t), presel(p), file(0), tree(0) {}
  };


 private:
  bool _verbose;
  vector<InputFileInfo> _inputFileInfo;
  void readInputFileInfo( Parameters* param, TString name );

  TString _outputDirectory;
  TString _outputPrefix;
  TString _treeName;

  TMVA::Types::EMVA _tmvaBookType;
  TString _tmvaBookName;
  TString _tmvaBookOptions;
  int _skipTrain;

  vector<FlavtagCategory> _categories;

  ClassDef(TrainMVA,1);
};

}

#endif
