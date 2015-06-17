// ReadMVA.h

#ifndef ReadMVA_h
#define ReadMVA_h 1

class TFile;
class TTree;

#include "lcfiplus.h"
#include "flavtag.h"

#include "TMVA/Types.h"
#include "TMVA/Reader.h"

namespace lcfiplus {

/**
	Lcfiplus algorithm for reading training data from TMVA.
	@author T. Tanabe, ICEPP, The University of Tokyo
	@version $Id$
 */
class ReadMVA : public Algorithm {
 public:
  ReadMVA() {}
  ~ReadMVA() {}

  void init(Parameters* param);
  void process();
  void end();

 protected:
  bool _verbose;

  TString _inputFileB;
  TString _inputFileC;
  TString _inputFileO;
  TString _treeNameB;
  TString _treeNameC;
  TString _treeNameO;
  TString _cutB;
  TString _cutC;
  TString _cutO;
  TString _outputDirectory;
  TString _outputPrefix;
  TString _treeName;

  TMVA::Types::EMVA _tmvaBookType;
  TString _tmvaBookName;
  TString _tmvaBookOptions;

  vector<FlavtagCategory> _categories;
  vector<TMVA::Reader*> _readers;

  ClassDef(ReadMVA,1);
};

}

#endif
