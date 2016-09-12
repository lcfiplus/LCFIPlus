#include "TrainMVA.h"

#include <string>
#include <assert.h>
#include "EventStore.h"
#include "LcfiInterface.h"
#include "JetFinder.h"
#include "TreeStorer.h"
#include "VertexFitterLCFI.h"
#include "VertexFinderTearDown.h"
#include "algoSigProb.h"
#include "algoEtc.h"
#include "flavtag.h"

#include "TROOT.h"
#include "TInterpreter.h"
#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TCut.h"
#include "TRandom3.h"
#include "TSystem.h"

#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/Config.h"

using namespace lcfiplus;
using namespace lcfiplus::algoSigProb;
using namespace lcfiplus::algoEtc;

void TrainMVA::readInputFileInfo( Parameters* param, TString name ) {
  TString tokenFile = "TrainMVA.InputRootFile";
  TString tokenTree = "TrainMVA.TreeName";
  TString tokenPresel = "TrainMVA.Preselection";
  TString file = param->get( tokenFile+name, string("") );

  if (file != "") {
    TString tree = param->get( tokenTree+name, string("") );
    TString presel = param->get( tokenPresel+name, string("") );
    _inputFileInfo.push_back( InputFileInfo( name, file, tree, presel ) );

    if (_verbose) {
      std::cout << "TrainMVA sample name=" << name << ", tree=" << tree << ", presel=" << presel << std::endl;
    }
  }
}

void TrainMVA::init(Parameters* param) {
  Algorithm::init(param);

  _verbose = param->get("TrainMVA.Verbose",true);

  readInputFileInfo( param, "B" );
  readInputFileInfo( param, "C" );
  readInputFileInfo( param, "O" );
  readInputFileInfo( param, "BB" );
  readInputFileInfo( param, "CC" );
  readInputFileInfo( param, "BC" );

  _skipTrain = param->get("TrainMVA.SkipTrain",int(0));

  // set output directory for weight files
  _outputDirectory = param->get("FlavorTag.WeightsDirectory",TString("lcfiweights"));
  cout << "WeightsDirectory set to: " << _outputDirectory << endl;

  //
  _outputPrefix = param->get("FlavorTag.WeightsPrefix",TString("zpole_v00"));
  cout << "WeightsPrefix set to: " << _outputPrefix << endl;

  // set TMVA method (e.g. BDT, MLP)
  TString bookTypeString = param->get("TrainMVA.BookType",TString("BDT"));

  if (bookTypeString == "BDT") {
    _tmvaBookType = TMVA::Types::kBDT;
  } else if (bookTypeString == "MLP") {
    _tmvaBookType = TMVA::Types::kMLP;
  } else {
    cout << "unknown TMVA type: " << bookTypeString << endl;
  }

  _tmvaBookName = param->get("FlavorTag.BookName",TString("bdt"));
  _tmvaBookOptions = param->get("TrainMVA.BookOptions",TString(""));

  if (_tmvaBookOptions == "") {
    if (_tmvaBookType == TMVA::Types::kBDT)
      // parameters updated to TMVA 4.2 parameter names
      _tmvaBookOptions = "!H:!V:NTrees=1000:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:NegWeightTreatment=IgnoreNegWeightsInTraining:BaggedSampleFraction=0.50:nCuts=20:MaxDepth=3";
    else if (_tmvaBookType == TMVA::Types::kMLP)
      _tmvaBookOptions = "!H:!V:NeuronType=tanh:NCycles=1000:HiddenLayers=N+5,5:TestRate=5:EstimatorType=MSE";
  }

  cout << "book type: " << _tmvaBookType << endl;
  cout << "book name: " << _tmvaBookName << endl;
  cout << "book opts: " << _tmvaBookOptions << endl;

  // process categories & variables
  for (int i=1; ; ++i) {
    FlavtagCategory c;

    stringstream catTag;
    catTag << "FlavorTag.CategoryDefinition" << i;
    c.definition = param->get(catTag.str().c_str(),string(""));
    if (c.definition == "") {
      cout << "definition for index " << i << " not found, skipping" << endl;
      break;
    }

    stringstream psTag;
    psTag << "FlavorTag.CategoryPreselection" << i;
    c.preselection = param->get(psTag.str().c_str(),string("1"));

    // assumes comma separated values
    stringstream varTag;
    varTag << "FlavorTag.CategoryVariables" << i;
    param->fetchArray( varTag.str().c_str(), c.vars );

    // read spectators
    stringstream specTag;
    specTag << "FlavorTag.CategorySpectators" << i;
    param->fetchArray( specTag.str().c_str(), c.spec );

    cout << "FlavorTag category: " << c.definition << endl;
    cout << "FlavorTag preselection: " << c.preselection << endl;
    for (unsigned int i=0; i<c.vars.size(); ++i)
      cout << "FlavorTag variable: " << c.vars[i] << endl;

    _categories.push_back(c);
  }

  if (_categories.size() == 0) {
    cout << "FlavorTag category definition was not found, aborting" << endl;
    exit(1);
  }
}

void TrainMVA::process() {
  // do nothing
}

void TrainMVA::end() {

  cout << "TrainMVA::end" << endl;

  for (vector<InputFileInfo>::iterator iter = _inputFileInfo.begin();
       iter != _inputFileInfo.end(); ++iter) {
    iter->file = new TFile(iter->fileName);
    if (iter->file->IsOpen() == false) throw Exception( "could not open file" );

    iter->tree = (TTree*) iter->file->Get(iter->treeName);
    if (iter->tree == 0) throw Exception( "could not find tree" );
  }


  gSystem->MakeDirectory(_outputDirectory);

  for (unsigned int icat=0; icat<_categories.size(); ++icat) {
    if (_skipTrain>0 && (int)icat < _skipTrain) {
      cout << "Skipping training for category #" << icat << endl;
      continue;
    }
    const FlavtagCategory& c = _categories[icat];

    // specify output directory
    TString prefix = _outputPrefix + "_c";
    prefix += icat;
    TString s = _outputDirectory + "/" + prefix + ".root";
    TFile* outputFile = new TFile(s,"RECREATE");

    TMVA::gConfig().GetIONames().fWeightFileDir = _outputDirectory;
    TMVA::Factory* factory = new TMVA::Factory(prefix,outputFile,
        "!V:!Silent:Color:!DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=multiclass" );
    //"!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=multiclass" );

    // define signal and background trees

    for (vector<InputFileInfo>::iterator iter = _inputFileInfo.begin();
         iter != _inputFileInfo.end(); ++iter) {
      factory->AddTree(
        iter->tree, TString("jet")+iter->name, 1.,
        TCut(iter->presel)+TCut(c.definition)+TCut(c.preselection) );
    }

    // add variables
    for (unsigned int iv=0; iv<c.vars.size(); ++iv) {
      factory->AddVariable( c.vars[iv], 'F' );
      if (_verbose) std::cout << "  - Adding variable '" << c.vars[iv] << "'" << std::endl;
    }

    // add spectators
    for (unsigned int is=0; is<c.spec.size(); ++is) {
      factory->AddSpectator( c.spec[is] );
      if (_verbose) std::cout << "  - Adding spectator '" << c.spec[is] << "'" << std::endl;
    }


    factory->PrepareTrainingAndTestTree( "", "SplitMode=Random:NormMode=NumEvents:!V" );
    factory->BookMethod( _tmvaBookType, _tmvaBookName, _tmvaBookOptions );
    factory->TrainAllMethods();
    factory->TestAllMethods();
    factory->EvaluateAllMethods();

    outputFile->Close();
    delete outputFile;
    delete factory;
  }

}
