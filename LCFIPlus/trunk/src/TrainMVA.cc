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
#include "Driver.h"
#include "Suehara.h"

#include "TROOT.h"
#include "TInterpreter.h"
#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TCut.h"
#include "TRandom3.h"
#include "TSystem.h"

using namespace flavtag;
using namespace flavtag::algoSigProb;
using namespace flavtag::algoEtc;

void TrainMVA::init(FlavtagParameters *param) {
	FlavtagAlgorithm::init(param);

	_inputFile = param->get("TrainNtupleFile",string("flavtag.root"));

	_verbose = param->get("TrainVerbose",true);

	_outputDirectory = param->get("TrainOutputDirectory",TString("flavtag"));
	_outputPrefix = param->get("TrainOutputPrefix",TString("BDT"));

	TString bookTypeString = param->get("TrainBookType",TString("BDT"));
	if (bookTypeString == "BDT") {
		_tmvaBookType = TMVA::Types::kBDT;
	} else if (bookTypeString == "MLP") {
		_tmvaBookType = TMVA::Types::kMLP;
	} else if (bookTypeString == "Fisher") {
		_tmvaBookType = TMVA::Types::kFisher;
	}

	_tmvaBookName = param->get("TrainBookName",TString("bdt"));

	//_tmvaBookOptions = "!H:!V:NTrees=800:nEventsMin=400:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning";

	_tmvaBookOptions = param->get("TrainBookOptions",TString("!H:!V:NTrees=1000:BoostType=Grad:Shrinkage=0.30:UseBaggedGrad:GradBaggingFraction=0.6:SeparationType=GiniIndex:nCuts=20:PruneMethod=CostComplexity:PruneStrength=50:NNodesMax=5"));
}

void TrainMVA::process() {
	// do nothing
}

void TrainMVA::end() {

	TFile* ifile = new TFile(_inputFile);
	TTree* tree = (TTree*)ifile->Get(_treeName);

	gSystem->MakeDirectory(_outputDirectory);

	for (unsigned int icat=0; icat<_categoryList.size(); ++icat) {
		TTree* copytree = tree->CopyTree(_categoryList[icat].cut);
		TString prefix = _outputPrefix+(ULong_t)icat;
		TString s = _outputDirectory + prefix + ".root";
		TFile* outputFile = new TFile(s,"RECREATE");
		TMVA::Factory* factory = new TMVA::Factory(prefix,outputFile,
			"!V:Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=multiclass" );
		//"!V:Silent:Color:DrawProgressBar:Transofrmations=I;D;P;G,D");

		// define signal and background trees
		for (unsigned int itype=0; itype<_typeList.size(); ++itype) {
			TCut cut( _categoryList[icat].cut );
			cut += _typeList[itype].cut;
			TTree* tmpTree = tree->CopyTree(cut);
			factory->AddTree(tmpTree, _typeList[itype].name);
		}

		// add variables
		const FlavtagCategory& c = _categoryList[icat];
		for (unsigned int iv=0; iv<c.vars.size(); ++iv) {
			factory->AddVariable( c.vars[iv], 'F' );
			if (_verbose) {
				std::cout << "  - Adding variable '" << c.vars[iv] << "'" << std::endl;
			}
		}

		// add spectators
		for (unsigned int is=0; is<c.spec.size(); ++is) {
			factory->AddSpectator( c.spec[is] );
			if (_verbose) {
				std::cout << "  - Adding spectator '" << c.spec[is] << "'" << std::endl;
			}
		}

		factory->BookMethod( _tmvaBookType, _tmvaBookName, _tmvaBookOptions );
		factory->TrainAllMethods();

		outputFile->Close();
		delete outputFile;
		delete factory;
	}
}
