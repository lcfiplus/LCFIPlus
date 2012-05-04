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

void TrainMVA::init(Parameters *param) {
	Algorithm::init(param);

	_inputFileB = param->get("TrainMVA.InputRootFileB",string("lcfiplusB.root"));
	_inputFileC = param->get("TrainMVA.InputRootFileC",string("lcfiplusC.root"));
	_inputFileO = param->get("TrainMVA.InputRootFileO",string("lcfiplusO.root"));

	_treeNameB = param->get("TrainMVA.TreeNameB",string("ntp"));
	_treeNameC = param->get("TrainMVA.TreeNameC",string("ntp"));
	_treeNameO = param->get("TrainMVA.TreeNameO",string("ntp"));

	_cutB = param->get("TrainMVA.PreSelectionB",string(""));
	_cutC = param->get("TrainMVA.PreSelectionC",string(""));
	_cutO = param->get("TrainMVA.PreSelectionO",string(""));

	_skipTrain = param->get("TrainMVA.SkipTrain",int(0));

	_verbose = param->get("TrainMVA.Verbose",true);

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
			_tmvaBookOptions = "!H:!V:NTrees=1000:BoostType=Grad:Shrinkage=0.10:UseBaggedGrad:GradBaggingFraction=0.50:nCuts=20:NNodesMax=8";
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
		c.preselection = param->get(psTag.str().c_str(),string(""));

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

	TFile* fileB = new TFile(_inputFileB);
	if (!fileB->IsOpen()) throw Exception( "could not open file" );

	TTree* treeB = (TTree*)fileB->Get(_treeNameB);
	if (treeB == 0) throw Exception( "could not find tree" );

	TFile* fileC = new TFile(_inputFileC);
	if (!fileC->IsOpen()) throw Exception( "could not open file" );

	TTree* treeC = (TTree*)fileC->Get(_treeNameC);
	if (treeC == 0) throw Exception( "could not find tree" );

	TFile* fileO = new TFile(_inputFileO);
	if (!fileO->IsOpen()) throw Exception( "could not open file" );

	TTree* treeO = (TTree*)fileO->Get(_treeNameO);
	if (treeO == 0) throw Exception( "could not find tree" );

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

		// define signal and background trees
		factory->AddTree( treeB, TString("jetB"), 1., TCut(c.definition) );
		factory->AddTree( treeC, TString("jetC"), 1., TCut(c.definition) );
		factory->AddTree( treeO, TString("jetO"), 1., TCut(c.definition) );

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
