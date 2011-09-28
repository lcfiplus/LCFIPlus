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


using namespace lcfiplus;
using namespace lcfiplus::algoSigProb;
using namespace lcfiplus::algoEtc;

void TrainMVA::init(LcfiplusParameters *param) {
	LcfiplusAlgorithm::init(param);

	_inputFileB = param->get("TrainNtupleFileB",string("lcfiplusB.root"));
	_inputFileC = param->get("TrainNtupleFileC",string("lcfiplusC.root"));
	_inputFileO = param->get("TrainNtupleFileO",string("lcfiplusO.root"));

	_treeNameB = param->get("TrainTreeNameB",string("ntp"));
	_treeNameC = param->get("TrainTreeNameC",string("ntp"));
	_treeNameO = param->get("TrainTreeNameO",string("ntp"));

	_cutB = param->get("TrainPreSelectionB",string(""));
	_cutC = param->get("TrainPreSelectionC",string(""));
	_cutO = param->get("TrainPreSelectionO",string(""));

	_verbose = param->get("TrainVerbose",true);

	_outputDirectory = param->get("TrainOutputDirectory",TString("lcfiplus"));
	_outputPrefix = param->get("TrainOutputPrefix",TString("BDT"));

	TString bookTypeString = param->get("TrainBookType",TString("BDT"));

	if (bookTypeString == "BDT") {
		_tmvaBookType = TMVA::Types::kBDT;
	} else if (bookTypeString == "MLP") {
		_tmvaBookType = TMVA::Types::kMLP;
	} else if (bookTypeString == "Fisher") {
		_tmvaBookType = TMVA::Types::kFisher;
	} else {
		cout << "unknown TMVA type: " << bookTypeString << endl;
	}

	_tmvaBookName = param->get("TrainBookName",TString("bdt"));

	_tmvaBookOptions = param->get("TrainBookOptions",TString("!H:!V:NTrees=1000:BoostType=Grad:Shrinkage=0.30:UseBaggedGrad:GradBaggingFraction=0.6:SeparationType=GiniIndex:nCuts=20:PruneMethod=CostComplexity:PruneStrength=50:NNodesMax=5"));


	// process categories & variables

	for (int i=1; ; ++i) {
		FlavtagCategory c;

		stringstream catTag;
		catTag << "FlavorTagCategoryDefinition" << i;
		c.definition = param->get(catTag.str().c_str(),string(""));
		if (c.definition == "") break;

		stringstream psTag;
		psTag << "FlavorTagCategoryPreselection" << i;
		c.preselection = param->get(psTag.str().c_str(),string(""));

		// assumes comma separated values
		stringstream varTag;
		varTag << "FlavorTagCategoryVariables" << i;
		string vars = param->get(varTag.str().c_str(),string(""));
		splitVars( vars, ',', c.vars );

		_categories.push_back(c);
	}

}

void TrainMVA::splitVars(const std::string& s, char c, std::vector<std::string>& v) {
	string::size_type i(0);
	string::size_type j( s.find(c) );
	while (j != string::npos) {
		v.push_back(s.substr(i, j-1));
		i = ++j;
		j = s.find(c,j);
		if (j==string::npos)
			v.push_back(s.substr(i,s.length()));
	}
}

void TrainMVA::process() {
	// do nothing
}

void TrainMVA::end() {

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

#if 0

	for (unsigned int icat=0; icat<_categories.size(); ++icat) {
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
#endif
}
