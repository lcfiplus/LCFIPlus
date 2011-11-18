#include "ReadMVA.h"

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

using namespace lcfiplus;
using namespace lcfiplus::algoSigProb;
using namespace lcfiplus::algoEtc;

void ReadMVA::init(Parameters *param) {
	Algorithm::init(param);

	// get FTManager for variable registration
	FTManager &mgr = FTManager::getInstance();
	mgr.setEval(true);

	// read output directory & file names
	_outputDirectory = param->get("TrainOutputDirectory",TString("lcfiweights"));
	_outputPrefix = param->get("TrainOutputPrefix",TString("zpole_v00"));
	_tmvaBookName = param->get("TrainBookName",TString("bdt"));

	// process categories & variables
	for (int i=1; ; ++i) {
		FlavtagCategory c;

		stringstream catTag;
		catTag << "FlavorTagCategoryDefinition" << i;
		c.definition = param->get(catTag.str().c_str(),string(""));
		if (c.definition == "") {
			cout << "definition for index " << i << " not found, skipping" << endl;
			break;
		}

		stringstream psTag;
		psTag << "FlavorTagCategoryPreselection" << i;
		c.preselection = param->get(psTag.str().c_str(),string(""));

		// assumes comma separated values
		stringstream varTag;
		varTag << "FlavorTagCategoryVariables" << i;
		param->fetchArray( varTag.str().c_str(), c.vars );

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

	for (unsigned int icat=0; icat<_categories.size(); ++icat) {
		const FlavtagCategory& c = _categories[icat];

		// specify output directory
		TString prefix = _outputPrefix+(ULong_t)icat;

		TMVA::Reader* reader = new TMVA::Reader( "Color:!Silent" );

		// add variables
		for (unsigned int iv=0; iv<c.vars.size(); ++iv) {
			reader->AddVariable( c.vars[iv], mgr.getVarAddress( c.vars[iv] ) );
			if (_verbose) std::cout << "  - Adding variable '" << c.vars[iv] << "'" << std::endl;
		}

		// add spectators
		for (unsigned int is=0; is<c.spec.size(); ++is) {
			reader->AddSpectator( c.spec[is], mgr.getVarAddress( c.spec[is] ) );
			if (_verbose) std::cout << "  - Adding spectator '" << c.spec[is] << "'" << std::endl;
		}

		TString wfile;
		wfile += _outputDirectory + "/" + _outputPrefix + "_c";
		wfile += icat;
		wfile += "_";
		wfile += _tmvaBookName;
		wfile += ".weights.xml";

		std::cout << "Opening weight file: " << wfile << endl;

		reader->BookMVA(_tmvaBookName, wfile);

		_readers.push_back(reader);

		mgr.addReader( reader, c );

		mgr.openTree();
	}
}

void ReadMVA::process() {
}

void ReadMVA::end() {
}
