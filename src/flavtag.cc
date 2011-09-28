#include "lcfiplus.h"
#include "flavtag.h"

#include "TTree.h"
#include "TFile.h"

namespace lcfiplus {

	double FTAlgo::getValue() {
		return _result;
	}

	void FTAlgo::setEvent(const Event* event) {
		_event = event;
	}

	void FTAlgo::setJet(const Jet* jet) {
		_jet = jet;
	}

	// singleton
	FTManager FTManager::_theInstance;

	FTManager::FTManager() : _file(0), _tree(0), _ntpName("ntp")
	{
	}

	void FTManager::add(FTAlgo* obj) {
		_algoList.push_back(obj);
	}

	void FTManager::fillTree() {
		if (_tree == 0) {
			cout << "FTManager: file not opened" << endl;
			exit(1);
		}
		_tree->Fill();
	}

	void FTManager::openFile(const char* filename) {
		if (_file) {
			cout << "FTManager: file already open" << endl;
			exit(1);
		}

		cout << "FTManager: opening file '" << filename << "' for output" << endl;
		_file = new TFile(filename,"RECREATE");
		cout << "FTManager: setting TTree name '" << _ntpName.c_str() << "'" << endl;
		_tree = new TTree(_ntpName.c_str(),"flavor tagging data");

		char buf[1024];
		for ( vector<FTAlgo*>::iterator iter = _algoList.begin(); iter != _algoList.end(); ++iter ) {
			FTAlgo* algo = *iter;
			snprintf( buf, 1000, "%s/D", algo->getName() );
			_tree->Branch( algo->getName(), algo->getValueAddress(), buf );

			cout << "Setting branch address for variable '" << algo->getName() << "'" << endl;
		}
	}

	void FTManager::closeFile() {
		if (_file) {
			cout << "FTManager: closing file" << endl;
			_tree->Write();
			_file->Close();
			_file = 0;
		} else {
			cout << "FTManager: no file was opened" << endl;
		}
	}

	void FTManager::process(const Event* event, const vector<Jet*>& jets) {
		for ( vector<FTAlgo*>::iterator iter = _algoList.begin(); iter != _algoList.end(); ++iter ) {
			FTAlgo* algo = *iter;
			algo->setEvent(event);
			algo->processEvent();
		}

		for (vector<Jet*>::const_iterator iter = jets.begin(); iter != jets.end(); ++iter) {
			const Jet* jet = *iter;
			for ( vector<FTAlgo*>::iterator iter2 = _algoList.begin(); iter2 != _algoList.end(); ++iter2 ) {
				FTAlgo* algo = *iter2;
				algo->setJet(jet);
				algo->process();
			}

			if (_file) fillTree();
		}
	}
}
