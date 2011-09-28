#include "MakeNtuple.h"

#include <assert.h>
#include "EventStore.h"
#include "lcfiplus.h"
#include "JetFinder.h"

#include "TROOT.h"
#include <string>
#include "flavtag.h"

using namespace lcfiplus;

namespace lcfiplus {
	void MakeNtuple::init(LcfiplusParameters *param) {
		LcfiplusAlgorithm::init(param);

		string outputFilename = param->get("TrainNtupleFile",string("flavtag.root"));
		cout << "FlavtagAlgorithm: Ntuple file set to " << outputFilename << endl;

		FTManager::getInstance().openFile(outputFilename.c_str());
	}

	void MakeNtuple::end() {
		FTManager::getInstance().closeFile();
	}
}
