// TrainMVA.h

#ifndef TrainMVA_h
#define TrainMVA_h 1

class TFile;
class TTree;

#include "interface.h"
#include "lcfiplus.h"
#include "flavtag.h"

#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TMVA/Types.h"

namespace lcfiplus{

	class TrainMVA : public LcfiplusAlgorithm {
		public:
			TrainMVA(){}
			virtual ~TrainMVA(){}

			void init(LcfiplusParameters *param);
			void process();
			void end();

		private:
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

			Int_t _tmvaBookType;
			TString _tmvaBookName;
			TString _tmvaBookOptions;

			vector<FlavtagCategory> _categories;

			ClassDef(TrainMVA,1);

			void splitVars(const std::string& s, char c, std::vector<std::string>& v);
	};

}

#endif
