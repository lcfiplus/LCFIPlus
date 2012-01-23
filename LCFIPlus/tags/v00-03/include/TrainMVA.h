// TrainMVA.h

#ifndef TrainMVA_h
#define TrainMVA_h 1

class TFile;
class TTree;

#include "lcfiplus.h"
#include "flavtag.h"

#include "TMVA/Types.h"

namespace lcfiplus{

	class TrainMVA : public Algorithm {
		public:
			TrainMVA(){}
			virtual ~TrainMVA(){}

			void init(Parameters *param);
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

			TMVA::Types::EMVA _tmvaBookType;
			TString _tmvaBookName;
			TString _tmvaBookOptions;

			vector<FlavtagCategory> _categories;

			ClassDef(TrainMVA,1);
	};

}

#endif
