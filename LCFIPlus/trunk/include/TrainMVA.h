// TrainMVA.h

#ifndef TrainMVA_h
#define TrainMVA_h 1

class TFile;
class TTree;

#include "interface.h"
#include "flavtag.h"
#include "FlavtagCategory.h"

#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TMVA/Types.h"

namespace flavtag{

	class TrainMVA : public FlavtagAlgorithm {
		public:
			TrainMVA(){}
			virtual ~TrainMVA(){}

			void init(FlavtagParameters *param);
			void process();
			void end();

		private:
			bool _verbose;

			TString _inputFile;
			TString _outputDirectory;
			TString _outputPrefix;
			TString _treeName;

			Int_t _tmvaBookType;
			TString _tmvaBookName;
			TString _tmvaBookOptions;

			vector<FlavtagCategory> _categoryList;
			vector<FlavtagType> _typeList;

			ClassDef(TrainMVA,1);
	};

}

#endif
