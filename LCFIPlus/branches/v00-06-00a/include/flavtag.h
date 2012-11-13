#ifndef flavtag_h
#define flavtag_h 1

#include "flavtag.h"
#include "JetFinder.h"

#include "TMVA/Reader.h"

class TFile;
class TTree;

namespace lcfiplus {

	// definition for flavor tag category
	struct FlavtagCategory {
		TString definition;
		TString preselection;
		std::vector<std::string> vars;
		std::vector<std::string> spec;
		void AddVariable(std::string s, char c) { vars.push_back(s); }
		void AddSpectator(std::string s) { spec.push_back(s); }
	};

	struct FlavtagType {
		TString name;
		TString cut;
	};


	// base class for algorithm to compute flavor tagging variables
	class FTAlgo {
		public:
			FTAlgo(string name) : _name(name) {}
			virtual ~FTAlgo() {}
			void setEvent(const Event* event, const Vertex *privtx);
			void setJet(const Jet* jet);
			void setNHitsJointProbD0(int value);
			void setNHitsJointProbZ0(int value);
			void setNHitsMostSignificantTrack(int value);
			float getValue();
			const string& getName() const { return _name; }
			float* getValueAddress() { return &_result; }

		protected:
			const Event* _event;
			const Vertex* _privtx;
			const Jet* _jet;
			int _nhitsJointProbD0;
			int _nhitsJointProbZ0;
			int _nhitsMostSignificantTrack;
			float _result;
			string _name;

		//////////////////////////////////////////
		// methods to be overloaded by subclass //
		//////////////////////////////////////////
		public:
			virtual void processEvent() {} // called once per event
			virtual void process() {} // called for each jet
	};

	// forward declaration for singleton
	class FTManager;

	class FTManager {
		private:
			static FTManager _theInstance;

		public:
			static FTManager& getInstance() { return _theInstance; }

			void add(FTAlgo* v);

			void fillTree();
			void openTree();
			void openFile(const char* filename);
			void closeFile();
			void process(const Event* event, const Vertex *privtx, int nhitsJointProbD0, int nhitsJointProbZ0, int nhitsMostSignificantTrack, JetVec & jets);

			float* getVarAddress(const string& varname);
			void setEval(bool seteval, bool exportAllVars) { _evaluate = seteval; _exportAllVars = exportAllVars;}

			void addReader(TMVA::Reader* reader, const FlavtagCategory& c);
			void setParamName(TString s) { _paramName = s; }

		private:
			FTManager();

			vector<FTAlgo*> _algoList;

			TFile* _file;
			TTree* _tree;
			string _ntpName;

			bool _evaluate;
			bool _exportAllVars;

			vector<TMVA::Reader*> _readers;
			vector<FlavtagCategory> _categories;
			TString _paramName;
	};

}

#endif

