#ifndef flavtag_h
#define flavtag_h 1

#include "flavtag.h"
#include "JetFinder.h"

class TFile;
class TTree;

namespace lcfiplus {

	// base class for algorithm to compute flavor tagging variables
	class FTAlgo {
		public:
			FTAlgo(string name) : _name(name) {}
			virtual ~FTAlgo() {}
			void setEvent(const Event* event);
			void setJet(const Jet* jet);
			double getValue();
			const char* getName() const { return _name.c_str(); }
			void* getValueAddress() { return &_result; }

		protected:
			const Event* _event;
			const Jet* _jet;
			double _result;
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
			void openFile(const char* filename);
			void closeFile();
			void process(const Event* event, const vector<Jet*>& jets);

		private:
			FTManager();

			vector<FTAlgo*> _algoList;

			TFile* _file;
			TTree* _tree;
			string _ntpName;
	};


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

}

#endif

