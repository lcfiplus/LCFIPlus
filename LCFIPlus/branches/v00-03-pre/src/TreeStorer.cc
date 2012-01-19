// TreeStorer.cc
#include "TreeStorer.h"
#include "lcfiplus.h"
#include "EventStore.h"

#include "TBranchElement.h"

namespace lcfiplus{

	// ctor/dtor
	TreeStorer::TreeStorer(const char *filename, const char *treename, int mode)
	{
		// open file
		_file = TFile::Open(filename, mode == mode_input ? "" : "recreate");
		if(!_file)throw(new Exception("Tree file cannot be opened."));

		// open tree
		if(mode == mode_input)
			_tree = dynamic_cast<TTree *>(_file->Get(treename));
		else
			_tree = new TTree(treename, treename);
		if(!_tree)throw(new Exception("Tree object cannot be obtained."));

		_mode = mode;
	}

	TreeStorer::~TreeStorer()
	{
		if(_file)
			delete _file; // tree is automatically deleted
	}

	// both for output and input
	void TreeStorer::Register(const char *name)
	{
		EventStore *store = Event::Instance();
		if(!_tree)throw(new Exception("TreeStorer::Register: Tree not initialized!"));

		if(_mode == mode_input){
			// obtain class information
			TBranchElement *br = dynamic_cast<TBranchElement*>(_tree->GetBranch(name));
			if(!br)throw(new Exception("TreeStorer::Register: Branch object cannot be obtained."));
			const char *classname = br->GetClassName();

			// secure the object
			void **bufptr = new (void *);
			*bufptr = store->RegisterObject(name, classname);
			if(!(*bufptr))throw(new Exception("TreeStorer::Register: EventStore::RegisterObject failed."));

			// set to the branch
			br->SetAddress(bufptr);
		}else{ // mode_output
			const char *classname = store->GetClassName(name);
			void ** bufptr = new (void *);
			*bufptr = store->GetObject(name);

			_tree->Branch(name, classname, bufptr);
		}
	}

	// for input only: register all branches
	void TreeStorer::RegisterAll()
	{
		if(!_tree)throw(new Exception("TreeStorer::RegisterAll: Tree not initialized!"));
		if(_mode != mode_input)throw(new Exception("TreeStorer::RegisterAll: cannot be called for output trees."));

		TObjArray *lbr = _tree->GetListOfBranches();
		for(int i=0;i<lbr->GetEntries();i++){
			const char *name = dynamic_cast<TNamed*>(lbr->At(i))->GetName();
			Register(name);
		}
	}

	// for output
	void TreeStorer::Fill()
	{
		if(_mode==mode_output && _tree)_tree->Fill();
	}
	void TreeStorer::Write()
	{
		if(_mode==mode_output && _tree)_tree->Write();
	}

	// for input
	void TreeStorer::GetEntry(int n)
	{
		if(_mode==mode_input && _tree)_tree->GetEntry(n);
	}
}
