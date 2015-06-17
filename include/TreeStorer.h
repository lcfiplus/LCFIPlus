// TreeStorer.h
#ifndef TreeStore_h
#define TreeStore_h 1

#include "TFile.h"
#include "TTree.h"

using namespace std;

namespace lcfiplus {

class TreeStorer {
 public:
  enum {mode_input, mode_output};

  TreeStorer(const char* filename, const char* treename, int mode);
  virtual ~TreeStorer();

  // both for output and input
  void Register(const char* name);

  // for input only: register all branches
  void RegisterAll();

  // output function
  void Fill();
  void Write();

  // for input
  void GetEntry(int n);

  // for manual operation
  TTree* GetTree() {
    return _tree;
  }
  int GetEntries()const {
    return _tree->GetEntries();
  }

 private:
  TFile* _file;
  TTree* _tree;

  int _mode;
};

}

#endif
