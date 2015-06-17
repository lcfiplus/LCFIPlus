#include "lcfiplus.h"
#include "flavtag.h"

#include "TTree.h"
#include "TFile.h"
#include "TTreeFormula.h"
#include "TH1F.h"


namespace lcfiplus {

float FTAlgo::getValue() {
  return _result;
}

void FTAlgo::setEvent(const Event* event, const Vertex* privtx) {
  _event = event;
  _privtx = privtx;
}

void FTAlgo::setJet(const Jet* jet) {
  _jet = jet;
}

void FTAlgo::setNHitsJointProbD0(int value) {
  _nhitsJointProbD0 = value;
}

void FTAlgo::setNHitsJointProbZ0(int value) {
  _nhitsJointProbZ0 = value;
}

void FTAlgo::setNHitsMostSignificantTrack(int value) {
  _nhitsMostSignificantTrack = value;
}

// singleton
FTManager FTManager::_theInstance;

FTManager::FTManager() : _file(0), _tree(0), _ntpName("ntp"), _evaluate(false), _paramName("lcfiplus") {
  _initialized = false;
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

void FTManager::openTree() {
  _tree = new TTree(_ntpName.c_str(),"flavor tagging data");
  char buf[1024];
  for ( vector<FTAlgo*>::iterator iter = _algoList.begin(); iter != _algoList.end(); ++iter ) {
    FTAlgo* algo = *iter;
    snprintf( buf, 1000, "%s/F", algo->getName().c_str() );
    _tree->Branch( algo->getName().c_str(), algo->getValueAddress(), buf );

    //cout << "Setting branch address for variable '" << algo->getName() << "'" << endl;
  }
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
    snprintf( buf, 1000, "%s/F", algo->getName().c_str() );
    _tree->Branch( algo->getName().c_str(), algo->getValueAddress(), buf );

    //cout << "Setting branch address for variable '" << algo->getName() << "'" << endl;
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

void FTManager::process(const Event* event, const Vertex* privtx, int nhitsJointProbD0, int nhitsJointProbZ0, int nhitsMostSignificantTrack, JetVec& jets) {
  for ( vector<FTAlgo*>::iterator iter = _algoList.begin(); iter != _algoList.end(); ++iter ) {
    FTAlgo* algo = *iter;
    algo->setEvent(event, privtx);
    algo->processEvent();
  }

  for (JetVecIte iter = jets.begin(); iter != jets.end(); ++iter) {
    const Jet* jet = *iter;
    for ( vector<FTAlgo*>::iterator iter2 = _algoList.begin(); iter2 != _algoList.end(); ++iter2 ) {
      FTAlgo* algo = *iter2;
      algo->setJet(jet);
      algo->setNHitsJointProbD0(nhitsJointProbD0);
      algo->setNHitsJointProbZ0(nhitsJointProbZ0);
      algo->setNHitsMostSignificantTrack(nhitsMostSignificantTrack);
      algo->process();
    }

    int category = -1;
    if (_evaluate) {
      double tags[10];
      int ncat = _categories.size();
      int ncls = _readers[0]->EvaluateMulticlass("bdt").size();

      for (int i=0; i<ncat; ++i) {
        TTreeFormula form( "form", _categories[i].definition, _tree );
        //cout << "category " << i << " evaluates to: " << form.EvalInstance() << endl;
        TTreeFormula presel( "presel", _categories[i].preselection, _tree );

        if ( form.EvalInstance() == 1 ) {
          category = i;
          if ( presel.EvalInstance() == 1 ) {
            for (int nc = 0; nc < ncls; nc++)
              tags[nc] = _readers[i]->EvaluateMulticlass("bdt")[nc];
          } else {
            memset(tags, 0, sizeof(tags));
          }
          //cout << "jet category is " << i << " [btag=" << btag << ", ctag=" << ctag << "]" << endl;
        }
      }

      Parameters param;
      param.add( "BTag", tags[0] );
      param.add( "CTag", tags[1] );
      param.add( "OTag", tags[2] );
      if (ncls > 3)param.add("BBTag", tags[3]);
      if (ncls > 4)param.add("CCTag", tags[4]);
      if (ncls > 5)param.add("BCTag", tags[5]);

      param.add( "Category", (double)category );

      // adding all variables
      if (_exportAllVars) {
        for ( vector<FTAlgo*>::iterator iter2 = _algoList.begin(); iter2 != _algoList.end(); ++iter2 ) {
          FTAlgo* algo = *iter2;
          TTreeFormula form("form", algo->getName().c_str(),_tree);
          //cout << algo->getName() << endl;
          param.add(algo->getName().c_str(), form.EvalInstance());
        }
      }
      jet->addParam( _paramName.Data(), param );

    }

    if (_file) fillTree();
  }
}

float* FTManager::getVarAddress(const string& varname) {
  std::vector<FTAlgo*>::iterator iter;
  for (iter = _algoList.begin(); iter != _algoList.end(); ++iter) {
    FTAlgo* algo = *iter;
    if ( varname == algo->getName().c_str() )
      return algo->getValueAddress();
  }

  throw (Exception("FTManager::getVarAddress(): variable name not found."));
}

void FTManager::addReader(TMVA::Reader* reader, const FlavtagCategory& c) {
  _readers.push_back(reader);
  _categories.push_back(c);
}

FtIPProbHolder::FtIPProbHolder(const char* d0probfile, const char* z0probfile) {
  TDirectory* dir = gDirectory;
  _fd0 = TFile::Open(d0probfile);
  if (!_fd0)throw(Exception("FlavorTag: D0 probability file open error!"));
  _fz0 = TFile::Open(z0probfile);
  if (!_fz0)throw(Exception("FlavorTag: Z0 probability file open error!"));
  dir->cd();

  _hd0[0] = dynamic_cast<TH1F*>(_fd0->Get("hb"));
  _hd0[1] = dynamic_cast<TH1F*>(_fd0->Get("hc"));
  _hd0[2] = dynamic_cast<TH1F*>(_fd0->Get("hq"));

  if (!(_hd0[0] && _hd0[1] && _hd0[2]))throw(Exception("FlavorTag: D0 probability histogram open error!"));

  _hd0p[0] = dynamic_cast<TH1F*>(_fd0->Get("hbp"));
  _hd0p[1] = dynamic_cast<TH1F*>(_fd0->Get("hcp"));
  _hd0p[2] = dynamic_cast<TH1F*>(_fd0->Get("hqp"));

  if (!(_hd0p[0] && _hd0p[1] && _hd0p[2]))throw(Exception("FlavorTag: D0 probability positive histogram open error!"));

  _hd0n[0] = dynamic_cast<TH1F*>(_fd0->Get("hbn"));
  _hd0n[1] = dynamic_cast<TH1F*>(_fd0->Get("hcn"));
  _hd0n[2] = dynamic_cast<TH1F*>(_fd0->Get("hqn"));

  if (!(_hd0n[0] && _hd0n[1] && _hd0n[2]))throw(Exception("FlavorTag: D0 probability negative histogram open error!"));

  _hd0ip[0] = dynamic_cast<TH1F*>(_fd0->Get("hbip"));
  _hd0ip[1] = dynamic_cast<TH1F*>(_fd0->Get("hcip"));
  _hd0ip[2] = dynamic_cast<TH1F*>(_fd0->Get("hqip"));

  double norm;
  norm = _hd0ip[2]->Integral(_hd0ip[2]->FindBin(-0.5), _hd0ip[2]->FindBin(0.5));
  _normd0ip[0] = _hd0ip[0]->Integral(_hd0ip[0]->FindBin(-0.5), _hd0ip[0]->FindBin(0.5)) / norm;
  _normd0ip[1] = _hd0ip[1]->Integral(_hd0ip[1]->FindBin(-0.5), _hd0ip[1]->FindBin(0.5)) / norm;

  cout << "normd0ip[0] = " << _normd0ip[0] << " normd0ip[1] = " << _normd0ip[1] << endl;

  _hd0jprob = dynamic_cast<TH1F*>(_fd0->Get("hjprob"));
  _hd0jprob2 = dynamic_cast<TH1F*>(_fd0->Get("hjprob2"));
  if (!(_hd0jprob && _hd0jprob2))cout << "Joint probability not included in the D0 probability file. jprobr2 is not usable." << endl;

  /////////////////////////////// Z0

  if (!(_hd0ip[0] && _hd0ip[1] && _hd0ip[2]))throw(Exception("FlavorTag: D0 probability near-IP histogram open error!"));

  _hz0[0] = dynamic_cast<TH1F*>(_fz0->Get("hb"));
  _hz0[1] = dynamic_cast<TH1F*>(_fz0->Get("hc"));
  _hz0[2] = dynamic_cast<TH1F*>(_fz0->Get("hq"));

  if (!(_hz0[0] && _hz0[1] && _hz0[2]))throw(Exception("FlavorTag: Z0 probability histogram open error!"));

  _hz0ip[0] = dynamic_cast<TH1F*>(_fz0->Get("hbip"));
  _hz0ip[1] = dynamic_cast<TH1F*>(_fz0->Get("hcip"));
  _hz0ip[2] = dynamic_cast<TH1F*>(_fz0->Get("hqip"));

  norm = _hz0ip[2]->Integral(_hz0ip[2]->FindBin(-0.5), _hz0ip[2]->FindBin(0.5));
  _normz0ip[0] = _hz0ip[0]->Integral(_hz0ip[0]->FindBin(-0.5), _hz0ip[0]->FindBin(0.5)) / norm;
  _normz0ip[1] = _hz0ip[1]->Integral(_hz0ip[1]->FindBin(-0.5), _hz0ip[1]->FindBin(0.5)) / norm;

  cout << "normz0ip[0] = " << _normz0ip[0] << " normz0ip[1] = " << _normz0ip[1] << endl;

  if (!(_hz0ip[0] && _hz0ip[1] && _hz0ip[2]))throw(Exception("FlavorTag: Z0 probability near-IP histogram open error!"));

  _hz0jprob = dynamic_cast<TH1F*>(_fz0->Get("hjprob"));
  _hz0jprob2 = dynamic_cast<TH1F*>(_fz0->Get("hjprob2"));
  if (!(_hz0jprob && _hz0jprob2))cout << "Joint probability not included in the Z0 probability file. jprobz2 is not usable." << endl;

}

FtIPProbHolder::~FtIPProbHolder() {
  if (_fd0)_fd0->Close();
  if (_fz0)_fz0->Close();
}

}
