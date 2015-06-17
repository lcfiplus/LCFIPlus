#include "TrackNtuple.h"

#include "EventStore.h"
#include "algoSigProb.h"
#include "algoEtc.h"

#include "TROOT.h"
#include "TFile.h"
#include "TNtuple.h"

using namespace lcfiplus;
using namespace lcfiplus::algoSigProb;
using namespace lcfiplus::algoEtc;

namespace lcfiplus {

void TrackNtuple::init(Parameters* param) {
  Algorithm::init(param);

  _primvtxcolname = param->get("PrimaryVertexCollectionName",string("PrimaryVertex"));
  _jetcolname = param->get("TrackNtuple.JetCollectionName",string("VertexJets"));
  Event::Instance()->setDefaultPrimaryVertex(_primvtxcolname.c_str()); // backward compatibility
  _hitcutJprob = param->get("TrackNtuple.HitsCountForJointProbability",int(4));

  string filename = param->get("TrackNtuple.RootFileName",string("tracks.root"));

  // file open
  _file = TFile::Open(filename.c_str(), "new");
  if (!_file)throw(Exception("TrackNtuple: ROOT file open error! Failed if the file already exists."));

  _tree = new TNtuple("tracks","d0/z0 of tracks","sd0:sd0sig:sz0:sz0sig:jprobcut");

}

void TrackNtuple::process() {

  Event* event = Event::Instance();
  //if (event->getTracks().size() == 0) return;

  const Vertex* privtx = event->getPrimaryVertex(_primvtxcolname.c_str());

  //TrackVec & tracks = event->getTracks();
  JetVec* jetsPtr(0);
  bool success = event->Get(_jetcolname.c_str(), jetsPtr);
  if (!success) {
    cout << "jets could not be found" << endl;
    return;
  }
  JetVec& jets = *jetsPtr;

  for (unsigned int nj = 0; nj < jets.size(); nj++) {
    const Jet* j = jets[nj];

    TrackVec tracks = j->getAllTracks(true);

    for (unsigned int n=0; n<tracks.size(); n++) {
      const Track* tr = tracks[n];

      float sd0 = signedD0(tr, j, privtx, true);
      float sd0sig = signedD0Significance(tr, j, privtx, true);
      float sz0 = signedZ0(tr, j, privtx, true);
      float sz0sig = signedZ0Significance(tr, j, privtx, true);

      bool hitCut = trackSelectionForFlavorTag(tr, _hitcutJprob);

      _tree->Fill(sd0,sd0sig,sz0,sz0sig,hitCut);
    }
  }
}

void TrackNtuple::end() {
  _file->Write();
}
}
