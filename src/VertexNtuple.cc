#include "VertexNtuple.h"

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

void VertexNtuple::init(Parameters* param) {
  Algorithm::init(param);

  _primvtxcolname = param->get("PrimaryVertexCollectionName",string("PrimaryVertex"));
  _jetcolname = param->get("VertexNtuple.JetCollectionName",string("VertexJets"));
  Event::Instance()->setDefaultPrimaryVertex(_primvtxcolname.c_str()); // backward compatibility

  string filename = param->get("VertexNtuple.RootFileName",string("vertices.root"));

  // file open
  _file = TFile::Open(filename.c_str(), "RECREATE");
  if (!_file)throw(Exception("VertexNtuple: ROOT file open error! Failed if the file already exists."));

  _tree = new TNtuple("vertices","","nvtxWithMultiTracks:nvtxWithSingleTrack:pdgOfSingleTrack:parentPdgOfSingleTrack");

}

void VertexNtuple::process() {

  Event* event = Event::Instance();

  //const Vertex* privtx = event->getPrimaryVertex(_primvtxcolname.c_str());

  JetVec* jetsPtr(0);
  bool success = event->Get(_jetcolname.c_str(), jetsPtr);
  if (!success) {
    cout << "jets could not be found" << endl;
    return;
  }
  JetVec& jets = *jetsPtr;

  for (unsigned int nj = 0; nj < jets.size(); nj++) {
    const Jet* j = jets[nj];

    VertexVec vertices = j->getVertices();

    int nvtxSingleTrack = 0;
    int nvtxMultiTracks = 0;
    int pdgOfSingleTrack = 0;
    int parentPdgOfSingleTrack = 0;
    for (unsigned int n=0; n<vertices.size(); n++) {
      const Vertex* vtx = vertices[n];
      if (vtx->getTracks().size() == 1) {
        nvtxSingleTrack++; 
        pdgOfSingleTrack = vtx->getTracks()[0]->getPDG();
        const MCParticle* mcp = vtx->getTracks()[0]->getMcp();
        if (mcp) {
          const MCParticle* parent = mcp->getParent();
          if (parent) {
            parentPdgOfSingleTrack = parent->getPDG();
          }
        }
      } else nvtxMultiTracks++;
    }

    _tree->Fill(nvtxMultiTracks,nvtxSingleTrack,pdgOfSingleTrack,parentPdgOfSingleTrack);

  }
}

void VertexNtuple::end() {
  _file->Write();
}
}
