// VertexFinderPerfect.cc

#include "lcfiplus.h"
#include "VertexFinderPerfect.h"

namespace lcfiplus {
void VertexFinderPerfect::findPerfectVertices(TrackVec& tracks, MCParticleVec& mcp, vector<MCVertex*>& selvtx, int minimumRecoTracks, double minimumDistance, bool print) {
  for (unsigned int n=0; n<selvtx.size(); n++)
    delete selvtx[n];

  selvtx.clear();

  vector<MCVertex*> vtx;

  // perfect vertex search
  for (unsigned int n = 0; n < mcp.size(); n++) {
    // skip unstable & neutral particles
    if (!mcp[n]->isStable())continue;
    if (mcp[n]->getCharge() == 0.0) continue;

    const TVector3& v = mcp[n]->getVertex();

    unsigned int nvtx;
    for (nvtx = 0; nvtx < vtx.size(); nvtx ++) {
      TVector3 dif = vtx[nvtx]->getPos() - v;
      if (dif.Mag() < minimumDistance)break;
    }
    if (nvtx == vtx.size()) {
      // new perfect vertex
      MCVertex* newvtx = new MCVertex;
      newvtx->setPos(v);
      newvtx->add(mcp[n]);

      const MCParticle* parent = mcp[n]->getSemiStableParent();
      if (parent)newvtx->setParent(parent);

      vtx.push_back(newvtx);
    } else {
      vtx[nvtx]->add(mcp[n]);

      const MCParticle* parent = mcp[n]->getSemiStableParent();
      if (parent && (!vtx[nvtx]->getParent() || vtx[nvtx]->getParent()->isParent(parent)))
        vtx[nvtx]->setParent(parent);
    }
  }
  // mcp scan end

  // associate reco particle
  for (unsigned int n = 0; n < tracks.size(); n++) {
    const MCParticle* mcp = tracks[n]->getMcp();
    if (!mcp)continue;

    for (unsigned int nvtx = 0; nvtx < vtx.size(); nvtx ++) {
      if (find(vtx[nvtx]->getDaughters().begin(), vtx[nvtx]->getDaughters().end(), mcp) != vtx[nvtx]->getDaughters().end()) {
        // vertex found
        vtx[nvtx]->add(tracks[n]);
        break;
      }
    }
  }

  // select vertices
  for (unsigned int nvtx = 0; nvtx < vtx.size(); nvtx ++) {
    if (vtx[nvtx]->getRecoTracks().size() >= (unsigned int)minimumRecoTracks) {
      selvtx.push_back(vtx[nvtx]);
    } else
      delete vtx[nvtx];
  }

  // print status
  if (print) {
    for (unsigned int nvtx = 0; nvtx < selvtx.size(); nvtx ++) {
      cout << "MCVertex found at ( " << selvtx[nvtx]->getPos().x() << ", " << selvtx[nvtx]->getPos().y() << ", " << selvtx[nvtx]->getPos().z()
           << "), # tracks = " << selvtx[nvtx]->getRecoTracks().size() << " ";
      if (nvtx>0) {
        for (unsigned int t=0; t<selvtx[nvtx]->getRecoTracks().size(); t++) {
          const MCParticle* mcp = selvtx[nvtx]->getRecoTracks()[t]->getMcp();
          const MCParticle* mcp2 = (mcp ? mcp->getSemiStableParent() : 0);
          cout << (mcp2 ? mcp2->getPDG() : 0) << " ";
        }
      }
      cout << endl;
    }
  }

}
}
