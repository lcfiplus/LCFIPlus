// VertexSelector.h

#ifndef VertexSelector_h
#define VertexSelector_h 1

#include "lcfiplus.h"
#include <vector>
#include <cmath>

using namespace std;
using namespace lcfiplus;

namespace lcfiplus {

class VertexSelectorConfig {
 public:
  bool rejectdist;
  double minpos;
  double maxpos;
  bool rejectdistnegative;
  bool rejectdistor;

  bool rejectk0;
  double k0width;
  double k0dirdot;

  bool rejectl0;
  double l0width;
  double l0dirdot;

  bool rejectconv;
  double convmass;
  double convdist;
  double convdirdot;

  VertexSelectorConfig() :
    rejectdist(false), minpos(0.), maxpos(1e+300),
    rejectdistnegative(false), rejectdistor(false),
    rejectk0(false), k0width(0.), k0dirdot(0.),
    rejectl0(false), l0width(0.), l0dirdot(0.),
    rejectconv(false), convmass(0.), convdist(0.), convdirdot(0.) {
  }

  void setV0Tight() { // for rejecting high purity v0 only
    rejectk0 = true;
    k0width = 0.005;
    k0dirdot = 0.999;
    rejectl0 = true;
    l0width = 0.005;
    l0dirdot = 0.99995;
    rejectconv = true;
    convmass = 0.005;
    convdist = 9;
    convdirdot = 0.99995;
  }
  void setV0Loose() { // for rejecting high purity v0 only
    rejectk0 = true;
    k0width = 0.010;
    k0dirdot = 0.999;
    rejectl0 = true;
    l0width = 0.010;
    l0dirdot = 0.999;
    rejectconv = true;
    convmass = 0.010;
    convdist = 9;
    convdirdot = 0.999;
  }
  void setNoV0Cut() {
    rejectk0 = false;
    rejectl0 = false;
    rejectconv = false;
  }

};

class VertexSelector {
 public:
  // const full version
  vector<const Vertex*> operator () (VertexVec& vertices, VertexSelectorConfig& config, vector<const Track*>& residualTracks, bool addTracks) {
    vector<const Vertex*> ret;

    for (unsigned int i=0; i<vertices.size(); i++) {
      if (vertices[i]->isPrimary()) {
        //cout << "Primary vertex found in Vertex Selector!!!!!!!" << endl;
        continue;
      }

      if (passesCut(vertices[i], config)) {
        ret.push_back(vertices[i]);
        if (!addTracks) {
          // TODO: possible use of STL algorithm set_difference - but need to sort
          for (TrackVecIte itt = vertices[i]->getTracks().begin(); itt != vertices[i]->getTracks().end(); itt++) {
            vector<const Track*>::iterator itt2 = remove_if(residualTracks.begin(), residualTracks.end(), bind2nd(equal_to<const Track*>(), *itt));
            residualTracks.erase(itt2, residualTracks.end());
          }
        }
      } else if (addTracks) {
        residualTracks.insert(residualTracks.end(),vertices[i]->getTracks().begin(), vertices[i]->getTracks().end());
      }
    }

    return ret;
  }
  // non-const full version
  vector<Vertex*> operator () (const vector<Vertex*>& vertices, VertexSelectorConfig& config, vector<const Track*>& residualTracks, bool addTracks) {
    vector<Vertex*> ret;

    for (unsigned int i=0; i<vertices.size(); i++) {
      if (passesCut(vertices[i], config)) {
        ret.push_back(vertices[i]);
        if (!addTracks) {
          // TODO: possible use of STL algorithm set_difference - but need to sort
          for (TrackVecIte itt = vertices[i]->getTracks().begin(); itt != vertices[i]->getTracks().end(); itt++) {
            vector<const Track*>::iterator itt2 = remove_if(residualTracks.begin(), residualTracks.end(), bind2nd(equal_to<const Track*>(), *itt));
            residualTracks.erase(itt2, residualTracks.end());
          }
        }
      } else if (addTracks) {
        residualTracks.insert(residualTracks.end(),vertices[i]->getTracks().begin(), vertices[i]->getTracks().end());
      }
    }

    return ret;
  }
  // const simple version
  vector<const Vertex*> operator () (VertexVec& vertices, VertexSelectorConfig& config) {
    vector<const Vertex*> ret;

    for (unsigned int i=0; i<vertices.size(); i++) {
      if (passesCut(vertices[i], config)) {
        ret.push_back(vertices[i]);
      }
    }

    return ret;
  }
  // non-const very simple version
  void operator () (vector<Vertex*>& vertices, VertexSelectorConfig& config) {

    for (vector<Vertex*>::iterator it = vertices.begin(); it!=vertices.end();) {
      if (!passesCut(*it, config)) {
        it = vertices.erase(it);
        continue;
      }
      it ++;
    }
  }

  // true if vtx is V0
  bool isV0(const Vertex* vtx, const VertexSelectorConfig& cfg, const Vertex* primary) {
    // two-track selection
    if ( vtx->getTracks().size() !=2 )
      return false;

    // charge selection
    if ( vtx->getTracks()[0]->getCharge() * vtx->getTracks()[1]->getCharge() != -1 )
      return false;

    // obtain momentum vectors at the V0 vertex (IMPORTANT)
    TVector3 mom1 = vtx->getTracks()[0]->momentumAtVertex(vtx);
    TVector3 mom2 = vtx->getTracks()[1]->momentumAtVertex(vtx);

    // V0 direction selection
    TVector3 posip;
    if (primary) posip = primary->getPos();
    double v0dirdot = (mom1+mom2).Unit().Dot( (vtx->getPos()-posip).Unit() );

    // K0 selection
    if (cfg.rejectk0) {
      TLorentzVector lvec1;
      TLorentzVector lvec2;
      lvec1.SetVectM( mom1, 0.1396 );
      lvec2.SetVectM( mom2, 0.1396 );
      double mass = (lvec1+lvec2).M();
      double k0mass = .498;
      if (v0dirdot > cfg.k0dirdot && fabs(mass-k0mass) < cfg.k0width )
        return true;
    }

    // Lambda selection
    if (cfg.rejectl0) {
      // assume higher momentum particle is the proton
      TLorentzVector protonForLambda;
      TLorentzVector pionForLambda;
      if (mom1.Mag() > mom2.Mag()) {
        protonForLambda.SetVectM( mom1, 0.9383 );
        pionForLambda.SetVectM( mom2, 0.1396 );
      } else {
        protonForLambda.SetVectM( mom2, 0.9383 );
        pionForLambda.SetVectM( mom1, 0.1396 );
      }

      double lambdaMass = (protonForLambda+pionForLambda).M();
      if (v0dirdot > cfg.l0dirdot && fabs(lambdaMass - 1.1157) < cfg.l0width)
        return true;
    }

    // check photon conversion; photon mass corrected geometrically
    if (cfg.rejectconv) {
      double ang1 = atan( vtx->getTracks()[0]->getTanLambda());
      double ang2 = atan( vtx->getTracks()[1]->getTanLambda());
      double convmassCorr = sqrt( vtx->getTracks()[0]->Vect().Mag()*vtx->getTracks()[1]->Vect().Mag()*(1-cos(ang1-ang2)) );
      if (v0dirdot > cfg.convdirdot && vtx->getPos().Mag()> cfg.convdist && convmassCorr < cfg.convmass)
        return true;
    }

    return false;
  }

  // true if passes cut
  bool passesCut(const Vertex* vtx, const VertexSelectorConfig& cfg, const Vertex* primary=0) {

    // position selection
    if (cfg.rejectdist) {
      bool cond = (vtx->getPos().Mag() > cfg.minpos && vtx->getPos().Mag() < cfg.maxpos);

      if (cfg.rejectdistnegative)cond = !cond;
      if (cfg.rejectdistor && cond)return true;
      if (!cfg.rejectdistor && !cond)return false;
    }

    // v0 selection
    if (cfg.rejectk0 || cfg.rejectl0 || cfg.rejectconv) {
      if (isV0(vtx,cfg,primary))return false;
//				if (isV0(vtx,cfg,primary)) {cout << "rejected by v0 selection" << endl;return false;}
    }

    return true;
  }

  //c-tor / d-tor
  VertexSelector() {}
  ~VertexSelector() {}
};
}

#endif //TrackSelector_h
