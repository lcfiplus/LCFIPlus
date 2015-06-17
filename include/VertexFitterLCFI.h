// VertexFitterLCFI.h

#ifndef VertexFitterLCFI_h
#define VertexFitterLCFI_h 1

#include "lcfiplus.h"
#include "LcfiInterface.h"

// ZVTOP includes
#include <zvtop/include/VertexFitterKalman.h>
#include <util/inc/util.h>

using namespace std;
using namespace lcfiplus;

namespace lcfiplus {

template<class Iterator>
class VertexFitterLCFI {
 public:
  Vertex* operator() (Iterator tracksBegin, Iterator tracksEnd, Vertex* pointConstraint = 0) {
    LcfiInterface lcfi;

    vector<vertex_lcfi::TrackState*> trackStates;
    map<vertex_lcfi::TrackState*, const Track*> trackMap;

    for (Iterator it = tracksBegin; it != tracksEnd; it++) {
      vertex_lcfi::Track* ptr = lcfi.lcfiTrack(*it);
      vertex_lcfi::TrackState* pstate = ptr->makeState();
      trackStates.push_back(pstate);
      trackMap[pstate] = *it;
    }

    vertex_lcfi::ZVTOP::VertexFitterKalman kalman;
//      void VertexFitterKalman::fitVertex(const std::vector<TrackState*> & Tracks,
//                      InteractionPoint* IP, Vector3 & Result,
//                      Matrix3x3 & ResultError,
//                      do uble & ChiSquaredOfFit,
//                      std::map<TrackState*,double> & ChiSquaredOfTrack,
//                      double & ChiSquaredOfIP);

    vertex_lcfi::util::Vector3 vtxpos;
    vertex_lcfi::util::Matrix3x3 vtxerr;
    double chi2fit;
    double chi2ip;
    map<vertex_lcfi::TrackState*, double> chi2map;

    vertex_lcfi::ZVTOP::InteractionPoint* ipConst = 0;
    if (pointConstraint) {
      vertex_lcfi::util::SymMatrix3x3 pcerr;
      vertex_lcfi::util::Vector3 pc;

      pc(0) = pointConstraint->getX();
      pc(1) = pointConstraint->getY();
      pc(2) = pointConstraint->getZ();

      pcerr(0,0) = pointConstraint->getCov()[Vertex::xx];
      pcerr(0,1) = pointConstraint->getCov()[Vertex::xy];
      pcerr(0,2) = pointConstraint->getCov()[Vertex::xz];
      pcerr(1,1) = pointConstraint->getCov()[Vertex::yy];
      pcerr(1,2) = pointConstraint->getCov()[Vertex::yz];
      pcerr(2,2) = pointConstraint->getCov()[Vertex::zz];

      ipConst = new vertex_lcfi::ZVTOP::InteractionPoint(pc,pcerr);
    }

    kalman.fitVertex(trackStates,ipConst, vtxpos, vtxerr, chi2fit, chi2map, chi2ip);
    if (ipConst)
      chi2fit -= ipConst->chi2(vtxpos);
    delete ipConst;

    double prob = vertex_lcfi::util::prob(chi2fit, trackStates.size() * 2 - 3);
    double cov[6];
    cov[0] = vtxerr(0,0);
    cov[1] = vtxerr(0,1);
    cov[2] = vtxerr(1,1);
    cov[3] = vtxerr(0,2);
    cov[4] = vtxerr(1,2);
    cov[5] = vtxerr(2,2);

    /*
    cout << "VertexFitterLCFI pos: " << scientific
    	<< vtxpos(0) << "  " << vtxpos(1) << "  " << vtxpos(2) << fixed << endl;
    cout << "VertexFitterLCFI cov: " << scientific
    	<< cov[0] << "  " << cov[1] << "  " << cov[2] << "  " << cov[3] << "  " << cov[4] << "  " << cov[5] << fixed << endl;
     */

    Vertex* vtx = new Vertex(chi2fit, prob,vtxpos(0), vtxpos(1), vtxpos(2),cov, false);
    for (unsigned int i=0; i<trackStates.size(); i++) {
      vtx->add(trackMap[trackStates[i]], chi2map[trackStates[i]]);
    }

    return vtx;
  }
};

typedef VertexFitterLCFI<vector<const Track*>::const_iterator > VertexFitterLCFI_V;
typedef VertexFitterLCFI<list<const Track*>::const_iterator > VertexFitterLCFI_L;
}

#endif
