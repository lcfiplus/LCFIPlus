#include <algorithm>
#include <algo/inc/pereventipfitter.h>

#include <inc/track.h>
#include <inc/trackstate.h>
#include <inc/event.h>
#include <inc/vertex.h>
#include <zvtop/include/candidatevertex.h>
#include <zvtop/include/VertexFitterKalman.h>
#include <util/inc/memorymanager.h>
#include <inc/decaychain.h>
#include <algo/inc/vertexmass.h>

#include "EventStore.h"
#include "JetFinder.h"
#include "LcfiInterface.h"

#include <stdlib.h>

namespace flavtag {

  LcfiInstance LcfiInstance::_instance;

  LcfiInstance::LcfiInstance() : _ipFitter(0), _zvres(0) {
  }

  LcfiInstance::~LcfiInstance() {
    if (_ipFitter) delete _ipFitter;
    if (_zvres) delete _zvres;
  }

  const vertex_lcfi::PerEventIPFitter* LcfiInstance::getIpFitter() {
    if (!_ipFitter) {
      _ipFitter = new vertex_lcfi::PerEventIPFitter();
      _ipFitter->setDoubleParameter("ProbThreshold",0.01);
      //_ipFitter->setDoubleParameter("ProbThreshold",0.005);
      //_ipFitter->setDoubleParameter("ProbThreshold",-1);
    }

    return _ipFitter;
  }

  vertex_lcfi::ZVRES* LcfiInstance::getZVRES() {
    if (!_zvres) {
      _zvres = new vertex_lcfi::ZVRES();
    }
    return _zvres;
  }

  LcfiInterface::LcfiInterface(Event* event, Vertex* primaryVertex) : debug(false), _primaryVertex(0) {
		if (primaryVertex) {
			_primaryVertex = lcfiVertex(primaryVertex,true); // it's important that this is set first
		}
    _event = lcfiEvent(event);
    vertex_lcfi::MemoryManager<vertex_lcfi::Event>::Event()->registerObject(_event);
		//debug = true;
  }

  LcfiInterface::~LcfiInterface() {
    vertex_lcfi::MetaMemoryManager::Event()->delAllObjects();
  }

  vertex_lcfi::Track*
  LcfiInterface::lcfiTrack(vertex_lcfi::Event* MyEvent, Track* track) const {
    HelixRep H;
    H.d0()	= track->getD0();
    H.z0()	= track->getZ0();
    H.phi()	= track->getPhi();
    H.invR()	= track->getOmega();
    H.tanLambda()  = track->getTanLambda();

    Vector3 Mom;
    Mom.x() = track->Px();
    Mom.y() = track->Py();
    Mom.z() = track->Pz();

    SymMatrix5x5 Cov;

    // Classic order
    const float* cov = track->getCovMatrix();
    Cov(0,0)=cov[0]; // d0d0
    Cov(3,0)=cov[6]; // d0z0
    Cov(3,3)=cov[9]; // z0z0

    Cov(1,0)=cov[1];
    Cov(1,1)=cov[2];
    Cov(0,2)=cov[3];
    Cov(1,2)=cov[4];
    Cov(2,2)=cov[5];
    Cov(1,3)=cov[7];

    Cov(2,3)=cov[8];
    Cov(0,4)=cov[10];
    Cov(1,4)=cov[11];
    Cov(2,4)=cov[12];
    Cov(3,4)=cov[13];
    Cov(4,4)=cov[14];

    double charge = H.invR() > 0 ? 1. : -1;

    vector<int> dummy;
    //for (int i=0; i<12; ++i) dummy.push_back(0);

    vertex_lcfi::Track* MyTrack = new vertex_lcfi::Track(MyEvent, H, Mom, charge, Cov, dummy, (void*)track);
    vertex_lcfi::MemoryManager<vertex_lcfi::Track>::Event()->registerObject(MyTrack);

    return MyTrack;
  }

  vertex_lcfi::Event* LcfiInterface::lcfiEvent(Event* event, vertex_lcfi::Vertex* ipVertex) const
  {
    vertex_lcfi::Event* ret;

		if (ipVertex) {
			ret = new vertex_lcfi::Event(ipVertex);
		} else {
			Vector3 ipPos;
			SymMatrix3x3 ipErr;
			ipPos.x() = 0.;
			ipPos.y() = 0.;
			ipPos.z() = 0.;
			ipErr(0,0) = 2.5e-04;
			ipErr(1,0) = 0.;
			ipErr(1,1) = 2.5e-04;
			ipErr(2,0) = 0.;
			ipErr(2,1) = 0.;
			ipErr(2,2) = 0.004;

			ret = new vertex_lcfi::Event(ipPos,ipErr);
		}

		if(event){
			if (debug) {
	    	printf("lcfiEvent: converting %d tracks\n",event->getTracks().size());
			}
	    for (vector<Track*>::const_iterator iter = event->getTracks().begin();
	        iter != event->getTracks().end(); ++iter) {
	      vertex_lcfi::Track* t = lcfiTrack(ret, *iter);
	      ret->addTrack(t);
	    }
		}

    return ret;
  }

  ////////

  Vertex*
  LcfiInterface::flavtagVertex(vertex_lcfi::Vertex* lcfiVertex) const {
    Vertex* vertex = new Vertex();

    vertex->_chi2 = lcfiVertex->chi2();
    vertex->_prob = lcfiVertex->probability();
    vertex->_x = lcfiVertex->position().x();
    vertex->_y = lcfiVertex->position().y();
    vertex->_z = lcfiVertex->position().z();
    vertex->_cov[0] = lcfiVertex->positionError()(0,0);
    vertex->_cov[1] = lcfiVertex->positionError()(1,0);
    vertex->_cov[2] = lcfiVertex->positionError()(1,1);
    vertex->_cov[3] = lcfiVertex->positionError()(2,0);
    vertex->_cov[4] = lcfiVertex->positionError()(2,1);
    vertex->_cov[5] = lcfiVertex->positionError()(2,2);

    for (vector<vertex_lcfi::Track*>::const_iterator iter = lcfiVertex->tracks().begin();
        iter < lcfiVertex->tracks().end(); ++iter) {
      vertex->add( (Track*) (*iter)->trackingNum() );
    }

		// make sure vertex contains at least two tracks
		if(debug) fprintf(stderr,"returning flavtagVertex (%.2e,%.2e,%.2e) with ntrks=%d\n",
				vertex->_x, vertex->_y, vertex->_z, vertex->getTracks().size());
		assert( vertex->getTracks().size() >= 2 );
		assert( vertex->getTracks().size() != 0 );
		if (vertex->getTracks().size() == 0) {
			fprintf(stderr,"error: vertex has no tracks!!\n");
			fprintf(stderr," lcfiVertex->tracks().size()=%d\n",lcfiVertex->tracks().size());
			exit(1);
		}

    return vertex;
  }

	vertex_lcfi::Vertex*
  LcfiInterface::lcfiVertex(Vertex* flavtagVertex, bool isPrimary) const {

		//Vertex(Event* Event, const std::vector<Track*> & Tracks, const Vector3 & Position, const SymMatrix3x3 & PosError, bool IsPrimary, double Chi2, double Probability, std::map<Track*,double> ChiTrack);

		vector<vertex_lcfi::Track*> tracks;
		Vector3 pos;
		SymMatrix3x3 cov;
		double chi2 = flavtagVertex->_chi2;
		double prob = flavtagVertex->_prob;

		pos.x() = flavtagVertex->_x;
		pos.y() = flavtagVertex->_y;
		pos.z() = flavtagVertex->_z;

		cov(0,0) = flavtagVertex->_cov[0];
		cov(1,0) = flavtagVertex->_cov[1];
		cov(1,1) = flavtagVertex->_cov[2];
		cov(2,0) = flavtagVertex->_cov[3];
		cov(2,1) = flavtagVertex->_cov[4];
		cov(2,2) = flavtagVertex->_cov[5];

    for (vector<Track*>::const_iterator iter = flavtagVertex->getTracks().begin();
        iter < flavtagVertex->getTracks().end(); ++iter) {
			tracks.push_back( lcfiTrack(*iter) );
    }

		vertex_lcfi::Vertex* vertex = new vertex_lcfi::Vertex(
				_event, tracks, pos, cov, isPrimary, chi2, prob);

    return vertex;
  }

  vector<Vertex*>
  LcfiInterface::flavtagVertices(vertex_lcfi::DecayChain* chain) const {
    vector<Vertex*> vtxList;

    // skip the first (primary) vertex
    assert(chain->vertices().size() > 0);
    for (vector<vertex_lcfi::Vertex*>::const_iterator iter = ++(chain->vertices().begin());
        iter < chain->vertices().end();++iter) {
      vtxList.push_back(flavtagVertex(*iter));
    }

    return vtxList;
  }

  bool LcfiInterface::passesCut(const Track* trk, const SecondaryVertexConfig& cfg) {
    if (fabs(trk->getD0()) > cfg.maxD0) {
      //printf("d0 cut fail: %f\n",trk->par[Track::d0]);
      return false;
    }
    if (trk->getCovMatrix()[tpar::d0d0] > cfg.maxD0Err) {
      //printf("d0d0 cut fail: %f\n",trk->cov[Track::d0d0]);
      return false;
    }
    if (fabs(trk->getZ0()) > cfg.maxZ0) {
      //printf("z0 cut fail: %f\n",trk->par[Track::z0]);
      return false;
    }
    if (trk->getCovMatrix()[tpar::z0z0] > cfg.maxZ0Err) {
      //printf("z0z0 cut fail: %f\n",trk->cov[Track::z0z0]);
      return false;
    }
    if (trk->Pt() < cfg.minPt) {
      return false;
    }
    /*
    if (trk->getRadiusOfInnermostHit() > cfg.maxInnermostHitRadius) {
      return false;
    }
    */

    // OR cuts
    if (trk->getTpcHits() >= cfg.minTpcHits) {
      return true;
    }
    if (trk->getFtdHits() >= cfg.minFtdHits) {
      return true;
    }
    if (trk->getVtxHits() >= cfg.minVtxHitsWithoutTpcFtd) {
      return true;
    }
    /*
    if (trk->getVtxHits() + trk->getFtdHits() >= cfg.minVtxPlusFtdHits) {
      return true;
    }
    */

    return false;
  }

  ///////////////////////
#include <marlin/Global.h>
#include <gear/BField.h>

  Vertex*
  LcfiInterface::findPrimaryVertex() {
		if (debug) fprintf(stderr,"findPrimaryVertex()\n");
    if (_primaryVertex) {
			if (debug) fprintf(stderr,"returning cached primary vertex\n");
      return flavtagVertex(_primaryVertex);
    }

    vertex_lcfi::Vertex* ipResult = LcfiInstance::getInstance().getIpFitter()->calculateFor(_event);
    _primaryVertex = ipResult;

    return flavtagVertex(ipResult);
  }

  /*
  void LcfiInterface::probMap( map<Track*,float>& probMap ) {
    vector<vertex_lcfi::TrackState*> TrackStates;
    for (vector<vertex_lcfi::Track*>::const_iterator iTrack = _event->tracks().begin();
        iTrack != _event->tracks().end(); ++iTrack)
    {
      TrackStates.push_back((*iTrack)->makeState());
    }

    vertex_lcfi::ZVTOP::VertexFitterKalman MyFitter;
    vertex_lcfi::ZVTOP::CandidateVertex CVertex(TrackStates, 0, &MyFitter);

    map<vertex_lcfi::TrackState*,double> tsProbMap;
    CVertex.mapProb(tsProbMap);

    for (map<vertex_lcfi::TrackState*,double>::iterator it = tsProbMap.begin(); it != tsProbMap.end(); ++it) {
      if (it->first) {
        probMap.insert( make_pair( (Track*)it->first->parentTrack()->trackingNum(), it->second ) );
      }
    }
  }
  */

  vector<Vertex*>
  LcfiInterface::findSecondaryVertices(Jet* jet, const SecondaryVertexConfig& cfg) {
    LcfiInstance::getInstance().getZVRES()->setDoubleParameter("TwoProngCut", cfg.TwoProngCut);
    LcfiInstance::getInstance().getZVRES()->setDoubleParameter("TrackTrimCut", cfg.TrackTrimCut);
		LcfiInstance::getInstance().getZVRES()->setDoubleParameter("ResolverCut",cfg.ResolverCut);

    assert(_primaryVertex != 0);

    vertex_lcfi::Jet* lcfiJet = new vertex_lcfi::Jet(_event, vector<vertex_lcfi::Track*>(),
        jet->E(),Vector3(jet->Px(), jet->Py(), jet->Pz()),0);
    vertex_lcfi::MemoryManager<vertex_lcfi::Jet>::Event()->registerObject(lcfiJet);

    const vector<Track*> tracks = jet->getTracks();
    const vector<Neutral*> neutrals = jet->getNeutrals();

    int usedTracks(0);
    for (vector<Track*>::const_iterator iter = tracks.begin(); iter != tracks.end(); ++iter) {
      if (passesCut(*iter,cfg)) {
        lcfiJet->addTrack( lcfiTrack(_event, *iter) );
        ++usedTracks;
      }
    }

		if (debug) {
			printf("finding secondary vertices for jet: (E=%.2e, px=%.2e, py=%.2e, pz=%.2e), ntrk=%d (%d passses cut)\n",
					jet->E(),jet->Px(),jet->Py(),jet->Pz(),tracks.size(),usedTracks);
		}

    // default _JetWeightingEnergyScaling = 5.0/40.0
    LcfiInstance::getInstance().getZVRES()->setDoubleParameter("Kalpha", 0.125 * lcfiJet->energy());

    vertex_lcfi::DecayChain* zvtopRes = LcfiInstance::getInstance().getZVRES()->calculateFor(lcfiJet);

    vector<Vertex*> vtxList = flavtagVertices(zvtopRes);
    return vtxList;
  }

  vector<Vertex*>
  LcfiInterface::forceZvtop(const Jet& jet) {
    LcfiInstance::getInstance().getZVRES()->setDoubleParameter("TwoProngCut", 1e10);
    LcfiInstance::getInstance().getZVRES()->setDoubleParameter("TrackTrimCut", 1e10);
		//LcfiInstance::getInstance().getZVRES()->setDoubleParameter("ResolverCut",1e10);
		LcfiInstance::getInstance().getZVRES()->setDoubleParameter("ResolverCut",0.6);

    vertex_lcfi::Jet* lcfiJet = new vertex_lcfi::Jet(_event, vector<vertex_lcfi::Track*>(),
        jet.E(),Vector3(jet.Px(), jet.Py(), jet.Pz()),0);
    vertex_lcfi::MemoryManager<vertex_lcfi::Jet>::Event()->registerObject(lcfiJet);

		const vector<Track*>& tracks = jet.getTracks();
    for (vector<Track*>::const_iterator iter = tracks.begin(); iter != tracks.end(); ++iter) {
			lcfiJet->addTrack( lcfiTrack(_event, *iter) );
    }

    //printf("lcfiJet tracks: %d\n", usedTracks);

    // default _JetWeightingEnergyScaling = 5.0/40.0
    //LcfiInstance::getInstance().getZVRES()->setDoubleParameter("Kalpha", 0.125 * lcfiJet->energy());
    LcfiInstance::getInstance().getZVRES()->setDoubleParameter("Kalpha", 0);

    vertex_lcfi::DecayChain* zvtopRes = LcfiInstance::getInstance().getZVRES()->calculateFor(lcfiJet);

    vector<Vertex*> vtxList = flavtagVertices(zvtopRes);
    return vtxList;
  }

	double LcfiInterface::getChi2TrackVtx(Vertex *vtx, Track *trk) const
	{
		vertex_lcfi::Track* ptr = lcfiTrack(trk);
		vertex_lcfi::TrackState *pstate = ptr->makeState();
		vertex_lcfi::TState state(pstate);

		// extract vertex
		//* v = [xyz], Cv=[Cxx,Cxy,Cyy,Cxz,Cyz,Czz]-covariance matrix
		double v[3], Cv[6];
		v[0] = vtx->getX(); v[1] = vtx->getY(); v[2] = vtx->getZ();
		for(int i=0;i<6;i++)
		{
			Cv[i] = vtx->getCov()[i];
		}

		vertex_lcfi::ZVTOP::VertexFitterKalman kal;
		double chi2 = kal.getDeviationFromVertex(&state, v, Cv);

		return chi2;
	}

  double LcfiInterface::vertexMassPtCorrection( Vertex* secondary, Vertex* primary, const TVector3& momentum, float sigmax ) const
	{
		vertex_lcfi::Vector3 mom;
		mom[0] = momentum.X();
		mom[1] = momentum.Y();
		mom[2] = momentum.Z();

		vertex_lcfi::Vertex* primary2 = lcfiVertex(primary,true);
		vertex_lcfi::Vertex* secondary2 = lcfiVertex(secondary);

		vertex_lcfi::VertexMass vm;

		//VertexMass::Ptcalc(  Vertex* IPVertex  , Vertex*  TheVertex, Vector3* momentum , float sigmax  ) const
		double pt = vm.Ptcalc(  primary2, secondary2, &mom, sigmax );
		return pt;
	}

}
