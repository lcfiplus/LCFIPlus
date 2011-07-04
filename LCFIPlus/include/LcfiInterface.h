#ifndef LcfiInterface_h
#define LcfiInterface_h 1

#include <inc/track.h>
#include <inc/trackstate.h>
#include <inc/event.h>
#include <inc/vertex.h>
#include <util/inc/memorymanager.h>
#include <inc/decaychain.h>
#include <algorithm>

#include <algo/inc/pereventipfitter.h>
#include <algo/inc/zvres.h>

#include "JetFinder.h"

using namespace std;

namespace flavtag {
  struct SecondaryVertexConfig {
    // cuts which are combined using the AND scheme
    float maxD0;
    float maxD0Err;
    float maxZ0;
    float maxZ0Err;
    float minPt;
    float maxInnermostHitRadius;
    // cuts which are combined using the OR scheme, then AND'd with the AND schemes above
    int minTpcHits;
    int minFtdHits;
    int minVtxHitsWithoutTpcFtd;
    int minVtxPlusFtdHits;
    // algorithm config
    float TwoProngCut;
    float TrackTrimCut;
		float ResolverCut;
    SecondaryVertexConfig() : TwoProngCut(10.), TrackTrimCut(10.), ResolverCut(0.6) {}
  };

  class LcfiInstance {
    private:
      LcfiInstance();
      ~LcfiInstance();

      static LcfiInstance _instance;

    public:
      static LcfiInstance& getInstance() { return _instance; };
      const vertex_lcfi::PerEventIPFitter* getIpFitter();
      vertex_lcfi::ZVRES* getZVRES();

    private:
      vertex_lcfi::PerEventIPFitter* _ipFitter;
      vertex_lcfi::ZVRES* _zvres;
  };

  class LcfiInterface {
    public:
			template<class Iterator> friend class VertexFitterLCFI;
 
     LcfiInterface(Event* event = NULL, Vertex* primaryVertex=0);
      ~LcfiInterface();

      Vertex* findPrimaryVertex();
      vector<Vertex*> findSecondaryVertices(Jet* jet, const SecondaryVertexConfig& cfg);
      //void probMap( map<Track*,float>& probMap );
			vector<Vertex*> forceZvtop(const Jet& jet);

			double getChi2TrackVtx(Vertex *vtx, Track *trk) const;
      bool passesCut(const Track* trk, const SecondaryVertexConfig& cfg);
			double vertexMassPtCorrection( Vertex* secondary, Vertex* primary, const TVector3& momentum, float sigmax ) const;
			bool debug;

    private:
      vertex_lcfi::Event* _event;
      vertex_lcfi::Vertex* _primaryVertex;

      vertex_lcfi::Event* lcfiEvent(Event* event, vertex_lcfi::Vertex* ipVertex=0) const;
      vertex_lcfi::Track* lcfiTrack(vertex_lcfi::Event* MyEvent, Track* track) const;
			vertex_lcfi::Vertex* lcfiVertex(Vertex* flavtagVertex, bool isPrimary=false) const;
      vertex_lcfi::Track* lcfiTrack(Track* track) const {return lcfiTrack(_event, track);}
      Vertex* flavtagVertex(vertex_lcfi::Vertex* MyLCFIVertex) const;
      vector<Vertex*> flavtagVertices(vertex_lcfi::DecayChain* chain) const;

  };

}

#endif
