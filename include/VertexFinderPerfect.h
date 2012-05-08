// VertexFinderPerfect.h

#ifndef VertexFinderPerfect_h
#define VertexFinderPerfect_h 1

#include "lcfiplus.h"
#include <list>
#include <vector>

using namespace std;

namespace lcfiplus {

	namespace VertexFinderPerfect{
		void findPerfectVertices(TrackVec &tracks, MCParticleVec &mcp, vector<MCVertex *> &selvtx, int minimumRecoTracks, double minimumDistance = 0.1);
	}
}

#endif
