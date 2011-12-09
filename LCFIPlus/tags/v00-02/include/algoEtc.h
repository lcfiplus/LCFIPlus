// algoEtc.h

#ifndef AlgoEtc_h
#define AlgoEtc_h 1

#include "lcfiplus.h"
#include "JetFinder.h"

namespace lcfiplus{
namespace algoEtc{

	// beam pseudo-tracks for primary vertex fitter
	void makeBeamTracks(Track *&t1, Track *&t2);
	void makeBeamVertex(Vertex *&vtx);

}}

#endif
