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
	void connectVerticesToJets(const JetVec &jets, const vector<Vertex *> &vtcs, vector<vector<Vertex *> > &jetVertices, vector<vector<const Track *> > &jetResidualTracks, const Vertex *ip = 0);
	vector<const Track *> extractTracks(VertexVec &vtx);
	double calcThrust( vector<TVector3>& list, TVector3 &taxis );

}}

#endif
