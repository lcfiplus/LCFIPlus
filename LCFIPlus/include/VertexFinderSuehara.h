// VertexFinderSuehara.h

#ifndef VertexFinderSuehara_h
#define VertexFinderSuehara_h 1

#include "flavtag.h"
#include "VertexFitterLCFI.h"
#include "VertexFitterSimple.h"
#include <list>
#include <vector>

#include "algoSigProb.h"

using namespace std;
using namespace flavtag::algoSigProb;

namespace flavtag {

	namespace VertexFinderSuehara{
		class SortTracksByIPSig{ // decending order
			public:
				SortTracksByIPSig(Vertex *v){priVertex = v;}
				bool operator() (const Track * p1, const Track * p2)
				{
//					double ipsig1 = sqrt(pow(trackD0Significance(p1, priVertex),2) + pow(trackZ0Significance(p1,priVertex),2));
//					double ipsig2 = sqrt(pow(trackD0Significance(p2, priVertex),2) + pow(trackZ0Significance(p2,priVertex),2));
					double ipsig1 = sqrt(pow(p1->getD0() / p1->getCovMatrix()[tpar::d0d0],2) + pow(p1->getZ0() / p1->getCovMatrix()[tpar::z0z0],2));
					double ipsig2 = sqrt(pow(p2->getD0() / p2->getCovMatrix()[tpar::d0d0],2) + pow(p2->getZ0() / p2->getCovMatrix()[tpar::z0z0],2));
					return ipsig1 > ipsig2;
				}

				Vertex *priVertex;
		};

		// find 2-track vertices
		// track, chi-square threshold, maximum invmass, minimum impact parameter significance of tracks (either)
		vector<flavtag::Vertex*> * findCandidates(const vector<Track *> &tracks, double chi2th, double massth, double ipsigth);

		// find one vertex with build-up method
		flavtag::Vertex* findOne(list<Track *> &tracks, double chi2th, double massth, bool removeTracks);
		// VertexFitterSimple version
		flavtag::Vertex* findOne2(list<Track *> &tracks, double chi2th, double massth, bool removeTracks);

		// compare functions
		bool VertexNearer(const Vertex *vtx1, const Vertex *vtx2);
		bool VertexProbLarger(const Vertex *vtx1, const Vertex *vtx2);

		// obtain vertex list
		void GetVertexList(list<Track *> &tracks, vector<Vertex *> &vtx, double chi2th, double massth, double posth = 0.3, double chi2orderinglimit = 1.0);

		// associating tracks to an existing vertex
		flavtag::Vertex * associateTracks(Vertex *vertex, list<Track *> &tracks, double chi2th, double massth, list<Track *> *residualTracks = 0);
		void associateIPTracks(vector<Vertex *> &vertices, double minimumdist = 0., int chi2mode = 0, double chi2ratio = 2.0);

		void buildUp(const vector<Track *> &tracks, vector<Vertex *> &vtx, double chi2thpri, double chi2thsec, double massth, double posth = 0.3, double chi2orderinglimit = 1.0, Vertex *ip = 0);
		void buildUpForJetClustering(const vector<Track *> &tracks, vector<Vertex *> &vtx);
	}

	//vector<flavtag::Vertex*> * findSueharaVertices(const Event& evt, const Jet& jet);
}

#endif
