// VertexFinderSuehara.h

#ifndef VertexFinderSuehara_h
#define VertexFinderSuehara_h 1

#include "lcfiplus.h"
#include "VertexFitterLCFI.h"
#include "VertexFitterSimple.h"
#include <list>
#include <vector>

#include "algoSigProb.h"

using namespace std;
using namespace lcfiplus::algoSigProb;

namespace lcfiplus {

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
		vector<lcfiplus::Vertex*> * findCandidates(TrackVec &tracks, double chi2th, double massth, double ipsigth);

		// find one vertex with build-up method
		lcfiplus::Vertex* findOne(list<const Track *> &tracks, double chi2th, double massth, bool removeTracks);
		// VertexFitterSimple version
		lcfiplus::Vertex* findOne2(list<const Track *> &tracks, double chi2th, double massth, bool removeTracks);

		// compare functions
		bool VertexNearer(const Vertex *vtx1, const Vertex *vtx2);
		bool VertexProbLarger(const Vertex *vtx1, const Vertex *vtx2);

		// obtain vertex list
		void GetVertexList(list<const Track *> &tracks, vector<Vertex *> &vtx, double chi2th, double massth, double posth = 0.3, double chi2orderinglimit = 1.0);

		// associating tracks to an existing vertex
		lcfiplus::Vertex * associateTracks(Vertex *vertex, list<const Track *> &tracks, double chi2th, double massth, list<const Track *> *residualTracks = 0);
		void associateIPTracks(vector<Vertex *> &vertices, double minimumdist = 0., int chi2mode = 0, double chi2ratio = 2.0);

		void buildUp(TrackVec &tracks, vector<Vertex *> &vtx, double chi2thpri, double chi2thsec, double massth, double posth = 0.3, double chi2orderinglimit = 1.0, Vertex *ip = 0);
		void buildUpForJetClustering(TrackVec &tracks, vector<Vertex *> &vtx);

		vector<Vertex *> makeSingleTrackVertices(Jet *jet, TrackVec &tracks, Vertex *ip, double minpos, double maxpos, double maxangle, double max_separation_per_pos);

	}

	//vector<lcfiplus::Vertex*> * findSueharaVertices(const Event& evt, const Jet& jet);
}

#endif
