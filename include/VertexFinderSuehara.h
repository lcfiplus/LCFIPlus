// VertexFinderSuehara.h

#ifndef VertexFinderSuehara_h
#define VertexFinderSuehara_h 1

#include "lcfiplus.h"
//#include "VertexFitterLCFI.h"
//#include "VertexFitterSimple.h"
#include <list>
#include <vector>
#include <functional>

#include "algoSigProb.h"
#include "VertexSelector.h"

using namespace std;
using namespace lcfiplus::algoSigProb;

namespace lcfiplus {

namespace VertexFinderSuehara {
class SortTracksByIPSig { // decending order
 public:
  SortTracksByIPSig(Vertex* v) {
    priVertex = v;
  }
  bool operator() (const Track* p1, const Track* p2) {
//					double ipsig1 = sqrt(pow(trackD0Significance(p1, priVertex),2) + pow(trackZ0Significance(p1,priVertex),2));
//					double ipsig2 = sqrt(pow(trackD0Significance(p2, priVertex),2) + pow(trackZ0Significance(p2,priVertex),2));
    double ipsig1 = sqrt(pow(p1->getD0() / p1->getCovMatrix()[tpar::d0d0],2) + pow(p1->getZ0() / p1->getCovMatrix()[tpar::z0z0],2));
    double ipsig2 = sqrt(pow(p2->getD0() / p2->getCovMatrix()[tpar::d0d0],2) + pow(p2->getZ0() / p2->getCovMatrix()[tpar::z0z0],2));
    return ipsig1 > ipsig2;
  }

  Vertex* priVertex;
};
/*
		// find 2-track vertices
		// track, chi-square threshold, maximum invmass, minimum impact parameter significance of tracks (either)
		vector<lcfiplus::Vertex*> * findCandidates(TrackVec &tracks, double chi2th, double massth, double ipsigth);

		// find one vertex with build-up method
		lcfiplus::Vertex* findOne(list<const Track *> &tracks, double chi2th, double massth, bool removeTracks);
		// VertexFitterSimple version
		lcfiplus::Vertex* findOne2(list<const Track *> &tracks, double chi2th, double massth, bool removeTracks);
*/
// compare functions
bool VertexNearer(const Vertex* vtx1, const Vertex* vtx2);
bool VertexProbLarger(const Vertex* vtx1, const Vertex* vtx2);

class VertexFinderSueharaConfig {
 public:
  // main parameters
  VertexSelectorConfig v0selTrack; // selector for tight v0 selection: used for track rejection
  VertexSelectorConfig v0selVertex;// selector for loose v0 selection: used for vertex rejection
  double chi2th;	// track chi2 threshold to accept vertices
  double chi2thV0SelTrack; // track chi2 threshold to reject tracks in v0 rejection
  double massth;	// maximum mass to accept vertices
  double chi2orderinglimit; // chi2 threshold to order distance rather than chi2 value

  // associateIPtracks parameters
  double minimumdistIP; // minimum distance to associate IP tracks
  double chi2ratioIP; 	// bias factor to associate IP rather than secondary

  // singletrackvertex parameters
  double minPosSingle;
  double maxPosSingle;
  double minEnergySingle;
  double maxAngleSingle;
  double maxSeparationPerPosSingle;
  double mind0SigSingle;
  double minz0SigSingle;

  //flg for AVF/chi2 algorithm
  bool avf;
  double temperature;   //parameter for avf

  //flg for BNess tagger fake rejection
  bool useBNess;
  double cutBNess;  //parameter for BNess
  double cutBNessE1;  //parameter for BNess

  // default values
  VertexFinderSueharaConfig() {
    v0selTrack.setV0Tight();
    v0selTrack.rejectdist = true;
    v0selTrack.rejectdistnegative = true;
    v0selTrack.rejectdistor = true;
    v0selTrack.minpos = 0.5;

    chi2thV0SelTrack = 2.;

    v0selVertex.setV0Loose();
    v0selVertex.rejectdist = true;
    v0selVertex.minpos = 0.3;

    chi2th = 9.;
    massth = 10.;
    chi2orderinglimit = 1.;

    minimumdistIP = 0.; // not used in default
    chi2ratioIP = 2.;		// biased twice for primary

    minPosSingle = 0.3;
    maxPosSingle = 30.;
    minEnergySingle = 1.;
    maxAngleSingle = 0.5;
    maxSeparationPerPosSingle = 0.1;
    mind0SigSingle = 5.;
    minz0SigSingle = 5.;

    avf = false;   //default: use chi2 algorithm
    temperature = 5.0;

    useBNess = false; //default: do not use BNess
    cutBNess = -0.80;
    cutBNessE1 = -0.15;
  }
};

// obtain vertex list
void GetVertexList(list<const Track*>& tracks, const Vertex* ip, vector<Vertex*>& vtx, vector<Vertex*>& v0vtx, VertexFinderSueharaConfig& cfg);

// associating tracks to an existing vertex
lcfiplus::Vertex* associateTracks(Vertex* vertex, const VertexVec& v0vtx, list<const Track*>& tracks, VertexFinderSueharaConfig& cfg, list<const Track*>* residualTracks = 0);
void associateIPTracks(vector<Vertex*>& vertices, Vertex* ip, VertexFinderSueharaConfig& cfg);
//using AVF method
void associateIPTracksAVF(vector<Vertex*>& vertices, Vertex* ip, VertexFinderSueharaConfig& cfg);

void buildUp(TrackVec& tracks, vector<Vertex*>& vtx, vector<Vertex*>& v0vtx, double chi2thpri, VertexFinderSueharaConfig& cfg, Vertex** ip = 0);
void buildUpForJetClustering(TrackVec& tracks, vector<Vertex*>& vtx);

vector<Vertex*> makeSingleTrackVertices
(VertexVec& vtcs, TrackVec& tracks, VertexVec& v0vtx, const Vertex* ip, VertexFinderSueharaConfig& cfg);
vector<Vertex*> makeSingleTrackVertices
(Jet* jet, TrackVec& tracks, VertexVec& v0vtx, const Vertex* ip, VertexFinderSueharaConfig& cfg);

void recombineVertices(vector<Vertex*>& vertices,vector<Vertex*>& singleVertices);

//new function for BNess tagger
void recombineVertices(vector<Vertex*>& vertices, vector<Vertex*>& singleVertices, VertexFinderSueharaConfig& cfg );

void optimizeTwoVertices(Vertex*& v1, Vertex*& v2, int nvr);

}

//vector<lcfiplus::Vertex*> * findSueharaVertices(const Event& evt, const Jet& jet);
}

#endif
