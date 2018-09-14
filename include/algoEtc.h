// algoEtc.h

#ifndef AlgoEtc_h
#define AlgoEtc_h 1

#include "lcfiplus.h"
#include "JetFinder.h"

namespace lcfiplus {
namespace algoEtc {

// beam pseudo-tracks for primary vertex fitter
void makeBeamTracks(Track*& t1, Track*& t2, bool smear=true);
void makeBeamVertex(Vertex*& vtx, bool smear=true);
void connectVerticesToJets(const JetVec& jets, const vector<Vertex*>& vtcs, vector<vector<Vertex*> >& jetVertices, vector<vector<const Track*> >& jetResidualTracks, const Vertex* ip = 0);
vector<const Track*> extractTracks(VertexVec& vtx);
double calcThrust( vector<TVector3>& list, TVector3& taxis );

bool SimpleSecMuonFinder(const Track* tr, double d0sigth, double z0sigth, double maxpos, double mudepmin,
                         double ecaldepmin, double ecaldepmax, double hcaldepmin, double hcaldepmax, double maxclusterpertrackenergy = 10.,
                         const Vertex* ip = 0);
bool SimpleSecElectronFinder(const Track* tr, double d0sigth, double z0sigth, double maxpos, double emin,
                             double minfracecal, double minecalpertrackenergy, double maxecalpertrackenergy,
                             const Vertex* ip = 0);

}
}

#endif
