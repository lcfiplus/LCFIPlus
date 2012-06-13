// algoSigProb.h

#ifndef AlgoSigProb_h
#define AlgoSigProb_h 1

#include "lcfiplus.h"
#include "JetFinder.h"

namespace lcfiplus{
namespace algoSigProb{

extern double trackD0Significance(const Track* trk, const Vertex* pri);
extern double trackZ0Significance(const Track* trk, const Vertex* pri);
extern double signedD0Significance(const Track* trk, const Jet* jet, const Vertex* pri, bool updateFlt=false);
extern double signedZ0Significance(const Track* trk, const Jet* jet, const Vertex* pri, bool updateFlt=false);
extern double signedD0(const Track* trk, const Jet* jet, const Vertex* pri, bool updateFlt=false);
extern double signedZ0(const Track* trk, const Jet* jet, const Vertex* pri, bool updateFlt=false);
extern void findMostSignificantTrack(const Jet* jet, const Vertex* pri, double sigVec[6]);

extern double prob1D(double sig, double maxsig, double* pars);
extern double trackProbD0(const Track* trk, const Vertex* pri);
extern double trackProbZ0(const Track* trk, const Vertex* pri);
extern double jointProbD0(const Jet* jet, const Vertex* pri, double maxd0sigcut = 1e+300);
extern double jointProbZ0(const Jet* jet, const Vertex* pri, double maxz0sigcut = 1e+300);

}}

#endif
