// algoSigProb.h

#ifndef AlgoSigProb_h
#define AlgoSigProb_h 1

#include "flavtag.h"
#include "JetFinder.h"

namespace flavtag{
namespace algoSigProb{

extern float trackD0Significance(const Track* trk, const Vertex* pri);
extern float trackZ0Significance(const Track* trk, const Vertex* pri);
extern float signedD0Significance(const Track* trk, const Jet* jet, const Vertex* pri, bool updateFlt=false);
extern float signedZ0Significance(const Track* trk, const Jet* jet, const Vertex* pri, bool updateFlt=false);
extern void findMostSignificantTrack(Jet* jet, Vertex* pri, float sigVec[6]);

extern float prob1D(float sig, float maxsig, float* pars);
extern float trackProbD0(Track* trk, Vertex* pri);
extern float trackProbZ0(Track* trk, Vertex* pri);
extern float jointProbD0(Jet* jet, Vertex* pri);
extern float jointProbZ0(Jet* jet, Vertex* pri);

}}

#endif
