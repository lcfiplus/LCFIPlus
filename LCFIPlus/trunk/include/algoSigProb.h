// algoSigProb.h

#ifndef AlgoSigProb_h
#define AlgoSigProb_h 1

#include "lcfiplus.h"
#include "JetFinder.h"

namespace lcfiplus{
namespace algoSigProb{

extern float trackD0Significance(const Track* trk, const Vertex* pri);
extern float trackZ0Significance(const Track* trk, const Vertex* pri);
extern float signedD0Significance(const Track* trk, const Jet* jet, const Vertex* pri, bool updateFlt=false);
extern float signedZ0Significance(const Track* trk, const Jet* jet, const Vertex* pri, bool updateFlt=false);
extern void findMostSignificantTrack(const Jet* jet, const Vertex* pri, float sigVec[6]);

extern float prob1D(float sig, float maxsig, float* pars);
extern float trackProbD0(const Track* trk, const Vertex* pri);
extern float trackProbZ0(const Track* trk, const Vertex* pri);
extern float jointProbD0(const Jet* jet, const Vertex* pri);
extern float jointProbZ0(const Jet* jet, const Vertex* pri);

}}

#endif
