#ifndef JetFinder_h
#define JetFinder_h 1

#include "lcfiplus.h"
#include "TObject.h"

using namespace std;
using namespace lcfiplus;

namespace lcfiplus {

  struct JetConfig {
    string algo;
    string algoY;
    int nJet;
    float Ycut;
    float coneR;
    float epsCut;
    string coreAlgo;
    float coreThreshold;
    float distCut;
    int nIteration;
    JetConfig() : nJet(0),Ycut(1),coneR(0),epsCut(0), coreThreshold(0),distCut(0),nIteration(5) {}
  };

	Jet* convertJetVertex(const Jet* jet);

  class JetFinder {
    public:
      static double funcDurham(Jet& jet1, Jet& jet2, double Evis2);
      static double funcJade(Jet& jet1, Jet& jet2, double Evis2);
      static double funcJadeE(Jet& jet1, Jet& jet2, double Evis2);
      static double funcDurhamCheat(Jet& jet1, Jet& jet2, double Evis2);
      static double funcDurhamVertex(Jet& jet1, Jet& jet2, double Evis2);

      JetFinder(const JetConfig& cfg);
      ~JetFinder() {};
      vector<Jet*> run(TrackVec &tracks);
      vector<Jet*> run(TrackVec &tracks, NeutralVec &neutrals, double *pymin = 0);
      vector<Jet*> run(TrackVec &tracks, NeutralVec &neutrals, VertexVec &vertices, double *pymin = 0, bool findmu = false);
      vector<Jet*> run(vector<Jet*> input, double *pymin = 0);

    private:
      double (*_Yfunc)(Jet& jet1, Jet& jet2, double Evis2);
      JetConfig _cfg;
  };

  class CheatedJetFinder {
    public:
      CheatedJetFinder(const JetConfig& cfg);
      ~CheatedJetFinder() {};
      vector<Jet*> run(Event* event);
      vector<MCParticle*> originalPartons( MCParticleVec & mcps );

    private:
      double (*_Yfunc)(Jet& jet1, Jet& jet2, double Evis2);
      JetConfig _cfg;
  };
}

#endif
