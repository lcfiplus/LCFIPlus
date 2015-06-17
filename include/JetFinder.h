#ifndef JetFinder_h
#define JetFinder_h 1

#include "lcfiplus.h"
#include "TObject.h"

using namespace std;
using namespace lcfiplus;

namespace lcfiplus {

/**
	Holds parameters for jet clustering algorithms.
	*/
struct JetConfig {
  string algo;
  string algoY;
  int nJet;
  double Ycut;
  int useBeamJets;
  double rParameter;
  double coneR;
  double epsCut;
  string coreAlgo;
  double coreThreshold;
  double distCut;
  int nIteration;
  double YaddVV;
  double YaddVL;
  double YaddLL;
  JetConfig() : nJet(0),Ycut(1),useBeamJets(0),rParameter(1.0),coneR(0),epsCut(0), coreThreshold(0),distCut(0),nIteration(5),YaddVV(100), YaddVL(100), YaddLL(100) {}
};

/**
	Converts a jet containing vertices in a tree structure
	into a jet containing all particles at the top level.

	The vertex information is stripped out.
	*/
Jet* convertJetVertex(const Jet* jet);

/**
	Finds jets using various jet clustering algorithms.

	@author T. Tanabe, ICEPP, The University of Tokyo
	@author T. Suehara, ICEPP, The University of Tokyo
	@version $Id$
	*/
class JetFinder {
 public:
  static double funcDurham(Jet& jet1, Jet& jet2, double Evis2, JetConfig& cfg);
  static double funcJade(Jet& jet1, Jet& jet2, double Evis2, JetConfig& cfg);
  static double funcJadeE(Jet& jet1, Jet& jet2, double Evis2, JetConfig& cfg);
  static double funcDurhamCheat(Jet& jet1, Jet& jet2, double Evis2, JetConfig& cfg);
  static double funcDurhamVertex(Jet& jet1, Jet& jet2, double Evis2, JetConfig& cfg);

  static double funcKt(Jet& jet1, Jet& jet2, double Evis2, JetConfig& cfg);
  static double funcKtVertex(Jet& jet1, Jet& jet2, double Evis2, JetConfig& cfg);

  /** Constructor.
  	@param[in] cfg specify the algorithm and the parameters
  	*/
  JetFinder(const JetConfig& cfg);
  /** Destructor. */
  ~JetFinder() {};
  /** Replace JetConfig.
  	@param[in] cfg specify the algorithm and the parameters
  	*/
  void Configure(const JetConfig& cfg);
  /**
  	Perform jet clustering with charged tracks only. */
  vector<Jet*> run(TrackVec& tracks);
  /**
  	Perform jet clustering using both charged tracks and neutral particles. */
  vector<Jet*> run(TrackVec& tracks, NeutralVec& neutrals, double* pymin = 0, int ynjetmax = 10);
  /**
  	Perform jet clustering using charged tracks and neutral particles
  	with the contraint that the given vertices are not merged.
  	*/
  vector<Jet*> run(TrackVec& tracks, NeutralVec& neutrals, VertexVec& vertices, double* pymin = 0, bool findmu = false, int ynjetmax = 10);
  /**
  	Lepton finder and vertex isolator/combiner, called from run() with vertices
  	*/
  vector<Jet*> prerun(TrackVec& tracks, NeutralVec& neutrals, VertexVec& vertices, bool findmu = false, int* pnVertexJets = 0);
  /**
  	Main function of jet clustering, called from other run() functions.
  	*/
  vector<Jet*> run(vector<Jet*> input, double* pymin = 0, int ynjetmax = 10);

 private:
  double (*_Yfunc)(Jet& jet1, Jet& jet2, double Evis2, JetConfig& cfg);
  JetConfig _cfg;
};

/*
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
*/
}

#endif
