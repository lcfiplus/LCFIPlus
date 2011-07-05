#ifndef JetFinder_h
#define JetFinder_h 1

#include "flavtag.h"
#include "TObject.h"

using namespace std;
using namespace flavtag;

namespace flavtag {

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

  class Jet : public TLorentzVector {
    public:
      static double funcDurham(Jet& jet1, Jet& jet2, double Evis2);
      static double funcJade(Jet& jet1, Jet& jet2, double Evis2);
      static double funcJadeE(Jet& jet1, Jet& jet2, double Evis2);
      static double funcDurhamCheat(Jet& jet1, Jet& jet2, double Evis2);
      static double funcDurhamVertex(Jet& jet1, Jet& jet2, double Evis2);

      // constructors
      Jet() : _id(-1) {};
      Jet(Track* trk);
      Jet(Neutral* neutral);
			Jet(Vertex *vtx) : _id(-1){add(vtx);}
      ~Jet() {};

			void setId(int id){_id = id;}
			int getId()const {return _id;}

      const vector<Track*>& getTracks() const { return _tracks; }
      const vector<Neutral*>& getNeutrals() const { return _neutrals; }
      const vector<Vertex*>& getVertices() const { return _vertices; }

      // methods
      void add(const Jet& jet);
      void add(Track *trk){_tracks.push_back(trk); *(TLorentzVector *)this += *(TLorentzVector *)trk;}
      void add(Neutral *neut){_neutrals.push_back(neut); *(TLorentzVector *)this += *(TLorentzVector *)neut;}
      void add(Vertex *vtx){
				_vertices.push_back(vtx);
				for(unsigned int i=0;i<vtx->getTracks().size();i++){
					*(TLorentzVector *)this += *(TLorentzVector *)(vtx->getTracks()[i]);
				}
			}
      //void recalculate();

      static int sort(const Jet* a, const Jet* b) {
        return (a->E() > b->E());
      }

			double sphericity() const;

    private:
      vector<Track*> _tracks;
      vector<Neutral*> _neutrals;
			vector<Vertex*> _vertices;

			int _id;

			ClassDefNV(Jet, 1);
  };

	Jet* convertJetVertex(Jet* jet);

  class JetFinder {
    public:
      JetFinder(const JetConfig& cfg);
      ~JetFinder() {};
      vector<Jet*> run(vector<Track*> tracks);
      vector<Jet*> run(vector<Track*> tracks, vector<Neutral*> neutrals, double *pymin = 0);
      vector<Jet*> run(vector<Track*> tracks, vector<Neutral*> neutrals, vector<Vertex *> vertices, double *pymin = 0, bool findmu = false);
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
      vector<MCParticle*> originalPartons( const vector<MCParticle*>& mcps );

    private:
      double (*_Yfunc)(Jet& jet1, Jet& jet2, double Evis2);
      JetConfig _cfg;
  };
}

#endif
