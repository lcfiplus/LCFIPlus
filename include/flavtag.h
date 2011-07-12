#ifndef flavtag_h
#define flavtag_h 1

//#ifdef lib_flavtag_EXPORTS
#define NO_EVE 1
//#endif

#include <vector>
#include <map>
#include <iostream>
#include <sstream>
#include <exception>
#include "HelixClass.h"
#include "Rtypes.h"
#include "TObject.h"
#include "TLorentzVector.h"
#include "TVector3.h"

using namespace std;

namespace flavtag {

  // used for: track parameters, track covariance matrix, calorimeter subsystems
  namespace tpar {
    enum par { d0=0, z0, ph, om, td, parN };
    enum cov { d0d0=0, d0ph, phph, d0om, phom, omom, d0z0,
      z0ph, z0om, z0z0, d0td, phtd, omtd, z0td, tdtd, covN };
    enum hit { VTX=0, FTD, SIT, TPC, SET, ETD, hitN };
    enum calo { ecal=0, hcal, yoke, lcal, lhcal, bcal, caloN };
  }

  template <class T> struct DeleteVector { bool operator() (T x) const { delete x; return true; } };

  // used for error rescaling of the track parameters
  // in bins of cosTheta and momentum
  struct ErrorRescale {
    int var;    // 0: d0, 1: z0
    float cmin; // min cosTheta
    float cmax; // max cosTheta
    float pmin; // min momentum
    float pmax; // max momentum
    float pull; // pull - scale factor
    float pullerr; // pull error
  };

//   class LorentzVector {
//     public:
// 			LorentzVector() : _px(0.),_py(0.),_pz(0.),_en(0.) {};
//       LorentzVector(float px, float py, float pz, float en);
//       ~LorentzVector() {}
// 
//       TVector3 getTVector3() const;
//       TLorentzVector getTLorentzVector() const;
// 
//       float getPx() const { return _px; }
//       float getPy() const { return _py; }
//       float getPz() const { return _pz; }
//       float getEnergy() const { return _en; }
// 
//       float getPt() const;
//       float getMomentum() const;
// 
//       void add(const LorentzVector& d);
// 
//     protected:
//       float _px;
//       float _py;
//       float _pz;
//       float _en;
//   };

  class Event;
  class Track;
  class Neutral;
  class MCParticle;
  class EventData;

//   class EventPointer {
//     public:
//       void setEvent(Event* event) { _event = event; }
//       Event* getEvent() const { return _event; }
//     protected:
//       Event* _event;
//   };

  struct TrackData {
    int id;
    int mcid;
    int pdg;
    // track parameter
    float par[tpar::parN];
    // covariance matrix
    float cov[tpar::covN];
    int nhits[tpar::hitN];
    float charge;
    float en;
    float px;
    float py;
    float pz;
    float calo[tpar::caloN];
    float rimh; // radius of innermost hit
    float chi2; // from track fit
    int ndf;

		ClassDefNV(TrackData, 1);
  };

//   struct NeutralData {
//     int id;
//     int mcid;
//     int pdg;
//     float en;
//     float px;
//     float py;
//     float pz;
//     // energy deposit in each calorimeter subsystem
//     // ecal, hcal, yoke, lcal, lhcal, bcal
//     float calo[tpar::caloN];
// 		ClassDefNV(NeutralData, 1);
//   };

  class Event {
    public:
      // constructors
/*      Event(const EventData& data);
      */
			Event(const Event& event);
			Event(const char *nameTracks = "Tracks", const char *nameNeutrals = "Neutrals", const char *namePFOs = "MCParticles");
//			Event(const vector<TrackData> &tracks, const vector<NeutralData> &neutrals, const vector<MCParticleData> &mcps) : _rescaleTableRead(false)
//			{Init(tracks, neutrals, mcps);}
      ~Event();
			void Init(const vector<Track *> *ptracks, const vector<Neutral *> *pneutrals, const vector<MCParticle *> *pmcps);

      const vector<Track*>& getTracks() const { return _tracks; }
      const vector<Neutral*>& getNeutrals() const { return _neutrals; }
      const vector<MCParticle*>& getMCParticles() const { return _mcps; }
      MCParticle* getMCParticle(int id) const;
      MCParticle* getMCParticle(const Track* trk) const;

      vector<MCParticle*> mcGetColorStrings() const;
      int mcNumberOfB() const;
      int mcNumberOfC() const;
      void rescaleErrors();
      void setTracks( vector<Track*> tracks ) { _tracks = tracks; }
      void deleteTracks();

    private:
      // members
      int id;
      vector<Track*> _tracks;
      vector<Neutral*> _neutrals;
      vector<MCParticle*> _mcps;
      vector<ErrorRescale> _errRescaleTable;
      void readErrorRescaleTable();
      bool _rescaleTableRead;

  };

  class Track : public TLorentzVector {//, protected TrackData {//, public EventPointer {
    public:
      // constructor
			Track(){}
//      Track(const TrackData& d);
      ~Track() {}

			int getId() const { return _id; }
			void setId(int id) {_id = id;}

			flavtag::MCParticle * getMcp() const { return _mcp; }
			void setMcp(flavtag::MCParticle *mcp) {_mcp = mcp;}

			int getPDG() const { return _pdg; }
			void setPDG(int pdg) { _pdg = pdg;}

			float getCharge() const { return _charge; }
			void setCharge(float charge){_charge = charge;}

      float getD0()         const { return _par[tpar::d0]; }
      float getZ0()         const { return _par[tpar::z0]; }
      float getPhi()        const { return _par[tpar::ph]; }
      float getOmega()      const { return _par[tpar::om]; }
      float getTanLambda()  const { return _par[tpar::td]; }
			void setHelix(float d0, float z0, float phi, float omega, float tanLambda)
			{
				_par[tpar::d0] = d0; _par[tpar::z0] = z0; _par[tpar::ph] = phi;
				_par[tpar::om] = omega; _par[tpar::td] = tanLambda;
			}
			void setHelix(float * par){memcpy(_par, par, sizeof(_par));}

      const float* getCovMatrix() const { return _cov; }
      void setCovMatrix(float* cov);

      int getVtxHits() const { return _nhits[tpar::VTX]; }
      int getFtdHits() const { return _nhits[tpar::FTD]; }
      int getSitHits() const { return _nhits[tpar::SIT]; }
      int getTpcHits() const { return _nhits[tpar::TPC]; }
      int getSetHits() const { return _nhits[tpar::SET]; }
      int getEtdHits() const { return _nhits[tpar::ETD]; }
			void setTrackHits(float * hits){memcpy(_nhits, hits, sizeof(_nhits));}

      const float* getCaloEdep() const { return _calo; }
			void setCaloEdep(float *calo) { memcpy(_calo, calo, sizeof(_calo));}

      float getRadiusOfInnermostHit() const { return _rimh; }
      void setRadiusOfInnermostHit(float rimh) { _rimh = rimh; }

      float getChi2() const { return _chi2; }
			void setChi2(float chi2){_chi2 = chi2;}

      int getNdf() const { return _ndf; }
			void setNdf(int ndf) {_ndf = ndf;}

      void setFlightLength(float flt) const{ _flt = flt; }
      float getFlightLength() const { return _flt; }

      float getX() const;
      float getY() const;
      float getZ() const;
			TVector3 getPos()const {return TVector3(getX(), getY(), getZ());}

    private:
      mutable float _flt;

			int _id;
			flavtag::MCParticle * _mcp;
			int _pdg;
			float _charge;

			// track parameter
			float _par[tpar::parN];
			// covariance matrix
			float _cov[tpar::covN];
			int _nhits[tpar::hitN];
			float _calo[tpar::caloN];

			float _rimh; // radius of innermost hit
			float _chi2; // from track fit
			int _ndf;

			ClassDef(flavtag::Track, 2);
	};

  class Neutral : public TLorentzVector{ // : public LorentzVector, protected NeutralData{//, public EventPointer {
    public:
      // constructor
//      Neutral(const NeutralData& d);
      Neutral() : _id(0), _pdg(0), _isV0(false), _mcp(0) {}
      ~Neutral() {}

			int getId() const { return _id; }
			void setId(int id) {_id = id;}
			flavtag::MCParticle * getMcp() const { return _mcp; }
			void setMcp(flavtag::MCParticle *mcp) {_mcp = mcp;}
			int getPDG() const { return _pdg; }
			void setPDG(int pdg) { _pdg = pdg;}

			const float * getCaloEdep() const {return _calo;}
			//TODO: range check
			void setCaloEdep(float *calo) { memcpy(_calo, calo, sizeof(_calo));}

			bool isV0() const { return _isV0; }
			void setV0(bool isV0=true) { _isV0 = isV0; }

		private:
			int _id;
			int _pdg;
			int _isV0;
			flavtag::MCParticle * _mcp;

			float _calo[tpar::caloN];

			ClassDef(flavtag::Neutral, 3);
	};

	class MCParticle : public TLorentzVector{

		public:

			//ctor/dtors ////////////////////////////////////////////////////////
			MCParticle(int id, int pdg, MCParticle *parent, float charge, const TLorentzVector &p, const TVector3 &v){Init(id, pdg, parent, charge, p, v);}
			MCParticle(){}
			~MCParticle(){}

			// initialization
			void Init(int id, int pdg, MCParticle *parent, float charge, const TLorentzVector &p, const TVector3 &v);

			// simple accessors ///////////////////////////////////////////////////
			int getId() const { return _id; }
			void setId(int id) {_id = id;}

			int getPDG() const { return _pdg; }
			void setPDG(int pdg) {_pdg = pdg;}

			float getCharge() const { return _charge; }
			void setCharge(float charge) { _charge = charge;}

			const TVector3 & getVertex() const { return _vertex;}
			void setVertex(const TVector3 &v){_vertex = v;}
			const TVector3 & getEndVertex();

			MCParticle* getParent() const {return _parent;}
			void setParent(MCParticle *parent){_parent = parent;}
			//vector<MCParticle*> getParents() const;

			const vector<flavtag::MCParticle*>& getDaughters() const { return _dau; }
			void addDaughter(MCParticle* mcp);

			bool isParent(MCParticle *mcp)const;

		// more intelligent accessors
      int getFlavor() const;
      MCParticle* getColorString();

      bool isStableTrack() const;
			bool isStable() const;
      bool isSemiStableB() const;
      bool isSemiStableC() const;
      bool isSemiStableS() const;
      bool isSemiStable() const;
      MCParticle* semileptonicDecay() const; // return MCParticle of lepton
      int getFlavorTagCategory() const;
      MCParticle* getSemiStableParent() const;

      // end points
      float decayDistance();
      MCParticle* findDauForDecay();
      float getEx();
      float getEy();
      float getEz();

      vector<flavtag::MCParticle*> promptTracks();

      // helix stuff
      float getD0();
      float getZ0();
      float getPhi();
      float getOmega();
      float getTanLambda();

    private:
			// basic properties
			int _id;
			int _pdg;
			float _charge;
			MCParticle *_parent;
			TVector3 _vertex;
      // list of daughters
      vector<flavtag::MCParticle*> _dau;

      MCParticle* _dauForDecay; //! cached object

      // helix stuff
      void makeHelix();
      bool _helixOK;						//!
      HelixClass _helix;				//!

		ClassDef(flavtag::MCParticle, 2);
  };

  class Vertex {
    friend class LcfiInterface;

    public:
      enum vtx { xx=0, xy, yy, xz, yz, zz };

      Vertex() : _chi2(0), _prob(0), _x(0), _y(0), _z(0) {}
			Vertex(float chi2, float prob, float x, float y, float z, const float cov[6])
				: _id(-1), _chi2(chi2), _prob(prob), _x(x), _y(y), _z(z)
			{
				memcpy(_cov, cov, sizeof(_cov));
			}

      ~Vertex() {};
      void add(Track* trk);
      void add(Track* trk,float chi2);
			void setId(int id){ _id = id; }

			int getId() const { return _id;}
      float getChi2() const { return _chi2; }
      float getProb() const { return _prob; }
      float getX() const { return _x; }
      float getY() const { return _y; }
      float getZ() const { return _z; }
			TVector3 getPos() const {return TVector3(_x, _y,_z);}
      const float* getCov() const { return _cov; }
      const vector<flavtag::Track*> & getTracks() const { return _tracks; }
			const map<flavtag::Track *, float> & getTracksChi2Map() const {return _chi2Tracks;}
			float getChi2Track(Track *tr){
				map<Track*,float>::iterator it = _chi2Tracks.find(tr);
				if(it != _chi2Tracks.end())return it->second;
				else return -1;
			}
			Track * getWorstTrack() const;

      float length(const Vertex* primary=0) const;
      float significance(const Vertex* primary) const;

			float getPparallel(const TVector3 &axis)const;
			float getVertexMass(const Vertex *daughter = 0, const TVector3 *paxis = 0, const double dmass = 1.87, double *ppt = 0, double *pp = 0)const;
			float getVertexAngle(const Vertex *daughter, const Vertex *primary = 0)const;

			MCParticle *getMcp()const;

			double dirdot(const Vertex* primary=0) const;
			bool passesV0selection(const Vertex* primary=0) const;
			TLorentzVector getFourMomentum() const;

    private:
			int _id; // necessary for LCIO relation; -1 if not used

      float _chi2;
      float _prob;
      float _x;
      float _y;
      float _z;
      float _cov[6];
      vector<Track*> _tracks;
			map<Track*, float> _chi2Tracks;

			ClassDefNV(Vertex, 1);
  };

	class MCVertex {
		public:
			MCVertex() : _parent(0), _recoVtx(0) {}
			void setParent(MCParticle* parent) { _parent = parent; }
			MCParticle* getParent() const  { return _parent; }
			void add(MCParticle* mcp) { _dauVec.push_back(mcp); }
			void add(Track* trk) { _recoTrks.push_back(trk); }
			void setRecoVertex(Vertex* vtx) { _recoVtx = vtx; }
			const vector<MCParticle*>& getDaughters() const { return _dauVec; }
			const vector<Track*>& getRecoTracks() const { return _recoTrks; }
			Vertex* getRecoVertex() const { return _recoVtx; }
		private:
			MCParticle* _parent;
			vector<MCParticle*> _dauVec;

			vector<Track*> _recoTrks;
			Vertex* _recoVtx;

			ClassDefNV(MCVertex, 1);
	};

  class TrackPocaXY {

    public:
      TrackPocaXY(const Track* trk, const Vertex* vtx);
      TrackPocaXY(const Track* trk, const TVector3 &v);
      bool success() { return _success; }
      double getFlightLength() { return _flt; }
      double getPoca() { return _poca; }

			double operator()(double* flt) const;

    private:
			void minimize();
      bool _success;
      const Track* _trk;
      double _flt;
      double _poca;
      TVector3 _v;

  };

	class Exception : std::exception{
		public:
			Exception(const char *message){_message = message;}
			virtual ~Exception()throw(){}

			virtual const char * what() const throw(){return _message.c_str();}
			void Print(){cerr << "flavtag_Exception: " << _message << endl;}

		private:
			string _message;
	};


	class Jet;
}

#endif
