#ifndef lcfiplus_h
#define lcfiplus_h 1

//#ifdef lib_lcfiplus_EXPORTS
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

#include "EventStore.h"

using namespace std;

namespace lcfiplus {

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

  class Event;
  class Track;
  class Neutral;
  class MCParticle;
	class Vertex;
	class Jet;


	class Exception : std::exception{
		public:
			Exception(const char *message){_message = message;}
			virtual ~Exception()throw(){}

			virtual const char * what() const throw(){return _message.c_str();}
			void Print(){cerr << "lcfiplus_Exception: " << _message << endl;}

		private:
			string _message;
	};

	// generic parameter class
	class LcfiplusParameters
	{
	public:
		LcfiplusParameters(bool as = true) : _allowstring(as){}
		~LcfiplusParameters(){
			// TODO: delete map objects!
		}

		// fetch for non-vector
		template<typename T> void fetch(const char *key, T &ret, T def = T())const
		{
			if(_map.find(key) == _map.end()){
				cout << "Parameter " << key << " not found." << endl;
				ret = def;
			}
			else if(_map.find(key)->second.first == typeid(T).name())ret = *(T *)(_map.find(key)->second.second);
			else if(_allowstring && _map.find(key)->second.first == typeid(string).name()){
				istringstream str(*(string *)_map.find(key)->second.second);
				str >> ret;
			}
			else if(_allowstring && _map.find(key)->second.first == typeid(vector<string>).name()){
				istringstream str((*(const vector<string> *)_map.find(key)->second.second)[0]);
				str >> ret;
			}
			else
				throw(Exception("Parameter type invalid."));
			return;
		}

		// fetch for vector
		template<typename T> void fetchArray(const char *key, vector<T> &ret)const
		{
			if(_map.find(key) == _map.end()){
				cout << "Parameter " << key << " not found." << endl;
			}
			else if(_map.find(key)->second.first == typeid(vector<T>).name())ret = *(vector<T> *)(_map.find(key)->second.second);
			else if(_map.find(key)->second.first == typeid(T).name())ret.push_back(*(T *)(_map.find(key)->second.second));
			else if(_allowstring && _map.find(key)->second.first == typeid(string).name()){
				ret.push_back(T());

				istringstream str(*(string *)_map.find(key)->second.second);
				str >> ret[0];
			}
			else if(_allowstring && _map.find(key)->second.first == typeid(vector<string>).name()){
				const vector<string> *svec = (const vector<string> *)_map.find(key)->second.second;
				ret.resize(svec->size());
				for(unsigned int n=0;n<svec->size();n++){
					istringstream str((*svec)[n]);
					str >> ret[n];
				}
			}
			else
				throw(Exception("Parameter type invalid."));
		}

	public:
		// for non-vector only (if string parameter)
		template<typename T> T get(const char *key, T def = T())const{
			T ret;
			fetch(key, ret, def);
			return ret;
		}

		template<typename T> std::vector<T> getVec(const char *key)const {
			std::vector<T> ret;
			fetchArray(key, ret);
			return ret;
		}

		bool exist(const char *key)const{return _map.find(key) != _map.end();}

		template<typename T> void add(const char *key, T &data){
			if(_map.find(key) != _map.end())throw(Exception("Double entry."));

			_map[key] = pair<string, void *>(typeid(T).name(), new T(data));
		}

	private:
		map<string, pair<string, void *> > _map;
		bool _allowstring;
	}; 

	class LcfiplusAlgorithm
	{
	public:
		LcfiplusAlgorithm(){_param = 0;}
		virtual ~LcfiplusAlgorithm(){}

		virtual void init(LcfiplusParameters *param){
			_param = param;
		}
		virtual void process() = 0;
		virtual void end(){}

	protected:
		LcfiplusParameters * GetParam()const{return _param;}
	private:
		LcfiplusParameters *_param;

		ClassDef(LcfiplusAlgorithm,1);
	};

	class Event : public EventStore {
		public:
			~Event();

			// singleton
			static Event * Instance();

			// standard collections retrievers
			// for other collections: use EventStore::Get()
			//const vector<Track*>& getTracks(const char *trackname = 0) const;
			const vector<Track*>& getTracks(const char *trackname = "Tracks") const;
			//const vector<Neutral*>& getNeutrals(const char *neutralname = 0) const;
			const vector<Neutral*>& getNeutrals(const char *neutralname = "Neutrals") const;
			//const vector<MCParticle*>& getMCParticles(const char *mcpname = 0) const;
			const vector<MCParticle*>& getMCParticles(const char *mcpname = "MCParticles") const;
			//const Vertex* getPrimaryVertex(const char *privtxname = 0) const;
			const Vertex* getPrimaryVertex(const char *privtxname = "PrimaryVertex") const;
			const vector<Vertex*>& getSecondaryVertices(const char *secvtxname = 0) const;
			const vector<Jet*>& getJets(const char *jetname = 0) const;

			// standard collection name setter/getter
			void setDefaultTracks(const char *name){_defTrackName = name;}
			void setDefaultNeutrals(const char *name){_defNeutralName = name;}
			void setDefaultMCParticles(const char *name){_defMcpName = name;}
			void setDefaultPrimaryVertex(const char *name){_defPriVtxName = name;}
			void setDefaultSecondaryVertices(const char *name){_defSecVtxName = name;}
			void setDefaultJets(const char *name){_defJetName = name;}

			const char * getDefaultTracks()const{return _defTrackName.c_str();}
			const char * getDefaultNeutrals()const{return _defNeutralName.c_str();}
			const char * getDefaultMCParticles()const{return _defMcpName.c_str();}
			const char * getDefaultPrimaryVertex()const{return _defPriVtxName.c_str();}
			const char * getDefaultSecondaryVertices()const{return _defSecVtxName.c_str();}
			const char * getDefaultJets()const{return _defJetName.c_str();}

			// utility functions for MCParticles
      MCParticle* getMCParticle(int id) const;
      MCParticle* getMCParticle(const Track* trk) const;

      vector<MCParticle*> mcGetColorStrings() const;
      int mcNumberOfB() const;
      int mcNumberOfC() const;

//      void rescaleErrors();

    private:

			// singleton
			static Event * _theInstance;

			// private constructor
			Event();

      // members
//      int id;
//      vector<ErrorRescale> _errRescaleTable;
//      void readErrorRescaleTable();
//      bool _rescaleTableRead;

			string _defTrackName;
			string _defNeutralName;
			string _defMcpName;
			string _defPriVtxName;
			string _defSecVtxName;
			string _defJetName;

  };

  class Track : public TLorentzVector {//, protected TrackData {//, public EventPointer {
    public:
      // constructor
			Track(){}
//      Track(const TrackData& d);
      ~Track() {}

			int getId() const { return _id; }
			void setId(int id) {_id = id;}

			lcfiplus::MCParticle * getMcp() const { return _mcp; }
			void setMcp(lcfiplus::MCParticle *mcp) {_mcp = mcp;}

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
			lcfiplus::MCParticle * _mcp;
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

			ClassDef(lcfiplus::Track, 2);
	};

  class Neutral : public TLorentzVector{ // : public LorentzVector, protected NeutralData{//, public EventPointer {
    public:
      // constructor
//      Neutral(const NeutralData& d);
      Neutral() : _id(0), _pdg(0), _isV0(false), _mcp(0) {}
      ~Neutral() {}

			int getId() const { return _id; }
			void setId(int id) {_id = id;}
			lcfiplus::MCParticle * getMcp() const { return _mcp; }
			void setMcp(lcfiplus::MCParticle *mcp) {_mcp = mcp;}
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
			lcfiplus::MCParticle * _mcp;

			float _calo[tpar::caloN];

			ClassDef(lcfiplus::Neutral, 3);
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

			const vector<lcfiplus::MCParticle*>& getDaughters() const { return _dau; }
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

      vector<lcfiplus::MCParticle*> promptTracks();

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
      vector<lcfiplus::MCParticle*> _dau;

      MCParticle* _dauForDecay; //! cached object

      // helix stuff
      void makeHelix();
      bool _helixOK;						//!
      HelixClass _helix;				//!

		ClassDef(lcfiplus::MCParticle, 2);
  };

  class Vertex {
    friend class LcfiInterface;

    public:
      enum vtx { xx=0, xy, yy, xz, yz, zz };

      Vertex() : _chi2(0), _prob(0), _x(0), _y(0), _z(0) {}
			Vertex(const float chi2, const float prob, const float x, const float y, const float z, const float cov[6])
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
      const vector<lcfiplus::Track*> & getTracks() const { return _tracks; }
			const map<lcfiplus::Track *, float> & getTracksChi2Map() const {return _chi2Tracks;}
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

  class Jet : public TLorentzVector {
    public:
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

			// parameter contrl
			void addParam(const char *paramname, LcfiplusParameters &param, bool forcereset = false){
				if(!forcereset && _params.count(paramname) == 1)throw(Exception("Jet::addParam: parameter of the specified name has been already registered."));
				_params[paramname] = param;
			}

			const LcfiplusParameters * getParam(const char *paramname)const{
				if(_params.count(paramname) == 0)throw(Exception("Jet::getParam: parameter of the specified name is not registered."));
				return &(_params.find(paramname)->second);
			}

    private:
      vector<Track*> _tracks;
      vector<Neutral*> _neutrals;
			vector<Vertex*> _vertices;

			map<string, LcfiplusParameters> _params;

			int _id;

			ClassDefNV(Jet, 1);
  };


}

#endif
