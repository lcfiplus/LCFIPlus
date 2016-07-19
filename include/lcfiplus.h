#ifndef lcfiplus_h
#define lcfiplus_h 1

#include <vector>
#include <map>
#include <iostream>
#include <sstream>
#include <exception>
#include "Rtypes.h"
#include "TObject.h"
#include "TLorentzVector.h"
#include "TVector3.h"

#include "EventStore.h"
#include "EVENT/Cluster.h"

#include <typeinfo>

using namespace std;

namespace lcfiplus {

// used for: track parameters, track covariance matrix, calorimeter subsystems
namespace tpar {
enum par { d0=0, z0, ph, om, td, parN };
enum cov { d0d0=0, d0ph, phph, d0om, phom, omom, d0z0,
           z0ph, z0om, z0z0, d0td, phtd, omtd, z0td, tdtd, covN
         };
enum hit { VTX=0, FTD, SIT, TPC, SET, ETD, hitN };
enum calo { ecal=0, hcal, yoke, lcal, lhcal, bcal, caloN };
}

template <class T> struct DeleteVector {
  bool operator() (T x) const {
    delete x;
    return true;
  }
};

template <class T> vector<const T*>* constVector(vector<T*>* ptr) {
  return reinterpret_cast<vector<const T*> *>(ptr);
}
template <class T> vector<const T*>& constVector(vector<T*>& ref) {
  return *reinterpret_cast<vector<const T*> *>(&ref);
}
template <class T> const vector<const T*>* constVector(const vector<T*>* ptr) {
  return reinterpret_cast<const vector<const T*> *>(ptr);
}
template <class T> const vector<const T*>& constVector(const vector<T*>& ref) {
  return *reinterpret_cast<const vector<const T*> *>(&ref);
}

// used for error rescaling of the track parameters
// in bins of cosTheta and momentum
struct ErrorRescale {
  int var;    // 0: d0, 1: z0
  double cmin; // min cosTheta
  double cmax; // max cosTheta
  double pmin; // min momentum
  double pmax; // max momentum
  double pull; // pull - scale factor
  double pullerr; // pull error
};

// basic types
class Event;
class Track;
class Neutral;
class MCParticle;
class MCColorSinglet;
class Vertex;
class Jet;

// const vectors
typedef const vector<const Track*> TrackVec;
typedef const vector<const Neutral*> NeutralVec;
typedef const vector<const MCParticle*> MCParticleVec;
typedef const vector<const MCColorSinglet*> MCColorSingletVec;
typedef const vector<const Vertex*> VertexVec;
typedef const vector<const Jet*> JetVec;

// const iterators
typedef vector<const Track*>::const_iterator TrackVecIte;
typedef vector<const Neutral*>::const_iterator NeutralVecIte;
typedef vector<const MCParticle*>::const_iterator MCParticleVecIte;
typedef vector<const MCColorSinglet*>::const_iterator MCColorSingletVecIte;
typedef vector<const Vertex*>::const_iterator VertexVecIte;
typedef vector<const Jet*>::const_iterator JetVecIte;

class Exception : std::exception {
 public:
  Exception(const char* message) {
    _message = message;
  }
  virtual ~Exception()throw() {}

  virtual const char* what() const throw() {
    return _message.c_str();
  }
  void Print() {
    cerr << "lcfiplus_Exception: " << _message << endl;
  }

 private:
  string _message;
};

// global parameters
class Globals {
 public:
  ~Globals();
  static Globals* Instance();
  void setBField(double bField) {
    _bField=bField;
  }
  double getBField()const {
    return _bField;
  }
  void setBeamSizeX(double beamSizeX) {
    _beamSizeX = beamSizeX;
  }
  double getBeamSizeX()const {
    return _beamSizeX;
  }
  void setBeamSizeY(double beamSizeY) {
    _beamSizeY = beamSizeY;
  }
  double getBeamSizeY()const {
    return _beamSizeY;
  }
  void setBeamSizeZ(double beamSizeZ) {
    _beamSizeZ = beamSizeZ;
  }
  double getBeamSizeZ()const {
    return _beamSizeZ;
  }

 private:
  // value of b field to be read from GEAR
  double _bField;
  double _beamSizeX;
  double _beamSizeY;
  double _beamSizeZ;
  static Globals* _theInstance;

  Globals();
};


// generic parameter class
class Parameters {
 public:
  Parameters(bool as = true) : _allowstring(as) {}

  // destructor for destroying objects in the map
  ~Parameters() {
    clear();
  }
  // copy operator/constructor
  Parameters(const Parameters& ref) {
    _allowstring = ref._allowstring;

    *this = ref;
  }
  Parameters& operator =(const Parameters& ref);  // in lcfiplus.cc

  // fetch for non-vector
  template<typename T> void fetch(const char* key, T& ret, T def = T())const {
    if (_map.find(key) == _map.end()) {
      cout << "Parameter " << key << " not found." << endl;
      ret = def;
    } else if (_map.find(key)->second.first == &typeid(T))ret = *(T*)(_map.find(key)->second.second);
    else if (_allowstring && _map.find(key)->second.first == &typeid(string)) {
      istringstream str(*(string*)_map.find(key)->second.second);
      str >> ret;
    } else if (_allowstring && _map.find(key)->second.first == &typeid(vector<string>)) {
      istringstream str((*(const vector<string>*)_map.find(key)->second.second)[0]);
      str >> ret;
    } else {
      cout << "Parameter type invalid: key = " << key << ", type = " << _map.find(key)->second.first->name() << " vs " << typeid(T).name() << endl;
      throw (Exception("Parameter type invalid."));
    }
    return;
  }

  // fetch for vector
  template<typename T> void fetchArray(const char* key, vector<T>& ret)const {
    if (_map.find(key) == _map.end()) {
      cout << "Parameter " << key << " not found." << endl;
    } else if (_map.find(key)->second.first == &typeid(vector<T>))ret = *(vector<T>*)(_map.find(key)->second.second);
    else if (_map.find(key)->second.first == &typeid(T))ret.push_back(*(T*)(_map.find(key)->second.second));
    else if (_allowstring && _map.find(key)->second.first == &typeid(string)) {
      ret.push_back(T());

      istringstream str(*(string*)_map.find(key)->second.second);
      str >> ret[0];
    } else if (_allowstring && _map.find(key)->second.first == &typeid(vector<string>)) {
      const vector<string>* svec = (const vector<string>*)_map.find(key)->second.second;
      ret.resize(svec->size());
      for (unsigned int n=0; n<svec->size(); n++) {
        istringstream str((*svec)[n]);
        str >> ret[n];
      }
    } else {
      cout << "Parameter type invalid: key = " << key << ", type = " << _map.find(key)->second.first->name() << " vs " << typeid(T).name() << endl;
      throw (Exception("Parameter type invalid."));
    }
  }

 public:
  // for non-vector only (if string parameter)
  template<typename T> T get(const char* key, T def = T())const {
    T ret;
    fetch(key, ret, def);
    return ret;
  }

  template<typename T> std::vector<T> getVec(const char* key)const {
    std::vector<T> ret;
    fetchArray(key, ret);
    return ret;
  }

  bool exist(const char* key)const {
    return _map.find(key) != _map.end();
  }

  template<typename T> void add(const char* key, T data) { // change T& data to T data for convenience
    if (_map.find(key) != _map.end())throw(Exception("Double entry."));

    _map[key] = pair<const type_info*, void*>(&typeid(T), new T(data));

  }
  void remove(const char* key) {
    remove(key, true);
  }
  void clear();

  template<typename T> void assign(const char* key, T data) {
    if (_map.find(key) == _map.end())throw(Exception("Parameters::assign(): the key has not been registered."));
    else if (_map.find(key)->second.first != &typeid(T))throw(Exception("Parameters::assign(): the value type is imcompatible."));

    // assign
    *static_cast<T*>(_map.find(key)->second.second) = data;
  }

  // direct accessors: avoiding copy to handle names
  const map<string, pair<const type_info*, void*> >& paramMap()const {
    return _map;
  }

 private:
  void remove(const string& key, bool delmap);

  template<typename T> void deleteData(const string& key) {
    T* p = static_cast<T*>(_map[key].second);
    delete p;
  }

  map<string, pair<const type_info*, void*> > _map;
  bool _allowstring;
};

class Algorithm {
 public:
  Algorithm() {
    _param = 0;
  }
  virtual ~Algorithm() {}

  virtual void init(Parameters* param) {
    _param = param;
  }
  virtual void process() = 0;
  virtual void end() {}

 protected:
  Parameters* GetParam()const {
    return _param;
  }
 private:
  Parameters* _param;

  ClassDef(Algorithm,1);
};

class Event : public EventStore {
 public:
  ~Event();

  // singleton
  static Event* Instance();

  // standard collections retrievers
  // for other collections: use EventStore::Get()
  const vector<const Track*>& getTracks(const char* trackname = 0) const;
  const vector<const Neutral*>& getNeutrals(const char* neutralname = 0) const;
  const vector<const MCParticle*>& getMCParticles(const char* mcpname = 0) const;
  const vector<const MCColorSinglet*>& getMCColorSinglets(const char* mcpname = 0) const;
  const Vertex* getPrimaryVertex(const char* privtxname = 0) const;
  const vector<const Vertex*>& getSecondaryVertices(const char* secvtxname = 0) const;
  const vector<const Jet*>& getJets(const char* jetname = 0) const;

  // standard collection name setter/getter
  void setDefaultTracks(const char* name) {
    _defTrackName = name;
  }
  void setDefaultNeutrals(const char* name) {
    _defNeutralName = name;
  }
  void setDefaultMCParticles(const char* name) {
    _defMcpName = name;
  }
  void setDefaultPrimaryVertex(const char* name) {
    _defPriVtxName = name;
  }
  void setDefaultSecondaryVertices(const char* name) {
    _defSecVtxName = name;
  }
  void setDefaultJets(const char* name) {
    _defJetName = name;
  }

  const char* getDefaultTracks()const {
    return _defTrackName.c_str();
  }
  const char* getDefaultNeutrals()const {
    return _defNeutralName.c_str();
  }
  const char* getDefaultMCParticles()const {
    return _defMcpName.c_str();
  }
  const char* getDefaultPrimaryVertex()const {
    return _defPriVtxName.c_str();
  }
  const char* getDefaultSecondaryVertices()const {
    return _defSecVtxName.c_str();
  }
  const char* getDefaultJets()const {
    return _defJetName.c_str();
  }

  // utility functions for MCParticles
  const MCParticle* getMCParticle(int id) const;
  const MCParticle* getMCParticle(const Track* trk) const;

  vector<const MCParticle*> mcGetColorStrings() const;
  int mcNumberOfB() const;
  int mcNumberOfC() const;
  vector<const MCParticle*> mcGetSemiStableBs() const;
  vector<const MCParticle*> mcGetSemiStableCs() const;
  vector<const MCParticle*> mcGetSemiStableBCs(bool separatebc) const;
  int mcFindParent(MCParticleVec& vec, const MCParticle* p) const;

//      void rescaleErrors();

 private:

  // singleton
  static Event* _theInstance;

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
  Track() {}
//      Track(const TrackData& d);
  ~Track() {}

  int getId() const {
    return _id;
  }
  void setId(int id) {
    _id = id;
  }

  const lcfiplus::MCParticle* getMcp() const {
    return _mcp;
  }
  void setMcp(const lcfiplus::MCParticle* mcp) {
    _mcp = mcp;
  }

  int getPDG() const {
    return _pdg;
  }
  void setPDG(int pdg) {
    _pdg = pdg;
  }

  double getCharge() const {
    return _charge;
  }
  void setCharge(double charge) {
    _charge = charge;
  }

  double getD0()         const {
    return _par[tpar::d0];
  }
  double getZ0()         const {
    return _par[tpar::z0];
  }
  double getPhi()        const {
    return _par[tpar::ph];
  }
  double getOmega()      const {
    return _par[tpar::om];
  }
  double getTanLambda()  const {
    return _par[tpar::td];
  }
  void setHelix(double d0, double z0, double phi, double omega, double tanLambda) {
    _par[tpar::d0] = d0;
    _par[tpar::z0] = z0;
    _par[tpar::ph] = phi;
    _par[tpar::om] = omega;
    _par[tpar::td] = tanLambda;
  }
  void setHelix(double* par) {
    memcpy(_par, par, sizeof(_par));
  }

  const double* getCovMatrix() const {
    return _cov;
  }
  void setCovMatrix(double* cov);

  int getVtxHits() const {
    return _nhits[tpar::VTX];
  }
  int getFtdHits() const {
    return _nhits[tpar::FTD];
  }
  int getSitHits() const {
    return _nhits[tpar::SIT];
  }
  int getTpcHits() const {
    return _nhits[tpar::TPC];
  }
  int getSetHits() const {
    return _nhits[tpar::SET];
  }
  int getEtdHits() const {
    return _nhits[tpar::ETD];
  }
  void setTrackHits(int* hits) {
    memcpy(_nhits, hits, sizeof(_nhits));
  }

  const double* getCaloEdep() const {
    return _calo;
  }
  void setCaloEdep(double* calo) {
    memcpy(_calo, calo, sizeof(_calo));
  }

  double getRadiusOfInnermostHit() const {
    return _rimh;
  }
  void setRadiusOfInnermostHit(double rimh) {
    _rimh = rimh;
  }

  double getChi2() const {
    return _chi2;
  }
  void setChi2(double chi2) {
    _chi2 = chi2;
  }

  int getNdf() const {
    return _ndf;
  }
  void setNdf(int ndf) {
    _ndf = ndf;
  }

  void setFlightLength(double flt) const {
    _flt = flt;
  }
  double getFlightLength() const {
    return _flt;
  }

  void setParticleIDProbability(string parName, double pidProbability) {
    _pidProbability[parName] = pidProbability;
  }

  double getParticleIDProbability(const char* parName) const {
    string partName = parName;
    double prob=0.0;
    if (_pidProbability.find(partName) != _pidProbability.end()){
      prob = _pidProbability.at(partName); // [] is not a const function
    }
    return prob;
  }

  //for track energy correction
  void setCorrEnergy(double mass) {
    _correnergy = sqrt(this->P()*this->P() + mass * mass);
  }
  
  double getCorrEnergy() const {
    return _correnergy;
  }

  void swapEnergy() {
    double tmpe = this->E();
    this->SetE(_correnergy);
    _correnergy = tmpe;
  }

  //for storing BNess
  void setBNess(double bness) {
    _bness = bness;
  }

  double getBNess() const {
    return _bness;
  }

  double getCNess() const {
    return _cness;
  }

  //for storing CNess
  void setCNess(double cness) {
    _cness = cness;
  }
  double getX() const;
  double getY() const;
  double getZ() const;
  TVector3 getPos()const {
    return TVector3(getX(), getY(), getZ());
  }
  TVector3 momentumAtVertex( const Vertex* vtx ) const;

 private:
  mutable double _flt;

  int _id;
  const lcfiplus::MCParticle* _mcp;
  int _pdg;
  double _charge;

  // track parameter
  double _par[tpar::parN];
  // covariance matrix
  double _cov[tpar::covN];
  int _nhits[tpar::hitN];
  double _calo[tpar::caloN];

  double _rimh; // radius of innermost hit
  double _chi2; // from track fit
  int _ndf;

  //ParticleID posterior probability
  map<string, double> _pidProbability;
  double _correnergy;

  //BNess
  double _bness, _cness;

  ClassDef(lcfiplus::Track, 2);
};

class Neutral : public TLorentzVector { // : public LorentzVector, protected NeutralData{//, public EventPointer {
 public:
  // constructor
//      Neutral(const NeutralData& d);
  Neutral() : _id(0), _pdg(0), _isV0(false), _mcp(0) {}
  ~Neutral() {}

  int getId() const {
    return _id;
  }
  void setId(int id) {
    _id = id;
  }
  const lcfiplus::MCParticle* getMcp() const {
    return _mcp;
  }
  void setMcp(const lcfiplus::MCParticle* mcp) {
    _mcp = mcp;
  }
  int getPDG() const {
    return _pdg;
  }
  void setPDG(int pdg) {
    _pdg = pdg;
  }

  const double* getCaloEdep() const {
    return _calo;
  }
  //TODO: range check
  void setCaloEdep(double* calo) {
    memcpy(_calo, calo, sizeof(_calo));
  }

  bool isV0() const {
    return _isV0;
  }
  void setV0(bool isV0=true) {
    _isV0 = isV0;
  }

  void setClusters(EVENT::ClusterVec clu) { _clstr = clu; }
  EVENT::ClusterVec getClusters() const { return _clstr; }
  
 private:
  int _id;
  int _pdg;
  int _isV0;
  const lcfiplus::MCParticle* _mcp;

  double _calo[tpar::caloN];

  EVENT::ClusterVec _clstr;

  ClassDef(lcfiplus::Neutral, 3);
};

class MCParticle : public TLorentzVector {

 public:

  //ctor/dtors ////////////////////////////////////////////////////////
  MCParticle(int id, int pdg, MCParticle* parent, double charge, const TLorentzVector& p, const TVector3& v) {
    Init(id, pdg, parent, charge, p, v);
  }
  MCParticle() {}
  ~MCParticle() {}

  // initialization: non-const parent is needed because of adding me to the daughter list
  void Init(int id, int pdg, MCParticle* parent, double charge, const TLorentzVector& p, const TVector3& v);

  // simple accessors ///////////////////////////////////////////////////
  int getId() const {
    return _id;
  }
  void setId(int id) {
    _id = id;
  }

  int getPDG() const {
    return _pdg;
  }
  void setPDG(int pdg) {
    _pdg = pdg;
  }

  double getCharge() const {
    return _charge;
  }
  void setCharge(double charge) {
    _charge = charge;
  }

  const TVector3& getVertex() const {
    return _vertex;
  }
  void setVertex(const TVector3& v) {
    _vertex = v;
  }
  const TVector3& getEndVertex()const;

  const MCParticle* getParent() const {
    return _parent;
  }
  void setParent(const MCParticle* parent) {
    _parent = parent;
  }
  //vector<MCParticle*> getParents() const;

  const vector<const MCParticle*>& getDaughters() const {
    return _dau;
  }
  void addDaughter(const MCParticle* mcp);

  bool isParent(const MCParticle* mcp)const;

  // more intelligent accessors
  int getFlavor() const;
  const MCParticle* getColorString()const;

  const MCColorSinglet* getColorSinglet(const vector<const MCColorSinglet*>* pcs)const;

  bool isStableTrack() const;
  bool isStable() const;
  bool isSemiStableB() const;
  bool isSemiStableC() const;
  bool isSemiStableS() const;
  bool isSemiStable() const;
  const MCParticle* semileptonicDecay() const; // return MCParticle of lepton
  int getFlavorTagCategory() const;
  const MCParticle* getSemiStableParent() const;
  const MCParticle* getSemiStableBParent() const;
  const MCParticle* getSemiStableCParent() const;

  // end points
  double decayDistance()const;
  const MCParticle* findDauForDecay()const;
  double getEx()const;
  double getEy()const;
  double getEz()const;

  vector<const lcfiplus::MCParticle*> promptTracks()const;

  // helix stuff
  double getD0()const;
  double getZ0()const;
  double getPhi()const;
  double getOmega()const;
  double getTanLambda()const;

 private:
  // basic properties
  int _id;
  int _pdg;
  double _charge;
  const MCParticle* _parent;
  TVector3 _vertex;
  // list of daughters
  vector<const lcfiplus::MCParticle*> _dau;

  mutable const MCParticle* _dauForDecay; //! cached object

  ClassDef(lcfiplus::MCParticle, 2);
};

class MCColorSinglet {
 public:
  const MCParticle* getMcp()const {
    return _mcp;
  }

  const MCParticle* _mcp;

  vector<const lcfiplus::MCParticle*> _initials;
  vector<const lcfiplus::MCParticle*> _finalstrings;
  vector<const lcfiplus::MCParticle*> _qqgluons;
  vector<const lcfiplus::MCParticle*> _finalcolorsinglets;
  vector<const lcfiplus::MCParticle*> _realparticles;

  ClassDefNV(lcfiplus::MCColorSinglet, 1);
};

class Vertex {
  friend class LcfiInterface;

 public:
  enum vtx { xx=0, xy, yy, xz, yz, zz };

  Vertex() : _chi2(0), _prob(0), _x(0), _y(0), _z(0), _isPrimary(false) {}
  Vertex(const double chi2, const double prob, const double x, const double y, const double z, const double cov[6], bool isPrimary)
    : _id(-1), _chi2(chi2), _prob(prob), _x(x), _y(y), _z(z), _isPrimary(isPrimary) {
    if (cov == 0) {
      memset(_cov,0, sizeof(_cov));
    } else {
      memcpy(_cov, cov, sizeof(_cov));
    }
  }

  // id is not copied
  Vertex(const Vertex& from) : _id(-1), _chi2(from._chi2), _prob(from._prob), _x(from._x), _y(from._y), _z(from._z), _isPrimary(from._isPrimary) {
    memcpy(_cov, from._cov, sizeof(_cov));
    _tracks = from._tracks;
    _chi2Tracks = from._chi2Tracks;
  }

  ~Vertex() {};
  void add(const Track* trk);
  void add(const Track* trk,double chi2);
  void setId(int id)const {
    _id = id;
  }

  int getId() const {
    return _id;
  }
  double getChi2() const {
    return _chi2;
  }
  double getProb() const {
    return _prob;
  }
  double getX() const {
    return _x;
  }
  double getY() const {
    return _y;
  }
  double getZ() const {
    return _z;
  }
  TVector3 getPos() const {
    return TVector3(_x, _y,_z);
  }
  const double* getCov() const {
    return _cov;
  }
  const vector<const Track*>& getTracks() const {
    return _tracks;
  }
  const map<const lcfiplus::Track*, double>& getTracksChi2Map() const {
    return _chi2Tracks;
  }
  double getChi2Track(const Track* tr)const {
    map<const Track*,double>::const_iterator it = _chi2Tracks.find(tr);
    if (it != _chi2Tracks.end())return it->second;
    else return -1;
  }
  const Track* getWorstTrack() const;
  double getChi2TrackFit(const Track* tr, int mode=1)const;

  double length(const Vertex* primary=0) const;
  double significance(const Vertex* primary) const;

  bool isPrimary()const {
    return _isPrimary;
  }
  void setPrimary(bool isPrimary) {
    _isPrimary = isPrimary;
  }

  double getPparallel(const TVector3& axis)const;
  double getVertexMass(const Vertex* daughter = 0, const TVector3* paxis = 0, const double dmass = 1.87, double* ppt = 0, double* pp = 0)const;
  double getVertexAngle(const Vertex* daughter, const Vertex* primary = 0)const;

  const MCParticle* getMcp()const;

  double dirdot(const Vertex* primary=0) const;
  bool passesV0selection(const Vertex* primary=0) const;
  TLorentzVector getFourMomentum() const;

  //for AVF
  void setVertexingName(string vtxname){ _vertexingname = vtxname;}
  string getVertexingName() const { return _vertexingname;}
  ////////////////

  //for vertex mass recovery
  double getRecoveredVertexMass() const {return _rvtxmass;}
  TLorentzVector getRecoveredFourMomentum() const {return _rvtxmom;}
  TLorentzVector getPi0sFourMomentum() const {return _pi0mom;}
  int getNPi0() const {return _npi0;}
  
  void setRecoveredVertexMass(double rvtxmass){_rvtxmass = rvtxmass;}
  void setRecoveredVertexMass() {_rvtxmass = _rvtxmom.M();}
  void setRecoveredFourMomentum(TLorentzVector rvtxmom) {_rvtxmom = rvtxmom;}
  void setPi0sFourMomentum(TLorentzVector pi0mom) {_pi0mom = pi0mom;}
  void setNPi0(int npi0) {_npi0 = npi0;}
  ////////////////

  void Print()const;
  static int dist_sort(const Vertex* a, const Vertex* b);
  bool covIsGood() const;

 private:
  mutable int _id; // necessary for LCIO relation; -1 if not used, mutable for modification in ConvertVertex()

  double _chi2;
  double _prob;
  double _x;
  double _y;
  double _z;
  double _cov[6];

  bool _isPrimary;

  vector<const lcfiplus::Track*> _tracks;
  map<const Track*, double> _chi2Tracks;


  //for AVF
  string _vertexingname;

  //for vertex mass recovery
  TLorentzVector _rvtxmom;
  TLorentzVector _pi0mom;
  double _rvtxmass;
  int _npi0;
  ////////////

  ClassDefNV(Vertex, 1);
};
class MCVertex {
 public:
  MCVertex() : _parent(0), _recoVtx(0) {}
  void setParent(const MCParticle* parent) {
    _parent = parent;
  }
  const MCParticle* getParent() const  {
    return _parent;
  }
  void add(const MCParticle* mcp) {
    _dauVec.push_back(mcp);
  }
  void add(const Track* trk) {
    _recoTrks.push_back(trk);
  }
  void setRecoVertex(const Vertex* vtx) {
    _recoVtx = vtx;
  }
  const vector<const MCParticle*>& getDaughters() const {
    return _dauVec;
  }
  const vector<const Track*>& getRecoTracks() const {
    return _recoTrks;
  }
  const Vertex* getRecoVertex() const {
    return _recoVtx;
  }

  const TVector3& getPos()const {
    return _pos;
  }
  void setPos(const TVector3& pos) {
    _pos = pos;
  }

 private:
  TVector3 _pos;
  const MCParticle* _parent;
  vector<const MCParticle*> _dauVec;

  vector<const Track*> _recoTrks;
  const Vertex* _recoVtx;

  ClassDefNV(MCVertex, 1);
};

class TrackPocaXY {

 public:
  TrackPocaXY(const Track* trk, const Vertex* vtx);
  TrackPocaXY(const Track* trk, const TVector3& v);
  bool success() {
    return _success;
  }
  double getFlightLength() {
    return _flt;
  }
  double getPoca() {
    return _poca;
  }

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
  Jet(const Track* trk);
  Jet(const Neutral* neutral);
  Jet(const Vertex* vtx) : _id(-1) {
    add(vtx);
  }
  Jet(const Jet& from, bool extractVertex = false);
  ~Jet() {};

  void setId(int id)const {
    _id = id;
  }
  int getId()const {
    return _id;
  }

  const vector<const Track*>& getTracks() const {
    return _tracks;
  }
  const vector<const Neutral*>& getNeutrals() const {
    return _neutrals;
  }
  const vector<const Vertex*>& getVertices() const {
    return _vertices;
  }

  /** returns list of vertices which are useful for flavor tagging. */
  vector<const Vertex*> getVerticesForFT() const;

  /** returns list of all tracks including those used to form vertices.
  	if withoutV0 is true, skip tracks which are identified as v0 (default: false)
   */
  vector<const Track*> getAllTracks(bool withoutV0=false) const;

  // methods
  void add(const Jet& jet);
  void add(const Track* trk) {
    _tracks.push_back(trk);
    *(TLorentzVector*)this += *(TLorentzVector*)trk;
  }
  void add(const Neutral* neut) {
    _neutrals.push_back(neut);
    *(TLorentzVector*)this += *(TLorentzVector*)neut;
  }
  void add(const Vertex* vtx, bool removeTracks = true) {
    _vertices.push_back(vtx);
    for (unsigned int i=0; i<vtx->getTracks().size(); i++) {
      *(TLorentzVector*)this += *(TLorentzVector*)(vtx->getTracks()[i]);

      if (removeTracks) {
        vector<const Track*>::iterator it;
        if ((it = find(_tracks.begin(), _tracks.end(), vtx->getTracks()[i])) != _tracks.end()) {
          *(TLorentzVector*)this -= *(TLorentzVector*)(*it);
          _tracks.erase(it);
        }
      }
    }
  }
  //void recalculate();

  static int sort_by_energy(const Jet* a, const Jet* b) {
    return (a->E() > b->E());
  }

  double sphericity() const;

  // parameter contrl
//			void addParam(const char *paramname, Parameters &param, bool forcereset = false){
  void addParam(const char* paramname, Parameters& param)const {
    if (_params.count(paramname) == 1)throw(Exception("Jet::addParam: parameter of the specified name has been already registered."));
    _params[paramname] = param;
  }

  const Parameters* getParam(const char* paramname)const {
    if (_params.count(paramname) == 0)throw(Exception("Jet::getParam: parameter of the specified name is not registered."));
    return &(_params.find(paramname)->second);
  }

  const map<string, Parameters>& params()const {
    return _params;
  }

  void recalcFourMomentum() {
    TLorentzVector v;
    vector<const Track*> tr = getAllTracks();
    for (unsigned int i=0; i<tr.size(); i++)
      v += (*tr[i]);
    for (unsigned int i=0; i<_neutrals.size(); i++)
      v += (*_neutrals[i]);

    *(TLorentzVector*)this = v;
  }

 private:
  vector<const Track*> _tracks;
  vector<const Neutral*> _neutrals;
  vector<const Vertex*> _vertices;

  mutable map<string, Parameters> _params;

  mutable int _id;

  ClassDefNV(Jet, 1);
};


}

#endif
