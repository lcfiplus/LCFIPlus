#include "lcfiplus.h"
#include "EventStore.h"

#include "TVector3.h"
#include "TLorentzVector.h"
#include "TMatrixFSym.h"
#include <TMatrixTSym.h>
#include "TMatrixDSym.h"
#include "HelixClass.h"
#include "TVectorD.h"
#include "TMethodCall.h"
#include <HelixClass.h>

#include "math.h"

#include "algo.h"
#include "geometry.h"

#include "VertexFinderSuehara.h"
#include "VertexFitterSimple.h"

namespace lcfiplus {

// singleton operations
Globals* Globals::_theInstance = NULL;
Globals* Globals::Instance() {
  if (_theInstance == NULL) _theInstance = new Globals;
  return _theInstance;
}
Globals::Globals() {}
Globals::~Globals() {
  _theInstance = 0;
}

// singleton operations
Event* Event::_theInstance = NULL;
Event* Event::Instance() {
  if (_theInstance == NULL) _theInstance = new Event;
  return _theInstance;
}

Event::Event() {}
Event::~Event() {
  _theInstance = 0;
}

TrackVec& Event::getTracks(const char* trackname) const {
  if (!trackname) trackname = _defTrackName.c_str();
  const char* classname = "vector<lcfiplus::Track*>";

//		if(!IsExist(trackname,classname))throw(Exception("Event::getTracks(): collection not found."));

  return *(TrackVec*&)GetObjectRef(trackname, classname);
}

NeutralVec& Event::getNeutrals(const char* neutralname) const {
  if (!neutralname) neutralname = _defNeutralName.c_str();
  const char* classname = "vector<lcfiplus::Neutral*>";

//		if(!IsExist(neutralname,classname))throw(Exception("Event::getNeutrals(): collection not found."));

  return *(NeutralVec*&)GetObjectRef(neutralname, classname);
}

MCParticleVec& Event::getMCParticles(const char* mcpname) const {
  if (!mcpname) mcpname = _defMcpName.c_str();
  const char* classname = "vector<lcfiplus::MCParticle*>";

//		if(!IsExist(mcpname,classname))throw(Exception("Event::getMCParticles(): collection not found."));

  return *(MCParticleVec*&)GetObjectRef(mcpname, classname);
}

MCColorSingletVec& Event::getMCColorSinglets(const char* mcpname) const {
  if (!mcpname) mcpname = _defMcpName.c_str();
  const char* classname = "vector<lcfiplus::MCColorSinglet*>";

//		if(!IsExist(mcpname,classname))throw(Exception("Event::getMCParticles(): collection not found."));

  return *(MCColorSingletVec*&)GetObjectRef(mcpname, classname);
}

const Vertex* Event::getPrimaryVertex(const char* privtxname) const {
  if (!privtxname) privtxname = _defPriVtxName.c_str();
  const char* classname = "vector<lcfiplus::Vertex*>";

//		if(!IsExist(privtxname,classname))throw(Exception("Event::getPrimaryVertex(): collection not found."));
  auto const& primaryVertices = (*(vector<const lcfiplus::Vertex*>*)GetObject(privtxname,classname));
  if(primaryVertices.size() > 0 ){
    return primaryVertices[0];
  } else {
    return nullptr;
  }
}

const vector<const Vertex*>& Event::getSecondaryVertices(const char* secvtxname) const {
  if (!secvtxname) secvtxname = _defSecVtxName.c_str();
  const char* classname = "vector<lcfiplus::Vertex*>";

//		if(!IsExist(secvtxname,classname))throw(Exception("Event::getSecondaryVertices(): collection not found."));

  return *(vector<const Vertex*>*&)GetObjectRef(secvtxname,classname);
}

const vector<const Jet*>& Event::getJets(const char* jetname) const {
  if (!jetname) jetname = _defJetName.c_str();
  const char* classname = "vector<lcfiplus::Jet*>";

//		if(!IsExist(jetname,classname))throw(Exception("Event::getJets(): collection not found."));

  return *(vector<const Jet*>*&)GetObjectRef(jetname,classname);
}


// MCParticle utilities
int Event::mcNumberOfB() const {
  int num(0);
  const vector<const MCParticle*>& mcps = getMCParticles();
  for (MCParticleVecIte it = mcps.begin(); it != mcps.end(); ++it) {
    const MCParticle* mcp = *it;
    if (mcp->isSemiStableB()) ++num;
  }
  return num;
}

int Event::mcNumberOfC() const {
  int num(0);
  MCParticleVec& mcps = getMCParticles();
  for (MCParticleVecIte it = mcps.begin(); it != mcps.end(); ++it) {
    const MCParticle* mcp = *it;
    if (mcp->isSemiStableC()) ++num;
  }
  return num;
}

const MCParticle* Event::getMCParticle(int id) const {
  if (id == 0) return 0;
  MCParticleVec& mcps = getMCParticles();
  for (MCParticleVecIte it = mcps.begin(); it != mcps.end(); ++it) {
    const MCParticle* mcp = *it;
    if (mcp->getId() == id) return mcp;
  }
  return 0;
}

const MCParticle* Event::getMCParticle(const Track* trk) const {
  return trk->getMcp();
}

vector<const MCParticle*> Event::mcGetColorStrings() const {
  vector<const MCParticle*> colorStrings;
  MCParticleVec& mcps = getMCParticles();
  for (MCParticleVecIte it = mcps.begin(); it != mcps.end(); ++it) {
    const MCParticle* cs = (*it)->getColorString();
    if (find( colorStrings.begin(), colorStrings.end(), cs ) == colorStrings.end() ) {
      colorStrings.push_back(cs);
    }
  }
  return colorStrings;
}

vector<const MCParticle*> Event::mcGetSemiStableBs() const {
  vector<const MCParticle*> ret;
  MCParticleVec& mcps = getMCParticles();
  for (MCParticleVecIte it = mcps.begin(); it != mcps.end(); ++it) {
    const MCParticle* ppart = (*it)->getSemiStableBParent();
    if (ppart && find(ret.begin(), ret.end(), ppart) == ret.end())ret.push_back(ppart);
  }

  return ret;
}

vector<const MCParticle*> Event::mcGetSemiStableCs() const {
  vector<const MCParticle*> ret;
  MCParticleVec& mcps = getMCParticles();
  for (MCParticleVecIte it = mcps.begin(); it != mcps.end(); ++it) {
    const MCParticle* ppart = (*it)->getSemiStableCParent();
    if (ppart && find(ret.begin(), ret.end(), ppart) == ret.end())ret.push_back(ppart);
  }

  return ret;
}

vector<const MCParticle*> Event::mcGetSemiStableBCs(bool separatebc) const {
  vector<const MCParticle*> ret;
  MCParticleVec& mcps = getMCParticles();
  for (MCParticleVecIte it = mcps.begin(); it != mcps.end(); ++it) {
    const MCParticle* ppartb = (*it)->getSemiStableBParent();
    const MCParticle* ppartc = (*it)->getSemiStableCParent();
    if (ppartb) {
      if (find(ret.begin(), ret.end(), ppartb) == ret.end())ret.push_back(ppartb);
    }
    if (ppartc && (separatebc || ppartb == 0)) {
      if (find(ret.begin(), ret.end(), ppartc) == ret.end())ret.push_back(ppartc);
    }
  }

  return ret;
}

int Event::mcFindParent(MCParticleVec& vec, const MCParticle* p) const {
  for (unsigned int i=0; i<vec.size(); i++) {
    if (p->isParent(vec[i]))return i;
  }
  return -1;
}

/*
  void Event::readErrorRescaleTable() {
    FILE* file = fopen("table.dat","r");
    if (file == 0) {
      printf("table.dat not found, error not rescaled\n");
      return;
    }
    ErrorRescale a;
    for (;;) {
      fscanf(file, "%d %f %f %f %f %f %f", &a.var, &a.cmin, &a.cmax, &a.pmin, &a.pmax, &a.pull, &a.pullerr);
      if (feof(file)) break;
      printf("%d %f %f %f %f %f %f\n", a.var, a.cmin, a.cmax, a.pmin, a.pmax, a.pull,a.pullerr);
      _errRescaleTable.push_back(a);
    }
    fclose(file);
    _rescaleTableRead = true;
  }


  // rescale track error
  void Event::rescaleErrors() {
    if (_rescaleTableRead == false) readErrorRescaleTable();

    printf("rescaleErrors() called\n");

    for (TrackVecIte iter = _tracks.begin(); iter != _tracks.end(); ++iter) {
      const Track* trk = *iter;
      double c = fabs(sin(atan(trk->getTanLambda())));
      double p = trk->P();
      double d0pull(1);
      double z0pull(1);
      for (vector<ErrorRescale>::const_iterator tabIter = _errRescaleTable.begin(); tabIter != _errRescaleTable.end(); ++tabIter) {
        const ErrorRescale& resc = *tabIter;
        printf(" testing cmin:cmax:pmin:pmax=%.2f:%.2f:%.2f:%.2f... ",
            resc.cmin, resc.cmax, resc.pmin, resc.pmax);
        if ( (resc.var == 0) && (resc.cmin <= c) && (c < resc.cmax) && (resc.pmin <= p) && (p < resc.pmax) ) {
          d0pull = resc.pull;
          printf("d0pull found %f\n",d0pull);
        } else {
          printf(" not found.\n");
        }
        if ( (resc.var == 1) && (resc.cmin <= c) && (c < resc.cmax) && (resc.pmin <= p) && (p < resc.pmax) ) {
          z0pull = resc.pull;
          printf("z0pull found %f\n",z0pull);
        }
      }

      const double* cov = trk->getCovMatrix();
      double newCov[tpar::covN];
      for (int i=0; i<tpar::covN; ++i) {
        newCov[i] = cov[i];
      }

      //d0d0=0, d0ph, phph, d0om, phom, omom, d0z0,
              //z0ph, z0om, z0z0, d0td, phtd, omtd, z0td, tdtd, covN

      newCov[tpar::d0d0] *= d0pull*d0pull;
      newCov[tpar::z0z0] *= z0pull*z0pull;
      newCov[tpar::d0z0] *= d0pull*z0pull;

      newCov[tpar::d0ph] *= d0pull;
      newCov[tpar::d0om] *= d0pull;
      newCov[tpar::d0td] *= d0pull;

      newCov[tpar::z0ph] *= z0pull;
      newCov[tpar::z0om] *= z0pull;
      newCov[tpar::z0td] *= z0pull;

      trk->setCovMatrix(newCov);
    }
  }
*/

//   Track::Track(const TrackData& data) :
//     LorentzVector(data.px, data.py, data.pz, data.en),
//     TrackData(data),
//     _flt(0) {
//   }

double Track::getX() const {
  return  ( -_par[tpar::d0] + 1./_par[tpar::om] )*sin( _par[tpar::ph] )
          -  ( sin( _par[tpar::ph] - _par[tpar::om]*_flt ) )/_par[tpar::om];
}

double Track::getY() const {
  return  (  _par[tpar::d0] - 1./_par[tpar::om] )*cos( _par[tpar::ph] )
          +  ( cos( _par[tpar::ph] - _par[tpar::om]*_flt ) )/_par[tpar::om];
}

double Track::getZ() const {
  return _par[tpar::z0] + _par[tpar::om]*_flt;
}

void Track::setCovMatrix(double* mycov) {
  for (int i=0; i<tpar::covN; ++i) {
    _cov[i] = mycov[i];
  }
}

//   Neutral::Neutral(const NeutralData& data) :
//     LorentzVector(data.px, data.py, data.pz, data.en), NeutralData(data) {
//   }

TVector3 Track::momentumAtVertex( const Vertex* vtx ) const {
  Helix hel(this);
  double tmin;
  hel.LogLikelihood( vtx->getPos(), tmin );
  TVector3 ret = hel.GetPosDerivT( tmin );
  ret.SetMag( Vect().Mag() );
  return ret;
}


void MCParticle::Init(int id,int pdg, MCParticle* parent, double charge, const TLorentzVector& p, const TVector3& v) {
  *(TLorentzVector*)this = p;
  _dauForDecay = 0;

  _id = id;
  _pdg = pdg;
  _parent = parent;
  _vertex = v;
  _charge = charge;

  // add myself to the parent particle
  if (parent)
    parent->addDaughter(this);
}

int MCParticle::getFlavor() const {
  int mypdg = abs(getPDG());
  if (mypdg < 300) return 0;
  if (mypdg > 100000) return 0;
  if (mypdg > 10000) mypdg %= 1000;
  if (mypdg > 1000) mypdg /= 10;
  if (mypdg > 500) return 5;
  if (mypdg > 400) return 4;
  return 0;
}

const MCParticle* MCParticle::getColorString() const {
  if (getPDG() == 92) return this;
  const MCParticle* parent = getParent();
  if (parent == 0) return parent;
  return parent->getColorString();
}

const MCColorSinglet* MCParticle::getColorSinglet(const vector<const MCColorSinglet*>* pcs)const {
  for (unsigned int n=0; n<pcs->size(); n++) {
    const MCParticle* mcp = (*pcs)[n]->getMcp();
    if (mcp == this)return (*pcs)[n];
  }
  if (getParent() ==0) return 0;

  // recursive
  return getParent()->getColorSinglet(pcs);
}

// checks if the MCParticle underwent semileptonic decay.
// only checks the immediate daughters.
// if it is found to have decayed semileptonically (e, mu, tau),
// the pointer to the lepton is returned; zero otherwise.
const MCParticle* MCParticle::semileptonicDecay() const {
  /*
      MCParticleVec & mcps = _event->getMCParticles();
      for (MCParticleVecIte it = mcps.begin(); it != mcps.end(); ++it) {
        const MCParticle* test = *it;
        int abspdg = abs(test->getPDG());
        if (abspdg == 11 || abspdg == 13 || abspdg == 15) {
          if (test->getParent() == this)
            return test;
        }
      }
      return 0;*/
  for (unsigned int i=0; i<_dau.size(); i++) {
    const MCParticle* test = _dau[i];
    int abspdg = abs(test->getPDG());
    if (abspdg == 11 || abspdg == 13 || abspdg == 15) {
      return test;
    }
  }

  return 0;
}

bool MCParticle::isStableTrack() const {
  switch (abs(_pdg)) {
  case 11:
  case 13:
  case 211:
  case 321:
  case 2212:
    return true;
  default:
    return false;
  }
}

bool MCParticle::isStable() const {
  switch (abs(_pdg)) {
  case 12: // nu_e
  case 14: // nu_mu
  case 16: // nu_tau
  case 22: // photon
  //case 111: // pi0
  //case 130: // K_L
  case 2112: // neutron
  case 1000022: // neutralino
    return true;
  default:
    return isStableTrack();
  }
}


bool MCParticle::isSemiStableB() const {
  if (getFlavor() != 5) return false;
  switch (abs(_pdg)) {
  case 511: // B0 - ctau = 458.7 um
  case 521: // B+ - ctau = 491.1 um
  case 531: // B_s0 - ctau = 441 um
  case 541: // B_c0 - ctau = 138 um
  case 5122: // Lambda_b0 - ctau = 415 um
  case 5332: // Omega_b- - ctau = unknown
  case 5132: // Xi_b- - ctau = 426 um
  case 5232: // Xi_b0 - ctau = 426 um
    return true;

  // unstable
  case 513: // B*+
  case 523: // B*0
  case 533: // B_s*0
  case 543: // B_c*0
  case 5112: // Sigma_b-
  case 5212: // Sigma_b0
  case 5222: // Sigma_b+
  case 5114: // Sigma_b*-
  case 5214: // Sigma_b*0
  case 5224: // Sigma_b*+
  case 5314: // Xi_b*-
  case 5324: // Xi_b*0
  case 5334: // Omega_b*-
  case 5312: // Sigma_b'-
  case 5322: // Sigma_b'0
  case 551: // eta_b
  case 553: // Upsilon(1S)
  case 515: // B_2^*0
  case 525: // B_2^*+
  case 535: // B_s2^*0
  case 10441: //
  case 10443: //
  case 10511: //
  case 10513: // B_1(L)^0
  case 10521: //
  case 10523: //
  case 10531: //
  case 10533: //
  case 20513: // B_1(H)^0
  case 20523: //
  case 20533: //
    return false;
  default:
    printf("unknown particle: %d\n",_pdg);
    return false;
  }
}

bool MCParticle::isSemiStableC() const {
  if (getFlavor() != 4) return false;
  switch (abs(_pdg)) {
  case 411: // D+
  case 421: // D0
  case 431: // D_s+
  case 4122: // Lambda_c+ - ctau = 59.9 um
  case 4132: // Xi_c0 - ctau = 33.6 um
  case 4232: // Xi_c+ - ctau = 132 um
  case 4332: // Omega_c0 - ctau = 21 um
    return true;
  case 413: // D*+
  case 415: // D_2*(2460)+
  case 423: // D*0
  case 425: // D_2*(2460)0
  case 433: // D_s*+
  case 435: // D_s**+
  case 441: // eta_c(1S)
  case 443: // J/psi(1S)
  case 4112: // Sigma_c0
  case 4212: // Sigma_c+
  case 4222: // Sigma_c++
  case 4322: // Sigma_c'+
  case 4312: // Sigma_c'0
  case 4114: // Sigma_c*0
  case 4214: // Sigma_c*+
  case 4224: // Sigma_c*++
  case 4314: // Xi_c*0
  case 4324: // Xi_c*+
  case 4334: // Omega_c*0
  case 10411: // D_0*(2400)+
  case 10413: // D_1(2420)+
  case 10421: // D_0*(2400)0
  case 10423: // wtf?
  case 10431: // D_s0*(2317)+
  case 10433: // D_s1(2536)+
  case 20413: // D_1(H)+
  case 20423: // D_1(2430)0
  case 20433: // D_s1(2460)+
  case 20443: // wtf is this?
    return false;
  default:
    printf("unknown particle: %d\n",_pdg);
    return false;
  }
}

bool MCParticle::isSemiStableS() const {
  if (getFlavor() > 3) return false;
  switch (abs(_pdg)) {
  //case 311: // K0
  case 130: // K_L
  case 310: // K_S
  case 321: // K+
  case 3122: // Lambda
  case 3222: // Sigma+
  case 3212: // Sigma0
  case 3112: // Sigma-
  case 3322: // Xi0
  case 3312: // Xi-
  case 3334: // Omega-
    return true;
  default:
    return false;
  }
}

bool MCParticle::isSemiStable() const {
  switch (abs(_pdg)) {
  // bottom
  case 511: // B0 - ctau = 458.7 um
  case 521: // B+ - ctau = 491.1 um
  case 531: // B_s0 - ctau = 441 um
  case 541: // B_c0 - ctau = 138 um
  case 5122: // Lambda_b0 - ctau = 415 um
  case 5132: // Xi_b- - ctau = 426 um
  case 5232: // Xi_b0 - ctau = 426 um
  case 5332: // Omega_b- - ctau = unknown

  // charm
  case 411: // D+
  case 421: // D0
  case 431: // D_s+
  case 4122: // Lambda_c+ - ctau = 59.9 um
  case 4132: // Xi_c0 - ctau = 33.6 um
  case 4232: // Xi_c+ - ctau = 132 um
  case 4332: // Omega_c0 - ctau = 21 um

  // strange
  //case 311: // K0
  case 130: // K_L
  case 310: // K_S
  case 321: // K+
  case 3122: // Lambda
  case 3222: // Sigma+
  case 3212: // Sigma0
  case 3112: // Sigma-
  case 3322: // Xi0
  case 3312: // Xi-
  case 3334: // Omega-

    return true;

  default:
    return false;
  }
}

const MCParticle* MCParticle::getSemiStableParent() const {
  const MCParticle* parent = getParent();
  if (parent == 0) return 0;
  /*
  printf(" current id: %d (pdg=%d), parent id: %d (pdg=%d)\n",
      id,getPDG(),parent->getId(),parent->getPDG() );
      */
  if (parent->isSemiStable()) return parent;
  if (parent->isStable()) return parent;
  return parent->getSemiStableParent();
}

const MCParticle* MCParticle::getSemiStableBParent() const {
  const MCParticle* parent = getParent();
  if (parent == 0) return 0;

  if (parent->isSemiStableB()) return parent;
  return parent->getSemiStableBParent();
}

const MCParticle* MCParticle::getSemiStableCParent() const {
  const MCParticle* parent = getParent();
  if (parent == 0) return 0;

  if (parent->isSemiStableC()) return parent;
  return parent->getSemiStableCParent();
}

// classify the MCParticle according to its position in the MC tree
//   1 = MCParticle originates from primary
//   2 = MCParticle is a result of a semi-stable b decay (e.g. B -> D l nu)
//   3 = MCParticle is a result of a semi-stable c decay (e.g. D -> K pi)
//   4 = MCParticle is a result of a semi-stable s decay (e.g. Ks -> pi+ pi-)
// iterates over parents to find the most immediate semi-stable particle
int MCParticle::getFlavorTagCategory() const {
  int tag(0);
  const MCParticle* ssp = getSemiStableParent();

  if (ssp == 0) {
    tag = 1;
  } else {
    switch (ssp->getFlavor()) {
    case 5:
      tag = 2; // bottom
      break;
    case 4:
      tag = 3; // charm;
      break;
    default:
      tag = 4;
      break;
    }
  }

  // fix tag by decay distance...
  if (tag > 1) {
    if (getVertex().Mag() < 1e-5) {
      tag = 1;
    }
  }

  return tag;
}

const TVector3& MCParticle::getEndVertex() const {
  findDauForDecay();
  return _dauForDecay->getVertex();
}

// find the decay distance by iterating over daughters
// and grand-daughters (and possibly more generations)
// to find the shortest distance where a daughter tracks splits off
double MCParticle::decayDistance() const {
  findDauForDecay();
  return ( (_dauForDecay->getVertex()-getVertex()).Mag());
}

const MCParticle* MCParticle::findDauForDecay() const {
  if (_dauForDecay) return _dauForDecay;

  for (unsigned int i=0; i<_dau.size(); ++i) {
    const MCParticle* dau = _dau[i];
    //printf("  checking dau: %d\n",dau->getPDG());
    double test = ( dau->getVertex()-getVertex() ).Mag();
    if (test > 0 && (dau->isStable() || dau->isSemiStable())) {
      //printf("  dau %d found \n",dau->getPDG());
      _dauForDecay = dau;
      return dau;
    }
  }

  // look at grand daughters if the immediate daughter is unstable
  // this will actually iterate until we find the first stable particle
  for (unsigned int i=0; i<_dau.size(); ++i) {
    const MCParticle* dau = _dau[i];
    //printf("  calling findDauForDecay for dau %d\n", dau->getPDG());
    const MCParticle* dautmp = dau->findDauForDecay();
    if (dautmp != dau) {
      _dauForDecay = dautmp;
      return _dauForDecay;
    }
  }

  _dauForDecay = this;
  return this;
}

double MCParticle::getEx() const {
  if (_dauForDecay == 0) findDauForDecay();
  return _dauForDecay->getVertex().x();
}

double MCParticle::getEy() const {
  if (_dauForDecay == 0) findDauForDecay();
  return _dauForDecay->getVertex().y();
}

double MCParticle::getEz() const {
  if (_dauForDecay == 0) findDauForDecay();
  return _dauForDecay->getVertex().z();
}

bool MCParticle::isParent(const MCParticle* mcp)const {
  if (getParent() == 0)return false;
  if (getParent() == mcp)return true;

  return getParent()->isParent(mcp);
}

// find all the tracks which originate promptly,
// including tracks from unstable daughters
vector<const MCParticle*> MCParticle::promptTracks() const {
  vector<const MCParticle*> ret;
  //printf("promptTracks called for particle: %d\n",getPDG());
  //printf("  has daughters: ");
  for (MCParticleVecIte it = _dau.begin(); it != _dau.end(); ++it) {
    //MCParticle* dau = *it;
    //printf("%d ",dau->getPDG());
  }
  //printf("\n");
  for (MCParticleVecIte it = _dau.begin(); it != _dau.end(); ++it) {
    const MCParticle* dau = *it;
    //printf("  - checking daughter: %d\n",dau->getPDG());
    if (dau->getCharge() != 0) {
      if ( dau->isStableTrack() ) {
        //printf("    [%d] is stable; added\n",dau->getPDG());
        ret.push_back(dau);
      } else if ( dau->isSemiStable() ) {
        //printf("    [%d] is semistable; ignored\n",dau->getPDG());
      } else {
        //printf("    [%d] not stable, looking at daughter\n",dau->getPDG());
        vector<const MCParticle*> more = dau->promptTracks();
        ret.insert(ret.end(), more.begin(), more.end());
      }
    } else {
      // treat neutral particles
      if ( dau->decayDistance() == decayDistance() ) {
        //printf("    [%d] is neutral and decays promptly, looking at daughter\n",dau->getPDG());
        vector<const MCParticle*> more = dau->promptTracks();
        ret.insert(ret.end(), more.begin(), more.end());
      } else {
        //printf("    [%d] is neutral and it flies; ignoring\n",dau->getPDG());
      }
    }
  }
  //printf("---- returning promptTracks for particle: %d\n",getPDG());
  return ret;
}

void MCParticle::addDaughter(const MCParticle* mcp) {
  _dau.push_back(mcp);
}

double MCParticle::getD0() const {
  float pos[3] = {
    static_cast<float>(getVertex().x()),
    static_cast<float>(getVertex().y()),
    static_cast<float>(getVertex().z())
  };
  float mom[3] = {
    static_cast<float>(Px()),
    static_cast<float>(Py()),
    static_cast<float>(Pz())
  };
  HelixClass h;
  h.Initialize_VP(pos,mom,getCharge(),Globals::Instance()->getBField());
  return h.getD0();
}
double MCParticle::getZ0() const {
  float pos[3] = {
    static_cast<float>(getVertex().x()),
    static_cast<float>(getVertex().y()),
    static_cast<float>(getVertex().z())
  };
  float mom[3] = {
    static_cast<float>(Px()),
    static_cast<float>(Py()),
    static_cast<float>(Pz())
  };
  HelixClass h;
  h.Initialize_VP(pos,mom,getCharge(),Globals::Instance()->getBField());
  return h.getZ0();
}
double MCParticle::getPhi() const {
  float pos[3] = {
    static_cast<float>(getVertex().x()),
    static_cast<float>(getVertex().y()),
    static_cast<float>(getVertex().z())
  };
  float mom[3] = {
    static_cast<float>(Px()),
    static_cast<float>(Py()),
    static_cast<float>(Pz())
  };
  HelixClass h;
  h.Initialize_VP(pos,mom,getCharge(),Globals::Instance()->getBField());
  return h.getPhi0();
}
double MCParticle::getOmega() const {
  float pos[3] = {
    static_cast<float>(getVertex().x()),
    static_cast<float>(getVertex().y()),
    static_cast<float>(getVertex().z())
  };
  float mom[3] = {
    static_cast<float>(Px()),
    static_cast<float>(Py()),
    static_cast<float>(Pz())
  };
  HelixClass h;
  h.Initialize_VP(pos,mom,getCharge(),Globals::Instance()->getBField());
  return h.getOmega();
}
double MCParticle::getTanLambda() const {
  float pos[3] = {
    static_cast<float>(getVertex().x()),
    static_cast<float>(getVertex().y()),
    static_cast<float>(getVertex().z())
  };
  float mom[3] = {
    static_cast<float>(Px()),
    static_cast<float>(Py()),
    static_cast<float>(Pz())
  };
  HelixClass h;
  h.Initialize_VP(pos,mom,getCharge(),Globals::Instance()->getBField());
  return h.getTanLambda();
}

void Vertex::add(const Track* trk) {
  _tracks.push_back(trk);
}
void Vertex::add(const Track* trk,double chi2) {
  _tracks.push_back(trk);
  _chi2Tracks[trk] = chi2;
}

const Track* Vertex::getWorstTrack() const {
  if (_chi2Tracks.size() ==0)return 0;

  map<const Track*, double>::const_iterator it;
  double chi2max = 0;
  const Track* worstTrack = 0;
  for (it = _chi2Tracks.begin(); it != _chi2Tracks.end(); it++) {
    if (chi2max < it->second) {
      worstTrack = it->first;
      chi2max = it->second;
    }
  }
  return worstTrack;
}

double Vertex::length(const Vertex* primary) const {
  TVector3 pos(_x,_y,_z);
  if (primary) {
    TVector3 ip(primary->_x,primary->_y,primary->_z);
    TVector3 dif(pos-ip);
    double dist = dif.Mag();
    return dist;
  } else {
    return pos.Mag();
  }
}

double Vertex::significance(const Vertex* primary) const {
  TVector3 pos(_x,_y,_z);
  TVector3 ip(primary->_x,primary->_y,primary->_z);
  TVector3 dif(pos-ip);
  //double dist = dif.Mag();

  TMatrixF dif1(1,3);
  TMatrixF dif2(3,1);
  dif1(0,0) = dif[0];
  dif1(0,1) = dif[1];
  dif1(0,2) = dif[2];

  dif2(0,0) = dif[0];
  dif2(1,0) = dif[1];
  dif2(2,0) = dif[2];

  TMatrixFSym mat(3);
  mat(0,0) = _cov[0];
  mat(0,1) = _cov[1];
  mat(1,1) = _cov[2];
  mat(0,2) = _cov[3];
  mat(1,2) = _cov[4];
  mat(2,2) = _cov[5];

  /*
  const double* pcov = primary->getCov();
  if (primary->covIsGood()) {
  	mat(0,0) += pcov[0];
  	mat(0,1) += pcov[1];
  	mat(1,1) += pcov[2];
  	mat(0,2) += pcov[3];
  	mat(1,2) += pcov[4];
  	mat(2,2) += pcov[5];
  }
  */

  double det;
  TMatrixFSym invmat = mat.Invert(&det);

  if (det<0) {
    cout << "Vertex::significance: matrix not invertible, using conservative errors" << endl;
    invmat(0,0) = 1e-2;
    invmat(1,1) = 1e-2;
    invmat(2,2) = 1e-2;
    invmat(0,1) = 0;
    invmat(0,2) = 0;
    invmat(1,2) = 0;
  }

  TMatrixF res1(3,1);
  res1.Mult(invmat,dif2);

  TMatrixF res2(1,1);
  res2.Mult(dif1,res1);

  /*
  printf("res1: %f,%f,%f\n",res1(0,0),res1(1,0),res1(2,0));
  printf("res2: %f\n",res2(0,0));
  */
  double ret = res2(0,0);
  if (ret<0) {
    cout << "Vertex::significance: significance is negative, forcing it to zero" << endl;
    ret = 0;
  }

  return sqrt(ret);
}

double Vertex::getPparallel(const TVector3& axis)const {
  TLorentzVector sum4v;
  for (unsigned int i=0; i<getTracks().size(); i++) {
    sum4v += *(TLorentzVector*)getTracks()[i];
  }
  return sum4v.Vect().Dot(axis.Unit());
}

double Vertex::getVertexMass(const Vertex* daughter, const TVector3* paxis, const double dmass, double* ppt, double* pp)const {
  TVector3 axis = (paxis ? *paxis : getPos()).Unit();

  TLorentzVector sum4v;
  if (daughter) {
    TVector3 vdau = daughter->getPos() - getPos();
    double gammabyboost = vdau.Mag() / 0.3 / 2;
    double pbyboost = gammabyboost * dmass;
    double pbytrack = daughter->getPparallel(vdau);

    cout << "getVertexMass: vdau = " << vdau.x() << ", " << vdau.y() << ", " << vdau.z();
    cout << ", gammabyboost = " << gammabyboost << ", pbytrack = " << pbytrack ;
    vdau.SetMag(max(pbyboost, pbytrack));
    sum4v.SetVectM(vdau, dmass);
    cout << ", p = " << sum4v.P() << endl;
  }

  for (unsigned int i=0; i<getTracks().size(); i++) {
    sum4v += *(TLorentzVector*)getTracks()[i];
    cout << "sum4v e = " << sum4v.E() << ", p = " << sum4v.P() << ", mass = " << sum4v.M() << endl;
  }

  // compensation vector (only size=energy)
  double pt = sum4v.Pt(axis);

  // vertex mass
  double e = sum4v.E() + pt;
  double p = sum4v.Vect().Dot(axis);
  if (ppt)*ppt = pt;
  if (pp)*pp = p;

  cout << "pt = " << pt << ", e  = " << e << ", p  = " << p << endl;
  return sqrt(e*e-p*p);
}

double Vertex::getVertexAngle(const Vertex* vdaughter, const Vertex* primary)const {
  TVector3 posip;
  if (primary)
    posip = primary->getPos();

  TVector3 posdaughter = vdaughter->getPos();
  TVector3 pos = getPos();

  // todo:

  return 0;
}

TLorentzVector Vertex::getFourMomentum() const {
  TLorentzVector ret;
  TrackVecIte iter;
  for (iter = _tracks.begin(); iter != _tracks.end(); ++iter) {
    ret += **iter;
  }
  return ret;
}

double Vertex::getChi2TrackFit(const Track* tr, int mode)const {
  VertexFitterSimple_V fitter;
  //cout << "getChi2TrackFit " << getCov()[0] << " " << getCov()[2] << " " << getCov()[5] << " " << tr->getCovMatrix()[tpar::d0d0] << " " << tr->getCovMatrix()[tpar::z0z0] << " " << tr->getCovMatrix()[tpar::omom] << " " << tr->getCovMatrix() [tpar::phph] << " " << tr->getCovMatrix()[tpar::tdtd] << " " << fitter.getChi2(this,tr) << endl;
  return fitter.getChi2(this,tr,mode);
}

TrackPocaXY::TrackPocaXY(const Track* trk, const Vertex* vtx) :
  _success(false), _trk(trk), _flt(0.), _poca(0.),
  _v(vtx->getX(),vtx->getY(),vtx->getZ()) {
  minimize();
}

TrackPocaXY::TrackPocaXY(const Track* trk, const TVector3& v) :
  _success(false), _trk(trk), _flt(0.), _poca(0.), _v(v) {
  minimize();
}

void TrackPocaXY::minimize() {
  // minimization by golden section
  // faster algorithm should be implemented after this is shown to work
  //double bdy = 0.5*3.1415926/_trk->getOmega();
  double bdy = 2.0*3.1415926/fabs(_trk->getOmega());
  _poca = golden<TrackPocaXY>(-bdy,0,bdy,this,1e-5,_flt);
  //printf("poca/flt %f/%f\n",_poca,_flt);
  _success = true;
}

double TrackPocaXY::operator()(double* flt) const {
  double d0 = _trk->getD0();
  double ph = _trk->getPhi();
  double om = _trk->getOmega();
  double omInv = 1./om;

  double x = (-d0+omInv)*sin(ph) - omInv*sin(ph-om*(*flt));
  double y = ( d0-omInv)*cos(ph) + omInv*cos(ph-om*(*flt));

  double x0 = _v.x();
  double y0 = _v.y();

  double distSq = pow(x-x0,2) + pow(y-y0,2);

  return distSq;
}

const MCParticle* Vertex::getMcp()const {
  map<const MCParticle*, int>::iterator it;
  map<const MCParticle*, int> mapmcp;
  for (unsigned int i=0; i<getTracks().size(); i++) {
    const MCParticle* mcp = getTracks()[i]->getMcp()->getSemiStableParent();
    if (!mcp)continue;
    if (mcp->getPDG() == 111) {
      const MCParticle* mcd = getTracks()[i]->getMcp();
      cout << "pi0 found!" << mcd->getPDG() << ", " << mcd->getParent()->getPDG() << ", " << mcd->getParent()->getParent()->getPDG() << endl;
    }
    it = mapmcp.find(mcp);
    if (it == mapmcp.end())
      mapmcp[mcp] = 1;
    else
      it->second = it->second + 1;
  }
  int nmcpmax = 0;
  const MCParticle* mcpmax = 0;
  for (it = mapmcp.begin(); it != mapmcp.end(); it++) {
    if (it->second > nmcpmax) {
      mcpmax = it->first;
      nmcpmax = it->second;
    }
    cout << it->second << ": " << it->first->getPDG() << ", ";
  }
  cout << endl;
  return mcpmax;
}

double Vertex::dirdot(const Vertex* primary) const {
  TLorentzVector sum;
  for (TrackVecIte iter = _tracks.begin(); iter != _tracks.end(); ++iter) {
    sum += **iter;
  }
  TVector3 posip;
  if (primary) posip = primary->getPos();
  return (getPos()-posip).Unit().Dot( sum.Vect().Unit() );
}

bool Vertex::passesV0selection(const Vertex* primary) const {
  // only look at 2-track vertices
  if (_tracks.size() != 2) return false;
  // require total zero charge
  if (_tracks[0]->getCharge()*_tracks[1]->getCharge() != -1) return false;

  // check K_S
  TLorentzVector sum(*_tracks[0]+*_tracks[1]);
  TVector3 posip;
  if (primary) posip = primary->getPos();
  double v0dirdot = dirdot(primary);
  double ksmass = sum.M();
  if (v0dirdot > 0.999 && fabs(ksmass-0.4976) < 0.015 )
    return true;

  // check lambda; assume higher momentum particle is the proton
  TLorentzVector protonForLambda;
  TLorentzVector pionForLambda;
  if (_tracks[0]->Vect().Mag() > _tracks[1]->Vect().Mag()) {
    protonForLambda.SetVectM( _tracks[0]->Vect(), 0.9383 );
    pionForLambda = *_tracks[1];
  } else {
    protonForLambda.SetVectM( _tracks[1]->Vect(), 0.9383 );
    pionForLambda = *_tracks[0];
  }
  double lambdaMass = (protonForLambda+pionForLambda).M();
  if (v0dirdot > 0.99995 && fabs(lambdaMass - 1.1157) < 0.02)
    return true;

  // check photon conversion; photon mass corrected geometrically
  double ang1 = atan(_tracks[0]->getTanLambda());
  double ang2 = atan(_tracks[1]->getTanLambda());
  double convmassCorr = sqrt( _tracks[0]->Vect().Mag()*_tracks[1]->Vect().Mag()*(1-cos(ang1-ang2)) );
  if (v0dirdot > 0.99995 && getPos().Mag()>9 && convmassCorr < 0.01)
    return true;

  // didn't pass the v0 selection; this is probably not a v0.
  return false;
}

void Vertex::Print()const {
  cout << "Vertex::Print(): vpos = ( " << _x << " " << _y << " " << _z << ") ";
  cout << ", chi2 = " << _chi2 << ", prob = " << _prob << (_isPrimary ? ", primary " : "" ) << endl;
  cout << "   cov matrix: ( " << _cov[xx] << " " << _cov[xy] << " " << _cov[xz] << ")" << endl;
  cout << "               ( " << _cov[xy] << " " << _cov[yy] << " " << _cov[yz] << ")" << endl;
  cout << "               ( " << _cov[xz] << " " << _cov[yz] << " " << _cov[zz] << ")" << endl;

  for (unsigned int n=0; n<_tracks.size(); n++) {
    const Track* tr = _tracks[n];
    const MCParticle* mcp = tr->getMcp();
    const MCParticle* pmcp = (mcp ? mcp->getSemiStableParent() : 0);
    const MCParticle* pmcpb = (mcp ? mcp->getSemiStableBParent() : 0);

    cout << "   [track] E = " << tr->E() << ", p = (" << tr->Px() << " " << tr->Py() << " " << tr->Pz() << "), chi2 = " << getChi2Track(tr) << ", ";
    cout << "MCP = " << (mcp ? mcp->getPDG() : 0) << ", parent = " << (pmcp ? pmcp->getPDG() : 0) << ", parent B = " << (pmcpb ? pmcpb->getPDG() : 0) << endl;
    /*
    const double* tcov = tr->getCovMatrix();
    cout << "   trk cov: "
    	<< "d0d0= " << tcov[tpar::d0d0] << ", "
    	<< "d0ph= " << tcov[tpar::d0ph] << ", "
    	<< "phph= " << tcov[tpar::phph] << ", "
    	<< "d0om= " << tcov[tpar::d0om] << ", " << endl
    	<< "phom= " << tcov[tpar::phom] << ", "
    	<< "omom= " << tcov[tpar::omom] << ", "
    	<< "d0z0= " << tcov[tpar::d0z0] << ", "
    	<< "z0ph= " << tcov[tpar::z0ph] << ", " << endl
    	<< "z0om= " << tcov[tpar::z0om] << ", "
    	<< "z0z0= " << tcov[tpar::z0z0] << ", "
    	<< "d0td= " << tcov[tpar::d0td] << ", "
    	<< "phtd= " << tcov[tpar::phtd] << ", " << endl
    	<< "omtd= " << tcov[tpar::omtd] << ", "
    	<< "z0td= " << tcov[tpar::z0td] << ", "
    	<< "tdtd= " << tcov[tpar::tdtd] << endl;
    	*/
  }
}

int Vertex::dist_sort(const Vertex* a, const Vertex* b) {
  double a1 = a->getX()*a->getX()+a->getY()*a->getY()+a->getZ()*a->getZ();
  double b1 = b->getX()*b->getX()+b->getY()*b->getY()+b->getZ()*b->getZ();
  return (a1 < b1);
}

bool Vertex::covIsGood() const {
  if (_cov[0] != _cov[0]) return false;
  if (_cov[1] != _cov[1]) return false;
  if (_cov[2] != _cov[2]) return false;
  if (_cov[3] != _cov[3]) return false;
  if (_cov[4] != _cov[4]) return false;
  if (_cov[5] != _cov[5]) return false;
  return true;
}

Jet::Jet(const Track* trk) : TLorentzVector(*trk), _id(-1) {
  _tracks.push_back(trk);
}

Jet::Jet(const Neutral* neut) : TLorentzVector(*neut), _id(-1) {
  _neutrals.push_back(neut);
}

/* combines this jet with another jet
 * energy and momentum are added
 * daughter lists are combined
 */
void Jet::add( const Jet& jet ) {
  *this += jet;
  for (TrackVecIte iter = jet.getTracks().begin(); iter != jet.getTracks().end(); ++iter) {
    _tracks.push_back(*iter);
  }
  for (NeutralVecIte iter = jet.getNeutrals().begin(); iter != jet.getNeutrals().end(); ++iter) {
    _neutrals.push_back(*iter);
  }
  for (VertexVecIte iter = jet.getVertices().begin(); iter != jet.getVertices().end(); ++iter) {
    _vertices.push_back(*iter);
  }
}


Jet::Jet(const Jet& from, bool extractVertex) : TLorentzVector(from), _id(-1) {
  _tracks = from.getTracks();
  _neutrals = from.getNeutrals();

  if (extractVertex) {
    for (VertexVecIte iter = from.getVertices().begin(); iter != from.getVertices().end(); ++iter) {
      _tracks.insert(_tracks.end(), (*iter)->getTracks().begin(), (*iter)->getTracks().end());
    }
  } else {
    _vertices = from.getVertices();
  }

  _params = from.params();
}

// calculate boosted sphericity in the jet's CoM frame
double Jet::sphericity() const {

  // assumes I'm a flat jet (no vertex structure)
  TVector3 jetBoost = BoostVector();
  TMatrixDSym sphMat(3);
  sphMat(0,0) = 0;
  sphMat(0,1) = 0;
  sphMat(0,2) = 0;
  sphMat(1,1) = 0;
  sphMat(1,2) = 0;
  sphMat(2,2) = 0;

  TrackVec& tracks = getAllTracks();
  for (TrackVecIte iter = tracks.begin(); iter != tracks.end(); ++iter) {
    TLorentzVector trkVec(**iter);
    trkVec.Boost(-jetBoost);
    sphMat(0,0) += trkVec.X()*trkVec.X();
    sphMat(0,1) += trkVec.X()*trkVec.Y();
    sphMat(0,2) += trkVec.X()*trkVec.Z();
    sphMat(1,1) += trkVec.Y()*trkVec.Y();
    sphMat(1,2) += trkVec.Y()*trkVec.Z();
    sphMat(2,2) += trkVec.Z()*trkVec.Z();
  }
  double norm = sphMat(0,0)+sphMat(1,1)+sphMat(2,2);
  sphMat *= 1.5/norm;
  TVectorD eig;
  sphMat.EigenVectors(eig);
  return eig(1)+eig(2);
}

vector<const Track*> Jet::getAllTracks(bool withoutV0) const {
  vector<const Track*> ret;

  vector<const Track*>::const_iterator iter;
  for (iter = getTracks().begin(); iter != getTracks().end(); ++iter) {
    ret.push_back( *iter );
  }

  vector<const Vertex*>::const_iterator iter2;
  for (iter2 = getVertices().begin(); iter2 != getVertices().end(); ++iter2) {
    const Vertex* vtx = *iter2;
    if (withoutV0 && vtx->passesV0selection(Event::Instance()->getPrimaryVertex())) continue;
    for (iter = vtx->getTracks().begin(); iter != vtx->getTracks().end(); ++iter) {
      ret.push_back( *iter );
    }
  }

  return ret;
}

vector<const Vertex*> Jet::getVerticesForFT() const {

  vector<const Vertex*> ret;
  vector<const Vertex*>::const_iterator iter;
  for (iter = getVertices().begin(); iter != getVertices().end(); ++iter) {
    const Vertex* vtx = *iter;
    if (vtx->getTracks().size() > 1)
      ret.push_back( vtx );
  }

  sort(ret.begin(),ret.end(),Vertex::dist_sort);

  return ret;
}


// under construction
void Parameters::remove(const string& key, bool delmap) {
  if (_map.find(key) == _map.end())throw(Exception("Parameters::remove: key not found."));
  else if (_map[key].first == &typeid(float))deleteData<float>(key);
  else if (_map[key].first == &typeid(double))deleteData<double>(key);
  else if (_map[key].first == &typeid(long))deleteData<long>(key);
  else if (_map[key].first == &typeid(int))deleteData<int>(key);
  else if (_map[key].first == &typeid(bool))deleteData<bool>(key);
  else {
    // looking for TClass
    TClass* cl = TClass::GetClass(*(_map[key].first));
    if (!cl)throw(Exception(TString::Format("Parameters::remove: definition of class %s not found in TClass table: cannot delete it.",_map[key].first->name())));

    cl->Destructor(_map[key].second);
  }

  if (delmap)
    _map.erase(key);
}

void Parameters::clear() {

  map<string, pair<const type_info*, void*> >::iterator it;
  for (it = _map.begin(); it != _map.end(); it++) {
    remove(it->first, false);
  }

  _map.clear();
}

Parameters& Parameters::operator =(const Parameters& ref) {

  clear();

  map<string, pair<const type_info*, void*> >::const_iterator it;
  for (it = ref._map.begin(); it != ref._map.end(); it++) {
    void* data = it->second.second;
    const char* key = it->first.c_str();
    const type_info* type = it->second.first;

    // built-in types
    if (type == &typeid(double))add(key,*(double*)data);
    else if (type == &typeid(float))add(key,*(float*)data);
    else if (type == &typeid(long))add(key,*(long*)data);
    else if (type == &typeid(int))add(key,*(int*)data);
    else if (type == &typeid(bool))add(key,*(bool*)data);
    else {
      // looking for TClass
      TClass* cl = TClass::GetClass(*type);
      if (!cl)throw(Exception(TString::Format("Parameters::operator=: definition of class %s not found in TClass table: cannot copy it.",type->name())));

      // looking for operator =
      TMethodCall mc;
      TString opetype = TString::Format("const %s &", cl->GetName());
      mc.InitWithPrototype(cl, "operator=", opetype);
      mc.SetParamPtrs(&data, 1);

      if (!mc.IsValid())
        throw(Exception(TString::Format("Parameters::operator=: class copy operator not found, class name = %s, parameter type = %s", type->name(), (const char*)opetype)));

      void* newdata = cl->New();
      mc.Execute(newdata);

      _map[it->first] = pair<const type_info*, void*>(type, newdata);
    }
  }
  return *this;
}

}


