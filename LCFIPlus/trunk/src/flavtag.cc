#include "flavtag.h"
#include "EventStore.h"

#include "TVector3.h"
#include "TLorentzVector.h"
#include "TMatrixFSym.h"
#include "HelixClass.h"

#include "math.h"

#include "algo.h"

namespace flavtag {
//   LorentzVector::LorentzVector(float px, float py, float pz, float en) :
//     _px(px), _py(py), _pz(pz), _en(en) {
//   }
// 
//   float LorentzVector::getPt() const {
//     return sqrt( _px*_px + _py*_py );
//   }
// 
//   float LorentzVector::getMomentum() const {
//     return sqrt( _px*_px + _py*_py + _pz*_pz );
//   }
// 
//   TVector3 LorentzVector::getTVector3() const {
//     TVector3 a(_px,_py,_pz);
//     return a;
//   }
// 
//   TLorentzVector LorentzVector::getTLorentzVector() const {
//     TLorentzVector a(_px,_py,_pz,_en);
//     return a;
//   }
// 
//   void LorentzVector::add(const LorentzVector& d) {
//     _px += d._px;
//     _py += d._py;
//     _pz += d._pz;
//     _en += d._en;
//   }

/*
  Event::Event(const EventData& data) : _rescaleTableRead(false) {
    for (vector<TrackData>::const_iterator it = data.tracks.begin(); it != data.tracks.end(); ++it) {
      Track* trk = new Track(*it);
      trk->setEvent(this);
      _tracks.push_back( trk );
    }
    for (vector<NeutralData>::const_iterator it = data.neutrals.begin(); it != data.neutrals.end(); ++it) {
      Neutral* neut = new Neutral(*it);
      neut->setEvent(this);
      _neutrals.push_back( neut );
    }
    for (vector<MCParticleData>::const_iterator it = data.mcps.begin(); it != data.mcps.end(); ++it) {
      MCParticle* mcp = new MCParticle(*it);
      mcp->setEvent(this);
      _mcps.push_back( mcp );

      /// build parent-daughter association here
      MCParticle* parent = mcp->getParent();
      if (parent) {
        parent->addDaughter(mcp);
      }

    }
  }
*/

  Event::Event(const Event& event) : _rescaleTableRead(false) {
    for (vector<Track*>::const_iterator it = event._tracks.begin(); it != event._tracks.end(); ++it) {
      Track* trk = new Track(**it);
//      trk->setEvent(this);
      _tracks.push_back( trk );
    }
    for (vector<Neutral*>::const_iterator it = event._neutrals.begin(); it != event._neutrals.end(); ++it) {
      Neutral* neut = new Neutral(**it);
//      neut->setEvent(this);
      _neutrals.push_back( neut );
    }
    for (vector<MCParticle*>::const_iterator it = event._mcps.begin(); it != event._mcps.end(); ++it) {
      MCParticle* mcp = new MCParticle(**it);
//      mcp->setEvent(this);
      _mcps.push_back( mcp );

			/*
				 no need to do this since presumably this has been
				 already done with the first call to Init()

      /// build parent-daughter association here
      MCParticle* parent = mcp->getParent();
      if (parent) {
        parent->addDaughter(mcp);
      }
			*/

    }
  }

	Event::Event(const char *nameTracks, const char *nameNeutrals, const char *nameMCPs) : _rescaleTableRead(false)
	{
		const vector<Track *> * pTracks;
		const vector<Neutral *> * pNeutrals;
		const vector<MCParticle *> *pMCPs;

		EventStore::Instance()->Get<Track>(nameTracks,pTracks);
		EventStore::Instance()->Get<Neutral>(nameNeutrals,pNeutrals);
		EventStore::Instance()->Get<MCParticle>(nameMCPs,pMCPs);

		if(!pTracks || !pNeutrals || !pMCPs )throw(new Exception("event constructor:Failed to get collections!"));

		Init(pTracks, pNeutrals, pMCPs);
	}

// copy mode
	void Event::Init(const vector<Track *> *ptracks, const vector<Neutral *> *pneutrals, const vector<MCParticle *> *pmcps)
	{
		_tracks.resize(ptracks->size());
		_neutrals.resize(pneutrals->size());
		_mcps.resize(pmcps->size());

		copy(ptracks->begin(), ptracks->end(), _tracks.begin());
		copy(pneutrals->begin(), pneutrals->end(), _neutrals.begin());
		copy(pmcps->begin(), pmcps->end(), _mcps.begin());
  }


  Event::~Event() {
//    for_each( _tracks.begin(), _tracks.end(), DeleteVector<Track*>() );
//    for_each( _neutrals.begin(), _neutrals.end(), DeleteVector<Neutral*>() );
//    for_each( _mcps.begin(), _mcps.end(), DeleteVector<MCParticle*>() );
  }

  int Event::mcNumberOfB() const {
    int num(0);
    for (vector<MCParticle*>::const_iterator it = _mcps.begin(); it != _mcps.end(); ++it) {
      MCParticle* mcp = *it;
      if (mcp->isSemiStableB()) ++num;
    }
    return num;
  }

  int Event::mcNumberOfC() const {
    int num(0);
    for (vector<MCParticle*>::const_iterator it = _mcps.begin(); it != _mcps.end(); ++it) {
      MCParticle* mcp = *it;
      if (mcp->isSemiStableC()) ++num;
    }
    return num;
  }

  MCParticle* Event::getMCParticle(int id) const {
    if (id == 0) return 0;
    for (vector<MCParticle*>::const_iterator it = _mcps.begin(); it != _mcps.end(); ++it) {
      MCParticle* mcp = *it;
      if (mcp->getId() == id) return mcp;
    }
    return 0;
  }

  MCParticle* Event::getMCParticle(const Track* trk) const {
    return trk->getMcp();
  }

  vector<MCParticle*> Event::mcGetColorStrings() const {
    vector<MCParticle*> colorStrings;
    for (vector<MCParticle*>::const_iterator it = _mcps.begin(); it != _mcps.end(); ++it) {
      MCParticle* cs = (*it)->getColorString();
      if (find( colorStrings.begin(), colorStrings.end(), cs ) == colorStrings.end() ) {
        colorStrings.push_back(cs);
      }
    }
    return colorStrings;
  }

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

    for (vector<Track*>::iterator iter = _tracks.begin(); iter != _tracks.end(); ++iter) {
      Track* trk = *iter;
      float c = fabs(sin(atan(trk->getTanLambda())));
      float p = trk->P();
      float d0pull(1);
      float z0pull(1);
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

      const float* cov = trk->getCovMatrix();
      float newCov[tpar::covN];
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

  void Event::deleteTracks() {
    int ntmp = _tracks.size();
    for_each( _tracks.begin(), _tracks.end(), DeleteVector<Track*>() );
    _tracks.clear();
    _tracks.resize(ntmp);
  }

//   Track::Track(const TrackData& data) :
//     LorentzVector(data.px, data.py, data.pz, data.en),
//     TrackData(data),
//     _flt(0) {
//   }

  float Track::getX() const {
    return  ( -_par[tpar::d0] + 1./_par[tpar::om] )*sin( _par[tpar::ph] )
      -  ( sin( _par[tpar::ph] - _par[tpar::om]*_flt ) )/_par[tpar::om];
  }

  float Track::getY() const {
    return  (  _par[tpar::d0] - 1./_par[tpar::om] )*cos( _par[tpar::ph] )
      +  ( cos( _par[tpar::ph] - _par[tpar::om]*_flt ) )/_par[tpar::om];
  }

  float Track::getZ() const {
    return _par[tpar::z0] + _par[tpar::om]*_flt;
  }

  void Track::setCovMatrix(float* mycov) {
    for (int i=0; i<tpar::covN; ++i) {
      _cov[i] = mycov[i];
    }
  }

//   Neutral::Neutral(const NeutralData& data) :
//     LorentzVector(data.px, data.py, data.pz, data.en), NeutralData(data) {
//   }

  void MCParticle::Init(int id,int pdg, MCParticle *parent, float charge, const TLorentzVector &p, const TVector3 &v)
	{
		*(TLorentzVector *)this = p;
		_dauForDecay = 0;
		_helixOK = false;

		_id = id;
		_pdg = pdg;
		_parent = parent;
		_vertex = v;
		_charge = charge;

		// add myself to the parent particle
		if(_parent)
			_parent->addDaughter(this);
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

  MCParticle* MCParticle::getColorString() {
    if (getPDG() == 92) return this;
    MCParticle* parent = getParent();
    if (parent == 0) return parent;
    return parent->getColorString();
  }

  // checks if the MCParticle underwent semileptonic decay.
  // only checks the immediate daughters.
  // if it is found to have decayed semileptonically (e, mu, tau),
  // the pointer to the lepton is returned; zero otherwise.
  MCParticle* MCParticle::semileptonicDecay() const {
/*    const vector<MCParticle*>& mcps = _event->getMCParticles();
    for (vector<MCParticle*>::const_iterator it = mcps.begin(); it != mcps.end(); ++it) {
      MCParticle* test = *it;
      int abspdg = abs(test->getPDG());
      if (abspdg == 11 || abspdg == 13 || abspdg == 15) {
        if (test->getParent() == this)
          return test;
      }
    }
    return 0;*/
		for(unsigned int i=0;i<_dau.size();i++){
			MCParticle *test = _dau[i];
      int abspdg = abs(test->getPDG());
      if (abspdg == 11 || abspdg == 13 || abspdg == 15) {
        return test;
      }
		}

		return 0;
  }

  bool MCParticle::isStableTrack() const {
    switch(abs(_pdg)) {
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
		switch(abs(_pdg)) {
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
    switch(abs(_pdg)) {
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
    switch(abs(_pdg)) {
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
    switch(abs(_pdg)) {
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
    switch(abs(_pdg)) {
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

  MCParticle* MCParticle::getSemiStableParent() const {
    MCParticle* parent = getParent();
    if (parent == 0) return 0;
    /*
    printf(" current id: %d (pdg=%d), parent id: %d (pdg=%d)\n",
        id,getPDG(),parent->getId(),parent->getPDG() );
        */
    if (parent->isSemiStable()) return parent;
    if (parent->isStable()) return parent;
    return parent->getSemiStableParent();
  }

  // classify the MCParticle according to its position in the MC tree
  //   1 = MCParticle originates from primary
  //   2 = MCParticle is a result of a semi-stable b decay (e.g. B -> D l nu)
  //   3 = MCParticle is a result of a semi-stable c decay (e.g. D -> K pi)
  //   4 = MCParticle is a result of a semi-stable s decay (e.g. Ks -> pi+ pi-)
  // iterates over parents to find the most immediate semi-stable particle
  int MCParticle::getFlavorTagCategory() const {
    int tag(0);
    MCParticle* ssp = getSemiStableParent();

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

	const TVector3& MCParticle::getEndVertex() {
		findDauForDecay();
		return _dauForDecay->getVertex();
	}

  // find the decay distance by iterating over daughters
  // and grand-daughters (and possibly more generations)
  // to find the shortest distance where a daughter tracks splits off
  float MCParticle::decayDistance() {
		findDauForDecay();
    return ( (_dauForDecay->getVertex()-getVertex()).Mag());
  }

  MCParticle* MCParticle::findDauForDecay() {
		if (_dauForDecay) return _dauForDecay;

		for (unsigned int i=0; i<_dau.size(); ++i) {
			MCParticle* dau = _dau[i];
			//printf("  checking dau: %d\n",dau->getPDG());
			float test = ( dau->getVertex()-getVertex() ).Mag();
			if (test > 0 && (dau->isStable() || dau->isSemiStable())) {
				//printf("  dau %d found \n",dau->getPDG());
				_dauForDecay = dau;
				return dau;
			}
		}

		// look at grand daughters if the immediate daughter is unstable
		// this will actually iterate until we find the first stable particle
		for (unsigned int i=0; i<_dau.size(); ++i) {
      MCParticle* dau = _dau[i];
			//printf("  calling findDauForDecay for dau %d\n", dau->getPDG());
			MCParticle* dautmp = dau->findDauForDecay();
			if (dautmp != dau) {
				_dauForDecay = dautmp;
				return _dauForDecay;
			}
    }

		_dauForDecay = this;
		return this;
  }

  float MCParticle::getEx() {
    if (_dauForDecay == 0) findDauForDecay();
    return _dauForDecay->getVertex().x();
  }

  float MCParticle::getEy() {
    if (_dauForDecay == 0) findDauForDecay();
    return _dauForDecay->getVertex().y();
  }

  float MCParticle::getEz() {
    if (_dauForDecay == 0) findDauForDecay();
    return _dauForDecay->getVertex().z();
  }

	bool MCParticle::isParent(MCParticle *mcp)const
	{
		if(getParent() == 0)return false;
		if(getParent() == mcp)return true;

		return getParent()->isParent(mcp);
	}

  // find all the tracks which originate promptly,
  // including tracks from unstable daughters
  vector<MCParticle*> MCParticle::promptTracks() {
    vector<MCParticle*> ret;
    //printf("promptTracks called for particle: %d\n",getPDG());
    //printf("  has daughters: ");
    for (vector<MCParticle*>::const_iterator it = _dau.begin(); it != _dau.end(); ++it) {
      //MCParticle* dau = *it;
      //printf("%d ",dau->getPDG());
    }
    //printf("\n");
    for (vector<MCParticle*>::const_iterator it = _dau.begin(); it != _dau.end(); ++it) {
      MCParticle* dau = *it;
      //printf("  - checking daughter: %d\n",dau->getPDG());
      if (dau->getCharge() != 0) {
        if ( dau->isStableTrack() ) {
          //printf("    [%d] is stable; added\n",dau->getPDG());
          ret.push_back(dau);
        } else if ( dau->isSemiStable() ) {
          //printf("    [%d] is semistable; ignored\n",dau->getPDG());
        } else {
          //printf("    [%d] not stable, looking at daughter\n",dau->getPDG());
          vector<MCParticle*> more = dau->promptTracks();
          ret.insert(ret.end(), more.begin(), more.end());
        }
      } else {
        // treat neutral particles
        if ( dau->decayDistance() == decayDistance() ) {
          //printf("    [%d] is neutral and decays promptly, looking at daughter\n",dau->getPDG());
          vector<MCParticle*> more = dau->promptTracks();
          ret.insert(ret.end(), more.begin(), more.end());
        } else {
          //printf("    [%d] is neutral and it flies; ignoring\n",dau->getPDG());
        }
      }
    }
    //printf("---- returning promptTracks for particle: %d\n",getPDG());
    return ret;
  }

  void MCParticle::addDaughter(MCParticle* mcp) {
    _dau.push_back(mcp);
  }

  void MCParticle::makeHelix() {
    if (_helixOK) return;
    float pos[3] = { getVertex().x(), getVertex().y(), getVertex().z() };
    float mom[3] = { Px(), Py(), Pz() };
    _helix.Initialize_VP(pos,mom,getCharge(),3.5);
    _helixOK = true;
  }

  float MCParticle::getD0() {
    makeHelix();
    return _helix.getD0();
  }
  float MCParticle::getZ0() {
    makeHelix();
    return _helix.getZ0();
  }
  float MCParticle::getPhi() {
    makeHelix();
    return _helix.getPhi0();
  }
  float MCParticle::getOmega() {
    makeHelix();
    return _helix.getOmega();
  }
  float MCParticle::getTanLambda() {
    makeHelix();
    return _helix.getTanLambda();
  }

  void Vertex::add(Track* trk) {
    _tracks.push_back(trk);
  }
  void Vertex::add(Track* trk,float chi2){
		_tracks.push_back(trk);
		_chi2Tracks[trk] = chi2;
	}

	Track * Vertex::getWorstTrack() const{
		if(_chi2Tracks.size() ==0)return 0;

		map<Track *, float>::const_iterator it;
		float chi2max = 0;
		Track * worstTrack = 0;
		for(it = _chi2Tracks.begin(); it != _chi2Tracks.end();it++){
			if(chi2max < it->second){
				worstTrack = it->first;
				chi2max = it->second;
			}
		}
		return worstTrack;
	}

  float Vertex::length(const Vertex* primary) const {
    TVector3 pos(_x,_y,_z);
		if (primary) {
			TVector3 ip(primary->_x,primary->_y,primary->_z);
			TVector3 dif(pos-ip);
			float dist = dif.Mag();
			return dist;
		} else {
			return pos.Mag();
		}
  }

  float Vertex::significance(const Vertex* primary) const {
    TVector3 pos(_x,_y,_z);
    TVector3 ip(primary->_x,primary->_y,primary->_z);
    TVector3 dif(pos-ip);
    //float dist = dif.Mag();

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

    double det;
    TMatrixFSym invmat = mat.Invert(&det);

    if (det<0) return -1;

    TMatrixF res1(3,1);
    res1.Mult(invmat,dif2);

    TMatrixF res2(1,1);
    res2.Mult(dif1,res1);

    /*
    printf("res1: %f,%f,%f\n",res1(0,0),res1(1,0),res1(2,0));
    printf("res2: %f\n",res2(0,0));
    */

    return sqrt(res2(0,0));
  }

	float Vertex::getPparallel(const TVector3 &axis)const{
		TLorentzVector sum4v;
		for(unsigned int i=0;i<getTracks().size();i++){
			sum4v += *(TLorentzVector *)getTracks()[i];
		}
		return sum4v.Vect().Dot(axis.Unit());
	}

  float Vertex::getVertexMass(const Vertex *daughter, const TVector3 *paxis, const double dmass, double *ppt, double *pp)const{
		TVector3 axis = (paxis ? *paxis : getPos()).Unit();

		TLorentzVector sum4v;
		if(daughter){
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

		for(unsigned int i=0;i<getTracks().size();i++){
			sum4v += *(TLorentzVector *)getTracks()[i];
			cout << "sum4v e = " << sum4v.E() << ", p = " << sum4v.P() << ", mass = " << sum4v.M() << endl;
		}

		// compensation vector (only size=energy)
		double pt = sum4v.Pt(axis);

		// vertex mass
		double e = sum4v.E() + pt;
		double p = sum4v.Vect().Dot(axis);
		if(ppt)*ppt = pt;
		if(pp)*pp = p;

		cout << "pt = " << pt << ", e  = " << e << ", p  = " << p << endl;
		return sqrt(e*e-p*p);
	}

	float Vertex::getVertexAngle(const Vertex *vdaughter, const Vertex *primary)const
	{
		TVector3 posip;
		if(primary)
			posip = primary->getPos();

		TVector3 posdaughter = vdaughter->getPos();
		TVector3 pos = getPos();

		// todo: 

		return 0;
	}

	TLorentzVector Vertex::getFourMomentum() const {
		TLorentzVector ret;
		vector<Track*>::const_iterator iter;
		for (iter = _tracks.begin(); iter != _tracks.end(); ++iter) {
			ret += **iter;
		}
		return ret;
	}

	TrackPocaXY::TrackPocaXY(const Track* trk, const Vertex* vtx) :
    _success(false), _trk(trk), _flt(0.), _poca(0.),
		_v(vtx->getX(),vtx->getY(),vtx->getZ()) {
			minimize();
	}

  TrackPocaXY::TrackPocaXY(const Track* trk, const TVector3 &v) :
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

MCParticle * Vertex::getMcp()const
{
	map<MCParticle *, int>::iterator it;
	map<MCParticle *, int> mapmcp;
	for(unsigned int i=0;i<getTracks().size();i++){
		MCParticle *mcp = getTracks()[i]->getMcp()->getSemiStableParent();
		if(!mcp)continue;
		if(mcp->getPDG() == 111){
		  MCParticle *mcd = getTracks()[i]->getMcp();
		  cout << "pi0 found!" << mcd->getPDG() << ", " << mcd->getParent()->getPDG() << ", " << mcd->getParent()->getParent()->getPDG() << endl;
		}
		it = mapmcp.find(mcp);
		if(it == mapmcp.end())
			mapmcp[mcp] = 1;
		else
		  it->second = it->second + 1;
	}
	int nmcpmax = 0;MCParticle *mcpmax = 0;
	for(it = mapmcp.begin(); it != mapmcp.end();it++){
		if(it->second > nmcpmax){mcpmax = it->first;nmcpmax = it->second;}
		cout << it->second << ": " << it->first->getPDG() << ", ";			
	}
	cout << endl;
	return mcpmax;	
}

double Vertex::dirdot(const Vertex* primary) const {
	TLorentzVector sum;
	for (vector<Track*>::const_iterator iter = _tracks.begin(); iter != _tracks.end(); ++iter) {
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

}
