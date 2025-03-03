// LCIOStorer.cc

// LCIO includes
#include "lcio.h"
#include "EVENT/LCRunHeader.h"
#include "UTIL/LCTOOLS.h"
#include "EVENT/LCCollection.h"

#include "IMPL/LCCollectionVec.h"
#include "IMPL/VertexImpl.h"
#include "IMPL/ReconstructedParticleImpl.h"

#include "EVENT/MCParticle.h"
#include <EVENT/Track.h>
#include <EVENT/Vertex.h>
#include "EVENT/LCIO.h"

#include "UTIL/LCRelationNavigator.h"
#include "UTIL/PIDHandler.h"
#include "UTIL/ILDConf.h"

// lcfiplus includes
#include "LCIOStorer.h"
#include "lcfiplus.h"
#include "EventStore.h"
#include "JetFinder.h"

#include <assert.h>
#include <math.h>
#include <fstream>

using namespace lcio;

namespace lcfiplus {

// ctor/dtor
LCIOStorer::LCIOStorer(const char* inputfile, const char* outputfile, bool autoread, bool autowrite, const char* outPrefix) {
  if (inputfile) {
    _reader = lcio::LCFactory::getInstance()->createLCReader();

    // if input file is not .slcio, treat as a list
    if (strstr(inputfile, ".slcio") != inputfile + strlen(inputfile) - strlen(".slcio")) {
      // open stream
      ifstream filelist(inputfile);
      vector<string>* svec = new vector<string>;
      copy(istream_iterator<string>(filelist), istream_iterator<string>(), back_inserter(*svec));

      cout << "Input " << svec->size()  << " LCIO files:" << endl;
      for (unsigned int i=0; i<svec->size(); i++)
        cout << "    " << (*svec)[i] << endl;

      _reader->open(*svec);
    } else {
      _reader->open(inputfile); // may throw exception...
      cout << "Input LCIO file open: " << inputfile << " successful." << endl;
    }
  } else
    _reader = NULL;

  if (outputfile) {
    _writer = lcio::LCFactory::getInstance()->createLCWriter();
    _writer->open(outputfile);
    cout << "Output LCIO file open: " << outputfile << " successful." << endl;
  } else
    _writer = NULL;

  // copy run header
  if (_reader && _writer) {
    /*			int ncp = 0;
    			lcio::LCRunHeader *header = 0;
    			while((header = _reader->readNextRunHeader())){
    				_writer->writeRunHeader(header);
    				ncp ++;
    			}
    			cout << "LCIOStorer: copy " << ncp << " run headers." << endl;
    */
    lcio::LCRunHeader* header = _reader->readNextRunHeader();
    if (header) {
      _writer->writeRunHeader(header);
      cout << "LCIOStorer: copied a run header." << endl;
    } else {
      cout << "LCIOStorer: the run header not found." << endl;
    }
  }

  _event = NULL;

  _autowrite = autowrite;
  _autoread = autoread;
  _savePrefix = (outPrefix ? outPrefix : "");


  // register to EventStore
  Event::Instance()->RegisterObserver(this);

  cout << "LCIOStorer: initialization finished." << endl;
}

LCIOStorer::~LCIOStorer() {
  if (_reader)_reader->close();
  if (_writer)_writer->close();
}

void LCIOStorer::InitMCPPFOCollections(const char* pfoColName, const char* mcColName, const char* mcpfoColName) {
  Event* event = Event::Instance();

  vector<Track*>* pTracks;
  vector<Neutral*>* pNeutrals;
  vector<MCParticle*>* pMCPs;
  vector<MCColorSinglet*>* pMCCSs;

  // pfo collection
  if (!event->IsExist(pfoColName)) {
    event->Register<Track>(pfoColName,pTracks);
    event->Register<Neutral>(pfoColName,pNeutrals);
    _importPFOCols[pfoColName] = make_pair(pTracks, pNeutrals);
  }

  // mcp collection
  if (!event->IsExist(mcColName)) {
    event->Register<MCParticle>(mcColName,pMCPs);
    _importMCPCols[mcColName] = pMCPs;

    event->Register<MCColorSinglet>(mcColName, pMCCSs);
    _importMCCSCols[mcColName] = pMCCSs;
  }

  // mcpfo relation
  if (_importMCPFOLinkCols.count(mcpfoColName) == 0) {
    _importMCPFOLinkCols[mcpfoColName] = make_pair(mcColName, pfoColName);
  }

  cout << "LCIOStorer initialized with MCParticle collection." << endl;
}

void LCIOStorer::InitPFOCollections(const char* pfoColName) {
  Event* event = Event::Instance();

  vector<Track*>* pTracks;
  vector<Neutral*>* pNeutrals;

  // pfo collection
  if (!event->IsExist(pfoColName)) {
    event->Register<Track>(pfoColName,pTracks);
    event->Register<Neutral>(pfoColName,pNeutrals);
    _importPFOCols[pfoColName] = make_pair(pTracks, pNeutrals);
  }

  cout << "LCIOStorer initialized without MCParticle collection." << endl;
}

void LCIOStorer::InitVertexCollection(const char* lcioName, const char* flavtagName, bool readnow) {
  cout << lcioName << " loading " << endl;
  vector<lcfiplus::Vertex*>* vtxcol = 0;
  Event::Instance()->Register<Vertex>(flavtagName,vtxcol);

  _importVertexCols[ lcioName ] = constVector(vtxcol);

  if (_event && readnow)
    ReadVertices(lcioName, constVector(vtxcol));
}

void LCIOStorer::InitJetCollection(const char* lcioName, const char* flavtagName, bool readnow, bool readvtx, const char* vtxname) {
  vector<lcfiplus::Jet*>* jetcol = 0;
  Event::Instance()->Register<Jet>(flavtagName,jetcol);

  _importJetCols[ lcioName ] = constVector(jetcol);
  if (readvtx) {
    string name = (vtxname ? vtxname : (string)lcioName + "_vtx");
    InitVertexCollection(name.c_str(), name.c_str(), readnow);
  }
  if (_event && readnow)
    ReadJets(lcioName, constVector(jetcol));
}

/*
	void LCIOStorer::InitVertexCollectionsAuto(lcio::LCEvent *evt)
	{
		const vector<string> * pcolnames = evt->getCollectionNames();
		for(unsigned int i=0;i<pcolnames->size();i++){
			lcio::LCCollection* colvtx = evt->getCollection((*pcolnames)[i]);
			const char *colname = (*pcolnames)[i].c_str();
			if((!Event::Instance()->IsExist(colname)) && colvtx->getTypeName() == lcio::LCIO::VERTEX){
				InitVertexCollection(colname, colname);
			}
		}
	}

	void LCIOStorer::InitJetCollectionsAuto(lcio::LCEvent *evt)
	{
		const vector<string> * pcolnames = evt->getCollectionNames();
		for(unsigned int i=0;i<pcolnames->size();i++){
			lcio::LCCollection* coljet = evt->getCollection((*pcolnames)[i]);
			const char *colname = (*pcolnames)[i].c_str();
			if((_pfoColName != colname) && (!Event::Instance()->IsExist(colname)) && coljet->getTypeName() == lcio::LCIO::RECONSTRUCTEDPARTICLE){
				InitJetCollection(colname, colname);
			}
		}
	}
*/

// to be used in standalone; not from LcfiplusProcessor
bool LCIOStorer::Next(bool /*autovertex*/, bool /*autojet*/) {
  lcio::LCEvent* evt = _reader->readNextEvent();
  if (!evt)return false;
  /*
  		if(autovertex)
  			InitVertexCollectionsAuto(evt);
  		if(autojet)
  			InitJetCollectionsAuto(evt);
  */
  SetEvent(evt);
  return true;
}

void LCIOStorer::AutoConvert() {
  EventStore* s = Event::Instance();
  const multimap<string, EventStore::StoredEntry>& emap = s->GetObjectMap();

  const vector<string>* lcionames = _event->getCollectionNames();

  map<string, EventStore::StoredEntry>::const_iterator it;
  for (it = emap.begin(); it != emap.end(); it++) {
    if (it->second.flag & EventStore::PERSIST) {
      // persist: to save
      string newname = _savePrefix + it->first;

      if (find(lcionames->begin(), lcionames->end(), newname) != lcionames->end()) {
        cout << "AutoConvert skipped because it is already included in lcio." << endl;
        continue; // skip collections that is already existed in lcio event
      }

      if (it->second.classname == "vector<lcfiplus::Vertex*>") {
        WriteVertices(it->first.c_str(), newname.c_str());
      } else if (it->second.classname == "vector<lcfiplus::Jet*>") {
        bool writeVertex = it->second.flag & EventStore::JET_WRITE_VERTEX;
        WriteJets(it->first.c_str(), newname.c_str(), writeVertex);
      } else {
        cout << "LCIOStorer::WriteEvent(): unknown conversion for lcio: class " << it->second.classname << ", ignored." << endl;
      }
    }
  }
}

void LCIOStorer::WriteEvent() {
  if (_writer && _event) {
    if (_autowrite) {
      AutoConvert();
    }
    _writer->writeEvent(_event);
  } else {
    cout << "LCIOStorer::WriteEvent(): event cannot be written." << endl;
  }
}

bool LCIOStorer::energy_sort_pfo(lcio::ReconstructedParticle* a, lcio::ReconstructedParticle* b) {
  assert(a != 0);
  assert(b != 0);
  return (a->getEnergy() > b->getEnergy());
}


void LCIOStorer::SetEvent(lcio::LCEvent* evt) {
  // buffer check
  _event = evt;

//		if(_importPFOCols.size() == 0)return;

  // clearing all collections
  Event::Instance()->ClearObjects();

  // clearing relations
  _trackLCIORel.clear();
  _neutralLCIORel.clear();
  _mcpLCIORel.clear();
  _vtxLCIORel.clear();
  _jetLCIORel.clear();

  _trackLCIORel2.clear();
  _neutralLCIORel2.clear();
  _mcpLCIORel2.clear();
  _vtxLCIORel2.clear();
  _jetLCIORel2.clear();

  // collections
//		lcio::LCCollection* colPFO = evt->getCollection(_pfoColName);
//		lcio::LCCollection* colMC = 0;
//		lcio::LCRelationNavigator* nav = 0;

  // convert MCParticle /////////////////////////////////////////////////////////////////////////
  map<string, vector<lcfiplus::MCParticle*> *>::iterator itMcpCol;
  for (itMcpCol = _importMCPCols.begin(); itMcpCol != _importMCPCols.end(); itMcpCol ++) {

    LCCollection* colMC = evt->getCollection(itMcpCol->first);

    // 160722 Separate parent-daughter connection from the main loop to avoid crash with overlaid particles which have opposite order
    int id;
    for (id = 0; id < colMC->getNumberOfElements(); ++id) {
      lcio::MCParticle* mcp = dynamic_cast<lcio::MCParticle*>( colMC->getElementAt(id) );
      assert(mcp != 0);
      //if (mcp->isBackscatter()) continue;

      // fix against weird 1e-308 numbers
      const double* vtmp = mcp->getVertex();
      double v[3] = { 0., 0., 0.};
      if (vtmp) {
        memcpy(v, vtmp, sizeof(v));
      }
      //printf("%d [%e,%e,%e]\n",mcp->getPDG(),v[0],v[1],v[2]);
      if (v[0]*v[0]+v[1]*v[1]+v[2]*v[2] < 1e-100) {
        v[0] = 0;
        v[1] = 0;
        v[2] = 0;
      }

      // convert into MCParticle: parent is empty at this point
      lcfiplus::MCParticle* mcpNew = new MCParticle(id, mcp->getPDG(), 0, mcp->getCharge(),
          TLorentzVector(TVector3(mcp->getMomentum()),mcp->getEnergy()), TVector3(v));

      // add to MCP list
      itMcpCol->second->push_back(mcpNew);

      // LCIO relation
      _mcpLCIORel[mcpNew] = mcp;
      _mcpLCIORel2[mcp] = mcpNew;
    }
    // second loop for parent daughter connection
    for (id = 0; id < colMC->getNumberOfElements(); ++id) {
      lcio::MCParticle* mcp = dynamic_cast<lcio::MCParticle*>( colMC->getElementAt(id) );
      assert(mcp != 0);
      //if (mcp->isBackscatter()) continue;

      // find parent id
      // TODO Really OK for single parent??
      lcfiplus::MCParticle* parent = 0;
      if (mcp->getParents().size()>0) {
        map<lcio::MCParticle*,lcfiplus::MCParticle*>::iterator iter = _mcpLCIORel2.find(mcp->getParents()[0]);

        if ( iter == _mcpLCIORel2.end() ) cout << "LCIOStorer::MCPconversion: parent not found in association map" << endl;
        else parent = iter->second;
      }

      lcfiplus::MCParticle *mcpNew = _mcpLCIORel2[mcp];
      mcpNew->setParent(parent);
      if (parent) parent->addDaughter(mcpNew);   

      // add daughters to other MCPs
      if (mcp->getParents().size()>1) {
        for (unsigned int i=1; i<mcp->getParents().size(); i++) {
          map<lcio::MCParticle*,lcfiplus::MCParticle*>::iterator iter = _mcpLCIORel2.find(mcp->getParents()[i]);
	  if ( iter == _mcpLCIORel2.end() ) cout << "LCIOStorer::MCPconversion: parent " << i << " not found in association map" << endl;
          else iter->second->addDaughter(mcpNew);
        }
      }
    }
    
    SetColorSinglets(*(itMcpCol->second), *(_importMCCSCols[itMcpCol->first]));
  }

  // convert PFO ////////////////////////////////////////////////////////////////////////////////
  map<string, pair<vector<lcfiplus::Track*> *, vector<lcfiplus::Neutral*> *> >::iterator itPfoCol;
  for (itPfoCol = _importPFOCols.begin(); itPfoCol != _importPFOCols.end(); itPfoCol ++) {

    LCCollection* colPFO = evt->getCollection(itPfoCol->first);
    PIDHandler PID(evt->getCollection(itPfoCol->first));

    // looking for LCRelation
    vector<lcio::LCRelationNavigator*> navs;
    map<string, pair<string,string> >::iterator itRelCol;
    for (itRelCol = _importMCPFOLinkCols.begin(); itRelCol != _importMCPFOLinkCols.end(); itRelCol ++) {
      if (itRelCol->second.second == itPfoCol->first) {
        navs.push_back(new LCRelationNavigator(evt->getCollection(itRelCol->first)));
      }
    }
    //cout << navs.size() << " relation found." << endl;

    // sorted PFO list
    vector<lcio::ReconstructedParticle*> pfo_list;
    for (int i=0; i<colPFO->getNumberOfElements(); ++i) {
      lcio::ReconstructedParticle* pfo = dynamic_cast<lcio::ReconstructedParticle*>( colPFO->getElementAt(i) );
      pfo_list.push_back(pfo);
    }
    // sort by energy order!
    sort(pfo_list.begin(),pfo_list.end(),energy_sort_pfo);

    unsigned int trkIdCounter(0);
    unsigned int neutIdCounter(0);
    // main loop for PFO
    for (unsigned int n=0; n<pfo_list.size(); ++n ) {
      lcio::ReconstructedParticle* pfo = pfo_list[n];
      lcio::MCParticle* mcp = NULL;
      lcfiplus::MCParticle* mcpf = NULL;

      // looking for MCParticle
      for (unsigned int n=0; n<navs.size(); n++) {
        if (navs[n]->getRelatedToObjects(pfo).size()) {
          mcp = dynamic_cast<lcio::MCParticle*>(navs[n]->getRelatedToObjects(pfo)[0]);  // TODO [0] OK?
          mcpf = _mcpLCIORel2[mcp];
          break;
        }
      }

      // find clusters
      vector<lcio::Cluster*> clusters = pfo->getClusters();
      double clusEnergy(0);
      double subE[6];
      for (int i=0; i<6; ++i) subE[i]=0;
      if (_readSubdetectorEnergies) {
        for (unsigned int iclus=0; iclus<clusters.size(); ++iclus) {
          clusEnergy += clusters[iclus]->getEnergy();
          for (int i=0; i<6; ++i) {
            subE[i] = clusters[iclus]->getSubdetectorEnergies()[i];
          }
        }
      }

//			cerr << (pfo->getCharge() ? "[Track]" : "[Neutral]") << pfo->getEnergy() << ", " << (unsigned int)pfo << endl;

      // convert Track ////////////////////////////////////////////////////////////////////////////////
      //pdg table
      map<int, double> pmass;
      pmass.insert(map<int, double>::value_type( 11, 0.000510998 ) );
      pmass.insert(map<int, double>::value_type( 13, 0.105658 ) );
      pmass.insert(map<int, double>::value_type( 211, 0.139570 ) );
      pmass.insert(map<int, double>::value_type( 321, 0.493677 ) );
      pmass.insert(map<int, double>::value_type( 2212, 0.938272 ) );
      
      if (pfo->getCharge() != 0) {
        assert(not pfo->getTracks().empty());

        lcfiplus::Track* track = new lcfiplus::Track;

        track->setId(trkIdCounter++); // start from 0. 110927 suehara
        track->setMcp(mcpf);

        track->setCharge(pfo->getCharge());

        track->SetE(pfo->getEnergy());
        TVector3 pTrack(pfo->getMomentum());
        track->SetVect(pTrack);

	//PIDs
	try{
	  int pidAlgoID = PID.getAlgorithmID(_pidAlgoName);
	  //pdg value
	  track->setPDG(PID.getParticleID(pfo,pidAlgoID).getPDG());
	  //posterior probabilities for each particle type hypothesis
	  int vecsize = PID.getParticleID(pfo,pidAlgoID).getParameters().size();
	  for(int i=0;i<vecsize;i++) 
	    track->setParticleIDProbability(PID.getParameterNames(pidAlgoID)[i],
					    (double)PID.getParticleID(pfo,pidAlgoID).getParameters()[i]);

	  //cal. corrected mass
	  track->setCorrEnergy(pmass[PID.getParticleID(pfo,pidAlgoID).getPDG()]);
	  //track->swapEnergy();  //really temporal need flag...
	}catch(UTIL::UnknownAlgorithm& e){
	}
	
	//tempolary
        // vector<ParticleID*> idvec = pfo->getParticleIDs();
        // if (idvec.size()>0) {
	//   track->setPDG(idvec[3]->getPDG());
        // }

        double pfoMom = pTrack.Mag();

        // use first track by default, but see immediately below...
        lcio::Track* trk = pfo->getTracks()[0];
        double delta(1e10);

        // now try to find if there's a better match (we use momentum),
        // because PandoraPFA can add a parent particle to the tracks container
        // after the first one for kinks and prongs, etc

        const double c = 2.99792458e8; // m*s^-1
        const double B = Globals::Instance()->getBField(); // Tesla
        const double mm2m = 1e-3;
        const double eV2GeV = 1e-9;
        const double eB = B*c*mm2m*eV2GeV;

        for (unsigned int i=0; i<pfo->getTracks().size(); ++i) {
          lcio::Track* testTrk = pfo->getTracks()[i];
          double om = testTrk->getOmega();
          double td = testTrk->getTanLambda();
          double cd = 1./sqrt(1+td*td);
          double pT = eB/fabs(om);
          double p = pT/cd;
          double testDelta = fabs(p-pfoMom);
          if (testDelta < delta) {
            delta = testDelta;
            trk = testTrk;
          }
        }

        // check against NaN's
        assert( trk->getD0() == trk->getD0() );
        assert( trk->getZ0() == trk->getZ0() );
        assert( trk->getPhi() == trk->getPhi() );
        assert( trk->getOmega() == trk->getOmega() );
        assert( trk->getTanLambda() == trk->getTanLambda() );

        // fill track parameters
        double par[lcfiplus::tpar::parN];
        par[lcfiplus::tpar::d0] = trk->getD0();
        par[lcfiplus::tpar::z0] = trk->getZ0();
        par[lcfiplus::tpar::ph] = trk->getPhi();
        par[lcfiplus::tpar::om] = trk->getOmega();
        par[lcfiplus::tpar::td] = trk->getTanLambda();
        track->setHelix(par);

        // ... and the covariance matrix
        const vector<float>& cov = trk->getCovMatrix();
        double fcov[lcfiplus::tpar::covN];
        for (int i=0; i<lcfiplus::tpar::covN; ++i) {
          fcov[i] = cov[i];
        }
        track->setCovMatrix(fcov);

        track->setChi2(trk->getChi2());
        track->setNdf(trk->getNdf());

        // store detector hit numbers
        int nhits[lcfiplus::tpar::hitN];
        if ( _trackHitOrdering == 0 ) {
          // ILD-LOI format, FTD and SIT switched w.r.t. ILD-DBD format
          // used also by SiD-DBD
          for (int i=0; i<lcfiplus::tpar::hitN; ++i) {
            nhits[i] = trk->getSubdetectorHitNumbers()[i];
          }
        } else if (_trackHitOrdering == 2) {
          // CLICdet format
          const vector<int>& vec = trk->getSubdetectorHitNumbers();
          int offset = 2; // 2=fit, 1=patrec
          nhits[tpar::VTX] = vec[2 * 1 - offset]; // vtx barrel
          nhits[tpar::FTD] = vec[2 * 2 - offset]; // vtx endcap
          nhits[tpar::SIT] = 0;
          // inner barrel and endcap, and outer barrel and endcap
          nhits[tpar::TPC] = vec[2 * 3 - offset] + vec[2 * 4 - offset] + vec[2 * 5 - offset] + vec[2 * 6 - offset];
          nhits[tpar::SET] = 0;
          nhits[tpar::ETD] = 0;
        } else {
          // ILD DBD format
          const vector<int>& vec = trk->getSubdetectorHitNumbers();
          int offset = 2; // 2=fit, 1=patrec
          nhits[tpar::VTX] = vec[2 * lcio::ILDDetID::VXD -offset];
          nhits[tpar::SIT] = vec[2 * lcio::ILDDetID::SIT -offset];
          nhits[tpar::FTD] = vec[2 * lcio::ILDDetID::FTD -offset];
          nhits[tpar::TPC] = vec[2 * lcio::ILDDetID::TPC -offset];
          nhits[tpar::SET] = vec[2 * lcio::ILDDetID::SET -offset];
          nhits[tpar::ETD] = vec[2 * lcio::ILDDetID::ETD -offset];
        }

        track->setTrackHits(nhits);
        track->setRadiusOfInnermostHit(trk->getRadiusOfInnermostHit());
        //printf("rimh = %f\n",trkData.rimh);

        track->setCaloEdep(subE);

        // register!
        itPfoCol->second.first->push_back(track);
        _trackLCIORel[track] = pfo;
        _trackLCIORel2[pfo] = track;
      } else {

        // convert neutrals ////////////////////////////////////////////////////////////////////////////////

        // assert(clusSize>0);
        // above is not necessarily true e.g. for a V0 which just has two tracks

        lcfiplus::Neutral* neut = new Neutral;
        neut->setId(neutIdCounter++); // start from 0. 110927 suehara
        neut->setMcp(mcpf);
        neut->setPDG(pfo->getType());
        /*
        vector<ParticleID*> idvec = pfo->getParticleIDs();
        if (idvec.size()>0) {
        	neut->setPDG(idvec[0]->getPDG());
        }
        */
        neut->SetE(pfo->getEnergy());
        neut->SetVect(TVector3(pfo->getMomentum()));

        // store calo
        neut->setCaloEdep(subE);

        itPfoCol->second.second->push_back(neut);
        _neutralLCIORel[neut] = pfo;
        _neutralLCIORel2[pfo] = neut;
	
	//set clustervec
	neut->setClusters(clusters);

        if (clusters.size()==0) { // v0
          int abspdg = abs(pfo->getType());

          if (abspdg == 310 || abspdg == 3122 || abspdg == 22) {
            neut->setV0();
            /*
            vector<lcio::Track*> trks = pfo->getTracks();
            assert(trks.size() == 2);
            MCParticle* v0dau = _mcpLCIORel2[trks[0]];
            MCParticle* par = v0dau->getParent();

            printf("ntrk %d, mcpf %d, mcp %d, par %d\n",
            		trks.size(),
            		mcpf ? mcpf->getPDG() : 666,
            		mcp ? mcp->getPDG() : 666,
            		par ? par->getPDG() : 666
            		);
             */
          } else {
            //printf("neutral pandora=%d, pfoid=%d\n",pfo->getType(),idvec[0]->getPDG());
            printf("neutral pandora=%d\n",pfo->getType());
          }
        }
      }

    }

    // clear navigators
    for (unsigned int inav=0; inav < navs.size(); inav ++) {
      delete navs[inav];
    }

  } // PFO end

  // vertex
  map<string, vector<const lcfiplus::Vertex*> *>::iterator itv;
  for (itv = _importVertexCols.begin(); itv != _importVertexCols.end(); itv++) {
    ReadVertices(itv->first.c_str(), itv->second);
  }

  // jet
  map<string, vector<const lcfiplus::Jet*> *>::iterator itj;
  for (itj = _importJetCols.begin(); itj != _importJetCols.end(); itj++) {
    ReadJets(itj->first.c_str(), itj->second);
  }

}

void LCIOStorer::ReadVertices(const char* vtxname, vector<const Vertex*>* lcficol) {
  lcio::LCCollection* colvtx;
  try {
    colvtx = _event->getCollection(vtxname);
  } catch (lcio::DataNotAvailableException& e) {
    cout << "LCIOStorer::ReadVertices(): vertex collection " << vtxname << " not found. Skip conversion." << endl;
    cout << e.what() << endl;
    return;
  }

  for (int n=0; n<colvtx->getNumberOfElements(); n++) {
    lcio::Vertex* lciovtx = dynamic_cast<lcio::Vertex*>(colvtx->getElementAt(n));

    double cov[6];
    for (int i=0; i<6; i++)cov[i] = lciovtx->getCovMatrix()[i];

    lcfiplus::Vertex* flavtx;
    if (_vtxLCIORel2[lciovtx] == 0) {
      flavtx = new lcfiplus::Vertex(lciovtx->getChi2(), lciovtx->getProbability(),
                                    lciovtx->getPosition()[0], lciovtx->getPosition()[1], lciovtx->getPosition()[2], cov, lciovtx->isPrimary());
      flavtx->setId(n);
      flavtx->setVertexingName(lciovtx->getAlgorithmType());

      lcio::ReconstructedParticle* lciorp = lciovtx->getAssociatedParticle();
      if (lciorp) {
        for (unsigned int ntr = 0; ntr < lciorp->getParticles().size(); ntr++) {
          lcio::ReconstructedParticle* lciotr = lciorp->getParticles()[ntr];
          if (lciotr->getParticles().size() == 1)lciotr = lciotr->getParticles()[0];
//						cerr << lciotr->getEnergy() << ", " << (unsigned int)lciotr << endl;

          // searching pfo
          if (_trackLCIORel2.count(lciotr)) {
            flavtx->add(_trackLCIORel2[lciotr]);
          } else {
            cerr << "LCIOStorer::SetEvent: Track associated to LCIO vertex is invalid!" << endl;
            continue;
          }
        }
      } else {
        if (_ignoreLackOfVertexRP == false) {
          cout << "LCIOStorer::SetEvent: no associated RP collection found for vertex collection '" << vtxname << "'. ";
          cout << "Set IgnoreLackOfVertexRP to 1 to proceed." << endl;
          throw (Exception("LCIOStorer::SetEvent: vertex RP collection not found"));
        }
      }

      _vtxLCIORel[flavtx] = lciovtx;
      _vtxLCIORel2[lciovtx] = flavtx;
      lcficol->push_back(flavtx);
    } else
      lcficol->push_back(const_cast<lcfiplus::Vertex*>(_vtxLCIORel2[lciovtx]));

  }
}

void LCIOStorer::ReadJets(const char* jetname, vector<const Jet*>* lcficol, const char* vtxrelname) {
//		cout << "Loading " << jetname << endl;
  lcio::LCCollection* coljet;

  try {
    coljet = _event->getCollection(jetname);
  } catch (lcio::DataNotAvailableException& e) {
    cout << "LCIOStorer::SetEvent(): jet collection " << jetname << "not found. Skip conversion." << endl;
    cout << e.what() << endl;
    return;
  }

  // initializing parameters
  lcio::PIDHandler pidh(coljet);

  // looking for vertices to be added

  LCRelationNavigator* nav = 0;
  try {
    string rname = (vtxrelname ? vtxrelname : string(jetname) + "_rel");
    lcio::LCCollection* colrel = _event->getCollection(rname);
    nav = new LCRelationNavigator(colrel);
  } catch (lcio::DataNotAvailableException& e) {
    // just ignore relation
  }

  for (int n=0; n<coljet->getNumberOfElements(); n++) {
    lcio::ReconstructedParticle* lciojet = dynamic_cast<lcio::ReconstructedParticle*>(coljet->getElementAt(n));

    lcfiplus::Jet* flajet;

    if (_jetLCIORel2[lciojet] == 0) {
      flajet = new lcfiplus::Jet;
      flajet->setId(n);

      for (unsigned int npart = 0; npart < lciojet->getParticles().size(); npart++) {
        lcio::ReconstructedParticle* lciorp = lciojet->getParticles()[npart];
        while (lciorp->getParticles().size()>0) {
          cout << "LCIO RP is nested. #P = " << lciorp->getParticles().size() << endl;
          lciorp = lciorp->getParticles()[0];
        }

        if (lciorp->getCharge() != 0) { //tracks
          // searching pfo
          if (_trackLCIORel2.count(lciorp)) {
            flajet->add(_trackLCIORel2[lciorp]);
          }
        } else { //neutrals
          // searching pfo
          if (_neutralLCIORel2.count(lciorp)) {
            flajet->add(_neutralLCIORel2[lciorp]);
          }
        }
      }

      // momentum/energy should be set at last since add() modifies the p/E
      flajet->SetPxPyPzE(lciojet->getMomentum()[0], lciojet->getMomentum()[1], lciojet->getMomentum()[2], lciojet->getEnergy());

      // vertex attachment
      if (nav) {
        const lcio::LCObjectVec& vtxlist = nav->getRelatedToObjects(lciojet);
        for (unsigned int n=0; n<vtxlist.size(); n++) {
          const lcfiplus::Vertex* vtx = _vtxLCIORel2[dynamic_cast<lcio::Vertex*>(vtxlist[n])];
          if (!vtx)throw("LCIOStorer::ReadJets: vertex object related to a jet is not found: please include vertex collection first.");
          flajet->add(vtx);
        }
      }

      // reading pids
      for (unsigned int na = 0; na < pidh.getAlgorithmIDs().size(); na++) {
        int algoid = pidh.getAlgorithmIDs()[na];
        string algoname = pidh.getAlgorithmName(algoid);
        const StringVec& paramnames = pidh.getParameterNames(algoid);

        const lcio::ParticleIDVec& pidv = lciojet->getParticleIDs();
        const lcio::ParticleID* pid = 0;
        for (unsigned int nv=0; nv<pidv.size(); nv++) {
          if (pidv[nv]->getAlgorithmType() == algoid) {
            pid = pidv[nv];
            break;
          }
        }
        if (pid) {
          int npm = pid->getParameters().size();

          lcfiplus::Parameters params;
          for (int np = 0; np < npm; np++) {
            params.add(paramnames[np].c_str(), (double)pid->getParameters()[np]);
          }

          flajet->addParam(algoname.c_str(), params);
        }
      }

      _jetLCIORel[flajet] = lciojet;
      _jetLCIORel2[lciojet] = flajet;
      lcficol->push_back(flajet);

    } else {
      lcficol->push_back(const_cast<lcfiplus::Jet*>(_jetLCIORel2[lciojet]));
    }
  }
  if (nav)delete nav;

}

void LCIOStorer::WriteVertices(const char* vertexName, const char* newName, const char* newRPName) {
  if (!_event) {
    cerr << "LCIOStorer::ConvertVertex: LCIO event has not been initialized." << endl;
    return ;
  }
  VertexVec* pvvtx;
  Event::Instance()->Get(vertexName, pvvtx);

  if (!pvvtx) {
    cerr << "LCIOStorer::ConvertVertex: failed to obtain vertices from the specified name: " << vertexName << endl;
    return ;
  }

  WriteVertices(pvvtx, (newName ? newName : vertexName), newRPName);
}

void LCIOStorer::WriteVertices(VertexVec* pvvtx, const char* newName, const char* newRPName) {
  // make collection
  lcio::LCCollectionVec* col = new lcio::LCCollectionVec(lcio::LCIO::VERTEX);
  lcio::LCCollectionVec* colRP = new lcio::LCCollectionVec(lcio::LCIO::RECONSTRUCTEDPARTICLE);
  for (unsigned int n=0; n<pvvtx->size(); n++) {
    // set ID to the lcfiplus::Vertex
    const lcfiplus::Vertex* flavtx = (*pvvtx)[n];
    if( not flavtx ) continue;
    flavtx->setId(n);

    // It seems that multiple entries in LCIO collection are not allowed: currently make independent entries
    /*
    			if(_vtxLCIORel[flavtx] != 0){
    				col->addElement(_vtxLCIORel[flavtx]);
    				continue;
    			}
    */
    // make new vertex
    lcio::VertexImpl* lciovtx = new lcio::VertexImpl;
    // make new recoparticle including vertex
    lcio::ReconstructedParticleImpl* lciorp = new lcio::ReconstructedParticleImpl;

    // associate particles
    TLorentzVector lv;
    double charge = 0.;
    for (unsigned int ntr = 0; ntr < flavtx->getTracks().size(); ntr++) {
      const lcfiplus::Track* flatr = flavtx->getTracks()[ntr];
      lcio::ReconstructedParticle* lciotr = _trackLCIORel[const_cast<lcfiplus::Track*>(flatr)];
      lv += (*flatr);
      charge += flatr->getCharge();
      lciorp->addParticle(lciotr);
      if (_updateVertexRPDaughters) {
        lcio::ReconstructedParticleImpl* lciotrimpl = dynamic_cast<lcio::ReconstructedParticleImpl*>(lciotr);
        if (lciotrimpl)
          lciotrimpl->setStartVertex(lciovtx);
        else
          throw (Exception(
                   "LCIOStorer: writing back associated vertices to the RP collection is failed."
                   "Please check if the input RP collection is writable."
                   "To switch off the writeback feature in the Marlin processor, specify set \"UpdateVertexRPDaughters\" to false."));
      }
    }
    double mom[3] = {lv.Px(), lv.Py(), lv.Pz()};
    float vpos[3] = {
      static_cast<float>(flavtx->getX()),
      static_cast<float>(flavtx->getY()),
      static_cast<float>(flavtx->getZ())
    };
    lciorp->setType(3); // 0: unknown 1: single 2:v0 3: compound 4:jet
    lciorp->setMomentum(mom);
    lciorp->setEnergy(lv.E());
    lciorp->setMass(lv.M());
    lciorp->setCharge(charge);
    lciorp->setReferencePoint(vpos);

    // ignore covmatrix of rp
    colRP->addElement(lciorp);

    // vertex initialization
    lciovtx->setPrimary(flavtx->isPrimary());
    if(!flavtx->isPrimary()){  //if secondary      
      if(flavtx->getVertexingName() == "AVF") lciovtx->setAlgorithmType("lcfiplusAVF");
      else lciovtx->setAlgorithmType("lcfiplus");
    }
    lciovtx->setPosition(vpos);
    // lcio vertex needs floats, not doubles
    // using {} to limit scope of temporary
    {
      float destCov[6];
      const double* sourceCov = flavtx->getCov();
      for (int i=0; i<6; ++i) {
        destCov[i] = static_cast<float>(sourceCov[i]);
      }
      lciovtx->setCovMatrix(destCov);
    }
    lciovtx->setChi2(flavtx->getChi2());
    lciovtx->setProbability(flavtx->getProb());
    lciovtx->setAssociatedParticle(lciorp);

    /*
    cout << scientific << "ConvertVertex position: "
    	<< vpos[0] << ","
    	<< vpos[1] << ","
    	<< vpos[2] << fixed << endl;

    const double* cov = flavtx->getCov();
    cout << scientific << "ConvertVertex covMatrix: "
    	<< cov[0] << ","
    	<< cov[1] << ","
    	<< cov[2] << ","
    	<< cov[3] << ","
    	<< cov[4] << ","
    	<< cov[5] << fixed << endl;

    const vector<double>& cov2 = lciovtx->getCovMatrix();
    cout << scientific << "ConvertVertex covMatrix (lcio): "
    	<< cov2[0] << ","
    	<< cov2[1] << ","
    	<< cov2[2] << ","
    	<< cov2[3] << ","
    	<< cov2[4] << ","
    	<< cov2[5] << fixed << endl;
     */

    // store relation
    _vtxLCIORel[flavtx] = lciovtx;
    _vtxLCIORel2[lciovtx] = flavtx;

    // register to collection
    col->addElement(lciovtx);
  }
  _event->addCollection(col, newName);
  _event->addCollection(colRP, newRPName ? newRPName : (const char*)TString::Format("%s_RP", newName));

  //cout << "ConvertVertex finished. # vertices = " << col->getNumberOfElements() << endl;

}

void LCIOStorer::WriteJets(const char* jetName, const char* newName, bool writeVertex, const char* vtxName, const char* relName) {
  if (!_event) {
    cerr << "LCIOStorer::ConvertJet: LCIO event has not been initialized." << endl;
    return ;
  }
  JetVec* pvjet;
  Event::Instance()->Get(jetName, pvjet);

  if (!pvjet) {
    cerr << "LCIOStorer::ConvertJet: failed to obtain jets from the specified name: " << jetName << endl;
    return ;
  }

//		_jetLCIORel.clear();
//		_jetLCIORel.resize(pvjet->size()); // map has no resize()

  // make collection
  lcio::LCCollectionVec* col = new lcio::LCCollectionVec(lcio::LCIO::RECONSTRUCTEDPARTICLE);
  lcio::LCRelationNavigator rel(lcio::LCIO::RECONSTRUCTEDPARTICLE, lcio::LCIO::VERTEX);

  // naming
  if (!newName)newName = jetName;
  string vn = string(newName) + "_vtx";
  if (!vtxName)vtxName = vn.c_str();
  string rn = string(newName) + "_rel";
  if (!relName)relName = rn.c_str();

  //TODO: check existence of collections

  // first reserve vertex collection
  if (writeVertex) {
    vector<const Vertex*> vv;
    for (unsigned int n=0; n<pvjet->size(); n++) {
      vv.insert(vv.end(), (*pvjet)[n]->getVertices().begin(),(*pvjet)[n]->getVertices().end());
    }
    WriteVertices(&vv, vtxName);
  }

  for (unsigned int n=0; n<pvjet->size(); n++) {
    // set ID to the lcfiplus::Vertex
    const lcfiplus::Jet* flajet = (*pvjet)[n];
    flajet->setId(n);

    lcio::ReconstructedParticle* lciojet2;

//			if(_jetLCIORel[flajet] == 0){ // new lciojet creation needed
    if (1) { // multiple entries are not allowed in LCIO: making duplicate objects
      // make new recoparticle including vertex
      lcio::ReconstructedParticleImpl* lciojet = new lcio::ReconstructedParticleImpl;

      // associate particles
      double charge = 0.;
      for (unsigned int ntr = 0; ntr < flajet->getTracks().size(); ntr++) {
        const lcfiplus::Track* flatr = flajet->getTracks()[ntr];
        lcio::ReconstructedParticle* lciotr = _trackLCIORel[const_cast<lcfiplus::Track*>(flatr)];
        charge += flatr->getCharge();
        lciojet->addParticle(lciotr);
        //cout << "LCIOStorer::ConvertJet: add track: id = " << flatr->getId() << ", energy = " << flatr->E() << flush;
        //cout << ", lcio energy = " << lciotr->getEnergy() << endl;
      }
      for (unsigned int nneut = 0; nneut < flajet->getNeutrals().size(); nneut++) {
        const lcfiplus::Neutral* flaneut = flajet->getNeutrals()[nneut];
        lcio::ReconstructedParticle* lcioneut = _neutralLCIORel[const_cast<lcfiplus::Neutral*>(flaneut)];
        lciojet->addParticle(lcioneut);
      }
      for (unsigned int nvtx = 0; nvtx < flajet->getVertices().size(); nvtx++) {
        const lcfiplus::Vertex* flavtx = flajet->getVertices()[nvtx];
        // first, extract all vertex tracks
        for (unsigned int ntr = 0; ntr < flavtx->getTracks().size(); ntr++) {
          const lcfiplus::Track* flatr = flavtx->getTracks()[ntr];
          lcio::ReconstructedParticle* lciotr = _trackLCIORel[const_cast<lcfiplus::Track*>(flatr)];
          charge += flatr->getCharge();
          lciojet->addParticle(lciotr);
        }
      }

      const TLorentzVector& lv = *flajet;
      double mom[3] = {lv.Px(), lv.Py(), lv.Pz()};
      lciojet->setType(4); // 0: unknown 1: single 2:v0 3: compound 4:jet
      lciojet->setMomentum(mom);
      lciojet->setEnergy(lv.E());
      lciojet->setMass(lv.M());
      lciojet->setCharge(charge);
      // ignore covmatrix of rp

      // add PID stuffs
      WriteAllPIDs(col, lciojet, flajet);

      // store relation
      _jetLCIORel[flajet] = lciojet;

      lciojet2 = lciojet;

    } else { // lciojet already exists
      cout << "LCIOStorer::WriteJets: jet is already registered in LCIO.";
      if (writeVertex) cout << " Just trying to create vertex relation.";
      cout << endl;

      lciojet2 = _jetLCIORel[flajet];
    }

    for (unsigned int nvtx = 0; nvtx < flajet->getVertices().size(); nvtx++) {
      const lcfiplus::Vertex* flavtx = flajet->getVertices()[nvtx];
      if (writeVertex) {
        // make relation jet - vtx
        lcio::Vertex* lciovtx = _vtxLCIORel[flavtx];
        rel.addRelation(lciojet2, lciovtx);
      }
    }

    // register to collection
    col->addElement(lciojet2);
  }
  _event->addCollection(col, newName);
  if (writeVertex)
    _event->addCollection(rel.createLCCollection(), relName);

  //cout << "ConvertJet finished. # jets = " << col->getNumberOfElements() << endl;
}
}

void LCIOStorer::WritePID(lcio::LCCollection* lciocol, lcio::ReconstructedParticle* lciojet, const lcfiplus::Jet* lcfijet, const char* paramname) {
  // init pid handler
  lcio::PIDHandler pidh(lciocol);

  // obtain Parameters
  const Parameters* lcfiparams = lcfijet->getParam(paramname);
  const map<string, pair<const type_info*, void*> >& paramMap = lcfiparams->paramMap();

  vector<string> outParamNames;
  vector<float> outParamData;

  map<string, pair<const type_info*, void*> >::const_iterator it;
  for (it = paramMap.begin(); it != paramMap.end(); it++) {
    string name = it->first;

    double d = lcfiparams->get(name.c_str(),(double)0.);

    outParamNames.push_back(name);
    outParamData.push_back(d);
  }

  // register parameters
  int algoID;
  try {
    algoID = pidh.getAlgorithmID(string(paramname));
  } catch (lcio::UnknownAlgorithm& e) {
    algoID = pidh.addAlgorithm(paramname, outParamNames);
  }

  //                         type, PDG, likelihood: as same as old LCFI
  pidh.setParticleID(lciojet, 42, 9999, 0.0, algoID, outParamData); 
}


void LCIOStorer::WriteAllPIDs(lcio::LCCollection* lciocol, lcio::ReconstructedParticle* lciojet, const lcfiplus::Jet* lcfijet) {
  const map<string, Parameters>& parammap = lcfijet->params();

  map<string, Parameters>::const_iterator it;
  for (it = parammap.begin(); it != parammap.end(); it++) {
    WritePID(lciocol, lciojet, lcfijet, it->first.c_str());
  }

  /*
  if(parammap.size())
  	cout << parammap.size() << " PID parameter sets written to jet." << endl;
  	//*/
}

void LCIOStorer::SetColorSinglets(vector<MCParticle*>& mcps, vector<MCColorSinglet*>& mccs) {
  int cspdg[] = {22, 23, 24, 25, 91, 92, 93, 94, -1}; // -1 is added to avoid warning
  int cspdg2[] = {22, 91, 92, 93, 94, -1}; // color singlet candidates which should check parents

  vector<MCParticle*> initialquarks;

  // looking for color singlets
  for (unsigned int n=0; n<mcps.size(); n++) {
    MCParticle* mcp = mcps[n];

    // ignore isolated particles
    if (mcp->getParent() == 0 && mcp->getDaughters().size() == 0)continue;

    int apdg = abs(mcp->getPDG());

    if (apdg <= 6 && ! mcp->getParent() && mcp->getDaughters().size() > 0)
      initialquarks.push_back(mcp);

    if (find(cspdg, cspdg + 8, apdg) != cspdg + 8) {
      // looking for preregistered parents
      bool parent = false;
      if (find(cspdg2, cspdg2 + 5, apdg) != cspdg2 + 5) {
        for (unsigned int ncs = 0; ncs < mccs.size(); ncs ++) {
          MCColorSinglet* pcs = mccs[ncs];
          if (mcp->isParent(pcs->_mcp)) {
            parent = true;
            break;
          }
        }

      }

      if (!parent) {
        // check daughters
        if (mcp->getDaughters().size() > 1) {
          bool qdau = false;
          // at least 2 daughters for color singlet
          for (unsigned int ndau = 0; ndau < mcp->getDaughters().size(); ndau++) {
            int dpdg = abs(mcp->getDaughters()[ndau]->getPDG());
            if (dpdg <= 6 || dpdg == 21) { // quark daughters
              qdau = true;
              break;
            }
          }
          if (qdau) {
            // all selections clear; add to mccs
            MCColorSinglet* newcs = new MCColorSinglet;
            newcs->_mcp = mcp;
            mccs.push_back(newcs);
          }
        }
      }

    }
  }

  // make contents of MCColorSinglet
  for (unsigned int n=0; n<mccs.size(); n++) {
    // obtain initial quarks
    MCColorSinglet* cs = mccs[n];
    const MCParticle* mcp = cs->_mcp;
    int apdg = abs(mcp->getPDG());

    // initials
    if (apdg > 90) { // intermediate particle
      // looking for initial particles to associate
      for (unsigned int q=0; q<initialquarks.size(); q++) {
        if (initialquarks[q]->getDaughters()[0] == mcp) {
          cs->_initials.push_back(initialquarks[q]);
        }
      }
    } else {
      for (unsigned int q=0; q<mcp->getDaughters().size(); q++) {
        if (abs(mcp->getDaughters()[q]->getPDG()) < 100)
          cs->_initials.push_back(mcp->getDaughters()[q]);
      }
    }
  }
  vector<const MCParticle*> finalcs;
  for (unsigned int n=0; n<mcps.size(); n++) {
    MCParticle* mcp = mcps[n];
    int apdg = abs(mcp->getPDG());

    // finalcolorsinglets & realparticles
    if (apdg > 90 && apdg<100 && mcp->getDaughters().size() > 0 && abs(mcp->getDaughters()[0]->getPDG()) > 100) {
      MCColorSinglet* cs = const_cast<MCColorSinglet*>(mcp->getColorSinglet(constVector(&mccs)));
      if (cs) {
        cs->_finalcolorsinglets.push_back(mcp);
        for (unsigned int i=0; i<mcp->getDaughters().size(); i++) {
          cs->_realparticles.push_back(mcp->getDaughters()[i]);
        }
        finalcs.push_back(mcp);
      }
    }
    // qqgluons
    if (apdg == 21 && mcp->getDaughters().size() > 1 && abs(mcp->getDaughters()[0]->getPDG())<=6) {
      MCColorSinglet* cs = const_cast<MCColorSinglet*>(mcp->getColorSinglet(constVector(&mccs)));
      if (cs)cs->_qqgluons.push_back(mcp);
    }
  }
  for (unsigned int n=0; n<mcps.size(); n++) {
    MCParticle* mcp = mcps[n];
    if (mcp->getDaughters().size()>0 && find(finalcs.begin(), finalcs.end(), mcp->getDaughters()[0]) != finalcs.end()) {
      MCColorSinglet* cs = const_cast<MCColorSinglet*>(mcp->getColorSinglet(constVector(&mccs)));
      if (cs)cs->_finalstrings.push_back(mcp);
    }
  }
  if (0) {
    for (unsigned int n=0; n<mccs.size(); n++) {
      MCColorSinglet* cs = mccs[n];
      const MCParticle* mcp = cs->_mcp;

      cout << "Color singlet found: pdg = " << mcp->getPDG() << ", E = " << mcp->E() << endl;
      cout << "# initial quarks = " << cs->_initials.size() << ", pdg = ";
      for (unsigned int q=0; q<cs->_initials.size(); q++)
        cout << cs->_initials[q]->getPDG() << " ";
      cout << endl;
      cout << "# qqgluons = " << cs->_qqgluons.size() << ", pdg = ";
      for (unsigned int q=0; q<cs->_qqgluons.size(); q++)
        cout << cs->_qqgluons[q]->getDaughters()[0]->getPDG() << " ";
      cout << endl;
      cout << "# finalstrings = " << cs->_finalstrings.size() << ", pdg = ";
      for (unsigned int q=0; q<cs->_finalstrings.size(); q++)
        cout << cs->_finalstrings[q]->getPDG() << " ";
      cout << endl;
      cout << "# finalcolorsinglets = " << cs->_finalcolorsinglets.size() << ", pdg = ";
      for (unsigned int q=0; q<cs->_finalcolorsinglets.size(); q++)
        cout << cs->_finalcolorsinglets[q]->getPDG() << " ";
      cout << endl;
      cout << "# realparticles = " << cs->_realparticles.size() << ", pdg = ";
      for (unsigned int q=0; q<cs->_realparticles.size(); q++)
        cout << cs->_realparticles[q]->getPDG() << " ";
      cout << endl;
    }
  }
}

void LCIOStorer::GetCallback(const char* name, const char* classname) {
  if (!_autoread)return;
  if (Event::Instance()->IsExist(name))return; // no-op since the collection has already been registered

  if (string("vector<lcfiplus::Vertex*>") == classname) {
    cout << "Loading LCIO vertex collection: " << name << endl;
    InitVertexCollection(name,name);
  } else if (string("vector<lcfiplus::Jet*>") == classname) {
    cout << "Loading LCIO jet collection: " << name << endl;
    InitJetCollection(name,name);
  }
}

