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

// lcfiplus includes
#include "LCIOStorer.h"
#include "lcfiplus.h"
#include "EventStore.h"
#include "JetFinder.h"

#include <assert.h>
#include <math.h>

//using namespace lcio;

namespace lcfiplus{

	// ctor/dtor
	LCIOStorer::LCIOStorer(const char *inputfile, const char *outputfile, bool autoconvert, const char *outPrefix)
	{
		if(inputfile){	
			_reader = lcio::LCFactory::getInstance()->createLCReader();
			_reader->open(inputfile); // may throw exception...
			cout << "Input LCIO file open: " << inputfile << " successful." << endl;
		}else
			_reader = NULL;

		if(outputfile){
			_writer = lcio::LCFactory::getInstance()->createLCWriter();
			_writer->open(outputfile);
			cout << "Output LCIO file open: " << outputfile << " successful." << endl;
		}else
			_writer = NULL;

		// copy run header
		if(_reader && _writer){
/*			int ncp = 0;
			lcio::LCRunHeader *header = 0;
			while((header = _reader->readNextRunHeader())){
				_writer->writeRunHeader(header);
				ncp ++;
			}
			cout << "LCIOStorer: copy " << ncp << " run headers." << endl;
*/
			lcio::LCRunHeader *header = _reader->readNextRunHeader();
			if(header){
				_writer->writeRunHeader(header);
				cout << "LCIOStorer: copied a run header." << endl;
			}else{
				cout << "LCIOStorer: the run header not found." << endl;
			}
		}

		_event = NULL;	

		_pTracks = NULL;
		_pNeutrals = NULL;
		_pMCPs = NULL;
		
		_autoconvert = autoconvert;
		_savePrefix = (outPrefix ? outPrefix : "");
		cout << "LCIOStorer: initialization finished." << endl;
	}

	LCIOStorer::~LCIOStorer()
	{
		if(_reader)_reader->close();
		if(_writer)_writer->close();
	}

	void LCIOStorer::InitCollections(const char *pfoColName, const char *mcColName, const char *mcpfoColName,
				const char *trackName, const char *neutralName, const char *mcpName)
	{
		/*
		cout << "LCIOStorer::InitCollections: registering collections" << endl;
		cout << "LCIOStorer::InitCollections:     pfo collection:   " << pfoColName << endl; 
		cout << "LCIOStorer::InitCollections:     MC collection:    " << mcColName << endl; 
		cout << "LCIOStorer::InitCollections:     mcpfo collection: " << mcpfoColName << endl; 
		cout << "LCIOStorer::InitCollections:     track name:       " << trackName << endl; 
		cout << "LCIOStorer::InitCollections:     neutral name:     " << neutralName << endl; 
		cout << "LCIOStorer::InitCollections:     mcp name:         " << mcpName << endl; 
		*/
		_pfoColName = pfoColName;
		_mcColName = mcColName;
		_mcpfoColName = mcpfoColName;

		Event *event = Event::Instance();
		event->Register<Track>(trackName,_pTracks);
		event->Register<Neutral>(neutralName,_pNeutrals);
		event->Register<MCParticle>(mcpName,_pMCPs);

		event->setDefaultTracks(trackName);
		event->setDefaultNeutrals(neutralName);
		event->setDefaultMCParticles(mcpName);

		if(!_pTracks || !_pNeutrals || !_pMCPs){
			cout << "[ERROR] LCIOStorer::InitCollections: failed to register LCIO collections in lcfiplus namespace." << endl;
			throw(Exception("LCIOStorer::InitCollections: failed to register LCIO collections in lcfiplus namespace."));
		}
	}
	
	void LCIOStorer::InitVertexCollection(const char *lcioName, const char *flavtagName)
	{
		vector<lcfiplus::Vertex *> *vtxcol = 0;
		Event::Instance()->Register<Vertex>(flavtagName,vtxcol);
		
		_importVertexCols[ lcioName ] = vtxcol;
	}
	
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

	void LCIOStorer::InitJetCollection(const char *lcioName, const char *flavtagName)
	{
		vector<lcfiplus::Jet *> *jetcol = 0;
		Event::Instance()->Register<Jet>(flavtagName,jetcol);
		
		_importJetCols[ lcioName ] = jetcol;
	}

	bool LCIOStorer::Next()
	{
		lcio::LCEvent *evt = _reader->readNextEvent();
		if(!evt)return false;
		SetEvent(evt);
		return true;
	}

	void LCIOStorer::AutoConvert()
	{
		EventStore *s = Event::Instance();
		const map<string, EventStore::StoredEntry> & emap = s->GetObjectMap();
		
		const vector<string> * lcionames = _event->getCollectionNames();

		map<string, EventStore::StoredEntry>::const_iterator it;
		for(it = emap.begin(); it != emap.end(); it++){
			if(it->second.flag & EventStore::PERSIST){
				// persist: to save
				string newname = _savePrefix + it->first;
				
				if(find(lcionames->begin(), lcionames->end(), newname) != lcionames->end())
				{
					cout << "AutoConvert skipped because it is already included in lcio." << endl;
					continue; // skip collections that is already existed in lcio event
				}
				
				if(it->second.classname == "vector<lcfiplus::Vertex*>"){
					ConvertVertex(it->first.c_str(), newname.c_str());
				}
				else if(it->second.classname == "vector<lcfiplus::Jet*>"){
					bool extractVertex = it->second.flag & EventStore::JET_EXTRACT_VERTEX;
					ConvertJet(it->first.c_str(), newname.c_str(), extractVertex);
				}else{
					cout << "LCIOStorer::WriteEvent(): unknown conversion for lcio: class " << it->second.classname << ", ignored." << endl;
				}
			}
		}
	}

	void LCIOStorer::WriteEvent()
	{
		if(_writer && _event){
			if(_autoconvert){
				AutoConvert();
			}
			_writer->writeEvent(_event);
		}
		else{
			cout << "LCIOStorer::WriteEvent(): event cannot be written." << endl;
		}
	}
	
	// helper functions for std::sort()
/*	bool LCIOStorer::energy_sort_trk(Track* a, Track* b) {
	  assert(a != 0);
	  assert(b != 0);
	  return (a->getOmega() < b->getOmega());
	}

	bool LCIOStorer::energy_sort_mc(MCParticle *a, MCParticle *b) {
	  assert(a != 0);
	  assert(b != 0);
	  // needs to be converted to float; otherwise
	  // energy_sort_mc(mcp,mcp) --- can return true...
	  float en1 = a->getEnergy();
	  float en2 = b->getEnergy();
	  return (en1 > en2);
	}
*/
	bool LCIOStorer::energy_sort_pfo(lcio::ReconstructedParticle *a, lcio::ReconstructedParticle *b) {
	  assert(a != 0);
	  assert(b != 0);
	  return (a->getEnergy() > b->getEnergy());
	}


	void LCIOStorer::SetEvent(lcio::LCEvent *evt)
	{
		// buffer check
//		if(_pTracks == NULL || _pNeutrals == NULL || _pMCPs == NULL)throw(new lcfiplus::Exception("LCIOStorer::SetEvent: Event buffer not set."));
		_event = evt;
		if(_pTracks == NULL || _pNeutrals == NULL || _pMCPs == NULL)return;

		// collections
		lcio::LCCollection* colMC = evt->getCollection(_mcColName);
		lcio::LCCollection* colPFO = evt->getCollection(_pfoColName);
		lcio::LCRelationNavigator* nav = new lcio::LCRelationNavigator(evt->getCollection(_mcpfoColName));
	
		// memory allocation
		Event::Instance()->ClearObjects();
		
		_pMCPs->reserve(colMC->getNumberOfElements());
		_pNeutrals->reserve(colPFO->getNumberOfElements());
		_pTracks->reserve(colPFO->getNumberOfElements());

		//cout << "PFO size: " << colPFO->getNumberOfElements() << endl;

		_mcpLCIORel.clear();
		_trackLCIORel.clear();
		_neutralLCIORel.clear();
		
		_mcpLCIORel.reserve(colMC->getNumberOfElements());
		_trackLCIORel.reserve(colPFO->getNumberOfElements());
		_neutralLCIORel.reserve(colPFO->getNumberOfElements());

		// Relation of MCParticle between LCIO and lcfiplus
		map<lcio::MCParticle *, lcfiplus::MCParticle *> mcpMap;
		
		int mcIdCounter(0);
		// convert MCParticle /////////////////////////////////////////////////////////////////////////
		for (mcIdCounter=0; mcIdCounter<colMC->getNumberOfElements(); ++mcIdCounter) {
			lcio::MCParticle* mcp = dynamic_cast<lcio::MCParticle*>( colMC->getElementAt(mcIdCounter) );
			assert(mcp != 0);
			//if (mcp->isBackscatter()) continue;
			
			// find parent id
			// TODO Really OK for single parent??
			lcfiplus::MCParticle *parent = 0;
			if (mcp->getParents().size()>0) {
				map<lcio::MCParticle*,lcfiplus::MCParticle *>::iterator iter = mcpMap.find(mcp->getParents()[0]);

				if ( iter == mcpMap.end() )	throw(new lcfiplus::Exception("parent not found in association map"));
			
				parent = iter->second;
			}

			// fix against weird 1e-308 numbers
			const double* vtmp = mcp->getVertex();
			double v[3] = { vtmp[0], vtmp[1], vtmp[2] };
			//printf("%d [%e,%e,%e]\n",mcp->getPDG(),v[0],v[1],v[2]);
			if (v[0]*v[0]+v[1]*v[1]+v[2]*v[2] < 1e-100) {
				v[0] = 0; v[1] = 0; v[2] = 0;
			}

			
			// convert into MCParticle
			lcfiplus::MCParticle *mcpNew = new MCParticle(mcIdCounter, mcp->getPDG(), parent, mcp->getCharge(),
					TLorentzVector(TVector3(mcp->getMomentum()),mcp->getEnergy()), TVector3(v));
			
			// add to MCP list
			_pMCPs->push_back(mcpNew);
			
			// MC - id relation
			mcpMap[mcp] = mcpNew;
			// LCIO relation
			_mcpLCIORel[mcIdCounter] = mcp;
		}
		
	
		// convert PFO ////////////////////////////////////////////////////////////////////////////////
		// sorted PFO list
		vector<lcio::ReconstructedParticle *> pfo_list;
		
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
			lcio::ReconstructedParticle *pfo = pfo_list[n];
			lcio::MCParticle *mcp = NULL;
			lcfiplus::MCParticle *mcpf = NULL;
			if(nav->getRelatedToObjects(pfo).size()){
				mcp = dynamic_cast<lcio::MCParticle *>(nav->getRelatedToObjects(pfo)[0]); // TODO [0] OK?
				mcpf = mcpMap[mcp];
			}
			
			// find clusters
			vector<lcio::Cluster*> clusters = pfo->getClusters();
			float clusEnergy(0);
			float subE[6];
			for (int i=0; i<6; ++i) subE[i]=0;
			for (unsigned int iclus=0; iclus<clusters.size(); ++iclus) {
				clusEnergy += clusters[iclus]->getEnergy();
				for (int i=0; i<6; ++i) {
					subE[i] = clusters[iclus]->getSubdetectorEnergies()[i];
				}
			}

//			cerr << (pfo->getCharge() ? "[Track]" : "[Neutral]") << pfo->getEnergy() << ", " << (unsigned int)pfo << endl;
			
			
			if (pfo->getCharge() != 0) {
				int trkSize = pfo->getTracks().size();
				assert(trkSize>0);

				lcfiplus::Track * track = new lcfiplus::Track;

				//				track->setId(++trkIdCounter); // start from 1...
				track->setId(trkIdCounter++); // start from 0. 110927 suehara
				track->setMcp(mcpf);
				track->setPDG(pfo->getType());
				/*
				vector<ParticleID*> idvec = pfo->getParticleIDs();
				if (idvec.size()>0) {
					track->setPDG(idvec[0]->getPDG());
				}
				*/
				track->setCharge(pfo->getCharge());

				track->SetE(pfo->getEnergy());
				TVector3 pTrack(pfo->getMomentum());
				track->SetVect(pTrack);
				
				double pfoMom = pTrack.Mag();

				// use first track by default, but see immediately below...
				lcio::Track* trk = pfo->getTracks()[0];
				double delta(1e10);

				// now try to find if there's a better match (we use momentum),
				// because PandoraPFA can add a parent particle to the tracks container
				// after the first one for kinks and prongs, etc
				
				const float c = 2.99792458e8; // m*s^-1
				const float B = 3.5;          // Tesla
				const float mm2m = 1e-3;
				const float eV2GeV = 1e-9;
				const float eB = B*c*mm2m*eV2GeV;

				for (unsigned int i=0; i<pfo->getTracks().size(); ++i) {
					lcio::Track* testTrk = pfo->getTracks()[i];
					float om = testTrk->getOmega();
					float td = testTrk->getTanLambda();
					float cd = 1./sqrt(1+td*td);
					float pT = eB/fabs(om);
					float p = pT/cd;
					float testDelta = fabs(p-pfoMom);
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
				float par[lcfiplus::tpar::parN];
				par[lcfiplus::tpar::d0] = trk->getD0();
				par[lcfiplus::tpar::z0] = trk->getZ0();
				par[lcfiplus::tpar::ph] = trk->getPhi();
				par[lcfiplus::tpar::om] = trk->getOmega();
				par[lcfiplus::tpar::td] = trk->getTanLambda();
				track->setHelix(par);
				
				// ... and the covariance matrix
				const vector<float>& cov = trk->getCovMatrix();
				float fcov[lcfiplus::tpar::covN];
				for (int i=0; i<lcfiplus::tpar::covN; ++i) {
					fcov[i] = cov[i];
				}
				track->setCovMatrix(fcov);

				track->setChi2(trk->getChi2());
				track->setNdf(trk->getNdf());

				// store detector hit numbers
				float nhits[lcfiplus::tpar::hitN];
				for (int i=0; i<lcfiplus::tpar::hitN; ++i) {
					nhits[i] = trk->getSubdetectorHitNumbers()[i];
				}
				track->setTrackHits(nhits);
				track->setRadiusOfInnermostHit(trk->getRadiusOfInnermostHit());
				//printf("rimh = %f\n",trkData.rimh);

				track->setCaloEdep(subE);
				
				// register!
				_pTracks->push_back(track);
				_trackLCIORel.push_back(pfo);
//				_trackLCIORel[trkIdCounter] = pfo;
			} else { // neutrals
				// assert(clusSize>0);
				// above is not necessarily true e.g. for a V0 which just has two tracks

				lcfiplus::Neutral *neut = new Neutral;
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
				
				_pNeutrals->push_back(neut);
				_neutralLCIORel.push_back(pfo);
//				_neutralLCIORel[neutIdCounter] = pfo;

				if (clusters.size()==0) { // v0
					int abspdg = abs(pfo->getType());

					if (abspdg == 310 || abspdg == 3122 || abspdg == 22) {
						neut->setV0();
						/*
						vector<lcio::Track*> trks = pfo->getTracks();
						assert(trks.size() == 2);
						MCParticle* v0dau = mcpMap[trks[0]];
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
		// vertex
		map<string, vector<lcfiplus::Vertex *> *>::iterator it;
		for(it = _importVertexCols.begin(); it != _importVertexCols.end(); it++){
			it->second->clear();
			lcio::LCCollection* colvtx;
			try{
				colvtx = evt->getCollection(it->first);
			}catch(lcio::DataNotAvailableException &e){
				cout << "LCIOStorer::SetEvent(): vertex collection not found. Skip conversion." << endl;
				cout << e.what() << endl;
				continue; 
			}		

			for(int n=0;n<colvtx->getNumberOfElements();n++){
				lcio::Vertex *lciovtx = dynamic_cast<lcio::Vertex*>(colvtx->getElementAt(n));
				
				float cov[6];
				for(int i=0;i<6;i++)cov[i] = lciovtx->getCovMatrix()[i];
				
				lcfiplus::Vertex *flavtx = new lcfiplus::Vertex(lciovtx->getChi2(), lciovtx->getProbability(), 
					lciovtx->getPosition()[0], lciovtx->getPosition()[1], lciovtx->getPosition()[2], cov);
				flavtx->setId(n);
				
				lcio::ReconstructedParticle *lciorp = lciovtx->getAssociatedParticle();
				if(lciorp){
					for(unsigned int ntr = 0; ntr < lciorp->getParticles().size(); ntr++){
						lcio::ReconstructedParticle *lciotr = lciorp->getParticles()[ntr];
						if(lciotr->getParticles().size() == 1)lciotr = lciotr->getParticles()[0];
//						cerr << lciotr->getEnergy() << ", " << (unsigned int)lciotr << endl;
						// searching pfo : todo: non-efficient!
						vector<lcio::ReconstructedParticle *>::iterator it2 = find(_trackLCIORel.begin(), _trackLCIORel.end(), lciotr);
						if(it2 == _trackLCIORel.end()){
							cerr << "LCIOStorer::SetEvent: Track associated to LCIO vertex is invalid!" << endl;
							continue;
						}
						flavtx->add((*_pTracks)[it2 - _trackLCIORel.begin()]);
					}
				}
				it->second->push_back(flavtx);
			}
		}
		
		// jet
		map<string, vector<lcfiplus::Jet *> *>::iterator itj;
		for(itj = _importJetCols.begin(); itj != _importJetCols.end(); itj++){
			itj->second->clear();
			lcio::LCCollection* coljet;
			try{
				coljet = evt->getCollection(itj->first);
			}catch(lcio::DataNotAvailableException &e){
				cout << "LCIOStorer::SetEvent(): jet collection not found. Skip conversion." << endl;
				cout << e.what() << endl;
				continue; 
			}		

			for(int n=0;n<coljet->getNumberOfElements();n++){
				lcio::ReconstructedParticle *lciojet = dynamic_cast<lcio::ReconstructedParticle*>(coljet->getElementAt(n));
				
				lcfiplus::Jet *flajet = new lcfiplus::Jet;
				flajet->setId(n);
				
				for(unsigned int npart = 0; npart < lciojet->getParticles().size(); npart++){
					lcio::ReconstructedParticle *lciorp = lciojet->getParticles()[npart];
					
					if(lciorp->getCharge() != 0){ //tracks
						// searching pfo : todo: non-efficient!
						vector<lcio::ReconstructedParticle *>::iterator it2 = find(_trackLCIORel.begin(), _trackLCIORel.end(), lciorp);
						if(it2 == _trackLCIORel.end()){
							cerr << "LCIOStorer::SetEvent: Track associated to LCIO jet is invalid!" << endl;
							continue;
						}
						flajet->add((*_pTracks)[it2 - _trackLCIORel.begin()]);
					}
					else{ //neutrals
						// searching pfo : todo: non-efficient!
						vector<lcio::ReconstructedParticle *>::iterator it2 = find(_neutralLCIORel.begin(), _neutralLCIORel.end(), lciorp);
						if(it2 == _neutralLCIORel.end()){
							cerr << "LCIOStorer::SetEvent: Neutral associated to LCIO jet is invalid!" << endl;
							continue;
						}
						flajet->add((*_pNeutrals)[it2 - _neutralLCIORel.begin()]);
					}
				}
				// momentum/energy should be set at last since add() modifies the p/E
				flajet->SetPxPyPzE(lciojet->getMomentum()[0], lciojet->getMomentum()[1], lciojet->getMomentum()[2], lciojet->getEnergy());
				itj->second->push_back(flajet);
			}
		}
	}

	void LCIOStorer::ConvertVertex(const char *vertexName, const char *newName, const char *newRPName)
	{
		if(!_event){
			cerr << "LCIOStorer::ConvertVertex: LCIO event has not been initialized." << endl;
			return ;
		}
		VertexVec *pvvtx;
		Event::Instance()->Get(vertexName, pvvtx);
	
		if(!pvvtx){
			cerr << "LCIOStorer::ConvertVertex: failed to obtain vertices from the specified name: " << vertexName << endl;
			return ;
		}
		
		_vtxLCIORel.clear();
//		_vtxLCIORel.resize(pvvtx->size()); // map has no resize()
		
		// make collection
		lcio::LCCollectionVec *col = new lcio::LCCollectionVec(lcio::LCIO::VERTEX);
		lcio::LCCollectionVec *colRP = new lcio::LCCollectionVec(lcio::LCIO::RECONSTRUCTEDPARTICLE);
		for(unsigned int n=0;n<pvvtx->size();n++){
			// set ID to the lcfiplus::Vertex
			const lcfiplus::Vertex *flavtx = (*pvvtx)[n];
			flavtx->setId(n);
			
			// make new vertex
			lcio::VertexImpl *lciovtx = new lcio::VertexImpl;
			// make new recoparticle including vertex
			lcio::ReconstructedParticleImpl *lciorp = new lcio::ReconstructedParticleImpl;
			
			// associate particles
			TLorentzVector lv;
			float charge = 0.;
			for(unsigned int ntr = 0; ntr < flavtx->getTracks().size(); ntr++){
				const lcfiplus::Track *flatr = flavtx->getTracks()[ntr];
				lcio::ReconstructedParticle *lciotr = _trackLCIORel[flatr->getId()];
				lv += (*flatr);
				charge += flatr->getCharge();
				lciorp->addParticle(lciotr);
			}
			float mom[3] = {lv.Px(), lv.Py(), lv.Pz()};
			float vpos[3] = {flavtx->getX(), flavtx->getY(), flavtx->getZ()};
			lciorp->setType(3); // 0: unknown 1: single 2:v0 3: compound 4:jet
			lciorp->setMomentum(mom);
			lciorp->setEnergy(lv.E());
			lciorp->setMass(lv.M());
			lciorp->setCharge(charge);
			lciorp->setReferencePoint(vpos);
			lciorp->setStartVertex(lciovtx);
			// ignore covmatrix of rp
			colRP->addElement(lciorp);
			
			// vertex initialization
			lciovtx->setPrimary(n==0); // TODO: too simple
			lciovtx->setAlgorithmType("lcfiplus");
			lciovtx->setPosition(vpos);
			lciovtx->setCovMatrix(flavtx->getCov());
			lciovtx->setChi2(flavtx->getChi2());
			lciovtx->setProbability(flavtx->getProb());
			lciovtx->setAssociatedParticle(lciorp);

			/*
			cout << scientific << "ConvertVertex position: "
				<< vpos[0] << ","
				<< vpos[1] << ","
				<< vpos[2] << fixed << endl;
			
			const float* cov = flavtx->getCov();
			cout << scientific << "ConvertVertex covMatrix: "
				<< cov[0] << ","
				<< cov[1] << ","
				<< cov[2] << ","
				<< cov[3] << ","
				<< cov[4] << ","
				<< cov[5] << fixed << endl;

			const vector<float>& cov2 = lciovtx->getCovMatrix();
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
			
			// register to collection
			col->addElement(lciovtx);
		}
		_event->addCollection(col, (newName ? newName : vertexName));
		_event->addCollection(colRP, newRPName ? newRPName : (const char *)TString::Format("%s_RP", newName ? newName : vertexName));
		
		//cout << "ConvertVertex finished. # vertices = " << col->getNumberOfElements() << endl;

	}

	void LCIOStorer::ConvertJet(const char *jetName, const char *newName, bool extractVertex)
	{
		if(!_event){
			cerr << "LCIOStorer::ConvertJet: LCIO event has not been initialized." << endl;
			return ;
		}
		JetVec *pvjet;
		Event::Instance()->Get(jetName, pvjet);
	
		if(!pvjet){
			cerr << "LCIOStorer::ConvertJet: failed to obtain jets from the specified name: " << jetName << endl;
			return ;
		}
		
		_jetLCIORel.clear();
//		_jetLCIORel.resize(pvjet->size()); // map has no resize()
		
		// make collection
		lcio::LCCollectionVec *col = new lcio::LCCollectionVec(lcio::LCIO::RECONSTRUCTEDPARTICLE);
		for(unsigned int n=0;n<pvjet->size();n++){
			// set ID to the lcfiplus::Vertex
			const lcfiplus::Jet *flajet = (*pvjet)[n];
			flajet->setId(n);
			
			// make new recoparticle including vertex
			lcio::ReconstructedParticleImpl *lciojet = new lcio::ReconstructedParticleImpl;
			
			// associate particles
			float charge = 0.;
			for(unsigned int ntr = 0; ntr < flajet->getTracks().size(); ntr++){
				const lcfiplus::Track *flatr = flajet->getTracks()[ntr];
				lcio::ReconstructedParticle *lciotr = _trackLCIORel[flatr->getId()];
				charge += flatr->getCharge();
				lciojet->addParticle(lciotr);
				//cout << "LCIOStorer::ConvertJet: add track: id = " << flatr->getId() << ", energy = " << flatr->E() << flush;
				//cout << ", lcio energy = " << lciotr->getEnergy() << endl;
			}
			for(unsigned int nneut = 0; nneut < flajet->getNeutrals().size(); nneut++){
				const lcfiplus::Neutral *flaneut = flajet->getNeutrals()[nneut];
				lcio::ReconstructedParticle *lcioneut = _neutralLCIORel[flaneut->getId()];
				lciojet->addParticle(lcioneut);
			}
			for(unsigned int nvtx = 0; nvtx < flajet->getVertices().size(); nvtx++){
				const lcfiplus::Vertex *flavtx = flajet->getVertices()[nvtx];
				if(!extractVertex && flavtx->getId() >= 0){ // valid ID
					lcio::ReconstructedParticle *lciovtx = _vtxLCIORel[flavtx]->getAssociatedParticle();
					charge += lciovtx->getCharge();
					lciojet->addParticle(lciovtx);
				}else{ // ID not available, add daughter particles
					for(unsigned int ntr = 0; ntr < flavtx->getTracks().size(); ntr++){
						const lcfiplus::Track *flatr = flavtx->getTracks()[ntr];
						lcio::ReconstructedParticle *lciotr = _trackLCIORel[flatr->getId()];
						charge += flatr->getCharge();
						lciojet->addParticle(lciotr);
					}
					//cout << "LCIOStorer::ConvertJet: add " << flavtx->getTracks().size() << "tracks in jet." << endl;
				}
			}
			
			const TLorentzVector &lv = *flajet;
			float mom[3] = {lv.Px(), lv.Py(), lv.Pz()};
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
			
			// register to collection
			col->addElement(lciojet);
		}
		_event->addCollection(col, (newName ? newName : jetName));
		
		cout << "ConvertJet finished. # jets = " << col->getNumberOfElements() << endl;
	}
}

void LCIOStorer::WritePID(lcio::LCCollection *lciocol, lcio::ReconstructedParticle *lciojet, const lcfiplus::Jet *lcfijet, const char *paramname)
{
	// init pid handler
	lcio::PIDHandler pidh(lciocol);

	// obtain Parameters
	const Parameters *lcfiparams = lcfijet->getParam(paramname);
	const map<string, pair<const type_info *, void *> > & paramMap = lcfiparams->paramMap();
	
	vector<string> outParamNames;
	vector<float> outParamData;

	map<string, pair<const type_info *, void *> >::const_iterator it;
	for(it = paramMap.begin(); it != paramMap.end(); it++){
		string name = it->first;

		double d = lcfiparams->get(name.c_str(),(double)0.);

		outParamNames.push_back(name);
		outParamData.push_back(d);
	}

	// register parameters
	int algoID;
	try{
		algoID = pidh.getAlgorithmID(string(paramname));
	}catch(lcio::UnknownAlgorithm &e){
		algoID = pidh.addAlgorithm(paramname, outParamNames);
	}

	//                         type, PDG, likelihood: as same as old LCFI
	pidh.setParticleID(lciojet, 42, 9999, 0.0, algoID, outParamData);
}


void LCIOStorer::WriteAllPIDs(lcio::LCCollection *lciocol, lcio::ReconstructedParticle *lciojet, const lcfiplus::Jet *lcfijet)
{
	const map<string, Parameters> & parammap = lcfijet->params();

	map<string, Parameters>::const_iterator it;
	for(it = parammap.begin(); it != parammap.end(); it++){
		WritePID(lciocol, lciojet, lcfijet, it->first.c_str());
	}

	if(parammap.size())
		cout << parammap.size() << " PID parameter sets written to jet." << endl;
}
