// LCIOStorer.cc
#include "LCIOStorer.h"
#include "flavtag.h"
#include "EventStore.h"
#include "JetFinder.h"

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

#include "UTIL/LCRelationNavigator.h"
//#include "UTIL/PIDHandler.h"

#include <assert.h>
#include <math.h>

//using namespace lcio;

namespace flavtag{

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
		_pfoColName = pfoColName;
		_mcColName = mcColName;
		_mcpfoColName = mcpfoColName;

		EventStore::Instance()->Register<Track>(trackName,_pTracks);
		EventStore::Instance()->Register<Neutral>(neutralName,_pNeutrals);
		EventStore::Instance()->Register<MCParticle>(mcpName,_pMCPs);
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
		EventStore *s = EventStore::Instance();
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
				
				if(it->second.classname == "vector<flavtag::Vertex*>"){
					ConvertVertex(it->first.c_str(), newname.c_str());
				}
				else if(it->second.classname == "vector<flavtag::Jet*>"){
					bool extractVertex = it->second.flag & EventStore::JET_EXTRACT_VERTEX;
					ConvertJet(it->first.c_str(), newname.c_str(), extractVertex);
				}else{
					cerr << "LCIOStorer::WriteEvent(): unknown conversion for lcio: class " << it->second.classname << ", ignored." << endl;
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
			cerr << "LCIOStorer::WriteEvent(): event cannot be written." << endl;
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
//		if(_pTracks == NULL || _pNeutrals == NULL || _pMCPs == NULL)throw(new flavtag::Exception("LCIOStorer::SetEvent: Event buffer not set."));
		_event = evt;
		if(_pTracks == NULL || _pNeutrals == NULL || _pMCPs == NULL)return;

		// collections
		lcio::LCCollection* colMC = evt->getCollection(_mcColName);
		lcio::LCCollection* colPFO = evt->getCollection(_pfoColName);
		lcio::LCRelationNavigator* nav = new lcio::LCRelationNavigator(evt->getCollection(_mcpfoColName));
		
		// memory allocation
		_pMCPs->clear();
		_pNeutrals->clear();
		_pTracks->clear();
		
		_pMCPs->reserve(colMC->getNumberOfElements());
		_pNeutrals->reserve(colPFO->getNumberOfElements());
		_pTracks->reserve(colPFO->getNumberOfElements());

		_mcpLCIORel.clear();
		_trackLCIORel.clear();
		_neutralLCIORel.clear();
		
		_mcpLCIORel.reserve(colMC->getNumberOfElements());
		_trackLCIORel.reserve(colPFO->getNumberOfElements());
		_neutralLCIORel.reserve(colPFO->getNumberOfElements());

		// Relation of MCParticle between LCIO and flavtag
		map<lcio::MCParticle *, flavtag::MCParticle *> mcpMap;
		
		int mcIdCounter(0);
		// convert MCParticle /////////////////////////////////////////////////////////////////////////
		for (mcIdCounter=0; mcIdCounter<colMC->getNumberOfElements(); ++mcIdCounter) {
			lcio::MCParticle* mcp = dynamic_cast<lcio::MCParticle*>( colMC->getElementAt(mcIdCounter) );
			assert(mcp != 0);
			//if (mcp->isBackscatter()) continue;
			
			// find parent id
			// TODO Really OK for single parent??
			flavtag::MCParticle *parent = 0;
			if (mcp->getParents().size()>0) {
				map<lcio::MCParticle*,flavtag::MCParticle *>::iterator iter = mcpMap.find(mcp->getParents()[0]);

				if ( iter == mcpMap.end() )	throw(new flavtag::Exception("parent not found in association map"));
			
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
			flavtag::MCParticle *mcpNew = new MCParticle(mcIdCounter, mcp->getPDG(), parent, mcp->getCharge(),
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
			flavtag::MCParticle *mcpf = NULL;
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

			if (pfo->getCharge() != 0) {
				int trkSize = pfo->getTracks().size();
				assert(trkSize>0);

				flavtag::Track * track = new flavtag::Track;

				track->setId(++trkIdCounter); // start from 1...
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
				float par[flavtag::tpar::parN];
				par[flavtag::tpar::d0] = trk->getD0();
				par[flavtag::tpar::z0] = trk->getZ0();
				par[flavtag::tpar::ph] = trk->getPhi();
				par[flavtag::tpar::om] = trk->getOmega();
				par[flavtag::tpar::td] = trk->getTanLambda();
				track->setHelix(par);
				
				// ... and the covariance matrix
				const vector<float>& cov = trk->getCovMatrix();
				float fcov[flavtag::tpar::covN];
				for (int i=0; i<flavtag::tpar::covN; ++i) {
					fcov[i] = cov[i];
				}
				track->setCovMatrix(fcov);

				track->setChi2(trk->getChi2());
				track->setNdf(trk->getNdf());

				// store detector hit numbers
				float nhits[flavtag::tpar::hitN];
				for (int i=0; i<flavtag::tpar::hitN; ++i) {
					nhits[i] = trk->getSubdetectorHitNumbers()[i];
				}
				track->setTrackHits(nhits);
				track->setRadiusOfInnermostHit(trk->getRadiusOfInnermostHit());
				//printf("rimh = %f\n",trkData.rimh);

				track->setCaloEdep(subE);
				
				// register!
				_pTracks->push_back(track);
				_trackLCIORel[trkIdCounter] = pfo;
			
			} else { // neutrals
				// assert(clusSize>0);
				// above is not necessarily true e.g. for a V0 which just has two tracks

				flavtag::Neutral *neut = new Neutral;
				neut->setId(++neutIdCounter);
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
				_neutralLCIORel[neutIdCounter] = pfo;

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
	}

	void LCIOStorer::ConvertVertex(const char *vertexName, const char *newName, const char *newRPName)
	{
		if(!_event){
			cerr << "LCIOStorer::ConvertVertex: LCIO event has not been initialized." << endl;
			return ;
		}
		const vector<Vertex *> *pvvtx;
		EventStore::Instance()->Get(vertexName, pvvtx);
	
		if(!pvvtx){
			cerr << "LCIOStorer::ConvertVertex: failed to obtain vertices from the specified name: " << vertexName << endl;
			return ;
		}
		
		_vtxLCIORel.clear();
		_vtxLCIORel.resize(pvvtx->size());
		
		// make collection
		lcio::LCCollectionVec *col = new lcio::LCCollectionVec(lcio::LCIO::VERTEX);
		lcio::LCCollectionVec *colRP = new lcio::LCCollectionVec(lcio::LCIO::RECONSTRUCTEDPARTICLE);
		for(unsigned int n=0;n<pvvtx->size();n++){
			// set ID to the flavtag::Vertex
			flavtag::Vertex *flavtx = (*pvvtx)[n];
			flavtx->setId(n);
			
			// make new vertex
			lcio::VertexImpl *lciovtx = new lcio::VertexImpl;
			// make new recoparticle including vertex
			lcio::ReconstructedParticleImpl *lciorp = new lcio::ReconstructedParticleImpl;
			
			// associate particles
			TLorentzVector lv;
			float charge = 0.;
			for(unsigned int ntr = 0; ntr < flavtx->getTracks().size(); ntr++){
				flavtag::Track *flatr = flavtx->getTracks()[ntr];
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
			lciovtx->setAlgorithmType("flavtag");
			lciovtx->setPosition(vpos);
			lciovtx->setCovMatrix(flavtx->getCov());
			lciovtx->setChi2(flavtx->getChi2());
			lciovtx->setProbability(flavtx->getProb());
			lciovtx->setAssociatedParticle(lciorp);
			
			// store relation
			_vtxLCIORel[n] = lciovtx;
			
			// register to collection
			col->addElement(lciovtx);
		}
		_event->addCollection(col, (newName ? newName : vertexName));
		_event->addCollection(colRP, newRPName ? newRPName : (const char *)TString::Format("%s_RP", newName ? newName : vertexName));
		
		cout << "ConvertVertex finished. # vertices = " << col->getNumberOfElements() << endl;

	}

	void LCIOStorer::ConvertJet(const char *jetName, const char *newName, bool extractVertex)
	{
		if(!_event){
			cerr << "LCIOStorer::ConvertJet: LCIO event has not been initialized." << endl;
			return ;
		}
		const vector<Jet *> *pvjet;
		EventStore::Instance()->Get(jetName, pvjet);
	
		if(!pvjet){
			cerr << "LCIOStorer::ConvertJet: failed to obtain jets from the specified name: " << jetName << endl;
			return ;
		}
		
		_jetLCIORel.clear();
		_jetLCIORel.resize(pvjet->size());
		
		// make collection
		lcio::LCCollectionVec *col = new lcio::LCCollectionVec(lcio::LCIO::RECONSTRUCTEDPARTICLE);
		for(unsigned int n=0;n<pvjet->size();n++){
			// set ID to the flavtag::Vertex
			flavtag::Jet *flajet = (*pvjet)[n];
			flajet->setId(n);
			
			// make new recoparticle including vertex
			lcio::ReconstructedParticleImpl *lciojet = new lcio::ReconstructedParticleImpl;
			
			// associate particles
			float charge = 0.;
			for(unsigned int ntr = 0; ntr < flajet->getTracks().size(); ntr++){
				flavtag::Track *flatr = flajet->getTracks()[ntr];
				lcio::ReconstructedParticle *lciotr = _trackLCIORel[flatr->getId()];
				charge += flatr->getCharge();
				lciojet->addParticle(lciotr);
			}
			for(unsigned int nneut = 0; nneut < flajet->getNeutrals().size(); nneut++){
				flavtag::Neutral *flaneut = flajet->getNeutrals()[nneut];
				lcio::ReconstructedParticle *lcioneut = _neutralLCIORel[flaneut->getId()];
				lciojet->addParticle(lcioneut);
			}
			for(unsigned int nvtx = 0; nvtx < flajet->getVertices().size(); nvtx++){
				flavtag::Vertex *flavtx = flajet->getVertices()[nvtx];
				if(!extractVertex && flavtx->getId() >= 0){ // valid ID
					lcio::ReconstructedParticle *lciovtx = _vtxLCIORel[flavtx->getId()]->getAssociatedParticle();
					charge += lciovtx->getCharge();
					lciojet->addParticle(lciovtx);
				}else{ // ID not available, add daughter particles
					for(unsigned int ntr = 0; ntr < flavtx->getTracks().size(); ntr++){
						flavtag::Track *flatr = flavtx->getTracks()[ntr];
						lcio::ReconstructedParticle *lciotr = _trackLCIORel[flatr->getId()];
						charge += flatr->getCharge();
						lciojet->addParticle(lciotr);
					}
				}
			}
			
			TLorentzVector &lv = *flajet;
			float mom[3] = {lv.Px(), lv.Py(), lv.Pz()};
			lciojet->setType(4); // 0: unknown 1: single 2:v0 3: compound 4:jet
			lciojet->setMomentum(mom);
			lciojet->setEnergy(lv.E());
			lciojet->setMass(lv.M());
			lciojet->setCharge(charge);
			// ignore covmatrix of rp
			
			// store relation
			_jetLCIORel[n] = lciojet;
			
			// register to collection
			col->addElement(lciojet);
		}
		_event->addCollection(col, (newName ? newName : jetName));
		
		cout << "ConvertJet finished. # jets = " << col->getNumberOfElements() << endl;
	}
}


