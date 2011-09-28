#include <assert.h>
#include <string>

#include "EventStore.h"
#include "LcfiInterface.h"
#include "JetFinder.h"
#include "TreeStorer.h"
#include "LCIOStorer.h"
#include "VertexFitterLCFI.h"
#include "VertexFinderTearDown.h"
#include "VertexFinderSuehara.h"

#include "TROOT.h"
#include "TInterpreter.h"
#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TCut.h"
#include "TRandom.h"
#include "TNtuple.h"

#include "Driver.h"
#include "Suehara.h"
#include "algoEtc.h"
#include "geometry.h"

using namespace lcfiplus;
using namespace lcfiplus::algoEtc;

		void lcioTest(const char *infile)
		{
			LCIOStorer ls(infile);
			ls.InitCollections("PandoraPFOs","MCParticlesSkimmed","RecoMCTruthLink","Tracks","Neutrals","MCParticles");
			ls.InitVertexCollection("ZVRESVertices_2Jets", "TestVertices");
			ls.InitJetCollection("Durham_2Jets", "TestJets");
			cout << "LCIO initialization successful." << endl;

			Event::Instance()->Print();

			int nev = 0;
			while(ls.Next()){
				const vector<Vertex *> *pvcol;
				bool b = Event::Instance()->Get("TestVertices", pvcol);

				if(!b){cerr << "Obtaining vertex collection failed." << endl;break;}

				cout << "Event # " << nev++ << ", # vertices = " << pvcol->size() << endl;

				for(unsigned int n=0; n < pvcol->size(); n++){
					Vertex *v = (*pvcol)[n];
					cout << "  Vertex # " << n << ", chi2 = " << v->getChi2() << ", prob = " << v->getProb();
					cout << ", pos = (" << v->getX() << ", " << v->getY() << ", " << v->getZ() << ")";
					cout << ", # tracks = " << v->getTracks().size() << endl;
				}

				const vector<Jet *> *pjcol;
				b = Event::Instance()->Get("TestJets", pjcol);
				if(!b){cerr << "Obtaining jet collection failed." << endl;break;}

				cout << "Event # " << nev << ", # jets = " << pjcol->size() << endl;

				for(unsigned int n=0;n < pjcol->size(); n++){
					Jet *j = (*pjcol)[n];
					cout << "  Jet # " << n << ", p = (" << j->Px() << ", " << j->Py() << ", " << j->Pz() << "), E = " << j->E();
					cout << ", # tracks = " << j->getTracks().size() << ", # neutrals = " << j->getNeutrals().size() << endl;
				}
			}

		}

		void treeTest(const char *treefile)
		{
			TreeStorer ts(treefile,"FlavTagTree",TreeStorer::mode_input);

			ts.RegisterAll();
			cout << "Tree registration finished." << endl;

			Event::Instance()->Print();
		}

		void lcioToTree(const char *infile, const char *outfile)
		{
			// initialization - input slcio
			LCIOStorer ls(infile);
			ls.InitCollections("PandoraPFOs","MCParticlesSkimmed","RecoMCTruthLink","Tracks","Neutrals","MCParticles");
			cout << "LCIO initialization successful." << endl;

			// initialization - output root
			TreeStorer ts(outfile,"FlavTagTree",TreeStorer::mode_output);
			ts.Register("Tracks");
			ts.Register("Neutrals");
			ts.Register("MCParticles");
			cout << "ROOT initialization successful." << endl;
			Event::Instance()->Print();

			while(ls.Next()){
				ts.Fill();
			}
			ts.Write();
		}

		void lcioToLcio(const char *infile, const char *outfile)
		{
			LCIOStorer ls(infile, outfile);
//			ls.InitCollections("PandoraPFOs","MCParticlesSkimmed","RecoMCTruthLink","Tracks","Neutrals","MCParticles");
			while(ls.Next()){
				ls.WriteEvent();
			}
		}

		void outVertex(const char *inputfile, const char *outputfile, int nStart, int nEnd, int bbhh, int inclVertex)
		{
			// initialization - input slcio
			LCIOStorer ls(inputfile, outputfile, true);
			ls.InitCollections("PandoraPFOs","MCParticlesSkimmed","RecoMCTruthLink","Tracks","Neutrals","MCParticles");
			cout << "LCIO initialization successful." << endl;

			vector<Vertex *> *pvertices = 0;
			Event::Instance()->Register<Vertex>("BuildUpVertices", pvertices,inclVertex ? EventStore::PERSIST : 0);
			vector<Jet *> *pjets = 0;
			Event::Instance()->Register<Jet>("SueharaJets", pjets,EventStore::PERSIST);

			int nev = -1;

			Event *event = Event::Instance();

			while(ls.Next()){
				nev ++;
				if(nev < nStart){continue;}
				if(nev >= nEnd)break;

				// check bbbbbb (reject H->WW etc.)
				if(bbhh){
					int hcount = 0;

					const vector<MCParticle *> *pmcp;
					Event::Instance()->Get("MCParticles", pmcp);
					const vector<MCParticle *> &mcps = *pmcp;

					for(unsigned int i=0;i<mcps.size();i++){
						if(mcps[i]->getPDG() != 25)continue;
	
						// higgs
						if(mcps[i]->getDaughters().size() != 2){
							cout << "ERR: # of higgs daughters = " << mcps[i]->getDaughters().size() << endl;
							break;
						}
						if(abs(mcps[i]->getDaughters()[0]->getPDG()) != 5)break;
						if(abs(mcps[i]->getDaughters()[1]->getPDG()) != 5)break;
						hcount ++;
					}
					if(hcount < 2)continue;
				}

				cout << "Event #" << nev << endl;

				pvertices->clear();
				pjets->clear();

				TrackSelectorConfig secVtxCfg;
				// original setup 110214
				secVtxCfg.maxD0 = 10;
				secVtxCfg.maxZ0 = 20;
				secVtxCfg.minPt = 0.1;
				secVtxCfg.maxInnermostHitRadius = 1e10;

				// 110214
				secVtxCfg.maxD0Err = 0.1;
				secVtxCfg.maxZ0Err = 0.1;

				// OR
				secVtxCfg.minTpcHits = 20;
				secVtxCfg.minFtdHits = 3;
				secVtxCfg.minVtxHitsWithoutTpcFtd = 3;
				secVtxCfg.minVtxPlusFtdHits = 0;
	
				JetConfig jetCfg;
				jetCfg.nJet = 6;
				JetFinder* jetFinder = new JetFinder(jetCfg);
				double ycut = 0;

//				LcfiInterface interface(event);

				// cut bad tracks
				const vector<Track *> &tracks = event->getTracks();
				vector<Track *> passedTracks = TrackSelector() (tracks, secVtxCfg);

// 				for(unsigned int i=0;i<tracks.size();i++){
// 					if(interface.passesCut(tracks[i], secVtxCfg))
// 						passedTracks.push_back(tracks[i]);
// 				}

				// build up vertexing
				VertexFinderSuehara::buildUp(passedTracks, *pvertices, 25, 9, 10);
				VertexFinderSuehara::associateIPTracks(*pvertices,1.);


				// select vertices
				vector<Vertex *> selectedVertices;
				vector<Track *> residualTracks = tracks;
				for(vector<Vertex *>::iterator it = pvertices->begin(); it != pvertices->end();it++){
					Vertex *v = *it;
					double mass = 0.;
					if(v->getTracks().size() == 2)
						mass = (*(TLorentzVector *)(v->getTracks()[0]) + *(TLorentzVector *)(v->getTracks()[1])).M();
					if((mass < .488 || mass > .508) && v->getPos().Mag() < 30 && v->getPos().Mag() > 0.3){
						selectedVertices.push_back(v);
						for(vector<Track *>::const_iterator itt = v->getTracks().begin(); itt != v->getTracks().end(); itt++){
							vector<Track *>::iterator itt2 = remove_if(residualTracks.begin(), residualTracks.end(), bind2nd(equal_to<Track *>(), *itt));
							residualTracks.erase(itt2, residualTracks.end());
						}
					}
				}
				*pjets = jetFinder->run(residualTracks, event->getNeutrals(), selectedVertices, &ycut);

//				if(inclVertex)
//					ls.ConvertVertex("BuildUpVertices");
//				ls.ConvertJet("SueharaJets");
				ls.WriteEvent();
			}
		}

		void testSuehara()
		{
			TreeStorer ts("share/zpole_v.root","FlavTagTree",TreeStorer::mode_input);
			ts.RegisterAll();
			Event::Instance()->Print();

			const vector<lcfiplus::Vertex *> *pvertices;
			Event::Instance()->Get("TearDownVertices",pvertices);

			for(int nev =0;nev< 100;nev++){
				ts.GetEntry(nev);
				cout << (*pvertices).size() << " " << (*pvertices)[0]->getChi2() << endl;
			}
		}

		void TearDownVertexing(const char *inputfile, const char *outputfile,
				double chi2th_primary, double chi2th_secondary, bool bonly, bool verbose)
		{
				LCIOStorer ls(inputfile);
				ls.InitCollections("PandoraPFOs","MCParticlesSkimmed","RecoMCTruthLink","Tracks","Neutrals","MCParticles");
				cout << "LCIO initialization successful." << endl;

//			TreeStorer ts(inputfile,"FlavTagTree",TreeStorer::mode_input);
//			ts.RegisterAll();

			// register vertex collection
			vector<lcfiplus::Vertex *> * pvertices;
			Event::Instance()->Register<lcfiplus::Vertex>("TearDownVertices", pvertices);

			Event::Instance()->Print();

			TreeStorer tsout(outputfile,"FlavTagTree",TreeStorer::mode_output);
			tsout.Register("Tracks");
			tsout.Register("Neutrals");
			tsout.Register("MCParticles");
			tsout.Register("TearDownVertices");

			Event *event = Event::Instance();

			int nev = -1;
			while(ls.Next()){
				nev ++;
//			for(int nev = 0;nev < ts.GetTree()->GetEntries();nev++){
//				ts.GetEntry(nev);

				if(verbose || nev%100 == 0)
					cout << "Event # " << nev << " # B is " << event->mcNumberOfB() << " ------------------------" << endl;

				if(bonly && event->mcNumberOfB()==0)continue;

				// Primary vertex by teardown
				vector<Track *> tracks = event->getTracks();
				// 		// IP constraint;
				// 		float iperr[6];
				// 		memset(iperr,0,sizeof(iperr));
				// 		iperr[Vertex::xx] = 640e-6; // 640 nm
				// 		iperr[Vertex::yy] = 5.7e-6; // 5.7 nm
				// 		iperr[Vertex::zz] = 64e-3;	// 64 um (crossing angle constraint)
				// 		float ippos[3];
				// 		ippos[0] = gRandom->Gaus(0,iperr[Vertex::xx]);
				// 		ippos[1] = gRandom->Gaus(0,iperr[Vertex::yy]);
				// 		ippos[2] = gRandom->Gaus(0,iperr[Vertex::zz]);
				// 		Vertex *ipConst = new Vertex(0,0,ippos[0], ippos[1], ippos[2], iperr);

				vector<Track*> beamTracks;
				beamTracks.resize(2);
				makeBeamTracks(beamTracks[0], beamTracks[1]);

				pvertices->clear();

//				lcfiplus::Vertex * vtx = lcfiplus::VertexFinderTearDown<vector>()(tracks, &beamTracks, chi2th_primary, 0);
				lcfiplus::Vertex * vtx = findPrimaryVertex(tracks,chi2th_primary);
				if(!vtx){
					cout << "Primary vertex not found! # tracks = " << tracks.size() << endl;
					continue;
				}
				pvertices->push_back(vtx);

				Track *worstTrack = vtx->getWorstTrack();
				if(verbose)
					cout << "Primary vertex: " << vtx->getTracks().size() << "/" << tracks.size() << ", " << vtx->getChi2Track(worstTrack) << endl;

				// jet clustering
				// TODO: y threshold
				JetConfig jetCfg;
				jetCfg.nJet = 2;
				JetFinder* jetFinder = new JetFinder(jetCfg);
				vector<Jet*> jets = jetFinder->run(event->getTracks(),event->getNeutrals());

				// jet loop
				for(int nj=0;nj<2;nj++){
					// copy tracks in the jet into a list
					const vector<Track *> &v = jets[nj]->getTracks();
					list<Track *> tracksInJet;
					tracksInJet.resize(v.size());
					copy(v.begin(), v.end(), tracksInJet.begin());

					// remove primary tracks
					for(unsigned int i=0;i<vtx->getTracks().size();i++){
						tracksInJet.remove(vtx->getTracks()[i]);
					}
					if(verbose)
						cout << "Jet #" << nj << " secondary tracks: " << tracksInJet.size() << "/" << v.size() << endl;

					while(tracksInJet.size() >= 2){
						// search secondary vertex
						list<Track *> residuals;
//						lcfiplus::Vertex *secvtx = lcfiplus::VertexFinderTearDown<list>()(tracksInJet,0, chi2th_secondary, &residuals);
						lcfiplus::Vertex *secvtx = lcfiplus::VertexFinderTearDown<list, VertexFitterSimple>()(tracksInJet,0, chi2th_secondary, &residuals);
						if(!secvtx)break;

						pvertices->push_back(secvtx);

						if(verbose)
							cout << "    Secondary vertex found! pos = (" << secvtx->getX() << "," << secvtx->getY() << "," << secvtx->getZ() << "), chi2 = "
								<< secvtx->getChi2() << endl;
						for(unsigned int i=0;i<secvtx->getTracks().size();i++){
							Track * tr = secvtx->getTracks()[i];
							if(verbose)
								cout << "        Track #" << i << ": p = (" << tr->Px() << "," << tr->Py() << "," << tr->Pz() << "), chi2 = "
									<< secvtx->getChi2Track(tr) << endl;
						}
						tracksInJet = residuals;
					}
				}

				tsout.Fill();

				/*		do{
							vtx = VertexFitterLCFI(tracks);
							Track *worstTrack = vtx->getWorstTrack();

							cout << "( " << vtx->getX() << ", " << vtx->getY() << ", " << vtx->getZ() << "), " << vtx->getProb() << ", " << vtx->getTracks().size() << ", ";
							cout << worstTrack->getD0() << ", " << worstTrack->getZ0() << ", ";
							cout << vtx->getChi2Track(worstTrack) << ", ";

							cout << event->getMCParticle(worstTrack)->getFlavorTagCategory() << ", ";
							cout << event->getMCParticle(worstTrack)->getVx() << ", " << event->getMCParticle(worstTrack)->getVy() << ", " << event->getMCParticle(worstTrack)->getVz() << endl;

				// todo: not suitable for vector
				for(vector<Track *>::iterator it = tracks.begin();it!=tracks.end();it++){
				if(*it == worstTrack){
				tracks.erase(it);
				break; // expects uniqueness
				}
				}

				}while(vtx->getTracks().size()>1);
				 */
			}
			tsout.Write();
		}

void checkMCTearDown(const char *inputfile, const char *outputfile, int nStart, int nEnd)
{
	LCIOStorer ls(inputfile);
	ls.InitCollections("PandoraPFOs","MCParticlesSkimmed","RecoMCTruthLink","Tracks","Neutrals","MCParticles");
	cout << "LCIO initialization successful." << endl;

	int nev=0;

// 		class VertexAnalysis{
// 			public:
// 				int event;
// 				double ipchi2[2];
// 				double vchi2;
// 				double trchi2[2];
//				int tag[2];
// 				int success;
// 				double truevtx[3][2];
// 				double recovtx[3];
// 
// 			ClassDefNV(VertexAnalysis,1);
// 		};
	VertexAnalysis *data = new VertexAnalysis;

// class MatchMCRecoVertex{
// 	public:
// 		int event;
// 		MCVertex *mcv;
// 		vector<Vertex *> recov;
// 
// 	ClassDefNV(MatchMCRecoVertex,1);
// };
	MatchMCRecoVertex *data2 = new MatchMCRecoVertex;

	TFile *pf = TFile::Open(outputfile,"recreate");
	TTree *tree = new TTree("tree","tree");
	tree->Branch("vana","VertexAnalysis",&data);
	TTree *tree2 = new TTree("tree2","tree2");
	tree2->Branch("mcreco","MatchMCRecoVertex",&data2);

	pf->Print();

	int summcv = 0;
	int summcrecov = 0;

	Event *event = Event::Instance();

	while(ls.Next()){
		if(nev < nStart){nev++;continue;}
		if(nev >= nEnd)break;


		// track check
 		cout << "Track check start: ###########################" << endl;
 		const vector<Track *> &tracks = event->getTracks();
 		cout << "# tracks = " << tracks.size() << endl;
 		for(unsigned int i=0;i<tracks.size();i++){
 			MCParticle *mc = tracks[i]->getMcp();
 			if(!mc){cout << "MC not found!, trp = " << tracks[i]->P() << endl;}
			else{
				cout << "MCP = " << mc->P() << ", McID = " << mc->getId() << ", RecoP = " << tracks[i]->P() << ", TrID = " << tracks[i]->getId() << endl;
			}
 		}
 		cout << "Track check end: ###########################" << endl;

		vector<MCVertex *> mcv = findMcVertex(event->getMCParticles());
		map<MCVertex*,int> table;
		matchMcVertex(*event, mcv, table);

		cout << "MCVertex size = " << mcv.size() << endl;
		for(unsigned int i=0;i<mcv.size();i++){
			MCVertex * pmcv = mcv[i];
			TVector3 v = pmcv->getParent()->getEndVertex();
			cout << v.x() << ", " << v.y() << ", " << v.z()
						<< " with " << pmcv->getDaughters().size() << " MC tracks";
			cout << ", (p,ct,pdg,id) = ";
			for(unsigned int j=0;j<pmcv->getDaughters().size();j++){
				cout << "(" << pmcv->getDaughters()[j]->P() << ", ";
//				TVector3 ve = pmcv->getDaughters()[j]->getEndVertex();
//				TVector3 vs = pmcv->getDaughters()[j]->getVertex();
//				cout << (ve-vs).Mag() << ")";
				cout << pmcv->getDaughters()[j]->CosTheta() << ", ";
				cout << pmcv->getDaughters()[j]->getPDG() << ", ";
				cout << pmcv->getDaughters()[j]->getId() << ")";
			}
			cout << ", " << pmcv->getRecoTracks().size() << " Reco tracks";
			for(unsigned int j=0;j<pmcv->getRecoTracks().size();j++){
				cout << "(" << pmcv->getRecoTracks()[j]->P() << ", ";
				cout << pmcv->getRecoTracks()[j]->CosTheta() << ", ";
				cout << pmcv->getRecoTracks()[j]->getId() << ")";
			}

			if(pmcv->getRecoTracks().size() >= 2){
				lcfiplus::Vertex *vtx = lcfiplus::VertexFinderTearDown<vector, VertexFitterSimple>()(pmcv->getRecoTracks(),0, 10000, NULL);
				if(vtx == NULL)cout << ", Vertex not found";
				else{
					cout << ", Vertex found at: " << vtx->getX() << ", " << vtx->getY() << ", " << vtx->getZ() << " with " << vtx->getTracks().size() << " tracks, ";
					cout << "chi2 = " << vtx->getChi2(); 
				}
			}
			cout << "." << endl;
		}


		data->event = nev;

		LcfiInterface lcfi;

//		vector<Vertex *> *pvtx = lcfiplus::VertexFinderSuehara::findCandidates(event->getTracks(), 10000, 100, 9);
		vector<Vertex *> vtx;
		vector<Vertex *> *pvtx = &vtx;
		lcfiplus::VertexFinderSuehara::buildUp(event->getTracks(), vtx, 25, 25, 100);
		cout << "Event #" << nev << ": " << pvtx->size() << " vertices found." << endl;
		for(unsigned int i=0;i<vtx.size();i++){
			cout << "position: (" << vtx[i]->getX() << ", " << vtx[i]->getY() << ", " << vtx[i]->getZ() << "), ";
			cout << "tracks: ";
			for(unsigned int j=0;j<vtx[i]->getTracks().size();j++){
				Track *tr = vtx[i]->getTracks()[j];
				cout << tr->getId() << " ( " << vtx[i]->getChi2Track(tr) << ", " << lcfi.getChi2TrackVtx(vtx[i], tr) << "), ";
			}
			cout << "chi2 = " << vtx[i]->getChi2();
			if(fabs(vtx[i]->getZ())>100){
				cout << endl << "err: (" << sqrt(vtx[i]->getCov()[Vertex::xx]) << ", "
																	<< sqrt(vtx[i]->getCov()[Vertex::yy]) << ", "
																	<< sqrt(vtx[i]->getCov()[Vertex::zz]) << ")" << endl;
				for(unsigned int j=0;j<vtx[i]->getTracks().size();j++){
					Track *tr = vtx[i]->getTracks()[j];
					Helix hel(tr);
					double tmin;
					TVector3 vpos = vtx[i]->getPos();
					hel.LogLikelihood(vpos, tmin);
					TVector3 trpos = hel.GetPos(tmin);
					cout << "Track #" << j << ", nearest pos: ( " << trpos.x() << ", " << trpos.y() << ", " << trpos.z() << ") with t = " << tmin << endl;
				}
			}
			cout << endl;
		}

		Vertex *priVertex = findPrimaryVertex(event->getTracks(), 25);
		cout << "Primary: (" << priVertex->getX() << ", " << priVertex->getY()  << ", " << priVertex->getZ() << ")";
		cout << ", err:(" << sqrt(priVertex->getCov()[Vertex::xx]) << ", " << sqrt(priVertex->getCov()[Vertex::yy]) << ", " << sqrt(priVertex->getCov()[Vertex::zz]) << ")" << endl;


		for(unsigned int i=0;i<pvtx->size();i++){
			Vertex *v = (*pvtx)[i];
			for(int j=0;j<2;j++){
				Track *tr = v->getTracks()[j];
				data->ipchi2[j] = lcfi.getChi2TrackVtx(priVertex, tr);
				data->trchi2[j] = v->getChi2Track(tr);
				data->tag[j] = tr->getMcp()->getFlavorTagCategory();
				MCParticle *p = tr->getMcp()->getSemiStableParent();
				if(p){
					data->truevtx[0][j] = p->getEndVertex().x();
					data->truevtx[1][j] = p->getEndVertex().y();
					data->truevtx[2][j] = p->getEndVertex().z();
				}else{
					data->truevtx[0][j] = data->truevtx[1][j] = data->truevtx[2][j] = 0;
				}
			}
			data->recovtx[0] = v->getX();
			data->recovtx[1] = v->getY();
			data->recovtx[2] = v->getZ();
			data->vchi2 = v->getChi2();

			MCParticle *p1 = v->getTracks()[0]->getMcp()->getSemiStableParent();
			MCParticle *p2 = v->getTracks()[1]->getMcp()->getSemiStableParent();

			bool sameVertex = false;
			list<MCParticle *>palist;

			MCParticle *p = p1;
			while(p){
				palist.push_back(p);
				p = p->getSemiStableParent();
			}
			p = p2;
			while(p){
				if(find(palist.begin(), palist.end(), p) != palist.end())sameVertex = true;
				p = p->getSemiStableParent();
			}

			data->success = (p1 == p2 ? 1 : sameVertex ? 2 : 0);

			tree->Fill();
		}

		cout << "BEGIN MatchMCReco ################################" << endl;
		// match MC-Reco vertex
		data2->event = nev;
		for(unsigned int i=0;i<mcv.size();i++){
			// remove mcv with recotrack < 2
			if(mcv[i]->getRecoTracks().size()<2)continue;
			TVector3 v = mcv[i]->getParent()->getEndVertex();
			if(v.Mag()<0.01 || v.Mag() > 100)continue;

			data2->recov.clear();
			data2->mcv = mcv[i];

			for(unsigned int j=0;j<pvtx->size();j++){
				Vertex *vtx = (*pvtx)[j];
				bool match = true;
				for(unsigned int k=0;k<vtx->getTracks().size();k++){
					if(find(mcv[i]->getRecoTracks().begin(), mcv[i]->getRecoTracks().end(), vtx->getTracks()[k]) == mcv[i]->getRecoTracks().end()){
						match = false;break;
					}
				}
				if(match){
					data2->recov.push_back(vtx);
					if(vtx->getChi2()>0){
						cout << "MCVertex: (" << v.x() << ", " << v.y() << ", " << v.z() << "), ";
						TrackPocaXY tp(vtx->getTracks()[0],v);
						cout << "RecoDist: " << sqrt(tp.getPoca()) << ", ";
						TrackPocaXY tp2(vtx->getTracks()[1],v);
						cout << sqrt(tp2.getPoca()) << ", chi2 = " << vtx->getChi2() << endl;
					}
				}
			}

			data2->recovsize = data2->recov.size();
			summcv ++;
			if(data2->recov.size())summcrecov ++;

			tree2->Fill();
		}
		cout << "END MatchMCReco ################################" << endl;

/*			Vertex *v = (*pvtx)[i];
			cout << v->getChi2() << ", (" << v->getX() << ", " << v->getY() << ", " << v->getZ() << "), ";
			cout << "#tracks = " << v->getTracks().size() << ", ";
			Track *tr1 = v->getTracks()[0];
			double ip1 = sqrt(pow(trackD0Significance(tr1, priVertex),2) + pow(trackZ0Significance(tr1,priVertex),2));
			Track *tr2 = v->getTracks()[1];
			double ip2 = sqrt(pow(trackD0Significance(tr2, priVertex),2) + pow(trackZ0Significance(tr2,priVertex),2));
			double mass = ((TLorentzVector)(*tr1) + (TLorentzVector)(*tr2)).M();
			cout << "M= " << mass << ", ";
			cout << tr1->getD0() << ", " << tr1->getZ0() << ", " << ip1 << ", (" << tr1->Px() << ", " << tr1->Py() << ", " << tr1->Pz() << "), ";
			cout << tr2->getD0() << ", " << tr2->getZ0() << ", " << ip2 << ", (" << tr2->Px() << ", " << tr2->Py() << ", " << tr2->Pz() << ")" << endl;;*/

		nev ++;
	}

	tree->Write();
	tree2->Write();

	cout << "SumMCV = " << summcv << endl;
	cout << "SumMCRecoV = " << summcrecov << endl;

	pf->Print();
}

void testSueharaVertex(const char *inputfile, const char *outputfile, int nStart, int nEnd)
{
	LCIOStorer ls(inputfile);
	ls.InitCollections("PandoraPFOs","MCParticlesSkimmed","RecoMCTruthLink","Tracks","Neutrals","MCParticles");
	cout << "LCIO initialization successful." << endl;

	int nev=0;
	while(ls.Next()){
		if(nev < nStart){nev++;continue;}
		if(nev >= nEnd)break;

		Event *event = Event::Instance();

		vector<Vertex *> *pvtx = lcfiplus::VertexFinderSuehara::findCandidates(event->getTracks(), 25, 10, 5);
		cout << "Event #" << nev << ": " << pvtx->size() << " vertices found." << endl;

		for(unsigned int i=0;i<pvtx->size();i++){
			Vertex *v = (*pvtx)[i];
			cout << v->getChi2() << ", " << v->getX() << ", " << v->getY() << ", " << v->getZ() << endl;
		}

		nev ++;
	}
}

Jet * JetMCMatch(vector<Jet *> &jets, MCParticle *mcp, vector<Track *> &assignedTracks, vector<Track *> &residualTracks)
{
	const vector<Track *> *pTracks;
	Event::Instance()->Get<Track>("Tracks",pTracks);

	vector<Track *> bTracks;

	vector<int> nTrackInJet;
	nTrackInJet.resize(jets.size());
	vector<int> nVertexTrackInJet;
	nVertexTrackInJet.resize(jets.size());

	int nvtx = 0;
	// get tracks
	for(unsigned int i=0;i<pTracks->size();i++){
		MCParticle *mcpc = (*pTracks)[i]->getMcp();

		if(mcpc==0)continue;
		if(mcpc->isParent(mcp)){
			bTracks.push_back((*pTracks)[i]);
			for(unsigned int j=0;j<jets.size();j++){
				for(unsigned int k=0;k<jets[j]->getVertices().size();k++){
					const vector<Track *> &vtr = jets[j]->getVertices()[k]->getTracks();
					if(find(vtr.begin(), vtr.end(), (*pTracks)[i]) != vtr.end()){
						// the track matched to this jet
						nTrackInJet[j] ++;
// 						if(nvtx == 0 && k > 0){cout << "CAUTION: vertices in the same jet might be from different semistables!" << endl;}
// 						if(nvtx > 0 && k == 0){cout << "CAUTION: vertices in different jets might be from the same samistable!" << endl;}
						nVertexTrackInJet[j] ++;
						nvtx ++;
					}
				}
				if(find(jets[j]->getTracks().begin(), jets[j]->getTracks().end(), (*pTracks)[i]) != jets[j]->getTracks().end()){
					// the track matched to this jet
					nTrackInJet[j] ++;
				}
			}
		}
	}

	int ntijMax = 0;
	int ntijMaxIndex = -1;

	int ntrsum = 0;
	int ntrvtxsum = 0;
	// determine best-match jet
	for(unsigned int j=0;j<jets.size();j++){
		if(ntijMax < nTrackInJet[j]){
			ntijMax = nTrackInJet[j];
			ntijMaxIndex = j;
		}
		ntrsum += nTrackInJet[j];
		ntrvtxsum += nVertexTrackInJet[j];
	}

	if(ntijMaxIndex == -1)
		return 0;

	vector<Track *> jetTracks = jets[ntijMaxIndex]->getTracks();
	for(unsigned int i=0;i<jets[ntijMaxIndex]->getVertices().size();i++){
		Vertex *vtx = jets[ntijMaxIndex]->getVertices()[i];
		jetTracks.insert(jetTracks.end(), vtx->getTracks().begin(), vtx->getTracks().end());
	}

	// obtain assignedtracks
	for(unsigned int i=0;i<bTracks.size();i++){
		if(find(jetTracks.begin(), jetTracks.end(), bTracks[i]) != jetTracks.end())
			assignedTracks.push_back(bTracks[i]);
		else
			residualTracks.push_back(bTracks[i]);
	}

	cout << "Assigned jet " << ntijMaxIndex << ", Vertex tracks: ";
	for(unsigned int j=0;j<jets.size();j++)cout << nVertexTrackInJet[j] << ",";
	cout << "/" << ntrvtxsum << ", all tracks: ";
	for(unsigned int j=0;j<jets.size();j++)cout << nTrackInJet[j] << ",";
	cout << "/" << bTracks.size() << ", PDG: " << mcp->getPDG() << endl;//"," ;
//	cout << " ( " << mcp->getVertex().x() << " " << mcp->getVertex().y() << " " << mcp->getVertex().z() << ")" << endl;

	return jets[ntijMaxIndex];
}

// vertex: vertex mode 0: no vertexing, 1: non-jet build-up, 2: jet build-up, 3: jet zvtop
void bbhhAnalysis3(const char *inputfile, const char *outputfile, int nStart, int nEnd, int vtx, int bbhh)
{
	const int njets = 6;


	LCIOStorer ls(inputfile);
	ls.InitCollections();

	TFile *pf = TFile::Open(outputfile,"recreate");

	TNtuple *ntJet2 = new TNtuple("ntJet2","ntJet2","nev:njet:nbjetmc:nvtx:nvtxjet:ngoodvtx:fracgoodvtxtrack:ycut:nbjet:fracgoodtrack");
	TNtuple *ntResidual = new TNtuple("ntResidual","ResidualTracks","nev:bid:btracks:mcvx:mcvy:mcvz:d0:d0err:z0:z0err:tre");

	TNtuple *nbJet = new TNtuple("nbJet", "number of b tracks in eachjet", "nev:nb1:nb2:nb3:nb4:nb5:nb6:nb11:nb12:nb13:nb14:nb15:nb16");

	int nev = -1;

	Event *event = Event::Instance();

	while(ls.Next()){
		nev ++;
		if(nev < nStart){continue;}
		if(nev >= nEnd)break;
		cout << "Event #" << nev << endl;

//		ts.GetEntry(nev);

		const vector<MCParticle *>& mcps = event->getMCParticles();
//		const vector<Track *>& tracks = event->getTracks();

		// check bbbbbb (reject H->WW etc.)
		if(bbhh){
			int hcount = 0;
			for(unsigned int i=0;i<mcps.size();i++){
				if(mcps[i]->getPDG() != 25)continue;
	
				// higgs
				if(mcps[i]->getDaughters().size() != 2){
					cout << "ERR: # of higgs daughters = " << mcps[i]->getDaughters().size() << endl;
					break;
				}
				if(abs(mcps[i]->getDaughters()[0]->getPDG()) != 5)break;
				if(abs(mcps[i]->getDaughters()[1]->getPDG()) != 5)break;
				hcount ++;
			}
			if(hcount < 2)continue;
		}

		// semistable B
		vector<MCParticle *> blist;
		for(unsigned int i=0;i<mcps.size();i++){
			if(mcps[i]->isSemiStableB()){
				blist.push_back(mcps[i]);
			}
		}
		cout << "Number of semistable B: " << blist.size() << endl;

		// vertex clustering

		TrackSelectorConfig secVtxCfg;
		// AND
		// original setup 110214
 		secVtxCfg.maxD0 = 10;
// 		secVtxCfg.maxD0Err = 0.25;
 		secVtxCfg.maxZ0 = 20;
// 		secVtxCfg.maxZ0Err = 1e10;
 		secVtxCfg.minPt = 0.1;
 		secVtxCfg.maxInnermostHitRadius = 1e10;

		// 110214
		secVtxCfg.maxD0Err = 0.1;
		secVtxCfg.maxZ0Err = 0.1;

		// OR
		secVtxCfg.minTpcHits = 20;
		secVtxCfg.minFtdHits = 3;
		secVtxCfg.minVtxHitsWithoutTpcFtd = 3;
		secVtxCfg.minVtxPlusFtdHits = 0;

		JetConfig jetCfg;
		jetCfg.nJet = njets;
		JetFinder* jetFinder = new JetFinder(jetCfg);
		double ycut = 0;

		// cut bad tracks
		const vector<Track *> &tracks = event->getTracks();
		vector<Track *> passedTracks = TrackSelector() (tracks, secVtxCfg);
		vector<Vertex *> vertices;

		LcfiInterface interface(event);

/*		for(unsigned int i=0;i<tracks.size();i++){
			if(interface.passesCut(tracks[i], secVtxCfg))
				passedTracks.push_back(tracks[i]);
		}*/
		cerr << "Track preselection: passed " << passedTracks.size() << "/" << tracks.size() << endl;

		// build up vertexing
		if(vtx == 1){
		  VertexFinderSuehara::buildUp(passedTracks, vertices, 25, 9, 10);
		  VertexFinderSuehara::associateIPTracks(vertices,1.);
		}

		if(vtx >= 2){
			vector<Jet*> prejets = jetFinder->run(tracks, event->getNeutrals(), &ycut);
			vector<Track *> passedTracksJet[njets];
			vector<Vertex *> verticesTmp[njets];

			if(vtx==2){
				Vertex *ip = findPrimaryVertex(tracks, 25);
				vertices.push_back(ip);

				for(unsigned int j=0;j<prejets.size();j++){
					passedTracksJet[j] = TrackSelector() (prejets[j]->getTracks(), secVtxCfg);
/*					for(unsigned int i=0;i<prejets[j]->getTracks().size();i++){
						if(interface.passesCut(prejets[j]->getTracks()[i],secVtxCfg))
							passedTracksJet[j].push_back(prejets[j]->getTracks()[i]);
					}*/
					VertexFinderSuehara::buildUp(passedTracksJet[j], verticesTmp[j], 25, 25,10,0.3, 1.0,ip);
				  VertexFinderSuehara::associateIPTracks(verticesTmp[j]);
					vertices.insert(vertices.end(), verticesTmp[j].begin()+1, verticesTmp[j].end());//exclude ip
				}
			}
			else if(vtx >= 3){
				SecondaryVertexConfig svc;
				svc.TrackQualityCuts = secVtxCfg;
				for(unsigned int j=0;j<prejets.size();j++){
					verticesTmp[j] = interface.findSecondaryVertices( prejets[j], svc );
					vertices.insert(vertices.end(), verticesTmp[j].begin(), verticesTmp[j].end());
				}
			}
		}

		// select vertices
		vector<Vertex *> selectedVertices;
		vector<Track *> residualTracks = tracks;
		for(vector<Vertex *>::iterator it = vertices.begin(); it != vertices.end();it++){
			Vertex *v = *it;
//			if(vtx==4 || (v->getTracks().size() >= 3 && v->getPos().Mag() < 20. && v->getPos().Mag() > 1.)
//					|| (v->getTracks().size() == 2 && v->getPos().Mag() < 20. && v->getPos().Mag() > 1. && v->getProb() > 0.1 )){
			double mass = 0.;
			if(v->getTracks().size() == 2)
				mass = (*(TLorentzVector *)(v->getTracks()[0]) + *(TLorentzVector *)(v->getTracks()[1])).M();
			if((mass < .488 || mass > .508) && v->getPos().Mag() < 30 && v->getPos().Mag() > 0.3){
				selectedVertices.push_back(v);
				for(vector<Track *>::const_iterator itt = v->getTracks().begin(); itt != v->getTracks().end(); itt++){
					vector<Track *>::iterator itt2 = remove_if(residualTracks.begin(), residualTracks.end(), bind2nd(equal_to<Track *>(), *itt));
					residualTracks.erase(itt2, residualTracks.end());
				}
			}
		}
		cout << "Selected # vertices = " << selectedVertices.size() << ", residual # tracks = " << residualTracks.size() << "/" << tracks.size() << endl;

//	TNtuple *ntResidual = new TNtuple("ntResidual","ResidualTracks","nev:bid:btracks:mcvx:mcvy:mcvz:d0:d0err:z0:z0err:tre");
		// calculate btracks
		int *btracks = new int[blist.size()];
		memset(btracks,0,sizeof(int)*blist.size());
		for(unsigned int i=0;i<tracks.size();i++){
			for(unsigned int k=0;k<blist.size();k++){
				if(tracks[i]->getMcp()->isParent(blist[k]))btracks[k] ++;
			}
		}

		for(unsigned int i=0;i<residualTracks.size();i++){
			Track *tr = residualTracks[i];
			MCParticle *mcp = tr->getMcp();
			TVector3 vtx = mcp->getVertex();
			for(unsigned int k=0;k<blist.size();k++){
				if(mcp->isParent(blist[k])) {
					ntResidual->Fill(nev,blist[k]->getId(), btracks[k], vtx.x(), vtx.y(), vtx.z(),
						tr->getD0(), sqrt(tr->getCovMatrix()[tpar::d0d0]), tr->getZ0(), sqrt(tr->getCovMatrix()[tpar::z0z0]), tr->E());
					cout << "Residual tracks: " << blist[k]->getId() << ",mc(" << vtx.x() << " " << vtx.y() << " " << vtx.z() << "),reco:(";
					cout << tr->E() << " " << tr->getD0() << " " << tr->getZ0() << ")" << endl;
					break;
				}
			}
		}

		// check vertex purity
		int goodvertex = 0;
		int goodtrack = 0;
		int allvtxtrack = 0;
		for(unsigned int i=0;i<selectedVertices.size();i++){
			int gtv = 0;
			for(unsigned int k=0;k<selectedVertices[i]->getTracks().size();k++){
				Track *tr = selectedVertices[i]->getTracks()[k];
				allvtxtrack ++;
				for(unsigned int ib=0;ib<blist.size();ib++){
					if(tr->getMcp()->isParent(blist[ib])){
						gtv ++;break;
					}
				}
			}
			if(gtv * 2 >= (int)selectedVertices[i]->getTracks().size())goodvertex ++;
			goodtrack += gtv;
		}
		double fracgoodvtxtrack = double(goodtrack)/double(allvtxtrack);

		// Jet clustering
		vector<Jet*> jets = jetFinder->run(residualTracks, event->getNeutrals(), selectedVertices, &ycut, true);
//		vector<Jet*> njets = jetFinder->run(tracks, event->getNeutrals(), &ymin);

		residualTracks.clear();
		vector<Track *> assignedTracks;

		map<Jet *, int > nbmap;
		int nvtxjet = 0;
		for(unsigned int nj=0;nj<jets.size();nj++){
			nbmap[jets[nj]] = 0;
			if(jets[nj]->getVertices().size()>0)nvtxjet ++;
		}

		// nbjet fill
		vector<int> nbjet0(max(int(jets.size()),6));
		vector<int> nbjet1(max(int(jets.size()),6));
		for(unsigned int nj=0;nj<jets.size();nj++){
			nbjet0[nj] = nbjet1[nj] = 0;

			for(unsigned int i=0;i<jets[nj]->getTracks().size();i++){
				Track *tr = jets[nj]->getTracks()[i];
				if(tr->getMcp() == 0)continue;
				if(tr->getMcp()->getSemiStableParent() == 0)continue;
				int pdg = tr->getMcp()->getSemiStableParent()->getPDG();
				if((abs(pdg)>400&&abs(pdg)<600) || (abs(pdg)>4000&&abs(pdg)<6000)){
					for(unsigned int k=0;k<blist.size();k++){
						if(tr->getMcp()->isParent(blist[k])){
							nbjet0[nj] ++;
							if(tr->E()>1.)nbjet1[nj] ++;
							break;
						}
					}
				}
			}
			for(unsigned int nv=0;nv<jets[nj]->getVertices().size();nv++){
				for(unsigned int i=0;i<jets[nj]->getVertices()[nv]->getTracks().size();i++){
					Track *tr = jets[nj]->getVertices()[nv]->getTracks()[i];
					if(tr->getMcp() == 0)continue;
					if(tr->getMcp()->getSemiStableParent() == 0)continue;
					int pdg = tr->getMcp()->getSemiStableParent()->getPDG();
					if((abs(pdg)>400&&abs(pdg)<600) || (abs(pdg)>4000&&abs(pdg)<6000)){
						for(unsigned int k=0;k<blist.size();k++){
							if(tr->getMcp()->isParent(blist[k])){
								nbjet0[nj] ++;
								if(tr->E()>1.)nbjet1[nj] ++;
								break;
							}
						}
					}
				}
			}
		}

		// sort by decending order
		sort(nbjet0.begin(), nbjet0.end(), greater<int>());
		sort(nbjet1.begin(), nbjet1.end(), greater<int>());
		nbJet->Fill(nev, nbjet0[0], nbjet0[1], nbjet0[2], nbjet0[3], nbjet0[4], nbjet0[5],
										 nbjet1[0], nbjet1[1], nbjet1[2], nbjet1[3], nbjet0[4], nbjet0[5]);

		for(unsigned int ib=0;ib<blist.size();ib++){
			vector<Track *> aTracks, rTracks;
//Jet * JetMCMatch(vector<Jet *> &jets, MCParticle *mcp, vector<Track *> &assignedTracks, vector<Track *> &residualTracks)
			Jet *jet = JetMCMatch(jets, blist[ib], aTracks, rTracks);
			if(jet){
				nbmap[jet] ++;
				assignedTracks.insert(assignedTracks.end(), aTracks.begin(), aTracks.end());
				residualTracks.insert(residualTracks.end(), rTracks.begin(), rTracks.end());
			}
		}

		int nbjet = 0;
		for(unsigned int nj=0;nj<jets.size();nj++){
			if(nbmap[jets[nj]] > 0)nbjet ++;
		}
//	TNtuple *ntJet2 = new TNtuple("ntJet2","ntJet2","nev:njet:nbjetmc:nvtx:nvtxjet:ngoodvtx:ngoodvtxtrack:ycut:nbjet:fracgoodtrack");
//	TNtuple *ntJet2 = new TNtuple("ntJet2","ntJet2","nev:njet:ycut:nbjet:fracgoodtrack");
		double fracgoodtrack = double(residualTracks.size()) / double(assignedTracks.size()+residualTracks.size());
		ntJet2->Fill(nev, jets.size(), blist.size(), selectedVertices.size(), nvtxjet, goodvertex, fracgoodvtxtrack, ycut, nbjet, fracgoodtrack);
	}

	ntJet2->Write();
	ntResidual->Write();
	nbJet->Write();
	pf->Close();
}


void bbhhAnalysis2(const char *inputfile, const char *outputfile, int nStart, int nEnd)
{
	LCIOStorer ls(inputfile);
	ls.InitCollections();

	TFile *pf = TFile::Open(outputfile,"recreate");

	TNtuple *ntJet2 = new TNtuple("ntJet2","ntJet2","nev:njet:ycut:nbjet:fracgoodtrack");
	TNtuple *ntJet22 = new TNtuple("ntJet22","ntJet22","nev:njet:ycut:nbjet:fracgoodtrack");

	int nev = -1;

	while(ls.Next()){
		nev ++;
		if(nev < nStart){continue;}
		if(nev >= nEnd)break;
		cout << "Event #" << nev << endl;

//		ts.GetEntry(nev);
		Event *event = Event::Instance();

		const vector<MCParticle *>& mcps = event->getMCParticles();
//		const vector<Track *>& tracks = event->getTracks();

		// check bbbbbb (reject H->WW etc.)
		int hcount = 0;
		for(unsigned int i=0;i<mcps.size();i++){
			if(mcps[i]->getPDG() != 25)continue;

			// higgs
			if(mcps[i]->getDaughters().size() != 2){
				cout << "ERR: # of higgs daughters = " << mcps[i]->getDaughters().size() << endl;
				break;
			}
			if(abs(mcps[i]->getDaughters()[0]->getPDG()) != 5)break;
			if(abs(mcps[i]->getDaughters()[1]->getPDG()) != 5)break;
			hcount ++;
		}
		if(hcount < 2)continue;

		// semistable B
		vector<MCParticle *> blist;
		for(unsigned int i=0;i<mcps.size();i++){
			if(mcps[i]->isSemiStableB()){
				blist.push_back(mcps[i]);
			}
		}
		cout << "Number of semistable B: " << blist.size() << endl;


		int prevnjet = 0;
		double prevycut = 0;
		int prevnbjet = 0;
		double prevfracgoodtrack = 0;

		// jet clustering
		for(int njet = 6;njet<=10;njet++){
			// normal 6-jet clustering
			JetConfig jetCfg;
			jetCfg.nJet = njet;
			JetFinder* jetFinder = new JetFinder(jetCfg);
			double ycut;

			vector<Jet*> jets = jetFinder->run(event->getTracks(),event->getNeutrals(),&ycut);

			vector<Track *> assignedTracks, residualTracks;
			map<Jet *, int > nbmap;
			for(unsigned int nj=0;nj<jets.size();nj++){
				nbmap[jets[nj]] = 0;
			}

			for(unsigned int ib=0;ib<blist.size();ib++){
				vector<Track *> aTracks, rTracks;
//Jet * JetMCMatch(vector<Jet *> &jets, MCParticle *mcp, vector<Track *> &assignedTracks, vector<Track *> &residualTracks)
				Jet *jet = JetMCMatch(jets, blist[ib], aTracks, rTracks);
				if(jet){
					nbmap[jet] ++;
					assignedTracks.insert(assignedTracks.end(), aTracks.begin(), aTracks.end());
					residualTracks.insert(residualTracks.end(), rTracks.begin(), rTracks.end());
				}
			}

			int nbjet = 0;
			for(unsigned int nj=0;nj<jets.size();nj++){
				if(nbmap[jets[nj]] > 0)nbjet ++;
			}

//	TNtuple *ntJet2 = new TNtuple("ntJet2","ntJet2","nev:njet:ycut:nbjet:fracgoodtrack");
			double fracgoodtrack = double(residualTracks.size()) / double(assignedTracks.size()+residualTracks.size());
			ntJet2->Fill(nev, jets.size(), ycut, nbjet, fracgoodtrack);

			if(jets.size()>6 && prevnbjet == 5 && nbjet == 6){
				ntJet22->Fill(nev, prevnjet, prevycut, prevnbjet, prevfracgoodtrack);
				ntJet22->Fill(nev, jets.size(), ycut, nbjet, fracgoodtrack);
			}

			prevnjet = jets.size();
			prevycut = ycut;
			prevnbjet = nbjet;
			prevfracgoodtrack = fracgoodtrack;
		}
	}
	ntJet2->Write();
	ntJet22->Write();
	pf->Print();
}

void bbhhAnalysis(const char *inputfile, const char *outputfile, int nStart, int nEnd)
{
/*	TreeStorer ts(inputfile,"FlavTagTree", TreeStorer::mode_input);
	ts.RegisterAll();
	cout << ts.GetEntries() << endl;

	for(int nev=0;nev < ts.GetEntries();nev++){*/

	bool zvtop = false;

	LCIOStorer ls(inputfile);
	ls.InitCollections();

	TFile *pf = TFile::Open(outputfile,"recreate");
	TNtuple *ntMissTracks = new TNtuple("ntMissTracks","ntMissTracks","px:py:pz:p:d0:z0");
	TNtuple *ntAllTracks = new TNtuple("ntAllTracks","ntAllTracks","px:py:pz:p:d0:z0");

	TNtuple *ntVertex = new TNtuple("ntVertex","ntVertex","nev:njet:nall:pjet:pall:vjettd:valltd:vjetzv:vallzv");
	TNtuple *ntJet = new TNtuple("ntJet","ntJet","nev:e:nvtx");

	int nev = -1;

	int sumvreco = 0, sumvall = 0;
	int sumvreco2 = 0, sumvall2 = 0;
	int sumtreco = 0, sumtall = 0;
	int sumdreco = 0, sumdall = 0;

	int sumee = 0;
	int summm = 0;
	int sumtt = 0;

	while(ls.Next()){
		nev ++;
		if(nev < nStart){continue;}
		if(nev >= nEnd)break;
		cout << "Event #" << nev << endl;

//		ts.GetEntry(nev);
		Event *event = Event::Instance();

		const vector<MCParticle *>& mcps = event->getMCParticles();
		const vector<Track *>& tracks = event->getTracks();
		// check bbbbbb (reject H->WW etc.)
		int hcount = 0;
		for(unsigned int i=0;i<mcps.size();i++){
			if(mcps[i]->getPDG() != 25)continue;

			// higgs
			if(mcps[i]->getDaughters().size() != 2){
				cout << "ERR: # of higgs daughters = " << mcps[i]->getDaughters().size() << endl;
				break;
			}
			if(abs(mcps[i]->getDaughters()[0]->getPDG()) != 5)break;
			if(abs(mcps[i]->getDaughters()[1]->getPDG()) != 5)break;
			hcount ++;
		}
		if(hcount < 2)continue;

		for(unsigned int i=0;i<mcps.size();i++){
			if(!mcps[i]->isSemiStableB())continue;

			for(unsigned int j=0;j<mcps[i]->getDaughters().size();j++){
				if(abs(mcps[i]->getDaughters()[j]->getPDG())==11){sumee++;break;}
				if(abs(mcps[i]->getDaughters()[j]->getPDG())==13){summm++;break;}
				if(abs(mcps[i]->getDaughters()[j]->getPDG())==15){sumtt++;break;}
			}
		}

		// HH->bbbb selected
// 		vector<MCVertex *> mcv = findMcVertex(mcps);
// 		map<MCVertex*,int> table;
// 		matchMcVertex(evt, mcv, table, false);

		// normal 6-jet clustering
		JetConfig jetCfg;
		jetCfg.nJet = 6;
		JetFinder* jetFinder = new JetFinder(jetCfg);
		vector<Jet*> jets = jetFinder->run(event->getTracks(),event->getNeutrals());

		// semistable B
		vector<MCParticle *> blist;
		for(unsigned int i=0;i<mcps.size();i++){
			if(mcps[i]->isSemiStableB()){
				blist.push_back(mcps[i]);
			}
		}
		cout << "Number of semistable B: " << blist.size() << endl;

		// jet-track match
		vector<pair<Track *, int> > *jtmatch = new vector<pair<Track *, int> >[blist.size()];
		for(unsigned int i=0;i<tracks.size();i++){
			MCParticle *mcp = tracks[i]->getMcp();
			int pcount = 0;
			for(unsigned int j=0;j<blist.size();j++){
				if(mcp->isParent(blist[j])){
					// mc-B found, search jet
					for(unsigned int k=0;k<6;k++){
						if(find(jets[k]->getTracks().begin(), jets[k]->getTracks().end(), tracks[i]) != jets[k]->getTracks().end()){
							jtmatch[j].push_back(pair<Track *, int>(tracks[i], k));
							break;
						}
					}
					pcount ++;
				}
			}
			if(pcount > 1)cout << "A track belonging to multiple Bs has been found!" << endl;
		}

		// jet-track match OK
		unsigned int nfrac[6];
		double pfrac[6];
		int nbjet[6] = {0,0,0,0,0,0};

		LcfiInterface lcfi;
		for(unsigned int i=0;i<blist.size();i++){ //i: MC B index
			memset(nfrac,0,sizeof(nfrac));
			memset(pfrac,0,sizeof(pfrac));

			unsigned int nfracmax = 0;
			double pfracmax = 0;
			double ptotal = 0;
			int nfracmaxindex = -1;

			for(unsigned int j=0;j<jtmatch[i].size();j++){ //j: jtmatch index
				pair<Track *, int> pair = jtmatch[i][j];
				int ijet = pair.second;
				nfrac[ijet] ++;
				pfrac[ijet] += pair.first->P();
				ptotal += pair.first->P();

				if(nfrac[ijet] > nfracmax){nfracmax = nfrac[ijet]; nfracmaxindex = ijet;}
				if(pfrac[ijet] > pfracmax)pfracmax = pfrac[ijet];
			}
			// summarize efrac/nfrac
			cout << "B #" << i << ", tracks " << nfracmax << "/" << jtmatch[i].size() << ", momentum " << pfracmax << "/" << ptotal;

			nbjet[nfracmaxindex]++;

			// run ZVTOP with b-origin tracks
			// formulate jets
			Jet jreco, jall;
			for(unsigned int j=0;j<jtmatch[i].size();j++){ //j: jtmatch index
				Track *tr = jtmatch[i][j].first;
				if(jtmatch[i][j].second == nfracmaxindex){ // in the main jet
					jreco.add(jtmatch[i][j].first);
				}
				else{
					ntMissTracks->Fill(tr->Px(), tr->Py(), tr->Pz(), tr->P(), tr->getD0(), tr->getZ0());
				}
				jall.add(jtmatch[i][j].first);
				ntAllTracks->Fill(tr->Px(), tr->Py(), tr->Pz(), tr->P(), tr->getD0(), tr->getZ0());
			}
			// run ZVTOP
			vector<Vertex *> *pvreco = findTearDownVertices(*event, jreco);
			vector<Vertex *> *pvall = findTearDownVertices(*event, jall);
			cout << ", recoVertex: " << pvreco->size() << "/" << pvall->size();// << endl;

			sumvreco += pvreco->size();
			sumvall += pvall->size();
			sumtreco += nfracmax;
			sumtall += jtmatch[i].size();

			if(nfracmax != jtmatch[i].size())sumdreco ++;
			sumdall ++;

			vector<Vertex *> vreco, vall;
			if(zvtop){
				vreco = lcfi.forceZvtop(jreco);
				vall = lcfi.forceZvtop(jall);

				cout << ", recoVertex: " << vreco.size() << "/" << vall.size() << endl;
			}
			else{
				Vertex *vtx = 0;
				list<Track *> trackList;
				trackList.resize(jreco.getTracks().size());
				copy(jreco.getTracks().begin(), jreco.getTracks().end(), trackList.begin());
				do{
					vtx = VertexFinderSuehara::findOne(trackList, 25, 10 , true);
					if(vtx)vreco.push_back(vtx);
				}while(vtx && trackList.size()>1);

				trackList.clear();
				trackList.resize(jall.getTracks().size());
				copy(jall.getTracks().begin(), jall.getTracks().end(), trackList.begin());
				do{
					vtx = VertexFinderSuehara::findOne(trackList, 25, 10 , true);
					if(vtx)vall.push_back(vtx);
				}while(vtx && trackList.size()>1);
				cout << ", recoVertex: " << vreco.size() << "/" << vall.size() << endl;
			}
			sumvreco2 += vreco.size();
			sumvall2 += vall.size();
	
//			TNtuple *ntVertex = new TNtuple("ntVertex","ntVertex","nev:njet:nall:pjet:pall:vjettd:valltd:vjetzv:vallzv");
			ntVertex->Fill(nev, nfracmax, jtmatch[i].size(), pfracmax, ptotal, pvreco->size(), pvall->size(), vreco.size(), vall.size());

		}

		for(int ij=0;ij<6;ij++){
			ntJet->Fill(nev, jets[ij]->E(), nbjet[ij]);
		}

		delete[] jtmatch;
	}

	ntMissTracks->Write();
	ntAllTracks->Write();
	ntVertex->Write();
	ntJet->Write();
	pf->Print();

	cout << "Sum # track Reco: " << sumtreco << ", Sum # track All: " << sumtall << endl;
	cout << "Sum # diffjet Reco: " << sumdreco << ", Sum # diffjet All: " << sumdall << endl;
	cout << "Sum # vertex TD Reco: " << sumvreco << ", Sum # vertex All: " << sumvall << endl;
	cout << "Sum # vertex BU Reco: " << sumvreco2 << ", Sum # vertex All: " << sumvall2 << endl;

	cout << "Sum e: " << sumee << ", Sum mu: " << summm << ", Sum tau: " << sumtt << endl;
}

void pointTest()
{
	lcfiplus::Point::SVector3 pos(0,0,0);
	lcfiplus::Point::SMatrixSym3 err;
	err(0,0) = 1;
	err(1,1) = 1;
	err(2,2) = 1;
	err(0,1) = -0.999;
	err(0,2) = 0;
	err(1,2) = 0;

	lcfiplus::Point pt(pos,err);
	TVector3 v(1,1,0);

	pt.LogLikelihood(v);
}

void helixTest()
{
	lcfiplus::Helix::SVector5 track(1000,0,0,0.001,1);
	lcfiplus::Helix::SVector5 track2(0,0,TMath::Pi()/2,0.01,1);
	lcfiplus::Helix::SMatrixSym5 err;

	lcfiplus::Helix hel(track,err,1);
	lcfiplus::Helix hel2(track2,err,-1);

	cout << "Track 1" << endl;
	for(int t=0;t<10;t++){
		TVector3 p = hel.GetPos(0.01*t);
		cout << "t: " << t << ", pos: (" << p.X() << ", " << p.Y() << ", " << p.Z() << ")" << endl;
	}
	cout << "Track 2" << endl;
	for(int t=0;t<10;t++){
		TVector3 p = hel2.GetPos(0.01*t);
		cout << "t: " << t << ", pos: (" << p.X() << ", " << p.Y() << ", " << p.Z() << ")" << endl;
	}

	cout << "Track 1-2" << endl;
	hel.ClosePoint(hel2);
	cout << "Track 2-1" << endl;
	hel2.ClosePoint(hel);
}

void helixVarianceTest()
{
	// d0, z0, phi0, omega, tandip
	lcfiplus::Helix::SVector5 track1( 1, 0.001,             0, -0.0002, -0.0001);
	lcfiplus::Helix::SVector5 track2(-1,-0.001, TMath::Pi()/2,  0.0001,  0.0001);
	lcfiplus::Helix::SMatrixSym5 err;
	err(0,0) = 10;
	err(1,1) = 1e-3;
	err(2,2) = 1e-3;
	err(3,3) = 1e-5;
	err(4,4) = 1e-2;

	lcfiplus::Helix hel1(track1,err,1);
	lcfiplus::Helix hel2(track2,err,-1);

	cout << "Track 1" << endl;
	for(int t=-10;t<10;t++){
		TVector3 p = hel1.GetPos(0.01*t);
		cout << "t: " << t << ", pos: (" << p.X() << ", " << p.Y() << ", " << p.Z() << ")" << endl;
	}
	cout << "Track 2" << endl;
	for(int t=-10;t<10;t++){
		TVector3 p = hel2.GetPos(0.01*t);
		cout << "t: " << t << ", pos: (" << p.X() << ", " << p.Y() << ", " << p.Z() << ")" << endl;
	}

	/*
	cout << "Track 1-2" << endl;
	hel1.ClosePoint(hel2);
	cout << "Track 2-1" << endl;
	hel2.ClosePoint(hel1);
	*/

	TVector3 v(1,1,0);
	double t=0;

	double var = hel1.LogLikelihood(v,t);
	printf("LL = %f @ t=%f\n",var,t);
}

double CompD0(Track *t1, Track *t2)
{
	return fabs(t1->getD0() / sqrt(t1->getCovMatrix()[tpar::d0d0])) > fabs(t2->getD0() / sqrt(t2->getCovMatrix()[tpar::d0d0]));
}

void simpleAnalysis(const char *inputfile, const char *outputfile, int nStart, int nEnd)
{
	LCIOStorer ls(inputfile);
	ls.InitCollections();

	int nev = -1;
	while(ls.Next()){
		nev ++;
		if(nev < nStart){continue;}
		if(nev >= nEnd)break;
		cout << "Event #" << nev << endl;

		Event *event = Event::Instance();
//		GeometryHandler *gh = GeometryHandler::Instance();

//		const vector<MCParticle *>& mcps = event->getMCParticles();
		vector<Track *> tracks = event->getTracks();

		// sort by D0
		sort(tracks.begin(),tracks.end(), CompD0);

		cout << "Track MC-pfoid relation -------------------------" << endl;
		// track mc-pfoid relation
		for(unsigned int i=0;i<tracks.size();i++){
			cout << tracks[i]->getPDG() << ", " << tracks[i]->getMcp()->getPDG() << ", " << tracks[i]->E() << endl;
		}
#if 0
		lcfiplus::Helix hel0(tracks[0]);
		for(unsigned int i=0;i<tracks.size();i++){
			lcfiplus::Helix hel(tracks[i]);

			cout << "Track #" << i << " ##########################" << endl;
			cout << "Track parameters: d0=" << tracks[i]->getD0() << ", z0=" << tracks[i]->getZ0()
						<< ", phi0=" << tracks[i]->getPhi() << ", omega=" << tracks[i]->getOmega() << ", tanlambda=" << tracks[i]->getTanLambda() << endl;
			const float *cov = tracks[i]->getCovMatrix();
			cout << "Track errors: " << endl;
			cout << cov[tpar::d0d0] << " " << cov[tpar::d0z0] << " " << cov[tpar::d0ph] << " " << cov[tpar::d0om] << " " << cov[tpar::d0td] << endl;
			cout << cov[tpar::d0z0] << " " << cov[tpar::z0z0] << " " << cov[tpar::z0ph] << " " << cov[tpar::z0om] << " " << cov[tpar::z0td] << endl;
			cout << cov[tpar::d0ph] << " " << cov[tpar::z0ph] << " " << cov[tpar::phph] << " " << cov[tpar::phom] << " " << cov[tpar::phtd] << endl;
			cout << cov[tpar::d0om] << " " << cov[tpar::z0om] << " " << cov[tpar::phom] << " " << cov[tpar::omom] << " " << cov[tpar::omtd] << endl;
			cout << cov[tpar::d0td] << " " << cov[tpar::z0td] << " " << cov[tpar::phtd] << " " << cov[tpar::omtd] << " " << cov[tpar::tdtd] << endl;

			TVector3 v1(0,0,0);
			TVector3 v2 = hel.GetPos(10);
			TVector3 v3 = hel.GetPos(20);
			v3.SetX(v3.x()+1);

			double d1 = hel.LogLikelihood(v1);
			double d2 = hel.LogLikelihood(v2);
			double d3 = hel.LogLikelihood(v3);
			cout << "Variances: " << d1 << ", " << d2 << ", " << d3 << endl;

			if(i>0){
				vector<PointBase *> vp;
				vector<Helix *> vh;
				vp.push_back(&hel0);
				vp.push_back(&hel);
				vh.push_back(&hel0);
				vh.push_back(&hel);
				TVector3 ip(0,0,0);

				gh->HelixPointFit(vh,0);
				gh->PointFit(vp,ip,0);
			}
/*
			for(int t=0;t<100;t++){
				Helix::SMatrixSym3 errXyz;
				Helix::SVector3 xyz;
				hel.GetPosErr(0.01 * t, xyz, errXyz);

				cout << "t: " << t << ", pos: (" << xyz(0) << ", " << xyz(1) << ", " << xyz(2) << ")" << endl;
				for(int j=0;j<3;j++){
					cout << errXyz(j,0) << " " << errXyz(j,1) << " " << errXyz(j,2) << endl;
				}
			}
*/
		}
#endif
	}
}

void VertexAnalysis110214(const char *inputfile, const char *outputfile, int nStart, int nEnd)
{
	LCIOStorer ls(inputfile);
	ls.InitCollections();

	TFile *pf = TFile::Open(outputfile,"recreate");

	TNtuple *ntVertex = new TNtuple("ntVertex","ntVertex","nev:nb:pdg:pdg2:good:ntr:pos:prob:sume:sump:mass:trm0:trm1:mce:mind0:minz0:mind0sig:minz0sig:maxmoverp");
	TNtuple *ntVcon = new TNtuple("ntVcon","ntVertexConnection","n1:n2:good1:good2:same:angle:e1:e2");

	int nev = -1;

	while(ls.Next()){
		nev ++;
		if(nev < nStart){continue;}
		if(nev >= nEnd)break;
		cout << "Event #" << nev << endl;

//		ts.GetEntry(nev);
		Event *event = Event::Instance();

		const vector<MCParticle *>& mcps = event->getMCParticles();
//		const vector<Track *>& tracks = event->getTracks();

		// check bbbbbb (reject H->WW etc.)
		int hcount = 0;
		for(unsigned int i=0;i<mcps.size();i++){
			if(mcps[i]->getPDG() != 25)continue;

			// higgs
			if(mcps[i]->getDaughters().size() != 2){
				cout << "ERR: # of higgs daughters = " << mcps[i]->getDaughters().size() << endl;
				break;
			}
			if(abs(mcps[i]->getDaughters()[0]->getPDG()) != 5)break;
			if(abs(mcps[i]->getDaughters()[1]->getPDG()) != 5)break;
			hcount ++;
		}
		if(hcount < 2)continue;

		// semistable B
		cout << "Semistable B: " << endl;
		vector<MCParticle *> blist;
		for(unsigned int i=0;i<mcps.size();i++){
			if(mcps[i]->isSemiStableB()){
				blist.push_back(mcps[i]);
				cout << mcps[i]->getId() << ", " << mcps[i]->getPDG() << ", " << mcps[i]->E() << ", (" << mcps[i]->getEx() << ", " << mcps[i]->getEy() << ", " << mcps[i]->getEz() << ") " << endl;
			}
		}

		// vertex clustering

		TrackSelectorConfig secVtxCfg;
/*		secVtxCfg.TwoProngCut = 10;
		secVtxCfg.TrackTrimCut = 10;
		secVtxCfg.ResolverCut = 0.6;*/
		//secVtxCfg.TwoProngCut = 10;
		//secVtxCfg.TrackTrimCut = 20;
		//secVtxCfg.ResolverCut = 20.;
		// AND
		// original setup 110214
 		secVtxCfg.maxD0 = 10;
// 		secVtxCfg.maxD0Err = 0.25;
 		secVtxCfg.maxZ0 = 20;
// 		secVtxCfg.maxZ0Err = 1e10;
 		secVtxCfg.minPt = 0.1;
 		secVtxCfg.maxInnermostHitRadius = 1e10;

		// 110214
		secVtxCfg.maxD0Err = 0.1;
		secVtxCfg.maxZ0Err = 0.1;

		// OR
		secVtxCfg.minTpcHits = 20;
		secVtxCfg.minFtdHits = 3;
		secVtxCfg.minVtxHitsWithoutTpcFtd = 3;
		secVtxCfg.minVtxPlusFtdHits = 0;

		// cut bad tracks
		const vector<Track *> &tracks = event->getTracks();
		vector<Track *> passedTracks = TrackSelector() (tracks, secVtxCfg);
		vector<Vertex *> vertices;

/*		for(unsigned int i=0;i<tracks.size();i++){
			if(interface.passesCut(tracks[i], secVtxCfg))
				passedTracks.push_back(tracks[i]);
		}*/
		cout << "Track preselection: passed " << passedTracks.size() << "/" << tracks.size() << endl;

		// passed tracks
		for(unsigned int k=0;k<blist.size();k++){
			MCParticle *mcp = blist[k];
			int nt=0, ntp = 0;
			for(unsigned int i=0;i<tracks.size();i++){
				if(tracks[i]->getMcp()->isParent(mcp))nt++;
			}
			for(unsigned int i=0;i<passedTracks.size();i++){
				if(passedTracks[i]->getMcp()->isParent(mcp)){
					ntp++;

					cout << "( " << passedTracks[i]->getMcp()->getPDG() << ", " << passedTracks[i]->P() << ", " << passedTracks[i]->getD0() << ", " << passedTracks[i]->getZ0() << "), ";
				}
			}
			cout << ntp << "/" << nt << endl;
		}

		// build up vertexing
		VertexFinderSuehara::buildUp(passedTracks, vertices, 25, 9, 10);
		VertexFinderSuehara::associateIPTracks(vertices,1.);

		cout << "Vertices: " << vertices.size() << endl;

		vector<int> nb(vertices.size());
		vector<bool> good(vertices.size());
		vector<double> energy(vertices.size());

		for(unsigned int i=1;i<vertices.size();i++){
			TVector3 pos = vertices[i]->getPos();
			cout << "( " << pos.x() << " " << pos.y() << " " << pos.z() << " )" << endl;
		}

		for(unsigned int i=1;i<vertices.size();i++){
			nb[i] = -1;
			int pdg = 0;
			int pdg2 = 0;
			good[i] = true;
			int ntr = vertices[i]->getTracks().size();
			double pos = vertices[i]->getPos().Mag();
			double prob = vertices[i]->getProb();
			TLorentzVector v;
			double mce = 0.;
			double mind0 = 100.;
			double minz0 = 100.;
			double mind0sig = 100.;
			double minz0sig = 100.;
			double maxmoverp = 0.;
			double trmass[2] = {0., 0.};

			for(unsigned int j=0;j<vertices[i]->getTracks().size();j++){
				v += *(vertices[i]->getTracks()[j]);
			}
			for(unsigned int j=0;j<vertices[i]->getTracks().size();j++){
				Track *tr = vertices[i]->getTracks()[j];
				if(fabs(tr->getD0()) < mind0)mind0 = fabs(tr->getD0());
				if(fabs(tr->getZ0()) < minz0)minz0 = fabs(tr->getZ0());
				double d0sig = fabs(tr->getD0() / sqrt(tr->getCovMatrix()[tpar::d0d0]));
				double z0sig = fabs(tr->getZ0() / sqrt(tr->getCovMatrix()[tpar::z0z0]));
				if(d0sig < mind0sig)mind0sig = d0sig;
				if(z0sig < minz0sig)minz0sig = z0sig;
				double difm = v.M() - (v - *tr).M();
				double difp = v.P() - (v - *tr).P();
				double moverp = difm/difp;
				if(maxmoverp < moverp)maxmoverp = moverp;
				if(j<2)trmass[j] = tr->M();

				MCParticle *mcp = tr->getMcp();

				if(pdg2 == 0 && mcp->getSemiStableParent())pdg2 = mcp->getSemiStableParent()->getPDG();

				for(unsigned int k=0;k<blist.size();k++){
					if(mcp->isParent(blist[k])){
						if(nb[i] == -1){
							nb[i] = k;
							pdg = blist[k]->getPDG();
							mce = blist[k]->E();
						}
						else if(nb[i] != (int)k){
							good[i] = false;
						}

						cout << blist[k]->getId() << " ";
						break;
					}
					if(k==blist.size()-1){
						good[i] = false;
						cout << "0 ";
					}
				}
			}
			cout << ", ";

//	TNtuple *ntVertex = new TNtuple("ntVertex","ntVertex","nev:nb:pdg:good:ntr:pos:prob:sume:sump:mass:mce:mind0:minz0:mind0sig:minz0sig");
			float f[] = {nev, nb[i], pdg, pdg2, good[i], ntr, pos, prob, v.E(), v.P(), v.M(), trmass[0], trmass[1], mce, mind0, minz0, mind0sig, minz0sig, maxmoverp};
			ntVertex->Fill(f);

			energy[i] = v.E();
		}
		cout << endl;

		for(unsigned int i=1;i<vertices.size();i++){
			Vertex *v = vertices[i];
			double mass = 0.;
			if(v->getTracks().size() == 2)
				mass = (*(TLorentzVector *)(v->getTracks()[0]) + *(TLorentzVector *)(v->getTracks()[1])).M();
			if((mass < .488 || mass > .508) && v->getPos().Mag() < 30 && v->getPos().Mag() > 0.3)
			{}else{cout << "Vertex " << i << " is thown away. mass = " << mass << ", dist = " << v->getPos().Mag() << endl;
			}
		}

		// vertex connection
		for(unsigned int i=1;i<vertices.size()-1;i++){
			for(unsigned int j=i+1;j<vertices.size();j++){
				// opening angle
				TVector3 pos1 = vertices[i]->getPos();
				TVector3 pos2 = vertices[j]->getPos();

				double angle = pos1.Angle(pos2);
				ntVcon->Fill(i,j, good[i], good[j], nb[i]==nb[j], angle, energy[i], energy[j]);
			}
		}
	}
	pf->Write();

}

void VertexAnalysis110215(const char *inputfile, const char *outputfile, int nStart, int nEnd, int hh)
{
	LCIOStorer ls(inputfile);
	ls.InitCollections();

	TFile *pf = TFile::Open(outputfile,"recreate");
	TNtuple *ntTrack = new TNtuple("ntTrack","ntTrack","pdg:x:y:z:gamma:ctau");
	TNtuple *ntLepton = new TNtuple("ntLepton","ntLepton","pdg:parpdg:parpdg2:tre:px:py:pz:d0:z0:d0err:z0err:calotot:ecal:hcal:yoke:chi2");

	int nev = -1;

	while(ls.Next()){
		nev ++;
		if(nev < nStart){continue;}
		if(nev >= nEnd)break;
		cout << "Event #" << nev << endl;

//		ts.GetEntry(nev);
		Event *event = Event::Instance();

		const vector<MCParticle *>& mcps = event->getMCParticles();
//		const vector<Track *>& tracks = event->getTracks();

		// check bbbbbb (reject H->WW etc.)
		int hcount = 0;
		for(unsigned int i=0;i<mcps.size();i++){
			if(mcps[i]->getPDG() != 25)continue;

			// higgs
			if(mcps[i]->getDaughters().size() != 2){
				cout << "ERR: # of higgs daughters = " << mcps[i]->getDaughters().size() << endl;
				break;
			}
			if(abs(mcps[i]->getDaughters()[0]->getPDG()) != 5)break;
			if(abs(mcps[i]->getDaughters()[1]->getPDG()) != 5)break;
			hcount ++;
		}
		if(hh && hcount < 2)continue;

		// semistable B
		cout << "Semistable B: " << endl;
		vector<MCParticle *> blist;
		for(unsigned int i=0;i<mcps.size();i++){
			if(mcps[i]->isSemiStableB()){
				blist.push_back(mcps[i]);
				cout << mcps[i]->getId() << ", " << mcps[i]->getPDG() << ", (" << mcps[i]->getEx() << ", " << mcps[i]->getEy() << ", " << mcps[i]->getEz() << ") " << endl;
			}
		}

		const vector<Track *> &tracks = event->getTracks();
		for(unsigned int i=0;i<tracks.size();i++){
			if(tracks[i]->getMcp() == 0)continue;

			const float *calo = tracks[i]->getCaloEdep();
			double calotot = 0.;
			for(int ncalo=0;ncalo<tpar::caloN;ncalo++)calotot += calo[ncalo];

			MCParticle *ssp = tracks[i]->getMcp()->getSemiStableParent();
			float f[] = {
				tracks[i]->getMcp()->getPDG(), tracks[i]->getMcp()->getParent()->getPDG(), (ssp ? ssp->getPDG() : 0),
				tracks[i]->E(), tracks[i]->Px(), tracks[i]->Py(), tracks[i]->Pz(), tracks[i]->getD0(), tracks[i]->getZ0(),
				sqrt(tracks[i]->getCovMatrix()[tpar::d0d0]), sqrt(tracks[i]->getCovMatrix()[tpar::z0z0]),
				calotot, calo[0],calo[1], calo[2], tracks[i]->getChi2() };

			ntLepton->Fill(f);

			for(unsigned int k=0;k<blist.size();k++){
				if(tracks[i]->getMcp()->isParent(blist[k])){
					MCParticle *parent = blist[k];
					MCParticle *mcptrack = tracks[i]->getMcp();
					TVector3 v = mcptrack->getVertex();
					double gamma = parent->E() / parent->M();
					ntTrack->Fill(parent->getPDG(), v.x(), v.y(), v.z(), gamma, v.Mag()/gamma);
				}
			}
		}
	}
	pf->Write();
	pf->Close();
}

void TrackDist110218(const char *inputfile, const char *outputfile, int nStart, int nEnd, int hh)
{
	LCIOStorer ls(inputfile);
	ls.InitCollections();

	TFile *pf = TFile::Open(outputfile,"recreate");
	TNtuple *ntVTcon = new TNtuple("ntVTcon","ntVertex-Track-connection","vpdg:ve:vmass:tpdg:te:angle:same1:same2");

	int nev = -1;

	while(ls.Next()){
		nev ++;
		if(nev < nStart){continue;}
		if(nev >= nEnd)break;
		cout << "Event #" << nev << endl;

//		ts.GetEntry(nev);
		Event *event = Event::Instance();

		const vector<MCParticle *>& mcps = event->getMCParticles();
		const vector<Track *>& tracks = event->getTracks();

		// check bbbbbb (reject H->WW etc.)
		int hcount = 0;
		for(unsigned int i=0;i<mcps.size();i++){
			if(mcps[i]->getPDG() != 25)continue;

			// higgs
			if(mcps[i]->getDaughters().size() != 2){
				cout << "ERR: # of higgs daughters = " << mcps[i]->getDaughters().size() << endl;
				break;
			}
			if(abs(mcps[i]->getDaughters()[0]->getPDG()) != 5)break;
			if(abs(mcps[i]->getDaughters()[1]->getPDG()) != 5)break;
			hcount ++;
		}
		if(hh && hcount < 2)continue;

		// semistable B
		cout << "Semistable B: " << endl;
		vector<MCParticle *> blist, bParentList, bPairList;
		for(unsigned int i=0;i<mcps.size();i++){
//			cout << mcps[i]->getPDG() << " " << (mcps[i]->getParent() ? mcps[i]->getParent()->getPDG() : 0) << endl;
			if(mcps[i]->isSemiStableB()){
				blist.push_back(mcps[i]);
				MCParticle *par = mcps[i];
				while(par && par->getPDG() != 92)par = par->getParent();
				bParentList.push_back(par);

				cout << mcps[i]->getId() << ", " << mcps[i]->getPDG() << ", (" << mcps[i]->getEx() << ", " << mcps[i]->getEy() << ", " << mcps[i]->getEz() << ") " << endl;
			}
		}
		if(blist.size()%2 == 1){cout << "ERR: odd number of semistable B: skip." << endl; continue;}

		// check b-list pair
		unsigned int k,kk;
		for(k=0;k<blist.size();k++){
			for(kk=0;kk<blist.size();kk++){
				if(k != kk && bParentList[k] == bParentList[kk]){
					bPairList.push_back(blist[kk]);
					break;
				}
			}
			if(kk == blist.size()){cout << "ERR: no pair found." << endl; break;}
		}
		if(k != blist.size()){continue;}

		for(unsigned int i=0;i<tracks.size();i++){
			MCParticle *mcptrack = tracks[i]->getMcp();
			if(!mcptrack)continue;

			for(k=0;k<blist.size();k++){
				double angle = blist[k]->Vect().Angle(mcptrack->Vect());
				bool same1 = mcptrack->isParent(bParentList[k]);
				double angle2 = bPairList[k]->Vect().Angle(mcptrack->Vect());

				bool same2 = same1 && (angle < angle2);
				ntVTcon->Fill(blist[k]->getPDG(), blist[k]->E(), blist[k]->getParent()->M(), mcptrack->getPDG(), mcptrack->E(), angle, same1, same2);
			}
		}
	}
	pf->Write();
	pf->Close();
}

vector<Jet *> SueharaJetClustering(Event &evt, int njets)
{
	TrackSelectorConfig secVtxCfg;
	// original setup 110214
	secVtxCfg.maxD0 = 10;
	secVtxCfg.maxZ0 = 20;
	secVtxCfg.minPt = 0.1;
	secVtxCfg.maxInnermostHitRadius = 1e10;

	// 110214
	secVtxCfg.maxD0Err = 0.1;
	secVtxCfg.maxZ0Err = 0.1;

	// OR
	secVtxCfg.minTpcHits = 20;
	secVtxCfg.minFtdHits = 3;
	secVtxCfg.minVtxHitsWithoutTpcFtd = 3;
	secVtxCfg.minVtxPlusFtdHits = 0;

	JetConfig jetCfg;
	jetCfg.nJet = njets;
	JetFinder* jetFinder = new JetFinder(jetCfg);
	double ycut = 0;

	Event *event = &evt;
	// cut bad tracks
	const vector<Track *> &tracks = event->getTracks();
	vector<Track *> passedTracks = TrackSelector() (tracks, secVtxCfg);
	vector<Vertex *> vertices;

// 	for(unsigned int i=0;i<tracks.size();i++){
// 		if(interface.passesCut(tracks[i], secVtxCfg))
// 			passedTracks.push_back(tracks[i]);
// 	}

	// build up vertexing
	VertexFinderSuehara::buildUp(passedTracks, vertices, 25,9, 10);
	VertexFinderSuehara::associateIPTracks(vertices,1.);

	// select vertices
	vector<Vertex *> selectedVertices;
	vector<Track *> residualTracks = tracks;
	for(vector<Vertex *>::iterator it = vertices.begin(); it != vertices.end();it++){
		Vertex *v = *it;
		double mass = 0.;
		if(v->getTracks().size() == 2)
			mass = (*(TLorentzVector *)(v->getTracks()[0]) + *(TLorentzVector *)(v->getTracks()[1])).M();
		if((mass < .488 || mass > .508) && v->getPos().Mag() < 30 && v->getPos().Mag() > 0.3){
			selectedVertices.push_back(v);
			for(vector<Track *>::const_iterator itt = v->getTracks().begin(); itt != v->getTracks().end(); itt++){
				vector<Track *>::iterator itt2 = remove_if(residualTracks.begin(), residualTracks.end(), bind2nd(equal_to<Track *>(), *itt));
				residualTracks.erase(itt2, residualTracks.end());
			}
		}
	}
	return jetFinder->run(residualTracks, event->getNeutrals(), selectedVertices, &ycut);
}
