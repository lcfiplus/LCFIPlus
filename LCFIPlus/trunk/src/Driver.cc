#include <assert.h>
#include "EventStore.h"
#include "LcfiInterface.h"
#include "JetFinder.h"
#include "TreeStorer.h"
#include "LCIOStorer.h"
#include "VertexFitterLCFI.h"
#include "VertexFinderTearDown.h"
#include "algoSigProb.h"
#include "algoEtc.h"
#include "Driver.h"

#include "TROOT.h"
#include "TInterpreter.h"
#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TCut.h"
#include <string>
#include "TRandom3.h"

// for event display
#include "TRootBrowser.h"
#include "TRint.h"
#include "TSystem.h"
#ifdef BUILD_EVE
#include "TEveManager.h"
#include "TEveBrowser.h"
#endif
#include "TGClient.h"
#include "TGFrame.h"
#include "TGButton.h"
#include "EventNavigator.h"

#ifdef BUILD_EVE
extern TEveManager* gEve;
extern TSystem* gSystem;
#endif

using namespace lcfiplus;
using namespace lcfiplus::algoSigProb;
using namespace lcfiplus::algoEtc;

TRandom3 _rand;

const double PION_MASS = 0.13957018;
const double PION_MASS2 = PION_MASS*PION_MASS;

#if 0
void ioTest(int argc, char* argv[]) {
	/*
		 EventCollection* store = new EventCollection();
		 store->openWrite("data.gz");
		 for (int i=0; i<10; ++i) {
		 Event* event = new Event();

		 for (int j=0; j<8; ++j) {
		 Track* track = new Track();
		 track->id = j;
		 track->par[Track::d0] = 0.11*j;
		 event->add(track);
		 }

		 for (int j=0; j<3; ++j) {
		 MCParticle* mcp = new MCParticle();
		 mcp->id = j;
		 mcp->px = 0.35*j;
		 event->add(mcp);
		 }

		 store->write(event);
		 }
		 store->closeWrite();
	 */

	char* input = "input.piu";
	char* output = "output.piu";

	if (argc>2) {
		input = argv[1];
		output = argv[2];
	}

	int nev(0);
	EventColRaw* store2 = new EventColRaw();
	store2->open(input);
	//store2->open("data.gz");
	store2->openWrite(output);
	EventData* eventData;
	while ( (eventData = store2->nextData()) ) {
		if ( (++nev % 100) == 0 ) {
			printf("nev = %d\n",nev);
		}
		store2->write(eventData);
		delete eventData;
	}
	store2->close();
	store2->closeWrite();
}
#endif


/* Associate MC particle with jets by closest angle.
	 If the jet has >=1 MC particle with heavy flavor,
	 that jet will be designated as a heavy jet.
 */
vector<int> findMcJetFlavor(vector<const Jet*> jets, vector<const MCParticle*> mcps) {
	vector<int> mcJetFlavor;
	vector<TVector3> jetDirs;

	for (JetVecIte iter = jets.begin(); iter != jets.end(); ++iter) {
		const Jet* jet = *iter;
		TVector3 v = jet->Vect();
		jetDirs.push_back(v.Unit());
		mcJetFlavor.push_back(1);
	}

	for (MCParticleVecIte iter = mcps.begin(); iter != mcps.end(); ++iter) {
		const MCParticle* mcp = *iter;
		int mcFlavor = mcp->getFlavor();
		if (mcFlavor < 4) continue;

		//printf("PDG code: %d (flavor=%d)\n", mcp->pdg, mcFlavor);
		TVector3 v = mcp->getVertex();
		TVector3 dir = v.Unit();

		double minAngle(1e10);
		unsigned int minIndex(0);

		for (unsigned int j=0; j<jetDirs.size(); ++j) {
			double angle = dir.Angle(jetDirs[j]);
			if (angle < minAngle) {
				minAngle = angle;
				minIndex = j;
			}
		}

		/*
			 if (minAngle > 0.5) {
			 printf("While matching MC particle to jets, a large angle deviation detected: %f\n", minAngle);
			 }
		 */

		if (mcFlavor > mcJetFlavor[minIndex])
			mcJetFlavor[minIndex] = mcFlavor;
	}

	return mcJetFlavor;
}

vector<MCVertex*> findMcVertex(vector<const MCParticle*> mcps) {
	bool debug(false);
	vector<MCVertex*> ret;

	// map to store daughter -> semi-stable parent association
	map<const MCParticle*,vector<const MCParticle*> > sspDau;

	for (MCParticleVecIte it = mcps.begin(); it != mcps.end(); ++it) {
		const MCParticle* mcp = *it;
		if (mcp->isStable() == false && mcp->isSemiStable() == false)
			continue;

		const MCParticle* ssp = mcp->getSemiStableParent();
		if (ssp == 0) continue;
		//printf("%d is stable with ssp %d\n",mcp->getPDG(),ssp->getPDG());

		// is the ssp very short-lived?
		float dif = ( ssp->getVertex()-mcp->getVertex() ).Mag();
		if (dif < 5e-3) {
			// if less than 5 microns, consider the vertex unresolvable
			//printf("less than 5e-3\n");
			const MCParticle* sspp = ssp->getSemiStableParent();
			if (sspp) ssp = sspp;
		}
		/*
			printf("sspp=%.3e, ssp=%.3e, dif=%.3e\n",
					sspp->getVertex().Mag(),
					ssp->getVertex().Mag(),
					dif);
		 */

		map<const MCParticle*,vector<const MCParticle*> >::iterator sspIter = sspDau.find(ssp);
		if ( sspIter == sspDau.end() ) {
			vector<const MCParticle*> a;
			a.push_back(mcp);
			sspDau.insert( make_pair(ssp, a) );
		} else {
			// find vector and push mcp
			sspDau[ssp].push_back(mcp);
		}
	}

	if (debug) {
		for (map<const MCParticle*,vector<const MCParticle*> >::iterator it = sspDau.begin(); it != sspDau.end(); ++it) {
			const MCParticle* ssp = it->first;
			vector<const MCParticle*>& dauVec = it->second;
			printf(" ssp=%d ndau=%d (start=%.3e, end=%.3e) \n",
					(int)ssp->getPDG(), (int)ssp->getDaughters().size(),
					ssp->getVertex().Mag(), ssp->getEndVertex().Mag());
			for (MCParticleVecIte dauIter = dauVec.begin(); dauIter != dauVec.end(); ++dauIter) {
				const MCParticle* dau = *dauIter;
				printf("   %d [p=%d] (start=%.3e, end=%.3e)\n", dau->getPDG(),
						dau->getParent()->getPDG(), dau->getVertex().Mag(),
						dau->getEndVertex().Mag() );
			}
		}
	}

	for (map<const MCParticle*,vector<const MCParticle*> >::iterator it = sspDau.begin(); it != sspDau.end(); ++it) {
		const MCParticle* ssp = it->first;
		if (ssp->getEndVertex().Mag() > 1000) continue;

		int ntrk(0);
		vector<const MCParticle*>& dauVec = it->second;
		for (MCParticleVecIte dauIter = dauVec.begin(); dauIter != dauVec.end(); ++dauIter) {
			const MCParticle* dau = *dauIter;
			if (dau->isStableTrack())
				++ntrk;
		}
		
		if (ntrk == 0) continue;

		MCVertex* v = new MCVertex;
		v->setParent(ssp);

		for (MCParticleVecIte dauIter = dauVec.begin(); dauIter != dauVec.end(); ++dauIter) {
			const MCParticle* dau = *dauIter;
			if (dau->isStableTrack()) {
				v->add(dau);
			}
		}
		ret.push_back(v);
	}

	if (debug) {
		for (vector<MCVertex*>::iterator it = ret.begin(); it != ret.end(); ++it) {
			MCVertex* v = *it;
			const MCParticle* ssp = v->getParent();
			vector<const MCParticle*> dauVec = v->getDaughters();

			printf(" ssp=%d ndau=%d (start=%.3e, end=%.3e) \n",
					(int)ssp->getPDG(), (int)ssp->getDaughters().size(),
					ssp->getVertex().Mag(), ssp->getEndVertex().Mag());
			for (vector<const MCParticle*>::iterator dauIter = dauVec.begin(); dauIter != dauVec.end(); ++dauIter) {
				const MCParticle* dau = *dauIter;
				if (dau->isStableTrack()) {
					printf("   %d [p=%d] (start=%.3e, end=%.3e)\n", (int)dau->getPDG(),
							(int)dau->getParent()->getPDG(), dau->getVertex().Mag(),
							dau->getEndVertex().Mag() );
				}
			}
		}
	}

	return ret;
}

const int SINGLE_MC_TRACK					= 1 << 0;
const int NOT_ENOUGH_RECOTRK			= 1 << 1;
const int TEARDOWN_FAILED					= 1 << 2;
const int TEARDOWN_LARGE_CHISQ		= 1 << 3;
const int ZVTOP_FAILED						= 1 << 4;
const int ZVTOP_LARGE_CHISQ				= 1 << 5;
const int TRACKS_NOT_IN_JETS			= 1 << 6;
const int TRACKS_IN_MULTIPLE_JETS	= 1 << 7;
const int RECO_MATCH							= 1 << 8;
const int V0_MATCH								= 1 << 9;

void matchMcVertex( const Event& evt, vector<MCVertex*>& vtxList, map<MCVertex*,int>& table, bool vertexing )
{
	bool debug(false);
	vector<const Track*> tracks = evt.getTracks();
	LcfiInterface lcfi;

	for (vector<MCVertex*>::const_iterator it = vtxList.begin(); it != vtxList.end(); ++it) {
		MCVertex* v = *it;
		table[v] = 0; // initialize flag

		const MCParticle* ssp = v->getParent();
		vector<const MCParticle*> dauVec = v->getDaughters();

		if (debug) {
			printf("================================================\n");
			printf("MC Vertex: pdg=%d ndau=%d decay=(%.3e,%.3e,%.3e) \n",
					(int)ssp->getPDG(), (int)dauVec.size(),
					ssp->getEndVertex().x(),
					ssp->getEndVertex().y(),
					ssp->getEndVertex().z()
					);
		}

		if (ssp->getDaughters().size() == 1) { table[v] |= SINGLE_MC_TRACK; }

		int nRecoTrk(0);
		vector<const Track*> trkList;

		for (MCParticleVecIte dauIter = dauVec.begin(); dauIter != dauVec.end(); ++dauIter) {
			const MCParticle* dau = *dauIter;
			const Track* trk(0);
			for (TrackVecIte trkIter = tracks.begin(); trkIter != tracks.end(); ++trkIter) {
//				if ( dau == evt.getMCParticle(*trkIter) ) {
				if ( dau == (*trkIter)->getMcp() ) {
					//printf(" match\n");
					if (trk != 0) {
						printf("double-match found: track1, track2 -> mcparticle\n");
					} else {
						trk = *trkIter;
						++nRecoTrk;
					}
				}
			}

			if (trk) {
				v->add(trk); // add recotrk to mcvertex
				trkList.push_back(trk);
			}

			if (debug) {
				printf("  %s dau=%4d", trk ? "      " : "*LOST*", dau->getPDG() );

				if (!trk) {
					float cosTheta = fabs( sin( atan( dau->getTanLambda() ) ) );
					printf("   **** track lost due to ");
					if (dau->Pt() < 0.2) {
						printf("low pt (%.4f)",dau->Pt());
					} else if (cosTheta > 0.95) {
						printf("being outside acceptance (cosTheta=%.4f)",cosTheta);
					} else {
						printf("unknown reasons");
					}
				}
				printf("\n");
			}
		}

		if (nRecoTrk <= 1) { table[v] |= NOT_ENOUGH_RECOTRK; }//printf("setting not enough reco: %d\n",table[v]);}

		if (vertexing && nRecoTrk > 1) {
			//printf("--------------- Running teardown vertexing...\n");
			Vertex* vtx = VertexFinderTearDown<vector>()(trkList, 0, 1e10, 0);

			if (!vtx) {
				//printf("*************** NOT FOUND\n");
				table[v] |= TEARDOWN_FAILED;
			}
			
			if (vtx) {
				if (vtx->getChi2() > 10) { table[v] |= TEARDOWN_LARGE_CHISQ; }

				if (debug) {
					printf("Teardown Vertex found:  (%.3e,%.3e,%.3e) chi2=%.3e ntrk=%d\n",
							vtx->getX(),vtx->getY(),vtx->getZ(),
							vtx->getChi2(),
							(int)vtx->getTracks().size()
							);
					for(unsigned int i=0;i<vtx->getTracks().size();i++){
						const Track * tr = vtx->getTracks()[i];
						printf("        Track #%d  pos=(%.3e,%.3e,%.3e), chi2=%e\n",
								i, tr->getX(), tr->getY(), tr->getZ(), vtx->getChi2Track(tr));
					}
				}
				delete vtx;
			}
			//printf("--------------- Running ZVTOP...\n");
			Jet jet;
			for (unsigned int i=0; i<trkList.size(); ++i) {
				jet.add( Jet(trkList[i]) );
			}
			//printf("jet e=%.3e\n",jet.E());
			vector<Vertex*> zvtopList = lcfi.forceZvtop(jet);
			if (zvtopList.size() == 0) {
				//printf("*************** NOT FOUND\n");
				table[v] |= ZVTOP_FAILED;
			}
			if (zvtopList.size() > 0) {
				Vertex* vtx = zvtopList[0];
				if (vtx->getChi2() > 10) { table[v] |= ZVTOP_LARGE_CHISQ; }

				if (debug) {
					printf("ZVTOP Vertex Found (%d): (%.3e,%.3e,%.3e) chi2=%.3e ntrk=%d\n",
							(int)zvtopList.size(),
							vtx->getX(),vtx->getY(),vtx->getZ(),
							vtx->getChi2(),
							(int)vtx->getTracks().size()
							);
					for(unsigned int i=0;i<vtx->getTracks().size();i++){
						const Track * tr = vtx->getTracks()[i];
						printf("        Track #%d  pos=(%.3e,%.3e,%.3e), chi2=%e\n",
								i, tr->getX(), tr->getY(), tr->getZ(), vtx->getChi2Track(tr) );
					}
				}
			}
		}
	}

	/*
	for (vector<Vertex*>::iterator rvtxIter = recovtx.begin(); rvtxIter != recovtx.end(); ++rvtxIter) {
		Vertex* recov = *rvtxIter;
		vector<Track*> dauList = recov->getTracks();
		Vertex * vtx = lcfiplus::VertexFinderTearDown<vector>()(tracksForPrimary, &beamTracks, 9.0, 0);
		printf("ndau = %d\n", dauList.size() );
	}
	*/
}

void matchMcVertexJet(
		const Event& evt,
		const vector<MCVertex*>& vtxList, map<MCVertex*,int>& table,
		const vector<const Jet*>& jets)
{
	for (vector<MCVertex*>::const_iterator it = vtxList.begin(); it != vtxList.end(); ++it) {
		MCVertex* v = *it;
		vector<const MCParticle*> dauVec = v->getDaughters();

		vector<const Track*> trkList;
		int nRecoTrk(0);
		TrackVec & tracks = evt.getTracks();

		for (MCParticleVecIte dauIter = dauVec.begin();
				dauIter != dauVec.end(); ++dauIter) {

			const MCParticle* dau = *dauIter;
			const Track* trk(0);
			for (TrackVecIte trkIter = tracks.begin(); trkIter != tracks.end(); ++trkIter) {
				if ( dau == evt.getMCParticle(*trkIter) ) {
					trk = *trkIter;
					break;
				}
			}
			if (trk) {
				trkList.push_back(trk);
				++nRecoTrk;
			}
		}

		if (nRecoTrk==0) continue;

		// check to see if the daughter tracks are split among jets
		int nJetMatch(0); // should be just one if there are no split
		for (JetVecIte jetIter = jets.begin(); jetIter != jets.end(); ++jetIter) {
			const Jet* j = *jetIter;
			vector<const Track*> jtrks = j->getTracks();
			for (TrackVecIte trkIter = trkList.begin(); trkIter != trkList.end(); ++trkIter) {
				if ( find( jtrks.begin(), jtrks.end(), *trkIter ) != jtrks.end() ) {
					++nJetMatch;
					break;
				}
			}
		}

		if (nJetMatch == 0) {
			table[v] |= TRACKS_NOT_IN_JETS;
		}
		if (nJetMatch > 1) {
			table[v] |= TRACKS_IN_MULTIPLE_JETS;
		}

		//printf("nJetMatch=%d\n",nJetMatch);
	}
}

void matchMcVertexReco(
		const Event& evt,
		const vector<MCVertex*>& vtxList, map<MCVertex*,int>& table,
		Vertex* vertex )
{
	const vector<const Track*>& tracks = vertex->getTracks();

	int nTrkMatchMAX(0);
	MCVertex* match(0);

	for (vector<MCVertex*>::const_iterator it = vtxList.begin(); it != vtxList.end(); ++it) {
		MCVertex* v = *it;
		vector<const Track*> recoTrksMatch = v->getRecoTracks();

		int nTrkMatch(0);

		for (TrackVecIte trkIter = recoTrksMatch.begin();
				trkIter != recoTrksMatch.end(); ++trkIter) {
			const Track* trk = *trkIter;
			if ( find( tracks.begin(), tracks.end(), trk ) != tracks.end() ) {
				++nTrkMatch;
			}
		}

		//if (nTrkMatch > 0) {
		if (nTrkMatch >= 2) {
			if (nTrkMatch > nTrkMatchMAX) {
				match = v;
			}
		}
	}

	if (match) {
		if (match->getRecoVertex() != 0) {
			printf(" vertex already matched!!\n");
		} else {
			match->setRecoVertex(vertex);
			table[match] |= RECO_MATCH;
		}
	}
}

void matchMcVertexRecoV0(
		const Event& evt,
		const vector<MCVertex*>& vtxList, map<MCVertex*,int>& table )
{
	const vector<const Neutral*>& neutrals = evt.getNeutrals();

	for (NeutralVecIte neutIter = neutrals.begin();
			neutIter != neutrals.end(); ++neutIter) {
		const Neutral* neut = *neutIter;
		if (neut->isV0()) {
			for (vector<MCVertex*>::const_iterator it = vtxList.begin(); it != vtxList.end(); ++it) {
				MCVertex* v = *it;
				const MCParticle* mcv0 = v->getParent();
				int vpdg = mcv0 ? mcv0->getPDG() : 0;
				if (mcv0 && vpdg == neut->getPDG()) {
					printf(" pfo %d, parent %d\n",neut->getPDG(), vpdg );
					printf(" reco (%.2e,%.2e,%.2e) mc (%.2e,%.2e,%.2e)\n",
							neut->Px(), neut->Py(), neut->Pz(),
							mcv0->Px(), mcv0->Py(), mcv0->Pz()
							);
					table[v] |= V0_MATCH;
				}
			}
		}
	}
}

vector<Vertex*> findAdditionalVertices(
		const Event& evt, const Jet* jet,
		const vector<Vertex*>& seedList,
		const Vertex* primaryVertex )
{
	vector<Vertex*> ret;

	vector<const Track*> useTracks;
	vector<const Track*> tracksToReject; // not require single

	TrackVec & ptrks = primaryVertex->getTracks();
	for (TrackVecIte trkIter = ptrks.begin(); trkIter != ptrks.end(); ++trkIter) {
		tracksToReject.push_back(*trkIter);
	}

	for (vector<Vertex *>::const_iterator it = seedList.begin(); it != seedList.end(); ++it) {
		Vertex* v = *it;
		ret.push_back(v);
		TrackVec & vtrks = v->getTracks();
		for (TrackVecIte vtrkIter = vtrks.begin(); vtrkIter != vtrks.end(); ++vtrkIter) {
			tracksToReject.push_back(*vtrkIter);
		}
	}

	TrackVec & jtrks = jet->getTracks();
	for (TrackVecIte jtrkIter = jtrks.begin(); jtrkIter != jtrks.end(); ++jtrkIter) {
		if ( find( tracksToReject.begin(), tracksToReject.end(), *jtrkIter ) == tracksToReject.end() ) {
			useTracks.push_back(*jtrkIter);
		}
	}

	Jet jet2;
	for (TrackVecIte utrkIter = useTracks.begin(); utrkIter != useTracks.end(); ++utrkIter) {
		jet2.add(*utrkIter);
	}

	vector<Vertex*>* extraVtxList = findTearDownVertices( evt, jet2 );

	for (vector<Vertex *>::const_iterator it = extraVtxList->begin(); it != extraVtxList->end(); ++it) {
		ret.push_back(*it);
	}

	delete extraVtxList;

	return ret;
}

vector<const Track*> findSingleTracks(const Event& evt, const Jet& jet, const vector<lcfiplus::Vertex*>& vtxList) {
	vector<const Track*> ret;
	vector<const Track*> cand;

	TrackVec &tracks = jet.getTracks();
	for (unsigned int i=0; i<tracks.size(); ++i) {
		const Track* trk = tracks[i];
		for (unsigned int j=0; j<vtxList.size(); ++j) {
			TrackVec & vtxTracks = vtxList[j]->getTracks();
			if (find(vtxTracks.begin(),vtxTracks.end(),trk) != vtxTracks.end()) {
				cand.push_back(trk);
			}
		}
	}

	//printf("%d\n",cand.size());

	return cand;
	return ret;
}

		void eventDisplay(const char* input, int start) {
#ifdef BUILD_EVE
			TRint* theApp = 0;
			theApp = new TRint("ROOT",0,0);

			TEveManager::Create();

			TEveBrowser* browser = gEve->GetBrowser();
			browser->StartEmbedding(TRootBrowser::kLeft);

			TGMainFrame* frmMain = new TGMainFrame(gClient->GetRoot(),600,400);
			frmMain->SetWindowName("XX GUI");
			frmMain->SetCleanup(kDeepCleanup);

			TGGroupFrame* frmEvent = new TGGroupFrame(frmMain, "Event Navigation", kHorizontalFrame);
			TGHorizontalFrame* hf = new TGHorizontalFrame(frmMain);

			EventNavigator* fh = new EventNavigator(input,start);
			//void* fh(0);

			TString icondir( Form("%s/icons/", gSystem->Getenv("ROOTSYS")) );
			TGPictureButton* b = 0;

			b = new TGPictureButton(hf, gClient->GetPicture(icondir + "GoBack.gif"));
			hf->AddFrame(b);
			b->Connect("Clicked()", "lcfiplus::EventNavigator", fh, "Bck()");

			b = new TGPictureButton(hf, gClient->GetPicture(icondir + "GoForward.gif"));
			hf->AddFrame(b);
			b->Connect("Clicked()", "lcfiplus::EventNavigator", fh, "Fwd()");

			frmEvent->AddFrame(hf);
			frmMain->AddFrame(frmEvent);

			frmMain->MapSubwindows();
			frmMain->Resize();
			frmMain->MapWindow();

			browser->StopEmbedding();
			browser->SetTabTitle("Navigation",0);


			gEve->FullRedraw3D(kTRUE);
			theApp->Run();
#else
			cerr << "BUILD_EVE not defined. event display is inactivated." << endl;
#endif
		}

int main(int argc, char* argv[]) {

	// if argc==1 abort
	// if argc==2 argv[1] macro is processed
	// if argc==3 argv[2] function in argv[1] file is processed
	// if argc==4 argv[2] function in argv[1] file is processed with param argv[3]
	// if argc>=5 abort 

	if(argc==1 || argc>=5){cout << "Usage: lcfiplus filename [ funcname [ params ] ]" << endl; return 1;}
	const char *macroname = argv[1];
	const char *funcname = (argc>2 ? argv[2] : "");
	const char *params = (argc>3 ? argv[3] : "");

	if(argc==2){
		cout << "Executing " << macroname << " ..." << endl;
		gROOT->Macro(macroname);
		//			gROOT->ProcessLine(TString::Format(".x %s",macroname));
	}
	else{
		if(strcmp(macroname, "0")){
			cout << "Loading " << macroname << " ..." << endl;
			gROOT->LoadMacro(macroname);
		}
		cout << "Calling " << funcname;
		if(argc==4)	cout << " with param " << params;
		cout << " ..." << endl;
		gROOT->GetInterpreter()->Execute(funcname, params);
	}

	return 0;
}
