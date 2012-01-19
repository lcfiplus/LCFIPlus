#include <assert.h>
#include "EventStore.h"
#include "LcfiInterface.h"
#include "JetFinder.h"
#include "TreeStorer.h"
#include "LCIOStorer.h"
#include "VertexFitterLCFI.h"
#include "VertexFinderTearDown.h"
#include "VertexFinderPerfect.h"
#include "algoSigProb.h"
#include "algoEtc.h"
#include "Driver.h"
#include "geometry.h"

#include "TROOT.h"
#include "TInterpreter.h"
#include "TFile.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TLorentzVector.h"
#include "TCut.h"
#include <string>
#include "TRandom3.h"

// for event display
#include "TRootBrowser.h"
#include "TRint.h"
#include "TSystem.h"
#include "TEveManager.h"
#include "TEveBrowser.h"
#include "TGClient.h"
#include "TGFrame.h"
#include "TGButton.h"
#include "EventNavigator.h"

//extern TEveManager* gEve;
//extern TSystem* gSystem;

using namespace lcfiplus;
using namespace lcfiplus::algoSigProb;
using namespace lcfiplus::algoEtc;

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


/////////////////////////////////////////////////////////////////////////

void testSuehara(const char *inputlist, const char *output)
{

	LCIOStorer ls(inputlist);
	ls.InitCollections("PandoraPFOs","MCParticlesSkimmed","RecoMCTruthLink","Tracks","Neutrals","MCParticles");
	ls.InitVertexCollection("BuildUpVertex", "BuildUpVertex");

	cout << "LCIO initialization successful." << endl;

	Event *event = Event::Instance();
	event->Print();

	TFile f(output,"RECREATE");
//	TNtuple *nt = new TNtuple("nt","nt","trackcategory:invertex:allsectracks:categorysectracks:d0:d0err:z0:z0err:e:pdg:mcvpos");
	TNtuple *nt = new TNtuple("nt","nt","trackcategory:categorysectracks:d0:d0err:z0:z0err:e:pdg:mcvpos:truecombi:linedistance:angle:trackpos:vtxpos:vtxtracks:vtxbtracks:vtxctracks");
//					nt->Fill(category, categorysectracks, d0,d0e,z0,z0e,e, mcp->getPDG(),mcvpos, truecombi, distance, angle, pos.Mag(), vtcs[nvtcs]->getPos());

	struct {
		float trackcategory;
		float categorysectracks;
		float d0;
		float d0err;
		float z0;
		float z0err;
		float e;
		float pdg;
		float mcvpos;
		float truecombi;
		float linedistance;
		float angle;
		float trackpos;
		float vtxpos;
		float vtxtracks;
		float vtxbtracks;
		float vtxctracks;
	}data;

	while(ls.Next()){
		MCParticleVec& mcps = event->getMCParticles();
/*
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
*/
		TrackVec &tracks = event->getTracks();
		VertexVec &vtcs = event->getSecondaryVertices("BuildUpVertex");
		MCParticleVec bs = event->mcGetSemiStableBs();
		MCParticleVec cs = event->mcGetSemiStableCs();

		cout << "We have " << bs.size() << " bs & " << cs.size() << " cs." << endl;
		cout << "# vertices = " << vtcs.size() << endl;
		for(unsigned int n=0;n<vtcs.size();n++){
			cout << "RecoVertices at (" << vtcs[n]->getX() << ", " << vtcs[n]->getY() << ", " << vtcs[n]->getZ() << "), # tracks = " << vtcs[n]->getTracks().size() << endl;
		}

		// perfect vertices
		vector<MCVertex *> mcvtcs;
		VertexFinderPerfect::findPerfectVertices(tracks, mcps, mcvtcs, 2);
		cout << "MCVertices ( >= 2 tracks) : " << mcvtcs.size() << endl;
		mcvtcs.clear();
		VertexFinderPerfect::findPerfectVertices(tracks, mcps, mcvtcs, 1);
		cout << "MCVertices ( >= 1 tracks) : " << mcvtcs.size() << endl;

		vector<int> btracks, bctracks;
		btracks.resize(bs.size());
		bctracks.resize(bs.size());

		// counting b/c tracks
		for(unsigned int n=0;n<tracks.size();n++){
			const MCParticle *mcp = tracks[n]->getMcp();

			int npar = event->mcFindParent(bs, mcp);
			if(npar < 0)continue; // non-bc track

			btracks[npar] ++;
			if(mcp->getSemiStableCParent())
				bctracks[npar] ++;
		}

		// fill ntuple
		for(unsigned int n=0;n<tracks.size();n++){
			const MCParticle *mcp = tracks[n]->getMcp();

			data.trackcategory = 0;
			if(mcp->getSemiStableCParent()) data.trackcategory = 2;
			else if(mcp->getSemiStableBParent()) data.trackcategory = 1;
			else if(mcp->getSemiStableParent()) data.trackcategory = 3;

			bool invertex = false;
			for(unsigned int nvtcs = 1; nvtcs < vtcs.size(); nvtcs ++){
				const Vertex *vtx = vtcs[nvtcs];

				if(find(vtx->getTracks().begin(), vtx->getTracks().end(), tracks[n]) != vtx->getTracks().end())
					invertex = true;
			}

			int nb = event->mcFindParent(bs, mcp);

/*
			if(!invertex && nb >= 0){
				TVector3 ip(0,0,0);
				VertexLine line(ip, bs[nb]->getEndVertex());
				Helix hel(tracks[n]);

				hel.ClosePoint(line);
				cout << "MCVertex: " << mcp->getVertex().X() << ", " << mcp->getVertex().Y() << ", " << mcp->getVertex().Z() << endl;
			}
*/
			int allsectracks = 0;
			data.categorysectracks = 0;

			if(nb >= 0){
				allsectracks = btracks[nb];
				if(data.trackcategory == 2)data.categorysectracks = bctracks[nb];
				else data.categorysectracks = btracks[nb] - bctracks[nb];
			}
			data.d0 = tracks[n]->getD0();
			data.d0err = sqrt(tracks[n]->getCovMatrix()[tpar::d0d0]);
			data.z0 = tracks[n]->getZ0();
			data.z0err = sqrt(tracks[n]->getCovMatrix()[tpar::z0z0]);
			data.e = tracks[n]->E();
			data.mcvpos = mcp->getVertex().Mag();
			data.pdg = mcp->getPDG();

//			nt->Fill(category, invertex, allsectracks, categorysectracks, d0,d0e,z0,z0e,e, mcp->getPDG(),mcvpos);
			if(!invertex){
				for(unsigned int nvtcs = 1; nvtcs < vtcs.size(); nvtcs ++){

					// obtain vertex tracks
					data.vtxtracks = vtcs[nvtcs]->getTracks().size();
					data.vtxbtracks = 0;
					data.vtxctracks = 0;
					data.truecombi = 0;
					for(unsigned int ntr = 0; ntr < data.vtxtracks; ntr ++)
					{
						const Track *tr = vtcs[nvtcs]->getTracks()[ntr];

						if(tr->getMcp()->getSemiStableCParent())data.vtxctracks ++;
						else if(tr->getMcp()->getSemiStableBParent())data.vtxbtracks ++;

						if(tr->getMcp()->isParent(bs[nb]))data.truecombi++;
					}

					VertexLine line(vtcs[0]->getPos(), vtcs[nvtcs]->getPos());
					Helix hel(tracks[n]);

					double linedist = 0;
					TVector3 pos = hel.ClosePoint(line, &linedist);
					data.linedistance = linedist;

					data.angle = vtcs[nvtcs]->getPos().Angle(tracks[n]->Vect());

					data.trackpos = pos.Mag();
					data.vtxpos = vtcs[nvtcs]->getPos().Mag();

					nt->Fill((float *)&data);
				}
			}
		}
	}
	f.Write();
}

#include "TrackSelector.h"

typedef ROOT::Math::SVector<double, 2> SVector2;
typedef ROOT::Math::SVector<double, 3> SVector3;
typedef ROOT::Math::SMatrix<double, 2,2,ROOT::Math::MatRepStd<double,2,2> > SMatrix2;
typedef ROOT::Math::SMatrix<double, 3,3,ROOT::Math::MatRepStd<double,3,3> > SMatrix3;
typedef ROOT::Math::SMatrix<double, 2,2,ROOT::Math::MatRepSym<double,2> > SMatrixSym2;
typedef ROOT::Math::SMatrix<double, 3,3,ROOT::Math::MatRepSym<double,3> > SMatrixSym3;
typedef ROOT::Math::SMatrix<double, 5,5,ROOT::Math::MatRepSym<double,5> > SMatrixSym5;
typedef ROOT::Math::SMatrix<double, 2,3,ROOT::Math::MatRepStd<double,2,3> > SMatrix23;
typedef ROOT::Math::SMatrix<double, 3,2,ROOT::Math::MatRepStd<double,3,2> > SMatrix32;

struct KalVtx {
	SVector3 pos;
	SMatrixSym3 cov;
};

struct KalTrk {
	SVector2 _dz;
	SVector3 _mom;
	SMatrixSym5 _cov;
	SMatrixSym2 _cov11;
	SMatrixSym3 _cov22;
	SMatrix23 _cov12;

	SMatrix32 _matGain;
	SVector2 _residual;

	KalTrk(const Track& trk) {
		_dz(0) = trk.getD0();
		_dz(1) = trk.getZ0();
		_mom(0) = trk.getPhi();
		_mom(1) = trk.getOmega();
		_mom(2) = trk.getTanLambda();
		_cov(0,0) = trk.getCovMatrix()[tpar::d0d0];
		_cov(0,1) = trk.getCovMatrix()[tpar::d0z0];
		_cov(0,2) = trk.getCovMatrix()[tpar::d0ph];
		_cov(0,3) = trk.getCovMatrix()[tpar::d0om];
		_cov(0,4) = trk.getCovMatrix()[tpar::d0td];
		_cov(1,1) = trk.getCovMatrix()[tpar::z0z0];
		_cov(1,2) = trk.getCovMatrix()[tpar::z0ph];
		_cov(1,3) = trk.getCovMatrix()[tpar::z0om];
		_cov(1,4) = trk.getCovMatrix()[tpar::z0td];
		_cov(2,2) = trk.getCovMatrix()[tpar::phph];
		_cov(2,3) = trk.getCovMatrix()[tpar::phom];
		_cov(2,4) = trk.getCovMatrix()[tpar::phtd];
		_cov(3,3) = trk.getCovMatrix()[tpar::omom];
		_cov(3,4) = trk.getCovMatrix()[tpar::omtd];
		_cov(4,4) = trk.getCovMatrix()[tpar::tdtd];

		for (int i=0; i<2; ++i) for (int j=i; j<2; ++j) _cov11(i,j) = _cov(i,j);
		for (int i=0; i<3; ++i) for (int j=i; j<3; ++j) _cov22(i,j) = _cov(i+2,j+2);
		for (int i=0; i<2; ++i) for (int j=i; j<3; ++j) _cov12(i,j) = _cov(i,j+2);
	}

	// computes & returns chi2 contribution
	double compute(const KalVtx& vtx) {
		SMatrix2 covRes;
		double chi2(0.);

		SVector2 dz2 = getClosestXY(vtx.pos);
		_residual = dz2-_dz; // residual

		/*
		printf(" trk dz: (%.4e,%.4e)\n", _dz(0),_dz(1));
		printf(" vtx dz: (%.4e,%.4e)\n", dz2(0),dz2(1));
		printf("----------\n");
		*/

		if ( fabs(_residual(0)) > 1e-3) printf(" res warning\n");
		//printf(" res: (%.3e,%.e)\n", res(0), res(1));

		double phi0 = _mom(0);
		double cphi0 = cos(phi0);
		double sphi0 = sin(phi0);

		SMatrix23 matA;
		matA(0,0) = -sphi0;
		matA(0,1) =  cphi0;
		matA(0,2) = 0;
		matA(1,0) = 0;
		matA(1,1) = 0;
		matA(1,2) = 1;

		SMatrix23 matB;
		matB(0,0) = -sphi0;
		matB(0,1) =  cphi0;
		matB(0,2) = 0;
		matB(1,0) = 0;
		matB(1,1) = 0;
		matB(1,2) = 1;

		// get linearized coefficient

		return chi2;
	}

	void update(KalVtx& vtx) {
		SVector3 delta = _matGain * _residual;
		vtx.pos += delta;

		//SMatrix

		vtx.cov = vtx.cov - _matGain * matMtransp;
	}

	// compute closest point in xy plane
	SVector2 getClosestXY( SVector3 pos ) {
		double x0 = pos(0);
		double y0 = pos(1);

		double d0 = _dz(0);
		double z0 = _dz(1);

		double phi0 = _mom(0);
		double cphi0 = cos(phi0);
		double sphi0 = sin(phi0);

		double om = _mom(1);
		double R=1./om;

		double tl = _mom(2);

		double xc = -(d0+R)*sphi0;
		double yc =  (d0+R)*cphi0;

		double tandphi = ( x0*cphi0+y0*sphi0 )/( x0*sphi0-y0*cphi0+d0+R );
		double dphi = atan( tandphi );
		double cphi = cos( phi0+dphi );
		double sphi = sin( phi0+dphi );

		double x2 = -( d0+R )*sphi0+R*sphi;
		double y2 =  ( d0+R )*cphi0-R*cphi;
		double z2 = z0 + dphi*R*tl;

		double x1 = -( d0+R )*sphi0+R*sphi0;
		double y1 =  ( d0+R )*cphi0-R*cphi0;
		double z1 = z0;

		int sign1 = ( pow(x2-xc,2) + pow(y2-yc,2) > R*R ) ? 1 : -1 ;
		int sign2 = ( pow(x0-xc,2) + pow(y0-yc,2) > R*R ) ? 1 : -1 ;
		int sign3 = om/fabs(om)>0 ? 1 : -1 ;
		
		/*
		printf(" dist criteria  = %d\n", sign1);
		printf(" dist criteria2 = %d\n", sign2);
		printf(" chrg criteria  = %d\n", sign3);
		printf(" xc=%.2e, yc=%.2e\n", xc, yc);
		printf(" x1=%.2e, y1=%.2e\n", x1, y1);
		printf(" x2=%.2e, y2=%.2e\n", x2, y2);
		printf(" x0=%.2e, y0=%.2e\n", x0, y0);
		printf(" (x2-xc) = %.2e, (y2-yc) = %.2e, discr = %.2e,  R = %.2e\n",
				x2-xc, y2-yc,
				sqrt( (x2-xc)*(x2-xc) + (y2-yc)*(y2-yc) ),
				R);
		 */

		double sign = sign2*sign3;

		SVector2 trkpos;
		trkpos(0) = sqrt( x2*x2+y2*y2 )*sign;
		trkpos(1) = z2;
		return trkpos;
	}
};

void testTomohiko() {
	LCIOStorer ls("input.slcio");
	ls.InitCollections("PandoraPFOs","MCParticlesSkimmed","RecoMCTruthLink","Tracks","Neutrals","MCParticles");

	LcfiInterface interface;

	Event *event = Event::Instance();
	event->Print();

	TrackSelectorConfig _secVtxCfg;
	_secVtxCfg.maxD0 = 20;
	_secVtxCfg.maxZ0 = 20;
	_secVtxCfg.maxInnermostHitRadius = 20;
	_secVtxCfg.minVtxPlusFtdHits = 5;

	int nEvt(0);

	while(ls.Next()){
		++nEvt; //if (nEvt>20) exit(0);

		TrackVec& tracks = event->getTracks();
		TrackVec passedTracks = TrackSelector() (tracks, _secVtxCfg);
		//printf("all tracks: %d, selected: %d\n",tracks.size(),passedTracks.size());
		Vertex *ip;
		makeBeamVertex(ip);

		/*
		Vertex * vtx =  VertexFinderTearDown<vector, VertexFitterSimple>()(tracks, 0, 25., 0, ip);
		//Vertex * vtx =  VertexFinderTearDown<vector>()(tracks, 0, 25., 0, ip);

		if (vtx) {
			printf("pri_vtx found: (%.2e,%.2e,%.2e)\n",vtx->getX(),vtx->getY(),vtx->getZ());
		} else {
			printf("pri_vtx not found\n");
		}
		*/

		int ntrk(0);

		KalVtx kvtx;
		kvtx.pos(0)=1e-2;
		kvtx.pos(1)=1e-2;
		kvtx.pos(0)=0;
		kvtx.pos(1)=0;

		vector<KalTrk> kvtrks;
		for (TrackVec::const_iterator iter = passedTracks.begin(); iter != passedTracks.end(); ++iter) {
			const Track* trk = *iter;

			/*
			Helix hel(trk);
			TVector3 vec = hel.GetPos(0);

			double d0 = trk->getD0();
			double z0 = trk->getZ0();
			double om = trk->getOmega();
			double tl = trk->getTanLambda();
			double R=1./om;
			double D=d0+R;

			double phi0 = trk->getPhi();
			double cphi0 = cos(phi0);
			double sphi0 = sin(phi0);

			double x=0;
			double y=0;
			double z=0;

			{
				double dphi = 1e-8;
				double cphi = cos(phi0+dphi);
				double sphi = sin(phi0+dphi);
				x = -(d0+R)*sphi0+R*sphi;
				y =  (d0+R)*cphi0-R*cphi;
				z =  z0 + dphi*R*tl;
				printf(" point closest to \n  (x,y,z)=(%.3e,%.3e,%.3e) is: \n", x,y,z);
			}

			double tandphi = (x*cphi0+y*sphi0)/(x*sphi0-y*cphi0+D);
			double dphi = atan( tandphi );
			double phi = phi0+dphi;
			double cphi = cos(phi0+dphi);
			double sphi = sin(phi0+dphi);

			double x2 = -(d0+R)*sphi0+R*sphi;
			double y2 =  (d0+R)*cphi0-R*cphi;
			double z2 = z0 + dphi*R*tl;
			printf("  (x,y,z)=(%.3e,%.3e,%.3e) with dphi=%.3e, phi=%.3e\n", x2,y2,z2, dphi, phi );

			if (++ntrk>5) exit(0);
			*/

			KalTrk kvtrk( **iter );
			kvtrks.push_back( kvtrk );
			kvtrk.compute( kvtx );
			kvtrk.update( kvtx );
		}

	}
	printf("processed %d events.\n",nEvt);
}
