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
#include "Suehara.h"

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
vector<int> findMcJetFlavor(vector<Jet*> jets, vector<MCParticle*> mcps) {
	vector<int> mcJetFlavor;
	vector<TVector3> jetDirs;

	for (vector<Jet*>::iterator iter = jets.begin(); iter != jets.end(); ++iter) {
		Jet* jet = *iter;
		TVector3 v = jet->Vect();
		jetDirs.push_back(v.Unit());
		mcJetFlavor.push_back(1);
	}

	for (vector<MCParticle*>::iterator iter = mcps.begin(); iter != mcps.end(); ++iter) {
		MCParticle* mcp = *iter;
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

vector<MCVertex*> findMcVertex(vector<MCParticle*> mcps) {
	bool debug(false);
	vector<MCVertex*> ret;

	// map to store daughter -> semi-stable parent association
	map<MCParticle*,vector<MCParticle*> > sspDau;

	for (vector<MCParticle*>::iterator it = mcps.begin(); it != mcps.end(); ++it) {
		MCParticle* mcp = *it;
		if (mcp->isStable() == false && mcp->isSemiStable() == false)
			continue;

		MCParticle* ssp = mcp->getSemiStableParent();
		if (ssp == 0) continue;
		//printf("%d is stable with ssp %d\n",mcp->getPDG(),ssp->getPDG());

		// is the ssp very short-lived?
		float dif = ( ssp->getVertex()-mcp->getVertex() ).Mag();
		if (dif < 5e-3) {
			// if less than 5 microns, consider the vertex unresolvable
			//printf("less than 5e-3\n");
			MCParticle* sspp = ssp->getSemiStableParent();
			if (sspp) ssp = sspp;
		}
		/*
			printf("sspp=%.3e, ssp=%.3e, dif=%.3e\n",
					sspp->getVertex().Mag(),
					ssp->getVertex().Mag(),
					dif);
		 */

		map<MCParticle*,vector<MCParticle*> >::iterator sspIter = sspDau.find(ssp);
		if ( sspIter == sspDau.end() ) {
			vector<MCParticle*> a;
			a.push_back(mcp);
			sspDau.insert( make_pair(ssp, a) );
		} else {
			// find vector and push mcp
			sspDau[ssp].push_back(mcp);
		}
	}

	if (debug) {
		for (map<MCParticle*,vector<MCParticle*> >::iterator it = sspDau.begin(); it != sspDau.end(); ++it) {
			MCParticle* ssp = it->first;
			vector<MCParticle*>& dauVec = it->second;
			printf(" ssp=%d ndau=%d (start=%.3e, end=%.3e) \n",
					(int)ssp->getPDG(), (int)ssp->getDaughters().size(),
					ssp->getVertex().Mag(), ssp->getEndVertex().Mag());
			for (vector<MCParticle*>::iterator dauIter = dauVec.begin(); dauIter != dauVec.end(); ++dauIter) {
				MCParticle* dau = *dauIter;
				printf("   %d [p=%d] (start=%.3e, end=%.3e)\n", dau->getPDG(),
						dau->getParent()->getPDG(), dau->getVertex().Mag(),
						dau->getEndVertex().Mag() );
			}
		}
	}

	for (map<MCParticle*,vector<MCParticle*> >::iterator it = sspDau.begin(); it != sspDau.end(); ++it) {
		MCParticle* ssp = it->first;
		if (ssp->getEndVertex().Mag() > 1000) continue;

		int ntrk(0);
		vector<MCParticle*>& dauVec = it->second;
		for (vector<MCParticle*>::iterator dauIter = dauVec.begin(); dauIter != dauVec.end(); ++dauIter) {
			MCParticle* dau = *dauIter;
			if (dau->isStableTrack())
				++ntrk;
		}
		
		if (ntrk == 0) continue;

		MCVertex* v = new MCVertex;
		v->setParent(ssp);

		for (vector<MCParticle*>::iterator dauIter = dauVec.begin(); dauIter != dauVec.end(); ++dauIter) {
			MCParticle* dau = *dauIter;
			if (dau->isStableTrack()) {
				v->add(dau);
			}
		}
		ret.push_back(v);
	}

	if (debug) {
		for (vector<MCVertex*>::iterator it = ret.begin(); it != ret.end(); ++it) {
			MCVertex* v = *it;
			MCParticle* ssp = v->getParent();
			vector<MCParticle*> dauVec = v->getDaughters();

			printf(" ssp=%d ndau=%d (start=%.3e, end=%.3e) \n",
					(int)ssp->getPDG(), (int)ssp->getDaughters().size(),
					ssp->getVertex().Mag(), ssp->getEndVertex().Mag());
			for (vector<MCParticle*>::iterator dauIter = dauVec.begin(); dauIter != dauVec.end(); ++dauIter) {
				MCParticle* dau = *dauIter;
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
	vector<Track*> tracks = evt.getTracks();
	LcfiInterface lcfi;

	for (vector<MCVertex*>::const_iterator it = vtxList.begin(); it != vtxList.end(); ++it) {
		MCVertex* v = *it;
		table[v] = 0; // initialize flag

		MCParticle* ssp = v->getParent();
		vector<MCParticle*> dauVec = v->getDaughters();

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
		vector<Track*> trkList;

		for (vector<MCParticle*>::iterator dauIter = dauVec.begin(); dauIter != dauVec.end(); ++dauIter) {
			MCParticle* dau = *dauIter;
			Track* trk(0);
			for (vector<Track*>::iterator trkIter = tracks.begin(); trkIter != tracks.end(); ++trkIter) {
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
						Track * tr = vtx->getTracks()[i];
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
						Track * tr = vtx->getTracks()[i];
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
		const vector<Jet*>& jets)
{
	for (vector<MCVertex*>::const_iterator it = vtxList.begin(); it != vtxList.end(); ++it) {
		MCVertex* v = *it;
		vector<MCParticle*> dauVec = v->getDaughters();

		vector<Track*> trkList;
		int nRecoTrk(0);
		const vector<Track*>& tracks = evt.getTracks();

		for (vector<MCParticle*>::iterator dauIter = dauVec.begin();
				dauIter != dauVec.end(); ++dauIter) {

			MCParticle* dau = *dauIter;
			Track* trk(0);
			for (vector<Track*>::const_iterator trkIter = tracks.begin(); trkIter != tracks.end(); ++trkIter) {
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
		for (vector<Jet*>::const_iterator jetIter = jets.begin(); jetIter != jets.end(); ++jetIter) {
			Jet* j = *jetIter;
			vector<Track*> jtrks = j->getTracks();
			for (vector<Track*>::iterator trkIter = trkList.begin(); trkIter != trkList.end(); ++trkIter) {
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
	const vector<Track*>& tracks = vertex->getTracks();

	int nTrkMatchMAX(0);
	MCVertex* match(0);

	for (vector<MCVertex*>::const_iterator it = vtxList.begin(); it != vtxList.end(); ++it) {
		MCVertex* v = *it;
		vector<Track*> recoTrksMatch = v->getRecoTracks();

		int nTrkMatch(0);

		for (vector<Track*>::const_iterator trkIter = recoTrksMatch.begin();
				trkIter != recoTrksMatch.end(); ++trkIter) {
			Track* trk = *trkIter;
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
	const vector<Neutral*>& neutrals = evt.getNeutrals();

	for (vector<Neutral*>::const_iterator neutIter = neutrals.begin();
			neutIter != neutrals.end(); ++neutIter) {
		Neutral* neut = *neutIter;
		if (neut->isV0()) {
			for (vector<MCVertex*>::const_iterator it = vtxList.begin(); it != vtxList.end(); ++it) {
				MCVertex* v = *it;
				MCParticle* mcv0 = v->getParent();
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

const int MAX_NJET = 10;
const int MAX_NTRK = 200;
const int MAX_NPTRK = 200;
const int MAX_NSTRK = 200;
const int MAX_NVTX = 50;

struct Data {
	TObjString* filename;
	int nevent;
	int mcnb;
	int mcnc;
	int njet;
	float pri_x;
	float pri_y;
	float pri_z;
	float pri_xx;
	float pri_xy;
	float pri_yy;
	float pri_xz;
	float pri_yz;
	float pri_zz;

	float ycut2;
	int jet_mcflav[MAX_NJET];
	int jet_mcnlep[MAX_NJET];
	float jet_mcbflight[MAX_NJET];
	int jet_mcntrk[MAX_NJET];
	int jet_mcntrkmatch[MAX_NJET];

	int jet_nvtx[MAX_NJET];
	int jet_nvtx2[MAX_NJET];
	int jet_nvtx3[MAX_NJET];
	int jet_nvtx4[MAX_NJET];
	int jet_ntrk[MAX_NJET];
	int jet_hasmu[MAX_NJET];
	int jet_hase[MAX_NJET];
	float jet_vtxlen[MAX_NJET];
	float jet_vtxsig[MAX_NJET];
	float jet_vtxmom[MAX_NJET];
	float jet_vtxmass[MAX_NJET];
	float jet_vtxmasspc[MAX_NJET];
	int jet_vtxmult[MAX_NJET];
	float jet_vtxprob[MAX_NJET];
	float jet_px[MAX_NJET];
	float jet_py[MAX_NJET];
	float jet_pz[MAX_NJET];
	float jet_e[MAX_NJET];
	int jet_nv0[MAX_NJET];
	int jet_bbaryon[MAX_NJET];

	float jet_trk1d0sig[MAX_NJET];
	float jet_trk2d0sig[MAX_NJET];
	float jet_trk1z0sig[MAX_NJET];
	float jet_trk2z0sig[MAX_NJET];
	float jet_trk1pt[MAX_NJET];
	float jet_trk2pt[MAX_NJET];
	float jet_jprobr[MAX_NJET];
	float jet_jprobz[MAX_NJET];
	int jet_nmcvtx[MAX_NJET];

	int nptrk;
	int ptrk_mctype[MAX_NPTRK];
	float ptrk_prob[MAX_NPTRK];

	int nstrk;
	int strk_mctype[MAX_NSTRK];

	int nvtx;
	int vtx_jetindex[MAX_NVTX];
	int vtx_ntrk[MAX_NVTX];
	float vtx_dist[MAX_NVTX];
	float vtx_distsig[MAX_NVTX];
	int vtx_mcmatch[MAX_NVTX];
	int vtx_mcvtxptr[MAX_NVTX];
	float vtx_chisq[MAX_NVTX];
	float vtx_prob[MAX_NVTX];
	int vtx_type[MAX_NVTX];
	float vtx_chrg[MAX_NVTX];
	float vtx_x[MAX_NVTX];
	float vtx_y[MAX_NVTX];
	float vtx_z[MAX_NVTX];
	float vtx_xx[MAX_NVTX];
	float vtx_xy[MAX_NVTX];
	float vtx_yy[MAX_NVTX];
	float vtx_xz[MAX_NVTX];
	float vtx_yz[MAX_NVTX];
	float vtx_zz[MAX_NVTX];
	float vtx_px[MAX_NVTX];
	float vtx_py[MAX_NVTX];
	float vtx_pz[MAX_NVTX];
	float vtx_e[MAX_NVTX];
	float vtx_mass[MAX_NVTX];
	int vtx_jptr[MAX_NVTX];

	int nmcvtx;
	int mcvtx_flag[MAX_NVTX];
	int mcvtx_ntrk[MAX_NVTX];
	int mcvtx_nreco[MAX_NVTX];
	float mcvtx_dist[MAX_NVTX];
	int mcvtx_jptr[MAX_NVTX];

	int ntrk;
	int trk_pri[MAX_NTRK];
	int trk_sec[MAX_NTRK];
	int trk_pdg[MAX_NTRK];
	int trk_mcpdg[MAX_NTRK];
	int trk_mcppdg[MAX_NTRK];
	int trk_mctype[MAX_NTRK];
	float trk_mcd0[MAX_NTRK];
	float trk_mcz0[MAX_NTRK];
	float trk_mcph[MAX_NTRK];
	float trk_mcom[MAX_NTRK];
	float trk_mctd[MAX_NTRK];
	float trk_mcpt[MAX_NTRK];
	float trk_mcmom[MAX_NTRK];
	float trk_mcen[MAX_NTRK];
	float trk_mcvx[MAX_NTRK];
	float trk_mcvy[MAX_NTRK];
	float trk_mcvz[MAX_NTRK];
	float trk_mcvdist[MAX_NTRK];
	float trk_priprob[MAX_NTRK];
	//enum hit { VTX=0, FTD, SIT, TPC, SET, ETD, hitN };
	int trk_nvtx[MAX_NTRK];
	int trk_nftd[MAX_NTRK];
	int trk_nsit[MAX_NTRK];
	int trk_ntpc[MAX_NTRK];
	int trk_nset[MAX_NTRK];
	int trk_netd[MAX_NTRK];
	float trk_rimh[MAX_NTRK];
	float trk_chi2[MAX_NTRK];
	int trk_ndf[MAX_NTRK];
	float trk_prob[MAX_NTRK];
	float trk_d0[MAX_NTRK];
	float trk_z0[MAX_NTRK];
	float trk_ph[MAX_NTRK];
	float trk_om[MAX_NTRK];
	float trk_td[MAX_NTRK];
	float trk_d0d0[MAX_NTRK];
	float trk_d0ph[MAX_NTRK];
	float trk_phph[MAX_NTRK];
	float trk_d0om[MAX_NTRK];
	float trk_phom[MAX_NTRK];
	float trk_omom[MAX_NTRK];
	float trk_d0z0[MAX_NTRK];
	float trk_z0ph[MAX_NTRK];
	float trk_z0om[MAX_NTRK];
	float trk_z0z0[MAX_NTRK];
	float trk_d0td[MAX_NTRK];
	float trk_phtd[MAX_NTRK];
	float trk_omtd[MAX_NTRK];
	float trk_z0td[MAX_NTRK];
	float trk_tdtd[MAX_NTRK];
	float trk_pt[MAX_NTRK];
	float trk_mom[MAX_NTRK];
	float trk_en[MAX_NTRK];
	float trk_px[MAX_NTRK];
	float trk_py[MAX_NTRK];
	float trk_pz[MAX_NTRK];
	float trk_ecal[MAX_NTRK];
	float trk_hcal[MAX_NTRK];
	float trk_yoke[MAX_NTRK];
	float trk_lcal[MAX_NTRK];
	float trk_lhcal[MAX_NTRK];
	float trk_bcal[MAX_NTRK];
	float trk_cal[MAX_NTRK];
	int trk_mclep[MAX_NTRK];
	float trk_d0sig[MAX_NTRK];
	float trk_z0sig[MAX_NTRK];
	int trk_secvtxflag[MAX_NTRK];
	int trk_vptr[MAX_NTRK];
	int trk_jptr[MAX_NTRK];
};

vector<Vertex*> findAdditionalVertices(
		const Event& evt, const Jet* jet,
		const vector<Vertex*>& seedList,
		const Vertex* primaryVertex )
{
	vector<Vertex*> ret;

	vector<Track*> useTracks;
	vector<Track*> tracksToReject; // not require single

	const vector<Track*>& ptrks = primaryVertex->getTracks();
	for (vector<Track*>::const_iterator trkIter = ptrks.begin(); trkIter != ptrks.end(); ++trkIter) {
		tracksToReject.push_back(*trkIter);
	}

	for (vector<Vertex *>::const_iterator it = seedList.begin(); it != seedList.end(); ++it) {
		Vertex* v = *it;
		ret.push_back(v);
		const vector<Track*>& vtrks = v->getTracks();
		for (vector<Track*>::const_iterator vtrkIter = vtrks.begin(); vtrkIter != vtrks.end(); ++vtrkIter) {
			tracksToReject.push_back(*vtrkIter);
		}
	}

	const vector<Track*>& jtrks = jet->getTracks();
	for (vector<Track*>::const_iterator jtrkIter = jtrks.begin(); jtrkIter != jtrks.end(); ++jtrkIter) {
		if ( find( tracksToReject.begin(), tracksToReject.end(), *jtrkIter ) == tracksToReject.end() ) {
			useTracks.push_back(*jtrkIter);
		}
	}

	Jet jet2;
	for (vector<Track*>::iterator utrkIter = useTracks.begin(); utrkIter != useTracks.end(); ++utrkIter) {
		jet2.add(*utrkIter);
	}

	vector<Vertex*>* extraVtxList = findTearDownVertices( evt, jet2 );

	for (vector<Vertex *>::const_iterator it = extraVtxList->begin(); it != extraVtxList->end(); ++it) {
		ret.push_back(*it);
	}

	delete extraVtxList;

	return ret;
}

vector<Track*> findSingleTracks(const Event& evt, const Jet& jet, const vector<lcfiplus::Vertex*>& vtxList) {
	vector<Track*> ret;
	vector<Track*> cand;

	const vector<Track *> &tracks = jet.getTracks();
	for (unsigned int i=0; i<tracks.size(); ++i) {
		Track* trk = tracks[i];
		for (unsigned int j=0; j<vtxList.size(); ++j) {
			const vector<Track*>& vtxTracks = vtxList[j]->getTracks();
			if (find(vtxTracks.begin(),vtxTracks.end(),trk) != vtxTracks.end()) {
				cand.push_back(trk);
			}
		}
	}

	//printf("%d\n",cand.size());

	return cand;
	return ret;
}


void processEvents(const char* input, const char* output, int nStart, int nEnd) {
	printf("processEvents\n");
	//bool rescaleError = true;
	bool rescaleError = false;


	struct DebugMsg {
		bool system;
		bool verbose;
		bool primary;
		bool secondary;
		bool ptrk;
		bool strk;
		bool jets;
		bool tracks;
	} debug;
	memset(&debug, 0, sizeof(debug));
	debug.system = true;
	//debug.verbose = true;
	//debug.primary = true;
	//debug.secondary = true;
	//debug.jets = true;
	//debug.ptrk = true;
	//debug.strk = true;

	/*
		 int nStart = 0;
		 int nEnd = 99;
		 char* input = "share/test2.slcio";
		 char* output = "share/lcfiplus.root";

		 if (argc>2) {
		 nStart = atoi(argv[1]);
		 nEnd = atoi(argv[2]);
		 }

		 if (argc>3) {
		 input = argv[3];
		 }

		 if (argc>4) {
		 output = argv[4];
		 }
	 */

	int nEvents = nEnd-nStart+1;

	if (debug.system) {
		printf("Using input file: %s\n",input);
		printf("Using output file: %s\n",output);
		printf("Processing %d events (%d -> %d)\n", nEvents, nStart, nEnd);
	}

	LCIOStorer ls(input);
	ls.InitCollections();

	char fileBuf[1024];
	strncpy(fileBuf, input, 1024);
	printf("input = %s\n",input);
	printf("fileBuf = %s\n",fileBuf);
	TObjString* fileBasename = new TObjString( basename(fileBuf) );

	JetConfig jetCfg;
	//jetCfg.nJet = 2;
	jetCfg.nJet = 6;
	//JetFinder* jetFinder = new JetFinder(jetCfg);
	//CheatedJetFinder* cheatedJetFinder = new CheatedJetFinder(jetCfg);

	//*
		 SecondaryVertexConfig secVtxCfg;
	// AND
	secVtxCfg.TrackQualityCuts.maxD0 = 10;
	secVtxCfg.TrackQualityCuts.maxD0Err = 0.25;
	secVtxCfg.TrackQualityCuts.maxZ0 = 20;
	secVtxCfg.TrackQualityCuts.maxZ0Err = 1e10;
	secVtxCfg.TrackQualityCuts.minPt = 0.1;
	secVtxCfg.TrackQualityCuts.maxInnermostHitRadius = 1e10;
	// OR
	secVtxCfg.TrackQualityCuts.minTpcHits = 20;
	secVtxCfg.TrackQualityCuts.minFtdHits = 3;
	secVtxCfg.TrackQualityCuts.minVtxHitsWithoutTpcFtd = 3;
	secVtxCfg.TrackQualityCuts.minVtxPlusFtdHits = 0;
	// */
	/*
	SecondaryVertexConfig secVtxCfg;
	// AND
	secVtxCfg.maxD0 = 100;
	secVtxCfg.maxD0Err = 1e10;
	secVtxCfg.maxZ0 = 200;
	secVtxCfg.maxZ0Err = 1e10;
	secVtxCfg.minPt = 0.1;
	secVtxCfg.maxInnermostHitRadius = 20;
	// OR
	secVtxCfg.minTpcHits = 0;
	secVtxCfg.minFtdHits = 0;
	secVtxCfg.minVtxHitsWithoutTpcFtd = 0;
	secVtxCfg.minVtxPlusFtdHits = 3;
	// */

	TFile* file = new TFile(output,"RECREATE");
	TTree* t = new TTree("ntp","events");
	Data d;

	t->Branch("filename", "TObjString", &fileBasename);
	t->Branch("nevent", &d.nevent, "nevent/I");
	t->Branch("mcnb", &d.mcnb, "mcnb/I");
	t->Branch("mcnc", &d.mcnc, "mcnc/I");

	t->Branch("pri_x", &d.pri_x, "pri_x/F");
	t->Branch("pri_y", &d.pri_y, "pri_y/F");
	t->Branch("pri_z", &d.pri_z, "pri_z/F");
	t->Branch("pri_xx", &d.pri_xx, "pri_xx/F");
	t->Branch("pri_xy", &d.pri_xy, "pri_xy/F");
	t->Branch("pri_yy", &d.pri_yy, "pri_yy/F");
	t->Branch("pri_xz", &d.pri_xz, "pri_xz/F");
	t->Branch("pri_yz", &d.pri_yz, "pri_yz/F");
	t->Branch("pri_zz", &d.pri_zz, "pri_zz/F");

	t->Branch("ycut2", &d.ycut2, "ycut2/F");

	t->Branch("njet", &d.njet, "njet/I");
	t->Branch("jet_mcflav", &d.jet_mcflav, "jet_mcflav[njet]/I");
	t->Branch("jet_mcnlep", &d.jet_mcnlep, "jet_mcnlep[njet]/I");
	t->Branch("jet_mcbflight", &d.jet_mcbflight, "jet_mcbflight[njet]/F");
	t->Branch("jet_mcntrk", &d.jet_mcntrk, "jet_mcntrk[njet]/I");
	t->Branch("jet_mcntrkmatch", &d.jet_mcntrkmatch, "jet_mcntrkmatch[njet]/I");
	t->Branch("jet_nvtx", &d.jet_nvtx, "jet_nvtx[njet]/I");
	t->Branch("jet_nvtx2", &d.jet_nvtx2, "jet_nvtx2[njet]/I");
	t->Branch("jet_nvtx3", &d.jet_nvtx3, "jet_nvtx3[njet]/I");
	t->Branch("jet_nvtx4", &d.jet_nvtx4, "jet_nvtx4[njet]/I");
	t->Branch("jet_ntrk", &d.jet_ntrk, "jet_ntrk[njet]/I");
	t->Branch("jet_hase", &d.jet_hase, "jet_hase[njet]/I");
	t->Branch("jet_hasmu", &d.jet_hasmu, "jet_hasmu[njet]/I");
	t->Branch("jet_vtxlen", &d.jet_vtxlen, "jet_vtxlen[njet]/F");
	t->Branch("jet_vtxsig", &d.jet_vtxsig, "jet_vtxsig[njet]/F");
	t->Branch("jet_vtxmom", &d.jet_vtxmom, "jet_vtxmom[njet]/F");
	t->Branch("jet_vtxmass", &d.jet_vtxmass, "jet_vtxmass[njet]/F");
	t->Branch("jet_vtxmasspc", &d.jet_vtxmasspc, "jet_vtxmasspc[njet]/F");
	t->Branch("jet_vtxmult", &d.jet_vtxmult, "jet_vtxmult[njet]/I");
	t->Branch("jet_vtxprob", &d.jet_vtxprob, "jet_vtxprob[njet]/F");
	t->Branch("jet_px", &d.jet_px, "jet_px[njet]/F");
	t->Branch("jet_py", &d.jet_py, "jet_py[njet]/F");
	t->Branch("jet_pz", &d.jet_pz, "jet_pz[njet]/F");
	t->Branch("jet_e", &d.jet_e, "jet_e[njet]/F");
	t->Branch("jet_bbaryon", &d.jet_bbaryon, "jet_bbaryon[njet]/I");
	t->Branch("jet_trk1d0sig", &d.jet_trk1d0sig, "jet_trk1d0sig[njet]/F");
	t->Branch("jet_trk2d0sig", &d.jet_trk2d0sig, "jet_trk2d0sig[njet]/F");
	t->Branch("jet_trk1z0sig", &d.jet_trk1z0sig, "jet_trk1z0sig[njet]/F");
	t->Branch("jet_trk2z0sig", &d.jet_trk2z0sig, "jet_trk2z0sig[njet]/F");
	t->Branch("jet_trk1pt", &d.jet_trk1pt, "jet_trk1pt[njet]/F");
	t->Branch("jet_trk2pt", &d.jet_trk2pt, "jet_trk2pt[njet]/F");
	t->Branch("jet_jprobr", &d.jet_jprobr, "jet_jprobr[njet]/F");
	t->Branch("jet_jprobz", &d.jet_jprobz, "jet_jprobz[njet]/F");
	t->Branch("jet_nv0", &d.jet_nv0, "jet_nv0[njet]/I");
	t->Branch("jet_nmcvtx", &d.jet_nmcvtx, "jet_nmcvtx[njet]/I");

	t->Branch("nptrk",&d.nptrk,"nptrk/I");
	t->Branch("ptrk_mctype",&d.ptrk_mctype,"ptrk_mctype[nptrk]/I");
	t->Branch("ptrk_prob",&d.ptrk_prob,"ptrk_prob[nptrk]/F");

	t->Branch("nstrk",&d.nstrk,"nstrk/I");
	t->Branch("strk_mctype",&d.strk_mctype,"strk_mctype[nstrk]/I");

	t->Branch("nvtx",&d.nvtx,"nvtx/I");
	t->Branch("vtx_jetindex",&d.vtx_jetindex,"vtx_jetindex[nvtx]/I");
	t->Branch("vtx_ntrk",&d.vtx_ntrk,"vtx_ntrk[nvtx]/I");
	t->Branch("vtx_dist",&d.vtx_dist,"vtx_dist[nvtx]/F");
	t->Branch("vtx_distsig",&d.vtx_distsig,"vtx_distsig[nvtx]/F");
	t->Branch("vtx_mcmatch",&d.vtx_mcmatch,"vtx_mcmatch[nvtx]/I");
	t->Branch("vtx_mcvtxptr",&d.vtx_mcvtxptr,"vtx_mcvtxptr[nvtx]/I");
	t->Branch("vtx_chisq",&d.vtx_chisq,"vtx_chisq[nvtx]/F");
	t->Branch("vtx_prob",&d.vtx_prob,"vtx_prob[nvtx]/F");
	t->Branch("vtx_type",&d.vtx_type,"vtx_type[nvtx]/I");
	t->Branch("vtx_chrg",&d.vtx_chrg,"vtx_chrg[nvtx]/F");
	t->Branch("vtx_x",&d.vtx_x,"vtx_x[nvtx]/F");
	t->Branch("vtx_y",&d.vtx_y,"vtx_y[nvtx]/F");
	t->Branch("vtx_z",&d.vtx_z,"vtx_z[nvtx]/F");
	t->Branch("vtx_xx",&d.vtx_xx,"vtx_xx[nvtx]/F");
	t->Branch("vtx_xy",&d.vtx_xy,"vtx_xy[nvtx]/F");
	t->Branch("vtx_yy",&d.vtx_yy,"vtx_yy[nvtx]/F");
	t->Branch("vtx_xz",&d.vtx_xz,"vtx_xz[nvtx]/F");
	t->Branch("vtx_yz",&d.vtx_yz,"vtx_yz[nvtx]/F");
	t->Branch("vtx_zz",&d.vtx_zz,"vtx_zz[nvtx]/F");
	t->Branch("vtx_px",&d.vtx_px,"vtx_px[nvtx]/F");
	t->Branch("vtx_py",&d.vtx_py,"vtx_py[nvtx]/F");
	t->Branch("vtx_pz",&d.vtx_pz,"vtx_pz[nvtx]/F");
	t->Branch("vtx_e",&d.vtx_e,"vtx_e[nvtx]/F");
	t->Branch("vtx_mass",&d.vtx_mass,"vtx_mass[nvtx]/F");
	t->Branch("vtx_jptr",&d.vtx_jptr,"vtx_jptr[nvtx]/I");

	t->Branch("nmcvtx",&d.nmcvtx,"nmcvtx/I");
	t->Branch("mcvtx_flag",&d.mcvtx_flag,"mcvtx_flag[nmcvtx]/I");
	t->Branch("mcvtx_ntrk",&d.mcvtx_ntrk,"mcvtx_ntrk[nmcvtx]/I");
	t->Branch("mcvtx_nreco",&d.mcvtx_nreco,"mcvtx_nreco[nmcvtx]/I");
	t->Branch("mcvtx_dist",&d.mcvtx_dist,"mcvtx_dist[nmcvtx]/F");
	t->Branch("mcvtx_jptr",&d.mcvtx_jptr,"mcvtx_jptr[nmcvtx]/I");

	t->Branch("ntrk",&d.ntrk,"ntrk/I");
	t->Branch("trk_pri",&d.trk_pri,"trk_pri[ntrk]/I");
	t->Branch("trk_sec",&d.trk_sec,"trk_sec[ntrk]/I");
	t->Branch("trk_pdg",&d.trk_pdg,"trk_pdg[ntrk]/I");
	t->Branch("trk_mcpdg",&d.trk_mcpdg,"trk_mcpdg[ntrk]/I");
	t->Branch("trk_mcppdg",&d.trk_mcppdg,"trk_mcppdg[ntrk]/I");
	t->Branch("trk_mctype",&d.trk_mctype,"trk_mctype[ntrk]/I");
	t->Branch("trk_mcd0",&d.trk_mcd0,"trk_mcd0[ntrk]/F");
	t->Branch("trk_mcz0",&d.trk_mcz0,"trk_mcz0[ntrk]/F");
	t->Branch("trk_mcph",&d.trk_mcph,"trk_mcph[ntrk]/F");
	t->Branch("trk_mcom",&d.trk_mcom,"trk_mcom[ntrk]/F");
	t->Branch("trk_mctd",&d.trk_mctd,"trk_mctd[ntrk]/F");
	t->Branch("trk_mcpt",&d.trk_mcpt,"trk_mcpt[ntrk]/F");
	t->Branch("trk_mcmom",&d.trk_mcmom,"trk_mcmom[ntrk]/F");
	t->Branch("trk_mcen",&d.trk_mcen,"trk_mcen[ntrk]/F");
	t->Branch("trk_mcvx",&d.trk_mcvx,"trk_mcvx[ntrk]/f");
	t->Branch("trk_mcvy",&d.trk_mcvy,"trk_mcvy[ntrk]/f");
	t->Branch("trk_mcvz",&d.trk_mcvz,"trk_mcvz[ntrk]/f");
	t->Branch("trk_mcvdist",&d.trk_mcvdist,"trk_mcvdist[ntrk]/f");
	t->Branch("trk_priprob",&d.trk_priprob,"trk_priprob[ntrk]/F");
	//enum hit { VTX=0, FTD, SIT, TPC, SET, ETD, hitN };
	t->Branch("trk_nvtx",&d.trk_nvtx,"trk_nvtx[ntrk]/I");
	t->Branch("trk_nftd",&d.trk_nftd,"trk_nftd[ntrk]/I");
	t->Branch("trk_nsit",&d.trk_nsit,"trk_nsit[ntrk]/I");
	t->Branch("trk_ntpc",&d.trk_ntpc,"trk_ntpc[ntrk]/I");
	t->Branch("trk_nset",&d.trk_nset,"trk_nset[ntrk]/I");
	t->Branch("trk_netd",&d.trk_netd,"trk_netd[ntrk]/I");
	t->Branch("trk_rimh",&d.trk_rimh,"trk_rimh[ntrk]/F");
	t->Branch("trk_chi2",&d.trk_chi2,"trk_chi2[ntrk]/F");
	t->Branch("trk_ndf",&d.trk_ndf,"trk_ndf[ntrk]/I");
	t->Branch("trk_prob",&d.trk_prob,"trk_prob[ntrk]/F");
	t->Branch("trk_d0",&d.trk_d0,"trk_d0[ntrk]/F");
	t->Branch("trk_z0",&d.trk_z0,"trk_z0[ntrk]/F");
	t->Branch("trk_ph",&d.trk_ph,"trk_ph[ntrk]/F");
	t->Branch("trk_om",&d.trk_om,"trk_om[ntrk]/F");
	t->Branch("trk_td",&d.trk_td,"trk_td[ntrk]/F");
	t->Branch("trk_d0d0",&d.trk_d0d0,"trk_d0d0[ntrk]/F");
	t->Branch("trk_d0ph",&d.trk_d0ph,"trk_d0ph[ntrk]/F");
	t->Branch("trk_phph",&d.trk_phph,"trk_phph[ntrk]/F");
	t->Branch("trk_d0om",&d.trk_d0om,"trk_d0om[ntrk]/F");
	t->Branch("trk_phom",&d.trk_phom,"trk_phom[ntrk]/F");
	t->Branch("trk_omom",&d.trk_omom,"trk_omom[ntrk]/F");
	t->Branch("trk_d0z0",&d.trk_d0z0,"trk_d0z0[ntrk]/F");
	t->Branch("trk_z0ph",&d.trk_z0ph,"trk_z0ph[ntrk]/F");
	t->Branch("trk_z0om",&d.trk_z0om,"trk_z0om[ntrk]/F");
	t->Branch("trk_z0z0",&d.trk_z0z0,"trk_z0z0[ntrk]/F");
	t->Branch("trk_d0td",&d.trk_d0td,"trk_d0td[ntrk]/F");
	t->Branch("trk_phtd",&d.trk_phtd,"trk_phtd[ntrk]/F");
	t->Branch("trk_omtd",&d.trk_omtd,"trk_omtd[ntrk]/F");
	t->Branch("trk_z0td",&d.trk_z0td,"trk_z0td[ntrk]/F");
	t->Branch("trk_tdtd",&d.trk_tdtd,"trk_tdtd[ntrk]/F");
	t->Branch("trk_pt",&d.trk_pt,"trk_pt[ntrk]/F");
	t->Branch("trk_mom",&d.trk_mom,"trk_mom[ntrk]/F");
	t->Branch("trk_en",&d.trk_en,"trk_en[ntrk]/F");
	t->Branch("trk_px",&d.trk_px,"trk_px[ntrk]/F");
	t->Branch("trk_py",&d.trk_py,"trk_py[ntrk]/F");
	t->Branch("trk_pz",&d.trk_pz,"trk_pz[ntrk]/F");
	t->Branch("trk_ecal",&d.trk_ecal,"trk_ecal[ntrk]/F");
	t->Branch("trk_hcal",&d.trk_hcal,"trk_hcal[ntrk]/F");
	t->Branch("trk_yoke",&d.trk_yoke,"trk_yoke[ntrk]/F");
	t->Branch("trk_lcal",&d.trk_lcal,"trk_lcal[ntrk]/F");
	t->Branch("trk_lhcal",&d.trk_lhcal,"trk_lhcal[ntrk]/F");
	t->Branch("trk_bcal",&d.trk_bcal,"trk_bcal[ntrk]/F");
	t->Branch("trk_cal",&d.trk_cal,"trk_cal[ntrk]/F");
	t->Branch("trk_mclep",&d.trk_mclep,"trk_mclep[ntrk]/I");
	t->Branch("trk_d0sig",&d.trk_d0sig,"trk_d0sig[ntrk]/F");
	t->Branch("trk_z0sig",&d.trk_z0sig,"trk_z0sig[ntrk]/F");
	t->Branch("trk_secvtxflag",&d.trk_secvtxflag,"trk_secvtxflag[ntrk]/I");
	t->Branch("trk_vptr",&d.trk_vptr,"trk_vptr[ntrk]/I");
	t->Branch("trk_jptr",&d.trk_jptr,"trk_jptr[ntrk]/I");

	// skip events until desired starting point
	for (int iEvent=0; iEvent<nStart; ++iEvent) {
		if (!ls.Next()) {
			if (debug.system) printf("reached end of collection while skipping event\n");
			exit(0);
		}
	}

	Event *event = Event::Instance();

	for (int iEvent=nStart; iEvent<=nEnd; ++iEvent) {
		if (debug.system) printf("evt: %d\n",iEvent);

		memset(&d, 0, sizeof(d));

		d.nevent = iEvent;

		if(!ls.Next()){
			if (debug.system) printf("end of collection\n");
			break;
		}

		//EventStore::Instance()->Print();
		//Event evt;	// default constructor initializes with default collection names (Tracks, Neutrals, MCParticles)

/*
		if (rescaleError) {
			event->rescaleErrors();
		}
*/
		/*
		///////////////////////////////////
		// check bbbbbb (reject H->WW etc.)
		///////////////////////////////////
		const vector<MCParticle*>& mcps = event->getMCParticles();
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

		// select true bb events only (reject bbg as well)
		d.mcnb = event->mcNumberOfB();
		d.mcnc = event->mcNumberOfC();

		//vector<Jet*> cheatedJets = cheatedJetFinder->run(event);
		//continue;

		if (debug.verbose) {
			printf("# of tracks/neutrals: %d/%d\n",(int)event->getTracks().size(),(int)event->getNeutrals().size());
			printf("number of B: %d, number of C: %d\n",d.mcnb,d.mcnc);
			printf("looking at tracks\n");
		}


		if (event->getTracks().size() == 0) {
			continue;
		}

		const vector<Track*>& tracks = event->getTracks();

		/////////////////////////////////////////////////////////////////////////////////
		// LCFI stuff: primary and secondary vertex finding
		/////////////////////////////////////////////////////////////////////////////////

		/*
		LcfiInterface* interface = new LcfiInterface(event);
		Vertex* primaryVertex = interface->findPrimaryVertex();
		if (debug.primary) {
			printf("IP: (%f,%f,%f)\n",primaryVertex->getX(),primaryVertex->getY(),primaryVertex->getZ());
		}
		*/

		// prepare tracks for primary vertex finder
		vector<Track*> tracksForPrimary;
		for (vector<Track*>::const_iterator it = tracks.begin(); it != tracks.end(); ++it) {
			Track* trk = *it;

			if (fabs(trk->getD0())<50 && fabs(trk->getZ0())<50
				&& trk->getVtxHits() + trk->getFtdHits() >= 3
				&& trk->getRadiusOfInnermostHit() <20)
			{
				tracksForPrimary.push_back( new Track(*trk) );
			}
		}

		Vertex* primaryVertex = findPrimaryVertex(tracksForPrimary,25);
		if (primaryVertex == 0) {
			fprintf(stderr,"primary vertex could not be found (ntrk=%d), skipping event\n",(int)tracksForPrimary.size());
			continue;
		}

		/// fill primary vertex
		d.pri_x = primaryVertex->getX();
		d.pri_y = primaryVertex->getY();
		d.pri_z = primaryVertex->getZ();
		d.pri_xx = primaryVertex->getCov()[Vertex::xx];
		d.pri_xy = primaryVertex->getCov()[Vertex::xy];
		d.pri_yy = primaryVertex->getCov()[Vertex::yy];
		d.pri_xz = primaryVertex->getCov()[Vertex::xz];
		d.pri_yz = primaryVertex->getCov()[Vertex::yz];
		d.pri_zz = primaryVertex->getCov()[Vertex::zz];


		LcfiInterface interface(event,primaryVertex);
		//vector<Track*> pTracks2 = primaryVertex->getTracks();
		//map<Track*,float> probMap;
		//interface->probMap(probMap);

		//delete interface;
		//event2->setTracks(vector<Track*>());
		//delete event2;

		if (debug.primary) {
			printf("IP: (%f,%f,%f) - beam constraint; ntrks+2=%d\n",primaryVertex->getX(),primaryVertex->getY(),primaryVertex->getZ(),
					(int)primaryVertex->getTracks().size());
			for(unsigned int i=0;i<primaryVertex->getTracks().size()-2;i++){
				Track* tr = primaryVertex->getTracks()[i];
				cout << "        Track #" << i << ": p = (" << tr->Px() << "," << tr->Py() << "," << tr->Pz() << "), chi2 = "
					<< primaryVertex->getChi2Track(tr)
					//<< "\n";
					<< ", mcFlavorType = " << event->getMCParticle(tr)->getFlavorTagCategory() << "\n";
			}
		}

#if 0
		// get rid of beam tracks
		tracksForPrimary.pop_back();
		tracksForPrimary.pop_back();

		lcfiplus::Vertex * vtx = lcfiplus::VertexFinderTearDown<vector>()(tracksForPrimary, &beamTracks, 9.0, 0);
		//lcfiplus::Vertex *secvtx = lcfiplus::VertexFinderTearDown<list>()(tracksInJet,0, chi2, &residuals);
		//lcfiplus::Vertex *secvtx = lcfiplus::VertexFinderTearDown<list>()(tracksInJet,0, chi2th_secondary, &residuals);
		if (debug.primary) {
			printf("IP: (%f,%f,%f) - teardown\n",vtx->getX(),vtx->getY(),vtx->getZ());
			//Track *worstTrack = vtx->getWorstTrack();
			cout << "Primary vertex: " << vtx->getTracks().size() << "/" << tracks.size() << "\n";
			//for(unsigned int i=0;i<vtx->getTracks().size();i++){
			for(unsigned int i=0;i<vtx->getTracks().size()-2;i++){
				Track* tr = vtx->getTracks()[i];
				cout << "        Track #" << i << ": "
					<< "d0 = " << tr->getD0() << ", z0 = " << tr->getZ0()
					<< " p = (" << tr->Px() << "," << tr->Py() << "," << tr->Pz() << "), chi2 = "
					<< vtx->getChi2Track(tr)
					//<< "\n";
					<< ", mcFlavorType = " << evt.getMCParticle(tr)->getFlavorTagCategory() << "\n";
			}
		}
#endif

		vector<Track*> pTracks = primaryVertex->getTracks();
		for (unsigned int iptrk = 0; iptrk < pTracks.size(); ++iptrk) {
			if (pTracks[iptrk]->getId() > 1000000) continue;
			//d.ptrk_prob[d.nptrk] = probMap[pTracks[iptrk]];
			MCParticle* mc = event->getMCParticle( pTracks[iptrk] );
			if (mc) {
				int type = mc->getFlavorTagCategory();
				d.ptrk_mctype[d.nptrk] = type;
				if (debug.ptrk) {
					MCParticle* parent = mc->getParent();
					MCParticle* gp = parent ? parent->getParent() : 0;
					MCParticle* mcssp = mc->getSemiStableParent();
					printf("primary track type: %d (pdg=%6d, ppdg=%6d, gppdg=%6d, ssp=%6d)\n",
							type,
							mc->getPDG(),
							parent ? parent->getPDG() : 0,
							gp ? gp->getPDG() : 0,
							mcssp ? mcssp->getPDG() : 0);
				}
			} else {
				if (debug.ptrk) printf("junk primary track\n");
			}
			++d.nptrk;
		}

		if (debug.verbose) {
			printf("making jets\n");
		}

		///////////////////////////////////////////
		// make jets
		///////////////////////////////////////////
		/*
		double ymin;
		vector<Jet*> jets = jetFinder->run(event->getTracks(),event->getNeutrals(),&ymin);
		d.ycut2 = ymin;
		*/

		// need to flatten vertex info 
		vector<Jet*> sueharaJets = SueharaJetClustering(*event,jetCfg.nJet);
		vector<Jet*> jets;
		for (vector<Jet*>::iterator it = sueharaJets.begin(); it != sueharaJets.end(); ++it) {
			jets.push_back( convertJetVertex(*it) );
		}

		//vector<Jet*> jets = jetFinder->run(event->getTracks(),event->getNeutrals());

		// mcJetFlavor as computed by the same algorithm by LCFI
		// groups MC particles by angles; not reliable when
		// when number of reco jets is not the same as the number of actual jets
		vector<int> mcJetFlavor = findMcJetFlavor(jets, event->getMCParticles());

		if (debug.jets) {
			printf("Jet flavor: %d %d\n", mcJetFlavor[0], mcJetFlavor[1]);
		}

		d.njet = jets.size();

		///////////////////////////////////////////
		// mc vtx
		///////////////////////////////////////////
		vector<MCVertex*> mcVtxList = findMcVertex(event->getMCParticles());
		//vector<MCVertex*> mcVtxList;
		map<MCVertex*,int> mcVtxFlagTable;
		matchMcVertex(*event, mcVtxList, mcVtxFlagTable);
		matchMcVertexJet(*event, mcVtxList, mcVtxFlagTable, jets);
		matchMcVertexRecoV0(*event, mcVtxList, mcVtxFlagTable);

		///////////////////////////////////////////
		// mc vtx end
		///////////////////////////////////////////


		vector<Vertex*> evtZvtopVtxList;
		vector<Vertex*> evtTearDownVtxList;
		vector<Vertex*> evtExtraVtxList;
		vector<Track*> singleTracksList;

		int ivtx(0); // event count

		//=========================================
		// jet loop
		//=========================================
		map<Vertex*,int> vertexJetMap;

		for (int iJet=0; iJet<jetCfg.nJet && iJet<10; ++iJet) {

			const vector<Track*>& tracks = jets[iJet]->getTracks();

			d.jet_mcflav[iJet] = mcJetFlavor[iJet];

			d.jet_px[iJet] = jets[iJet]->Px();
			d.jet_py[iJet] = jets[iJet]->Py();
			d.jet_pz[iJet] = jets[iJet]->Pz();
			d.jet_e[iJet] = jets[iJet]->E();

			d.jet_bbaryon[iJet] = 0;

			for (unsigned int k=0; k<tracks.size(); ++k) {
				Track* trk = tracks[k];
				MCParticle* mcp = event->getMCParticle(trk);
				MCParticle* ssp = mcp->getSemiStableParent();
				if (ssp) {
					int abspdg = abs(ssp->getPDG());
					if (abspdg >= 5000 && abspdg < 6000) {
						d.jet_bbaryon[iJet] = 1;
					}
				}
			}

			map<MCParticle*,int> sspFreq;
			for (unsigned int k=0; k<tracks.size(); ++k) {
				Track* trk = tracks[k];
				MCParticle* mcp = event->getMCParticle(trk);
				if (mcp) {
					// iterate to find the first semistable particle
					// (first from the primary vertex)
					MCParticle* ssp(0);
					while ( (mcp = mcp->getSemiStableParent() ) ) {
						ssp = mcp;
					}
					if (ssp) ++sspFreq[ssp];
					//printf("trk[%d]: mcp=%p\n",k,(void*)ssp);
				}
			}

			MCParticle* ssb(0); // semi-stable B hadron

			map<MCParticle*,int>::const_iterator iter;
			for (iter = sspFreq.begin(); iter != sspFreq.end(); ++iter) {
				MCParticle* mcp = iter->first;
				//printf("%p (%d): %d\n",(void*)mcp,mcp ? mcp->getFlavor() : -1,iter->second);
				if ( mcp && mcp->getFlavor() == 5 ) {
					ssb = mcp;
					break;
				}
			}

			if (ssb) {
				d.jet_mcbflight[iJet] = ssb->decayDistance();
				if (d.jet_mcbflight[iJet] < 1e-2) {
					printf("small b flight: %e (PDG=%d)\n",d.jet_mcbflight[iJet], ssb->getPDG() );
				}

				vector<MCParticle*> mcprompt = ssb->promptTracks();
				if (debug.verbose) {
					for (unsigned int i=0; i<mcprompt.size(); ++i) {
						printf(" prompt trk [%d] pdg=%d\n",i,mcprompt[i]->getPDG());
					}
				}

				d.jet_mcntrk[iJet] = mcprompt.size();
				int nmatch(0);
				for (vector<Track*>::const_iterator it = tracks.begin(); it != tracks.end(); ++it) {
					MCParticle* mcp = event->getMCParticle(*it);
					if (mcp) {
						vector<MCParticle*>::iterator found = find( mcprompt.begin(), mcprompt.end(), mcp );
						if (found != mcprompt.end()) ++nmatch;
					}
				}
				d.jet_mcntrkmatch[iJet] = nmatch;
			}

			// secondary track loop

			d.jet_ntrk[iJet] = tracks.size();

			// lepton id from pandora...
			for (unsigned int iTrk=0; iTrk<tracks.size(); ++iTrk) {
				const Track* trk = tracks[iTrk];
				if ( abs(trk->getPDG()) == 11 ) {
					d.jet_hase[iJet] = 1;
				} else if ( abs(trk->getPDG()) == 13 ) {
					d.jet_hasmu[iJet] = 1;
				}

				MCParticle* mcp = event->getMCParticle(trk);
				if (mcp) {
					int abspdg = abs(mcp->getPDG());
					if (abspdg == 11 || abspdg == 13) {
						++d.jet_mcnlep[iJet];
					}
				}
			}

			// first find trk with best and second best significance in
			// the r-phi plane (following LCFI recipe)
			float sigVec[6];
			findMostSignificantTrack(jets[iJet],primaryVertex,sigVec);
			d.jet_trk1d0sig[iJet] = sigVec[0];
			d.jet_trk2d0sig[iJet] = sigVec[1];
			d.jet_trk1z0sig[iJet] = sigVec[2];
			d.jet_trk2z0sig[iJet] = sigVec[3];
			d.jet_trk1pt[iJet] = sigVec[4];
			d.jet_trk2pt[iJet] = sigVec[5];

			d.jet_jprobr[iJet] = jointProbD0(jets[iJet],primaryVertex);
			d.jet_jprobz[iJet] = jointProbZ0(jets[iJet],primaryVertex);

			assert(d.jet_jprobr[iJet] == d.jet_jprobr[iJet]);
			assert(d.jet_jprobz[iJet] == d.jet_jprobz[iJet]);

			if (debug.tracks) {
				printf("[%d] %f:%f:%f:%f\n",iJet,sigVec[0],sigVec[1],sigVec[2],sigVec[3]);
			}

			double e = jets[iJet]->E();
			TVector3 mom = jets[iJet]->Vect();
			if (debug.jets) {
				printf("(%d) [%-.2f, %-.2f, %-.2f, %-.2f]\n",iJet,e,mom.X(),mom.Y(),mom.Z());
			}

			////////////////////////////////////////////////////////////
			// preparing jet contents for secondary vertex finding
			////////////////////////////////////////////////////////////

			Jet jet;
			const vector<Neutral*>& jetNeutrals = jets[iJet]->getNeutrals();
			for (vector<Neutral*>::const_iterator it = jetNeutrals.begin(); it != jetNeutrals.end(); ++it) {
				jet.add( Jet(*it) );
			}
			const vector<Track*>& jetTracks = jets[iJet]->getTracks();
			for (vector<Track*>::const_iterator it = jetTracks.begin(); it != jetTracks.end(); ++it) {
				Track* trk = *it;
				const vector<Track*>& pTracks = primaryVertex->getTracks();
				if ( find(pTracks.begin(), pTracks.end(), trk) != pTracks.end() ) continue;
				//if ( probMap[trk] > 0.9 ) continue;
				//if ( fabs(trk->getD0()) < 0.1 && fabs(trk->getZ0()) < 0.1 ) continue;
				//if ( trk->getVtxHits()>0 ) continue;
				jet.add( Jet(trk) );
			}

			//vector<Vertex*> zvtopVtxList = interface->findSecondaryVertices( &jet, secVtxCfg );

			vector<Vertex*> zvtopVtxList = interface.findSecondaryVertices( jets[iJet], secVtxCfg );
			vector<Vertex*>* tearDownVtxList = findTearDownVertices( *event, jet );
			vector<Vertex*> moreVtxList1 = findAdditionalVertices( *event, jets[iJet], zvtopVtxList, primaryVertex );

			//vector<Vertex*> zvtopVtxList;
			//vector<Vertex*>* tearDownVtxList = new vector<Vertex*>;
			//vector<Vertex*> moreVtxList1;

			if (debug.strk) {
				printf(" zvtop.size=%d, morevtx.size=%d\n", (int)zvtopVtxList.size(), (int)moreVtxList1.size() );
			}

			if (debug.strk) {
				for (unsigned int j=0; j<zvtopVtxList.size(); ++j) {
					double x = zvtopVtxList[j]->getX();
					double y = zvtopVtxList[j]->getY();
					double z = zvtopVtxList[j]->getZ();
					int ntrk = zvtopVtxList[j]->getTracks().size();
					printf("  zvtop vertex [%2d] %d tracks, pos=(%6.3e,%6.3e,%6.3e)\n",j,ntrk,x,y,z);
				}

				for (unsigned int j=0; j<tearDownVtxList->size(); ++j) {
					double x = (*tearDownVtxList)[j]->getX();
					double y = (*tearDownVtxList)[j]->getY();
					double z = (*tearDownVtxList)[j]->getZ();
					int ntrk = (*tearDownVtxList)[j]->getTracks().size();
					printf("  teardown vertex [%2d] %d tracks, pos=(%6.3e,%6.3e,%6.3e)\n",j,ntrk,x,y,z);
				}
			}

			//
			////const vector<Vertex*>& vtxList = *tearDownVtxList;

			//const vector<Vertex*>& vtxList = zvtopVtxList;
			//const vector<Vertex*>& vtxList2 = *tearDownVtxList;

			const vector<Vertex*>& vtxList = moreVtxList1;
			//const vector<Vertex*>& vtxList2 = zvtopVtxList;

			//const vector<Vertex*>& vtxList = zvtopVtxList;
			//const vector<Vertex*>& vtxList2 = moreVtxList1;

			//const vector<Vertex*>& vtxList = sueharaJets[iJet]->getVertices();

			//////////////////////////////////////////
			// match reco-mc vertex
			for (unsigned int i=0; i<vtxList.size(); ++i) {
				matchMcVertexReco(*event, mcVtxList, mcVtxFlagTable, vtxList[i]);
			}
			//////////////////////////////////////////

			vector<Track*> singleTracks = findSingleTracks( *event, jet, vtxList );

			evtZvtopVtxList.insert( evtZvtopVtxList.end(), zvtopVtxList.begin(), zvtopVtxList.end() );
			evtTearDownVtxList.insert( evtTearDownVtxList.end(), tearDownVtxList->begin(), tearDownVtxList->end() );
			evtExtraVtxList.insert( evtExtraVtxList.end(), moreVtxList1.begin(), moreVtxList1.end() );
			singleTracksList.insert( singleTracksList.end(), singleTracks.begin(), singleTracks.end() );

			double bestSig(-1);
			double bestLen(-1);
			int vtxmult(0);
			TLorentzVector vtxmomSum;
			double oneMinusProb(1.);

			// the vertex loop
			for (unsigned int j=0; j<vtxList.size(); ++j) {
				double len = vtxList[j]->length(primaryVertex);
				double sig = vtxList[j]->significance(primaryVertex);
				if (sig > bestSig) {
					bestSig = sig;
					bestLen = len;
				}

				const vector<Track*>& vtxTracks = vtxList[j]->getTracks();
				vtxmult += vtxTracks.size();

				for (unsigned int k=0; k<vtxTracks.size(); ++k) {
					Track* trk = vtxTracks[k];
          TLorentzVector vtxmom = *trk;
          vtxmomSum += vtxmom;
				}

				// if multiple vertices are found, combine probability by
				// 1-pALL = (1-p1)*(1-p2)*(1-p3)*...
				int ndf = vtxTracks.size()*2 - 3;
				assert(ndf>0);
				double prob = TMath::Prob(vtxList[j]->getChi2(),ndf);
				oneMinusProb *= (1-prob);

				d.vtx_jetindex[ivtx] = iJet;
				++ivtx;
			}

			/*
			int nvtxDecSig(0);
			for (unsigned int j=0; j<vtxList.size(); ++j) {
				double sig = vtxList[j]->significance(primaryVertex);
				if (sig>3) ++nvtxDecSig;
			}

			int nvtxDecSig2(0);
			for (unsigned int j=0; j<vtxList2.size(); ++j) {
				double sig = vtxList2[j]->significance(primaryVertex);
				if (sig>3) ++nvtxDecSig2;
			}
			*/

			d.jet_nvtx[iJet] = vtxList.size();
			/*
			d.jet_nvtx2[iJet] = vtxList2.size();
			d.jet_nvtx3[iJet] = nvtxDecSig;
			d.jet_nvtx4[iJet] = nvtxDecSig2;
			*/
			d.jet_vtxlen[iJet] = bestLen;
			d.jet_vtxsig[iJet] = bestSig;
			d.jet_vtxmult[iJet] = vtxmult;
			d.jet_vtxmom[iJet] = vtxmomSum.Vect().Mag();
			d.jet_vtxmass[iJet] = vtxmomSum.M();
			// should really order them by distance
			if ( vtxList.size() > 0 ) {
				float pt = interface.vertexMassPtCorrection(vtxList[0],primaryVertex,vtxmomSum.Vect(),2);
				float vm = vtxmomSum.M();
				d.jet_vtxmasspc[iJet] = sqrt( vm*vm+pt*pt ) + pt;
			}
			d.jet_vtxprob[iJet] = 1-oneMinusProb;

			int nv0(0);
			const vector<Neutral*>& neuts = jets[iJet]->getNeutrals();
			for (unsigned int ineu=0; ineu<neuts.size(); ++ineu) {
				if (neuts[ineu]->isV0()) {
					++nv0;
				}
			}
			d.jet_nv0[iJet] = nv0;

			for (unsigned int j=0; j<vtxList.size(); ++j) {
				vertexJetMap.insert( make_pair(vtxList[j],iJet) );
			}

			delete tearDownVtxList;
		}

		// match mc vertex by angle... this is dangerous if jet clustering is bad!
		// should only be used for 2-jet samples
		map<int,int> mcVertexJetMap;
		for (unsigned int imcv=0; imcv<mcVtxList.size(); ++imcv) {
			MCVertex* mcv = mcVtxList[imcv];
			MCParticle* mcp = mcv->getParent();
			TVector3 vtxvec( mcp->getEx(), mcp->getEy(), mcp->getEz() );
			if (vtxvec.Mag() == 0) {
				fprintf(stderr,"mc vertex has zero\n");
				continue;
			}
			vtxvec.SetMag(1);

			float maxval(-1);
			int index=0;

			for (int iJet=0; iJet<jetCfg.nJet && iJet<10; ++iJet) {
				Jet* jet = jets[iJet];
				TVector3 jetvec( jet->Px(), jet->Py(), jet->Pz() );
				jetvec.SetMag(1);

				float val = vtxvec.Dot(jetvec);
				if (val > maxval) {
					maxval = val;
					index = iJet;
				}
			}
			mcVertexJetMap.insert( make_pair( imcv, index ) );
		}
		for (int iJet=0; iJet<jetCfg.nJet && iJet<10; ++iJet) {
			int nmcvtx(0);
			for (unsigned int imcv=0; imcv<mcVtxList.size(); ++imcv) {
				if ( mcVertexJetMap.find( imcv )->second == iJet ) {
					++nmcvtx;
				}
			}
			d.jet_nmcvtx[iJet] = nmcvtx;
		}

		/////////////////////////////////////////////////////////////////////////////////
		// jet loop end
		/////////////////////////////////////////////////////////////////////////////////

		const vector<Vertex*>& evtVtxList = evtExtraVtxList;

		/////////////////////////////////////////////////////////////////////////////////
		// the vertex loop
		/////////////////////////////////////////////////////////////////////////////////
		for (unsigned int ivtx=0; ivtx<evtVtxList.size(); ++ivtx) {
			const vector<Track*>& vtxTracks = evtVtxList[ivtx]->getTracks();

			for (unsigned int istrk=0; istrk<vtxTracks.size(); ++istrk) {
				const Track* trk = vtxTracks[istrk];
				MCParticle* mc = event->getMCParticle( trk );
				if (mc) {
					d.strk_mctype[d.nstrk] = mc->getFlavorTagCategory();
				} else {
					if (debug.strk) {
						printf("junk secondary track\n");
					}
				}
				++d.nstrk;
			}

			int mcmatch(0);
			int mcvtxptr(-1);
			for (unsigned int imcv=0; imcv<mcVtxList.size(); ++imcv) {
				MCVertex* mcv = mcVtxList[imcv];
				if (mcv->getRecoVertex() && mcv->getRecoVertex() == evtVtxList[ivtx]) {
					mcmatch = 1;
					mcvtxptr = imcv;
					break;
				}
			}

			d.vtx_ntrk[ivtx] = evtVtxList[ivtx]->getTracks().size();
			d.vtx_dist[ivtx] = evtVtxList[ivtx]->length();
			d.vtx_distsig[ivtx] = evtVtxList[ivtx]->significance(primaryVertex);
			d.vtx_mcmatch[ivtx] = mcmatch;
			d.vtx_mcvtxptr[ivtx] = mcvtxptr;
			d.vtx_chisq[ivtx] = evtVtxList[ivtx]->getChi2();
			int ndf = vtxTracks.size()*2 - 3;
			double prob = TMath::Prob(evtVtxList[ivtx]->getChi2(),ndf);
			d.vtx_prob[ivtx] = prob;

			int type(0);
			if (find( evtZvtopVtxList.begin(), evtZvtopVtxList.end(), evtVtxList[ivtx] ) == evtZvtopVtxList.end() ) {
				type = 1;
			}
			d.vtx_type[ivtx] = type;

			float totChrg(0);
			for (unsigned int itrk=0; itrk != evtVtxList[ivtx]->getTracks().size(); ++itrk) {
				Track* t = evtVtxList[ivtx]->getTracks()[itrk];
				totChrg += t->getCharge();
			}
			d.vtx_chrg[ivtx] = totChrg;

			TVector3 vec;
			float px(0),py(0),pz(0),e(0);
			for (unsigned int itrk=0; itrk != evtVtxList[ivtx]->getTracks().size(); ++itrk) {
				Track* t = evtVtxList[ivtx]->getTracks()[itrk];
				float px0 = t->Px();
				float py0 = t->Py();
				float pz0 = t->Pz();
				px += px0;
				py += py0;
				pz += pz0;
				e += sqrt( px0*px0+py0*py0+pz0*pz0 + PION_MASS2 );
			}

			d.vtx_x[ivtx] = evtVtxList[ivtx]->getX();
			d.vtx_y[ivtx] = evtVtxList[ivtx]->getY();
			d.vtx_z[ivtx] = evtVtxList[ivtx]->getZ();
			d.vtx_xx[ivtx] = evtVtxList[ivtx]->getCov()[Vertex::xx];
			d.vtx_xy[ivtx] = evtVtxList[ivtx]->getCov()[Vertex::xy];
			d.vtx_yy[ivtx] = evtVtxList[ivtx]->getCov()[Vertex::yy];
			d.vtx_xz[ivtx] = evtVtxList[ivtx]->getCov()[Vertex::xz];
			d.vtx_yz[ivtx] = evtVtxList[ivtx]->getCov()[Vertex::yz];
			d.vtx_zz[ivtx] = evtVtxList[ivtx]->getCov()[Vertex::zz];
			d.vtx_px[ivtx] = px;
			d.vtx_py[ivtx] = py;
			d.vtx_pz[ivtx] = pz;
			d.vtx_e[ivtx] = e;
			d.vtx_mass[ivtx] = sqrt( e*e-px*px-py*py-pz*pz );

			int jptr = -1;
			map<Vertex*,int>::iterator iterVertexJetMap = vertexJetMap.find( evtVtxList[ivtx] );
			if (iterVertexJetMap != vertexJetMap.end() ) {
				jptr = iterVertexJetMap->second;
			}
			d.vtx_jptr[ivtx] = jptr;
		}
		/////////////////////////////////////////////////////////////////////////////////
		// the vertex end
		/////////////////////////////////////////////////////////////////////////////////


		/////////////////////////////////////////////////////////////////////////////////
		// fill mc vertex info
		/////////////////////////////////////////////////////////////////////////////////
		d.nmcvtx = mcVtxList.size();
		for (int i=0; i<d.nmcvtx; ++i) {
			MCVertex* mcv = mcVtxList[i];
			d.mcvtx_ntrk[i] = mcv->getDaughters().size();
			d.mcvtx_nreco[i] = mcv->getRecoTracks().size();
			d.mcvtx_flag[i] = mcVtxFlagTable[mcv];
			d.mcvtx_dist[i] = mcv->getParent()->getEndVertex().Mag();
			d.mcvtx_jptr[i] = mcVertexJetMap.find(i)->second;
		}


		/////////////////////////////////////////////////////////////////////////////////
		// fill tracks info into ntuple
		/////////////////////////////////////////////////////////////////////////////////
		d.ntrk = tracks.size();
		for (unsigned int itrk=0; itrk < tracks.size(); ++itrk) {
			int mclep(0);

			MCParticle* mc = event->getMCParticle( tracks[itrk] );
			if (mc) {
				d.trk_mcpdg[itrk] = mc->getPDG();
				MCParticle* mcparent = mc->getParent();
				if (mcparent) d.trk_mcppdg[itrk] = mcparent->getPDG();
				switch (abs(mc->getPDG())) {
					case 11: mclep = 2; break;
					case 13: mclep = 3; break;
					default: mclep = 1; break;
				}
				d.trk_mclep[itrk] = mclep;
				d.trk_mctype[itrk] = mc->getFlavorTagCategory();

        d.trk_mcmom[itrk] = mc->P();
        d.trk_mcpt[itrk] = mc->Pt();
        d.trk_mcen[itrk] = mc->E();

        d.trk_mcvx[itrk] = mc->getVertex().x();
        d.trk_mcvy[itrk] = mc->getVertex().y();
        d.trk_mcvz[itrk] = mc->getVertex().z();
        d.trk_mcvdist[itrk] = mc->getVertex().Mag();

				d.trk_mcd0[itrk] = mc->getD0();
				d.trk_mcz0[itrk] = mc->getZ0();
				d.trk_mcph[itrk] = mc->getPhi();
				d.trk_mcom[itrk] = mc->getOmega();
				d.trk_mctd[itrk] = mc->getTanLambda();
			}

			//d.trk_priprob[itrk] = probMap[tracks[itrk]];

			d.trk_pri[itrk] = 0;
			vector<Track*> pTracks = primaryVertex->getTracks();
			for (unsigned int iptrk = 0; iptrk < pTracks.size(); ++iptrk) {
				if (tracks[itrk]->getId() == pTracks[iptrk]->getId()) {
					d.trk_pri[itrk] = 1;
					break;
				}
			}
			d.trk_sec[itrk] = 0;

			d.trk_rimh[itrk] = tracks[itrk]->getRadiusOfInnermostHit();
			d.trk_chi2[itrk] = tracks[itrk]->getChi2();
			d.trk_ndf[itrk] = tracks[itrk]->getNdf();
			d.trk_prob[itrk] = TMath::Prob( d.trk_chi2[itrk], d.trk_ndf[itrk] );

			d.trk_d0[itrk] = tracks[itrk]->getD0();
			d.trk_z0[itrk] = tracks[itrk]->getZ0();
			d.trk_ph[itrk] = tracks[itrk]->getPhi();
			d.trk_om[itrk] = tracks[itrk]->getOmega();
			d.trk_td[itrk] = tracks[itrk]->getTanLambda();

			const float* cov = tracks[itrk]->getCovMatrix();
			d.trk_d0d0[itrk] = cov[tpar::d0d0];
			d.trk_d0ph[itrk] = cov[tpar::d0ph];
			d.trk_phph[itrk] = cov[tpar::phph];
			d.trk_d0om[itrk] = cov[tpar::d0om];
			d.trk_phom[itrk] = cov[tpar::phom];
			d.trk_omom[itrk] = cov[tpar::omom];
			d.trk_d0z0[itrk] = cov[tpar::d0z0];
			d.trk_z0ph[itrk] = cov[tpar::z0ph];
			d.trk_z0om[itrk] = cov[tpar::z0om];
			d.trk_z0z0[itrk] = cov[tpar::z0z0];
			d.trk_d0td[itrk] = cov[tpar::d0td];
			d.trk_phtd[itrk] = cov[tpar::phtd];
			d.trk_omtd[itrk] = cov[tpar::omtd];
			d.trk_z0td[itrk] = cov[tpar::z0td];
			d.trk_tdtd[itrk] = cov[tpar::tdtd];

			//enum hit { VTX=0, FTD, SIT, TPC, SET, ETD, hitN };
			d.trk_nvtx[itrk] = tracks[itrk]->getVtxHits();
			d.trk_nftd[itrk] = tracks[itrk]->getFtdHits();
			d.trk_nsit[itrk] = tracks[itrk]->getSitHits();
			d.trk_ntpc[itrk] = tracks[itrk]->getTpcHits();
			d.trk_nset[itrk] = tracks[itrk]->getSetHits();
			d.trk_netd[itrk] = tracks[itrk]->getEtdHits();

			/*
				 enum cov { d0d0=0, d0ph, phph, d0om, phom, omom, d0z0,
				 z0ph, z0om, z0z0, d0td, phtd, omtd, z0td, tdtd, covN };
			 */

			d.trk_pdg[itrk] = tracks[itrk]->getPDG();
      d.trk_pt[itrk] = tracks[itrk]->Pt();
      d.trk_mom[itrk] = tracks[itrk]->P();
      d.trk_en[itrk] = tracks[itrk]->E();
      d.trk_px[itrk] = tracks[itrk]->Px();
      d.trk_py[itrk] = tracks[itrk]->Py();
      d.trk_pz[itrk] = tracks[itrk]->Pz();

      d.trk_ecal[itrk] = tracks[itrk]->getCaloEdep()[tpar::ecal];
      d.trk_hcal[itrk] = tracks[itrk]->getCaloEdep()[tpar::hcal];
      d.trk_yoke[itrk] = tracks[itrk]->getCaloEdep()[tpar::yoke];
      d.trk_lcal[itrk] = tracks[itrk]->getCaloEdep()[tpar::lcal];
      d.trk_lhcal[itrk] = tracks[itrk]->getCaloEdep()[tpar::lhcal];
      d.trk_bcal[itrk] = tracks[itrk]->getCaloEdep()[tpar::bcal];
      float caloTot(0);
      for (int ical=0; ical<tpar::caloN; ++ical) {
        caloTot += tracks[itrk]->getCaloEdep()[ical];
      }
      d.trk_cal[itrk] = caloTot;

			d.trk_d0sig[itrk] = trackD0Significance(tracks[itrk],primaryVertex);
			d.trk_z0sig[itrk] = trackZ0Significance(tracks[itrk],primaryVertex);


			// match track with secondary vertices
			vector<Vertex*>::iterator it;
			bool zvtopMatch(false);
			bool tearDownMatch(false);
			bool singleTrackMatch(false);
			bool extraMatch(false);

			for (it = evtZvtopVtxList.begin();
					!zvtopMatch &&
					it != evtZvtopVtxList.end(); ++it) {
				Vertex* vtx = *it;
				vector<Track*> vtrks = vtx->getTracks();
				for (vector<Track*>::iterator it = vtrks.begin();
						!zvtopMatch &&
						it != vtrks.end(); ++it)
				{
					Track* vtrk = *it;
					if (tracks[itrk] == vtrk) {
						zvtopMatch = true;
					}
				}
			}
			for (it = evtTearDownVtxList.begin();
					!tearDownMatch &&
					it != evtTearDownVtxList.end(); ++it) {
				Vertex* vtx = *it;
				vector<Track*> vtrks = vtx->getTracks();
				for (vector<Track*>::iterator it = vtrks.begin();
						!tearDownMatch &&
						it != vtrks.end(); ++it)
				{
					Track* vtrk = *it;
					if (tracks[itrk] == vtrk) {
						tearDownMatch = true;
					}
				}
			}
			for (it = evtExtraVtxList.begin();
					!extraMatch &&
					it != evtExtraVtxList.end(); ++it) {
				Vertex* vtx = *it;
				vector<Track*> vtrks = vtx->getTracks();
				for (vector<Track*>::iterator it = vtrks.begin();
						!extraMatch &&
						it != vtrks.end(); ++it)
				{
					Track* vtrk = *it;
					if (tracks[itrk] == vtrk) {
						extraMatch = true;
					}
				}
			}

			for (vector<Track*>::iterator it = singleTracksList.begin();
					it != singleTracksList.end(); ++it) {
				if (tracks[itrk] == *it) {
					singleTrackMatch = true;
					break;
				}
			}

			if (zvtopMatch)    d.trk_secvtxflag[itrk] |= 1;
			if (tearDownMatch) d.trk_secvtxflag[itrk] |= 2;
			if (extraMatch) d.trk_secvtxflag[itrk] |= 4;
			if (singleTrackMatch) d.trk_secvtxflag[itrk] |= 8;


			int vptr = -1;
			for (unsigned int ivtx=0; ivtx<evtVtxList.size(); ++ivtx) {
				const vector<Track*> tList = evtVtxList[ivtx]->getTracks();
				if ( find( tList.begin(), tList.end(), tracks[itrk] ) != tList.end() ) {
					vptr = ivtx;
					break;
				}
			}
			d.trk_vptr[itrk] = vptr;

			int jptr = -1;
			for (unsigned int ijet=0; ijet<jets.size(); ++ijet) {
				const vector<Track*> tList = jets[ijet]->getTracks();
				if ( find( tList.begin(), tList.end(), tracks[itrk] ) != tList.end() ) {
					jptr = ijet;
					break;
				}
			}
			d.trk_jptr[itrk] = jptr;
		}

		d.nvtx = ivtx;

		t->Fill();

		// clean up

		while (!jets.empty()) {
			delete jets.back();
			jets.pop_back();
		}
		while (!sueharaJets.empty()) {
			delete sueharaJets.back();
			sueharaJets.pop_back();
		}


		delete primaryVertex;
	}

	file->Write();
	file->Close();

	//delete jetFinder;
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

void lcfiplus_default()
{
	//	lcioToTree();
	processEvents("/data5/soft/samples/grid/users/walsh/zpole/samples/reconstructed/ILD_00/01-06/DST-01-06_M-06-07-p01_zpole-RAL_90GeV_Z_to_qq_ILD_00_LCPhys_0001.slcio","test.root",0,1000);
}

int main(int argc, char* argv[]) {
	if (argc == 5) {
		processEvents(argv[1],argv[2],atoi(argv[3]),atoi(argv[4]));
		return 0;
	}

	// if argc==1 the default command is executed
	// if argc==2 argv[1] macro is processed
	// if argc==3 argv[2] function in argv[1] file is processed
	// if argc==4 argv[2] function in argv[1] file is processed with param argv[3]
	// if argc>=5 abort 

	if(argc>=5){cout << "usage: lcfiplus [ filename [ funcname [ params ] ] ]" << endl; return 1;}
	if(argc==1){lcfiplus_default();return 0;}
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
