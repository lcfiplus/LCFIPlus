#include <iostream>
#include <math.h>
#include "EventNavigator.h"
#include "EventStore.h"
#include "lcfiplus.h"
#include "JetFinder.h"
#include "LcfiInterface.h"
#include "LCIOStorer.h"
#include "VertexFinderTearDown.h"
#include "VertexFinderSuehara.h"
#include "algoEtc.h"

#include "TError.h"
#include "TEveArrow.h"
#include "TEveBoxSet.h"
#include "TEveBrowser.h"
#include "TEveCompound.h"
#include "TEveGeoNode.h"
#include "TEveManager.h"
#include "TEvePointSet.h"
#include "TEveRGBAPalette.h"
#include "TEveStraightLineSet.h"
#include "TEveTrack.h"
#include "TEveTrackPropagator.h"
#include "TEveTrans.h"
#include "TEveVSDStructs.h"

using namespace lcfiplus;
using namespace lcfiplus::algoEtc;

namespace lcfiplus {

#ifdef BUILD_EVE

  EventNavigator::EventNavigator(const char* input, int start) : _start(start) {
		_ls = new LCIOStorer(input);
		_ls->InitCollections();
		for(int i=0;i<start;i++)if(!_ls->Next()){printf("end of collection before start\n");return;}
  }

  EventNavigator::~EventNavigator() {
		delete _ls;
	}

  void EventNavigator::Fwd() {
    std::cout << "Fwd()\n";

		if (!_ls->Next()) {
      printf("end of collection\n");
			return;
		}

    int nB(0);

		// default constructor initializes with default collection names (Tracks, Neutrals, MCParticles)
		Event* event = Event::Instance();

		/*
    while (nB != 2) {
      nB = event->mcNumberOfB();
			delete event;
			if (!_ls->Next()) {
				printf("end of collection\n");
				return;
			}
			event = new Event;
    }
		*/

    if (event) {
      //event->rescaleErrors();
      vector<Track*> tracks = event->getTracks();
      printf(" nB=%d, ntrks=%d\n",nB,(int)tracks.size());
      drawEvent(event);
    } else {
      printf("end of collection\n");
    }
  }

  void EventNavigator::Bck() {
    std::cout << "Bck() - currently not implemented :(\n";
  }

  void EventNavigator::drawEvent(Event* event) {
    const float MAX_R(1000);
    const float MAX_Z(1000);
    //const float MAX_R(100);
    //const float MAX_Z(100);
    //const float MAX_R(10);
    //const float MAX_Z(10);

		const int njet = 2;

    vector<MCParticle*> mcps = event->getMCParticles();
    vector<Track*> tracks = event->getTracks();

    JetConfig jetCfg;
    jetCfg.nJet = njet;
    JetFinder jetFinder(jetCfg);
    vector<Jet*> jets = jetFinder.run(event->getTracks(),event->getNeutrals());

		TrackSelectorConfig secVtxCfg;

		// AND
		secVtxCfg.maxD0 = 10;
		secVtxCfg.maxD0Err = .1;
		secVtxCfg.maxZ0 = 20;
		secVtxCfg.maxZ0Err = .1;
		secVtxCfg.minPt = 0.1;
		secVtxCfg.maxInnermostHitRadius = 1e10;
		// OR
		secVtxCfg.minTpcHits = 20;
		secVtxCfg.minFtdHits = 3;
		secVtxCfg.minVtxHitsWithoutTpcFtd = 3;
		secVtxCfg.minVtxPlusFtdHits = 0;

		/////////////////////////////////////////////////////////////////////////////////
		// LCFI stuff: primary and secondary vertex finding
		/////////////////////////////////////////////////////////////////////////////////
/*
		//Vertex* primaryVertex = interface.findPrimaryVertex();

		// prepare tracks for primary vertex finder
		vector<Track*> tracksForPrimary;
		for (vector<Track*>::const_iterator it = tracks.begin(); it != tracks.end(); ++it) {
			Track* trk = *it;

			//tracksForPrimary.push_back( trk );
			//if (fabs(trk->getD0())<50 && fabs(trk->getZ0())<50) {
			if (fabs(trk->getD0())<10 && fabs(trk->getZ0())<10) {
				tracksForPrimary.push_back( trk );
			}
		}

		vector<Track*> beamTracks;
		beamTracks.resize(2);
		makeBeamTracks(beamTracks[0], beamTracks[1]);

		tracksForPrimary.push_back( beamTracks[0] );
		tracksForPrimary.push_back( beamTracks[1] );

		Event* event2 = new Event(*event);
		event2->deleteTracks();
		event2->setTracks( tracksForPrimary );

		LcfiInterface interface2(event2);
		Vertex* primaryVertex2 = interface2.findPrimaryVertex();
*/

    vector<Track*> vtxTracks;
    //vector<Vertex*> vtxList[njet];

    //vector<Vertex*> tearDownVtxList[njet];
    vector<Track*> tearDownVtxTracks;

		vector<Vertex*> buildUpVtxList;
    vector<Track*> vtxTracksBU;
 
		// buildup vertex
		vector<Track *> passedTracks = TrackSelector() (tracks, secVtxCfg);
//		for(unsigned int i=0;i<tracks.size();i++){
//			if(interface.passesCut(tracks[i], secVtxCfg))
//				passedTracks.push_back(tracks[i]);
//		}
//		VertexFinderSuehara::buildUp(tracks, buildUpVtxList, 25, 10);
		VertexFinderSuehara::buildUp(passedTracks, buildUpVtxList, 25, 9, 10,0.3,1.0);
		VertexFinderSuehara::associateIPTracks(buildUpVtxList,0,0,2);

		cout << "Buildup: " << buildUpVtxList.size() << " vertices found." << endl;
		for(unsigned int i=0;i<buildUpVtxList.size();i++){
			cout << "position: (" << buildUpVtxList[i]->getX() << ", " << buildUpVtxList[i]->getY() << ", " << buildUpVtxList[i]->getZ() << "), ";
			cout << "tracks: ";
			for(unsigned int j=0;j<buildUpVtxList[i]->getTracks().size();j++){
				Track *tr = buildUpVtxList[i]->getTracks()[j];
				cout << tr->getId() << ": " << tr->getMcp()->getId() << " : " << buildUpVtxList[i]->getChi2Track(tr) << ", ";

				if(i>0)
					vtxTracksBU.push_back(tr);
			}
			cout << ", chi2 = " << buildUpVtxList[i]->getChi2();
			cout << endl;

			if(i>0&&i<buildUpVtxList.size()-1){
				for(unsigned int j=i+1; j<buildUpVtxList.size();j++){
					bool near = (buildUpVtxList[j]->getPos().Mag() > buildUpVtxList[i]->getPos().Mag());

					int nnear = (near ? i : j);
					int nfar = (near ? j : i);

					cout << "Vertex combination " << nnear << "&" << nfar << ": ";
					cout << "pmass: " << buildUpVtxList[nnear]->getVertexMass(buildUpVtxList[nfar]);
					TVector3 daxis = (buildUpVtxList[nfar]->getPos() - buildUpVtxList[nnear]->getPos()).Unit();
					cout << ", dmass: " << buildUpVtxList[nfar]->getVertexMass(0,&daxis) << endl;
				}
			}
		}
/*
    for (unsigned int iJet=0; iJet<jets.size(); ++iJet) {

			// zvtop
      vtxList[iJet] = interface.findSecondaryVertices( jets[iJet], secVtxCfg );
      printf(" jet %d has %d vertices (zvtop)\n",iJet+1,(int)vtxList[iJet].size() );

      for (unsigned int iVtx=0; iVtx<vtxList[iJet].size(); ++iVtx) {
        const vector<Track*>& trks = vtxList[iJet][iVtx]->getTracks();
        vtxTracks.insert( vtxTracks.end(), trks.begin(), trks.end() );
      }

      for (unsigned int iVtx=0; iVtx<vtxList[iJet].size(); ++iVtx) {
				Vertex* vtx = vtxList[iJet][iVtx];
				printf("pos=(%.3e,%.3e,%.3e) chi2=%.3e\n",
						vtx->getX(),vtx->getY(),vtx->getZ(),
						vtx->getChi2()
						);
				for(unsigned int i=0;i<vtx->getTracks().size();i++){
					Track * tr = vtx->getTracks()[i];
					printf("        Track #%d  pos=(%.3e,%.3e,%.3e), chi2=%e\n",
							i, tr->getX(), tr->getY(), tr->getZ(), vtx->getChi2Track(tr) );
				}
			}


			////////////////////////////////////////////////////////////
			// preparing jet contents for secondary vertex finding      ////////////////////////////////////////////////////////////

			Jet jet;
			const vector<Neutral*>& jetNeutrals = jets[iJet]->getNeutrals();
			for (vector<Neutral*>::const_iterator it = jetNeutrals.begin(); it != jetNeutrals.end(); ++it)
			{  
				jet.add( Jet(*it) );
			}
			const vector<Track*>& jetTracks = jets[iJet]->getTracks();
			for (vector<Track*>::const_iterator it = jetTracks.begin(); it != jetTracks.end(); ++it) {
				Track* trk = *it;
				const vector<Track*>& pTracks = primaryVertex2->getTracks();
				if ( find(pTracks.begin(), pTracks.end(), trk) != pTracks.end() ) continue;
				jet.add( Jet(trk) );
			}

			vector<Vertex*>* tmp = findTearDownVertices( *event, jet );

      tearDownVtxList[iJet] = *tmp;
      printf(" jet %d has %d vertices (teardown)\n",iJet+1,(int)tearDownVtxList[iJet].size() );

      for (unsigned int iVtx=0; iVtx<tearDownVtxList[iJet].size(); ++iVtx) {
        const vector<Track*>& trks = tearDownVtxList[iJet][iVtx]->getTracks();
        tearDownVtxTracks.insert( tearDownVtxTracks.end(), trks.begin(), trks.end() );
      }

    }
*/

    static TEveCompound* eveMCParticles(0);
    static TEveCompound* eveTracks(0);
    static TEveCompound* eveVertices(0);
    static TEveCompound* eveTDVertices(0);
    static TEveCompound* eveBUVertices(0);
    if (eveMCParticles != 0) {
      eveMCParticles->DestroyElements();
      eveMCParticles->Destroy();
    }
    if (eveTracks != 0) {
      eveTracks->DestroyElements();
      eveTracks->Destroy();
    }
    if (eveVertices != 0) {
      eveVertices->DestroyElements();
      eveVertices->Destroy();
    }
    if (eveTDVertices != 0) {
      eveTDVertices->DestroyElements();
      eveTDVertices->Destroy();
    }
    if (eveBUVertices != 0) {
      eveBUVertices->DestroyElements();
      eveBUVertices->Destroy();
    }

    eveMCParticles = new TEveCompound();
    eveMCParticles->SetMainColor(kBlue);
    eveMCParticles->SetName("MCParticles");
    eveMCParticles->OpenCompound();

    eveTracks = new TEveCompound();
    eveTracks->SetMainColor(kRed);
    eveTracks->SetName("Tracks");
    eveTracks->OpenCompound();
/*
    eveVertices = new TEveCompound();
    eveVertices->SetMainColor(kRed);
    eveVertices->SetName("Vertices");
    eveVertices->OpenCompound();

    eveTDVertices = new TEveCompound();
    eveTDVertices->SetMainColor(kGreen);
    eveTDVertices->SetName("TearDownVertices");
    eveTDVertices->OpenCompound();
*/
    eveBUVertices = new TEveCompound();
    eveBUVertices->SetMainColor(kCyan);
    eveBUVertices->SetName("BuildUpVertices");
    eveBUVertices->OpenCompound();

    TEveTrackPropagator* propsetMCParticle = new TEveTrackPropagator();
    TEveTrackPropagator* propsetTrack = new TEveTrackPropagator();

    TEveCompound* compoundMCParticles   = new TEveCompound("MCParticles", "All charged MCParticles");
    compoundMCParticles->SetMainColor(kBlue);

    TEveCompound* compoundTracks   = new TEveCompound("Tracks", "All charged tracks");
    compoundTracks->SetMainColor(kRed);
/*
    TEveCompound* compoundVertices   = new TEveCompound("Vertices", "All charged tracks");
    compoundVertices->SetMainColor(kGreen);

    TEveCompound* compoundTDVertices   = new TEveCompound("TearDownVertices", "All charged tracks");
    compoundTDVertices->SetMainColor(kGreen);
*/
    TEveCompound* compoundBUVertices   = new TEveCompound("BuildUpVertices", "All charged tracks");
    compoundBUVertices->SetMainColor(kGreen);

    int nMCParticle =  mcps.size();
    int nTrack = tracks.size();
/*
    TEvePointSet* vbox0 = new TEvePointSet();
    vbox0->SetMarkerColor(33);
    vbox0->SetMarkerStyle(3);
    vbox0->SetMarkerSize(3);
		vbox0->SetNextPoint(0,0,0);
    compoundVertices->AddElement(vbox0);

    for (unsigned int iJet=0; iJet<jets.size(); ++iJet) {
      for (unsigned int iVtx=0; iVtx<vtxList[iJet].size(); ++iVtx) {
        Vertex* vtx = vtxList[iJet][iVtx];
        float x0 = vtx->getX();
        float y0 = vtx->getY();
        float z0 = vtx->getZ();
        //if (x0*x0+y0*y0+z0*z0 > MAX_R*MAX_R+MAX_Z*MAX_Z) continue;
				TEvePointSet* vbox = new TEvePointSet();
				vbox->SetMarkerColor(kMagenta);
				vbox->SetMarkerStyle(4);
				vbox->SetMarkerSize(2);
        vbox->SetNextPoint(x0,y0,z0);
				vbox->SetTitle(Form("ZVTOP Vertex No.=%d\n"
							"Ntrk=%d, Chi2=%f\n"
							"(Vx, Vy, Vz) = (%.3f, %.3f, %.3f)",
							iVtx,
							(int)vtx->getTracks().size(), vtx->getChi2(),
							x0, y0, z0));
				compoundVertices->AddElement(vbox);
      }
    }

    for (unsigned int iJet=0; iJet<jets.size(); ++iJet) {
      for (unsigned int iVtx=0; iVtx<tearDownVtxList[iJet].size(); ++iVtx) {
        Vertex* vtx = tearDownVtxList[iJet][iVtx];
        float x0 = vtx->getX();
        float y0 = vtx->getY();
        float z0 = vtx->getZ();
        //if (x0*x0+y0*y0+z0*z0 > MAX_R*MAX_R+MAX_Z*MAX_Z) continue;
				TEvePointSet* vbox2 = new TEvePointSet();
				vbox2->SetMarkerColor(kGreen);
				vbox2->SetMarkerStyle(4);
				vbox2->SetMarkerSize(2);
        vbox2->SetNextPoint(x0,y0,z0);
				vbox2->SetTitle(Form("Teardown Vertex No.=%d\n"
							"Ntrk=%d, Chi2=%f\n"
							"(Vx, Vy, Vz) = (%.3f, %.3f, %.3f)",
							iVtx,
							(int)vtx->getTracks().size(), vtx->getChi2(),
							x0, y0, z0));
				compoundTDVertices->AddElement(vbox2);
      }
    }
*/
    for (unsigned int iVtx=0; iVtx<buildUpVtxList.size(); ++iVtx) {
      Vertex* vtx = buildUpVtxList[iVtx];
      float x0 = vtx->getX();
      float y0 = vtx->getY();
      float z0 = vtx->getZ();
      //if (x0*x0+y0*y0+z0*z0 > MAX_R*MAX_R+MAX_Z*MAX_Z) continue;
			TEvePointSet* vbox2 = new TEvePointSet();
			vbox2->SetMarkerColor(vtx->getTracks().size() == 2 ? kWhite : kCyan);
			vbox2->SetMarkerStyle(4);
			vbox2->SetMarkerSize(2);
      vbox2->SetNextPoint(x0,y0,z0);
			vbox2->SetTitle(Form("Buildup Vertex No.=%d\n"
						"Ntrk=%d, Chi2=%f\n"
						"(Vx, Vy, Vz) = (%.3f, %.3f, %.3f)",
						iVtx,
						(int)vtx->getTracks().size(), vtx->getChi2(),
						x0, y0, z0));
			compoundBUVertices->AddElement(vbox2);
    }

/*
    compoundVertices->CloseCompound();
    eveVertices->AddElement(compoundVertices);
    eveVertices->CloseCompound();

    compoundTDVertices->CloseCompound();
    eveTDVertices->AddElement(compoundTDVertices);
    eveTDVertices->CloseCompound();
*/    compoundBUVertices->CloseCompound();
    eveBUVertices->AddElement(compoundBUVertices);
    eveBUVertices->CloseCompound();

    for(int i=0; i<nMCParticle; i++) {
      MCParticle* part = mcps[i];
      // draw charged MCParticles only
      if (part->getCharge() == 0) continue;
			// draw stable particles only
      if ( !part->isStableTrack() ) continue;
      // skip primary particles
      //if (part->getFlavorTagCategory() == 1) continue;

      float px=part->Px();
      float py=part->Py();
      float pz=part->Pz();
      int pdg=part->getPDG();
      float Vx=part->getVertex().x();
      float Vy=part->getVertex().y();
      float Vz=part->getVertex().z();
      float Ex=part->getEx();
      float Ey=part->getEy();
      float Ez=part->getEz();
      int charge=(int)part->getCharge();
      //mass=part->getMass();
      float energy=part->E();
      float pt = part->Pt();
      //ParentNum=part->getParents().size();
      //DaughterNum=part->getDaughters().size();
      TEveVector Vtx(Vx, Vy, Vz);
      TEveVector End(Ex, Ey, Ez);
			//if ( !part->isStableTrack() && (Vtx-End).Mag() < MAX_R*MAX_R+MAX_Z*MAX_Z ) { continue; }

      TEveRecTrack* chargedTrack = new TEveRecTrack();
      chargedTrack->fV.Set(Vtx);
      chargedTrack->fP.Set(px, py, pz);
      chargedTrack->fSign = (int) charge;

      propsetMCParticle->SetMagFieldObj(new TEveMagFieldConst(0., 0., -3.5));
      //propsetMCParticle->SetMagFieldObj(new TEveMagFieldDuo(350, -3.5, 2.0));

      //Test for Charged Tracks
      propsetMCParticle->SetName("Track propagator for charged particles");
      propsetMCParticle->SetMaxR(MAX_R);
      propsetMCParticle->SetMaxZ(MAX_Z);
      propsetMCParticle->SetMaxOrbs(3.0);
      //propsetMCParticle->SetDelta(0.01); // Step
      propsetMCParticle->SetDelta(0.001); // Step
      // propsetMCParticle->SetRnrCluster2Ds(kTRUE); // Show Clusters 2D
      propsetMCParticle->RefPMAtt().SetMarkerColor(kYellow);
      propsetMCParticle->RefPMAtt().SetMarkerStyle(kCircle);
      propsetMCParticle->RefPMAtt().SetMarkerSize(1.0);
      //End of Test

      TEveTrack* track = new TEveTrack(chargedTrack, propsetMCParticle);
      /*
      if ( (Vtx-End).Mag() > 0 ) {
        TEvePathMark* pm3 = new TEvePathMark(TEvePathMark::kDecay);
        pm3->fV.Set(End);
        track->AddPathMark(*pm3);
      }
      */

      track->SetName(Form("Track %d", i));    // i = tracknum
      track->SetLineColor(kRed);
      track->SetLineWidth(1);
      //track->SetLineStyle(7);
      track->SetLineStyle(1);
			/*
			printf(" pdg=%d / ssp=%d / cat=%d / dist=%.3e\n",
					part->getPDG(),
					part->getSemiStableParent() ? part->getSemiStableParent()->getPDG() : 666,
					part->getFlavorTagCategory(),
					part->getVertex().Mag()
					);
			 */
      switch (part->getFlavorTagCategory()) {
        case 1: //p
          track->SetLineColor(kBlue);
          break;
        case 2: //b
          track->SetLineColor(kCyan);
          break;
        case 3: //c
          track->SetLineColor(kYellow);
          break;
        case 4: //s
          track->SetLineColor(kGreen);
          break;
        default:
          track->SetLineColor(kWhite);
          break;
      }

      /*
      if (Vtx.Mag() > 0) {
        track->SetLineColor(kGreen);
      }
      */

      track->SetSmooth(kTRUE);
      track->SetTitle(Form("MCParticle No.=%d\n""Charge=%d, pdg=%d\n"
            "Energy=%f, Pt=%f\n"
            "(Vx, Vy, Vz) = (%.3f, %.3f, %.3f)\n"
            "(Ex, Ey, Ez) = (%.3f, %.3f, %.3f)\n"
            "(Px, Py, Pz) = (%.3f, %.3f, %.3f)",
            i, charge, pdg, energy, pt,
            Vx, Vy, Vz, Ex, Ey, Ez, px, py, pz));
      compoundMCParticles->AddElement(track);
    }

    compoundMCParticles->CloseCompound();
    eveMCParticles->AddElement(compoundMCParticles);
    eveMCParticles->CloseCompound();

    for(int i=0; i<nTrack; i++) {
      Track* part = tracks[i];
      if (part->getCharge() == 0) continue;

      float px=part->Px();
      float py=part->Py();
      float pz=part->Pz();
      int pdg=part->getPDG();
			// x0,y0,z0
			double d0 = part->getD0();
			double z0 = part->getZ0();
			double phi0 = part->getPhi();
			float Vx = -d0 * sin(phi0);
			float Vy = d0 * cos(phi0);
			float Vz = z0;
/*      float Vx=part->getX();
      float Vy=part->getY();
      float Vz=part->getZ();*/
      //float Ex=part->getEx();
      //float Ey=part->getEy();
      //float Ez=part->getEz();
      int charge=(int)part->getCharge();
      //mass=part->getMass();
      float energy=part->E();
      float pt = part->Pt();
      //ParentNum=part->getParents().size();
      //DaughterNum=part->getDaughters().size();
      TEveVector Vtx(Vx, Vy, Vz);
      //TEveVector End(Ex, Ey, Ez);

      TEveRecTrack* chargedTrack = new TEveRecTrack();
      chargedTrack->fV.Set(Vtx);
      chargedTrack->fP.Set(px, py, pz);
      chargedTrack->fSign = (int) charge;

      propsetTrack->SetMagFieldObj(new TEveMagFieldConst(0., 0., -3.5));
      //propsetTrack->SetMagFieldObj(new TEveMagFieldDuo(350, -3.5, 2.0));

      //Test for Charged Tracks
      propsetTrack->SetName("Track propagator for charged particles");
      propsetTrack->SetMaxR(MAX_R);
      propsetTrack->SetMaxZ(MAX_Z);
      propsetTrack->SetMaxOrbs(2.0);
      //propsetTrack->SetDelta(0.01); // Step
      propsetTrack->SetDelta(0.001); // Step
      // propsetTrack->SetRnrCluster2Ds(kTRUE); // Show Clusters 2D
      propsetTrack->RefPMAtt().SetMarkerColor(kYellow);
      propsetTrack->RefPMAtt().SetMarkerStyle(kCircle);
      propsetTrack->RefPMAtt().SetMarkerSize(1.0);
      //End of Test

      TEveTrack* track = new TEveTrack(chargedTrack, propsetTrack);
      /*
      if ( (Vtx-End).Mag() > 0 ) {
        TEvePathMark* pm3 = new TEvePathMark(TEvePathMark::kDecay);
        pm3->fV.Set(End);
        track->AddPathMark(*pm3);
      }
      */

      track->SetName(Form("Track %d", i));    // i = tracknum
      if ( find( vtxTracksBU.begin(), vtxTracksBU.end(), part )
          != vtxTracksBU.end() ) {
        track->SetLineColor(kGreen);
        track->SetLineWidth(2);
      } else if ( find( vtxTracks.begin(), vtxTracks.end(), part )
          != vtxTracks.end() ) {
        track->SetLineColor(kRed);
        track->SetLineWidth(2);
      } else {
        track->SetLineColor(kBlue);
        track->SetLineWidth(2);
      }
      track->SetLineStyle(7);

      track->SetSmooth(kTRUE);
      track->SetTitle(Form("Track No.=%d\n""Charge=%d, pdg=%d\n"
            "Energy=%f, Pt=%f\n"
            "(Vx, Vy, Vz) = (%.3f, %.3f, %.3f)\n"
            //"(Ex, Ey, Ez) = (%.3f, %.3f, %.3f)\n"
            "(Px, Py, Pz) = (%.3f, %.3f, %.3f)",
            i, charge, pdg, energy, pt,
            Vx, Vy, Vz, /*Ex, Ey, Ez,*/ px, py, pz));
      compoundTracks->AddElement(track);
    }

    compoundTracks->CloseCompound();
    eveTracks->AddElement(compoundTracks);
    eveTracks->CloseCompound();

    std::cout<<"End of Building MC Tracks"<<endl;

    gEve->AddElement(eveBUVertices);
//    gEve->AddElement(eveVertices);
//    gEve->AddElement(eveTDVertices);
    gEve->AddElement(eveMCParticles);
    gEve->AddElement(eveTracks);
    gEve->FullRedraw3D(kTRUE);
  }
#else
	//void EventNavigator::Fwd(){}
	//void EventNavigator::Bck(){}
	//void EventNavigator::drawEvent(Event* event){}

#endif // BUILD_EVE
}
