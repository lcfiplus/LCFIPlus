#include "MakeNtuple.h"

#include <assert.h>
#include "EventStore.h"
#include "LcfiInterface.h"
#include "JetFinder.h"
#include "TreeStorer.h"
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

using namespace flavtag;
using namespace flavtag::algoSigProb;
using namespace flavtag::algoEtc;

const double PION_MASS = 0.13957018;
const double PION_MASS2 = PION_MASS*PION_MASS;

void MakeNtuple::init(FlavtagParameters *param) {
	FlavtagAlgorithm::init(param);

	string outputFilename = param->get("TrainNtupleFile",string("flavtag.root"));

	FlavtagData& d = _data;

	_file = new TFile(outputFilename.c_str(),"RECREATE");
	_t = new TTree("ntp","events");

	_t->Branch("flavor", &d.flavor, "flavor/D");
	_t->Branch("nvtx", &d.nvtx, "nvtx/D");
	_t->Branch("vtxlen1", &d.vtxlen1, "vtxlen1/D");
	_t->Branch("vtxsig1", &d.vtxsig1, "vtxsig1/D");
	_t->Branch("vtxlen2", &d.vtxlen2, "vtxlen2/D");
	_t->Branch("vtxsig2", &d.vtxsig2, "vtxsig2/D");
	_t->Branch("vtxlen12", &d.vtxlen12, "vtxlen12/D");
	_t->Branch("vtxsig12", &d.vtxsig12, "vtxsig12/D");
	_t->Branch("vtxmom", &d.vtxmom, "vtxmom/D");
	_t->Branch("vtxmass", &d.vtxmass, "vtxmass/D");
	_t->Branch("vtxmasspc", &d.vtxmasspc, "vtxmasspc/D");
	_t->Branch("vtxmult", &d.vtxmult, "vtxmult/D");
	_t->Branch("vtxmult1", &d.vtxmult1, "vtxmult1/D");
	_t->Branch("vtxmult2", &d.vtxmult2, "vtxmult2/D");
	_t->Branch("vtxmom1", &d.vtxmom1, "vtxmom1/D");
	_t->Branch("vtxmom2", &d.vtxmom2, "vtxmom2/D");
	_t->Branch("vtxmass1", &d.vtxmass1, "vtxmass1/D");
	_t->Branch("vtxmass2", &d.vtxmass2, "vtxmass2/D");
	_t->Branch("vtxprob", &d.vtxprob, "vtxprob/D");
	_t->Branch("px", &d.px, "px/D");
	_t->Branch("py", &d.py, "py/D");
	_t->Branch("pz", &d.pz, "pz/D");
	_t->Branch("e", &d.e, "e/D");
	_t->Branch("trk1d0sig", &d.trk1d0sig, "trk1d0sig/D");
	_t->Branch("trk2d0sig", &d.trk2d0sig, "trk2d0sig/D");
	_t->Branch("trk1z0sig", &d.trk1z0sig, "trk1z0sig/D");
	_t->Branch("trk2z0sig", &d.trk2z0sig, "trk2z0sig/D");
	_t->Branch("trk1pt", &d.trk1pt, "trk1pt/D");
	_t->Branch("trk2pt", &d.trk2pt, "trk2pt/D");
	_t->Branch("jprobr", &d.jprobr, "jprobr/D");
	_t->Branch("jprobz", &d.jprobz, "jprobz/D");
	_t->Branch("sphericity", &d.sphericity, "sphericity/D");
}

void MakeNtuple::process() {
	///////////////////////////////////////////
	// initialization start
	///////////////////////////////////////////
	Event evt;
	if (evt.getTracks().size() == 0) return;
	FlavtagData& d = _data;
	memset(&d, 0, sizeof(d));
	JetConfig jetCfg;
	jetCfg.nJet = 2;
	SecondaryVertexConfig secVtxCfg;
	// AND
	secVtxCfg.maxD0 = 10;
	secVtxCfg.maxD0Err = 0.25;
	secVtxCfg.maxZ0 = 20;
	secVtxCfg.maxZ0Err = 1e10;
	secVtxCfg.minPt = 0.1;
	secVtxCfg.maxInnermostHitRadius = 1e10;
	// OR
	secVtxCfg.minTpcHits = 20;
	secVtxCfg.minFtdHits = 3;
	secVtxCfg.minVtxHitsWithoutTpcFtd = 3;
	secVtxCfg.minVtxPlusFtdHits = 0;
	///////////////////////////////////////////
	// initialization end
	///////////////////////////////////////////

	// fill jet type based on the number of semistable b & c baryons
	d.mcnb = evt.mcNumberOfB();
	d.mcnc = evt.mcNumberOfC();
	if (d.mcnb == 2)
		d.flavor = 5;
	else if (d.mcnb == 0 && d.mcnc == 2)
		d.flavor = 4;
	else if (d.mcnb == 0 && d.mcnc == 0)
		d.flavor = 3;
	else {
		d.flavor = 0;
	}

	///////////////////////////////////////////
	// primary vertex finding
	///////////////////////////////////////////
	const vector<Track*>& tracks = evt.getTracks();
	vector<Track*> tracksForPrimary;
	for (vector<Track*>::const_iterator it = tracks.begin(); it != tracks.end(); ++it) {
		Track* trk = *it;
		if (fabs(trk->getD0())<50 && fabs(trk->getZ0())<50) {
			if (trk->getVtxHits() + trk->getFtdHits() >= 3 && trk->getRadiusOfInnermostHit() <20) {
				tracksForPrimary.push_back( trk );
			}
		}
	}
	Vertex* primaryVertex = findPrimaryVertex(tracksForPrimary,25);
	if (primaryVertex == 0) {
		fprintf(stderr,"primary vertex could not be found (ntrk=%d), skipping event\n",(int)tracksForPrimary.size());
		return;
	}

	///////////////////////////////////////////
	// LCFIVertex interface
	///////////////////////////////////////////
	LcfiInterface interface(&evt,primaryVertex);

	///////////////////////////////////////////
	// perform jet clustering
	///////////////////////////////////////////
	vector<Jet*> sueharaJets = SueharaJetClustering(evt,jetCfg.nJet);
	// create "basic" jet without vertex info
	vector<Jet*> jets;
	for (vector<Jet*>::iterator it = sueharaJets.begin(); it != sueharaJets.end(); ++it) {
		jets.push_back( convertJetVertex(*it) );
	}

	///////////////////////////////////////////
	// fill jet info
	///////////////////////////////////////////
	for (unsigned int iJet=0; iJet<jets.size(); ++iJet) {

		d.px = jets[iJet]->Px();
		d.py = jets[iJet]->Py();
		d.pz = jets[iJet]->Pz();
		d.e  = jets[iJet]->E();

		// find trk with best (and 2nd best) significance in the r-phi plane
		float sigVec[6];
		findMostSignificantTrack(jets[iJet],primaryVertex,sigVec);
		d.trk1d0sig = sigVec[0];
		d.trk2d0sig = sigVec[1];
		d.trk1z0sig = sigVec[2];
		d.trk2z0sig = sigVec[3];
		d.trk1pt = sigVec[4];
		d.trk2pt = sigVec[5];

		// fill joint probabilities
		d.jprobr = jointProbD0(jets[iJet],primaryVertex);
		d.jprobz = jointProbZ0(jets[iJet],primaryVertex);

		// fill boosted sphericity
		d.sphericity= jets[iJet]->sphericity();

		////////////////////////////////////////////////////////////
		// look at vertex information
		////////////////////////////////////////////////////////////
		int nvtxAccepted(0);
		int vtxmult(0);
		TLorentzVector vtxmomSum;
		double oneMinusProb(1.);

		int ivtx1(-1);
		int ivtx2(-1);
		double bestLen1(-1);
		double bestLen2(-1);

		// the vertex loop
		const vector<Vertex*>& vtxList = sueharaJets[iJet]->getVertices();
		for (unsigned int j=0; j<vtxList.size(); ++j) {

			// apply v0 veto
			if (vtxList[j]->passesV0selection(primaryVertex)) continue;

			// keep track of vertex with smallest (2nd smallest) decay distance
			double len = vtxList[j]->length(primaryVertex);
			if (len < bestLen1) {
				bestLen2 = bestLen1;
				bestLen1 = len;
				ivtx2 = ivtx1;
				ivtx1 = j;
			} else if (len < bestLen2) {
				bestLen2 = len;
				ivtx2 = j;
			}

			// total vertex multiplicity
			const vector<Track*>& vtxTracks = vtxList[j]->getTracks();
			vtxmult += vtxTracks.size();

			for (unsigned int k=0; k<vtxTracks.size(); ++k) {
				Track* trk = vtxTracks[k];
				TLorentzVector vtxmom = *trk;
				vtxmomSum += vtxmom;
			}

			// for multiple vertices, combine probability by
			// 1-pALL = (1-p1)*(1-p2)*(1-p3)*...
			int ndf = vtxTracks.size()*2-3;
			double prob = TMath::Prob(vtxList[j]->getChi2(),ndf);
			oneMinusProb *= (1-prob);

			// count the number of accepted vertices
			++nvtxAccepted;
		}

		d.nvtx= nvtxAccepted;

		if (ivtx1 != -1) {
			d.vtxlen1 = vtxList[ivtx1]->length(primaryVertex);
			d.vtxsig1 = vtxList[ivtx1]->significance(primaryVertex);
			d.vtxdirdot1 = vtxList[ivtx1]->dirdot(primaryVertex);
			d.vtxmult1 = vtxList[ivtx1]->getTracks().size();

			TLorentzVector vtx1p4 = vtxList[ivtx1]->getFourMomentum();
			d.vtxmom1 = vtx1p4.Vect().Mag();
			d.vtxmass1 = vtx1p4.M();

			double pt = interface.vertexMassPtCorrection(vtxList[ivtx1],primaryVertex,vtxmomSum.Vect(),2);
			double vm = vtxmomSum.M();
			d.vtxmasspc = sqrt( vm*vm+pt*pt ) + pt;
		}

		if (ivtx2 != -1) {
			d.vtxlen2 = vtxList[ivtx2]->length(primaryVertex);
			d.vtxsig2 = vtxList[ivtx2]->significance(primaryVertex);
			d.vtxdirdot2 = vtxList[ivtx2]->dirdot(primaryVertex);
			d.vtxmult2 = vtxList[ivtx2]->getTracks().size();

			TLorentzVector vtx2p4 = vtxList[ivtx2]->getFourMomentum();
			d.vtxmom2 = vtx2p4.Vect().Mag();
			d.vtxmass2 = vtx2p4.M();

			d.vtxlen12 = vtxList[ivtx2]->length(vtxList[ivtx1]);
			d.vtxsig12 = vtxList[ivtx2]->significance(vtxList[ivtx1]);
			d.vtxdirdot12 = vtxList[ivtx2]->dirdot(vtxList[ivtx1]);
		}

		d.vtxmult = vtxmult;
		d.vtxmom  = vtxmomSum.Vect().Mag();
		d.vtxmass = vtxmomSum.M();
		d.vtxprob = 1-oneMinusProb;

		_t->Fill();
	}

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

void MakeNtuple::end() {
	_file->Write();
	_file->Close();
}
