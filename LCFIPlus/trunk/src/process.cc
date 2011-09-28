#include <assert.h>
#include <string>

#include "EventStore.h"
#include "TrackSelector.h"
#include "JetFinder.h"
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

#include "process.h"
#include "geometry.h"

using namespace lcfiplus;
using namespace lcfiplus::algoEtc;

namespace lcfiplus{

	void PrimaryVertexFinder::init(LcfiplusParameters *param){
		LcfiplusAlgorithm::init(param);

		_vertex = 0;

		string vcolname = param->get("PrimaryVertexCollectionName",string("PrimaryVertex"));
		Event::Instance()->Register(vcolname.c_str(), _vertex, EventStore::PERSIST);

		// default setting
		Event::Instance()->setDefaultPrimaryVertex(vcolname.c_str());

		// parameters...(should be exported)
		_secVtxCfg = new TrackSelectorConfig;

		/*
		_secVtxCfg->maxD0 = 50.;
		_secVtxCfg->maxZ0 = 50.;
		_secVtxCfg->maxInnermostHitRadius = 20.;
		_secVtxCfg->minVtxPlusFtdHits = 3;
		*/

		_secVtxCfg->maxD0 = 20.;
		_secVtxCfg->maxZ0 = 20.;
		_secVtxCfg->maxInnermostHitRadius = 20.;
		_secVtxCfg->minVtxPlusFtdHits = 5;

		_chi2th = 25.;
	}

	void PrimaryVertexFinder::process(){
		Event *event = Event::Instance();

		// clearing old vertices
		if (_vertex->size()>0){delete (*_vertex)[0]; _vertex->clear();}

		// cut bad tracks
		const vector<Track *> &tracks = event->getTracks();
		vector<Track *> passedTracks = TrackSelector() (tracks, *_secVtxCfg);

		// primary vertex finder
		Vertex * vtx = findPrimaryVertex(passedTracks,_chi2th);
		_vertex->push_back(vtx);
	}

	void PrimaryVertexFinder::end() {
		delete _secVtxCfg;
	}


	void BuildUpVertex::init(LcfiplusParameters *param){
		LcfiplusAlgorithm::init(param);

		// register collection
//		string vcolname = (*param)["VertexCollectionName"];
		string vcolname = param->get("VertexCollectionName",string("BuildUpVertex"));
		Event::Instance()->Register(vcolname.c_str(), _vertices, EventStore::PERSIST);

		// configuration
		// track cut
		_secVtxCfg = new TrackSelectorConfig; // todo: memory leak

		_secVtxCfg->maxD0 = param->get("TrackCut.MaxD0", 10.);
		_secVtxCfg->maxZ0 = param->get("TrackCut.MaxZ0", 20.);
		_secVtxCfg->minPt = param->get("TrackCut.MinPt", 0.1);
		_secVtxCfg->maxInnermostHitRadius = 1e10;

		_secVtxCfg->maxD0Err = param->get("TrackCut.MaxD0Err", .1);
		_secVtxCfg->maxZ0Err = param->get("TrackCut.MaxZ0Err", .1);

		_secVtxCfg->minTpcHits = param->get("TrackCut.MinTpcHits", 20);
		_secVtxCfg->minFtdHits = param->get("TrackCut.MinFtdHits", 3);
		_secVtxCfg->minVtxHitsWithoutTpcFtd = param->get("TrackCut.MinVtxHits", 3);
		_secVtxCfg->minVtxPlusFtdHits = param->get("TrackCut.MinVtxFtdHits", 0);

		// buildup parameters
		_chi2thpri = param->get("BuildUp.PrimaryChi2Threshold", 25.);
		_chi2thsec = param->get("BuildUp.SecondaryChi2Threshold", 9.);
		_massth = param->get("BuildUp.MassThreshold", 10.);
		_posth = param->get("BuildUp.MinimumDistIP", 0.3);
		_chi2orderinglimit = param->get("BuildUp.MaximumChi2ForDistOrder", 1.0);

		_doassoc = param->get("AssocIPTracks.DoAssoc", 1);
		_minimumdist = param->get("AssocIPTracks.MinimumDist", 0.);
		_chi2ratio = param->get("AssocIPTracks.Chi2RatioSecToPri", 2.0);
	}

	void BuildUpVertex::process(){
		Event *event = Event::Instance();

		// clearing old vertices
		for(unsigned int n=0;n<_vertices->size();n++)
			delete (*_vertices)[n];
		_vertices->clear();

		// cut bad tracks
		const vector<Track *> &tracks = event->getTracks();
		vector<Track *> passedTracks = TrackSelector() (tracks, *_secVtxCfg);

		// build up vertexing
		VertexFinderSuehara::buildUp(passedTracks, *_vertices, _chi2thpri, _chi2thsec, _massth, _posth, _chi2orderinglimit);
		if(_doassoc)
			VertexFinderSuehara::associateIPTracks(*_vertices,_minimumdist, 0, _chi2ratio);
	}

	void JetClustering::init(LcfiplusParameters *param){
		LcfiplusAlgorithm::init(param);

		string vcolname = param->get("VertexCollectionName",string("BuildUpVertex"));
		Event::Instance()->Get(vcolname.c_str(), _vertices);
		string jcolname = param->get("JetCollectionName",string("Jets"));
		Event::Instance()->Register(jcolname.c_str(), _jets, EventStore::PERSIST | EventStore::JET_EXTRACT_VERTEX);

		_njets = param->get("NJetsRequested",int(6));
		_ycut = param->get("YCut", double(0));

		_useMuonID = param->get("UseMuonID", int(1));

		_vsMinDist = param->get("VertexSelectionMinimumDistance", double(0.3));
		_vsMaxDist = param->get("VertexSelectionMaximumDistance", double(30.));
		_vsK0MassWidth = param->get("VertexSelectionK0MassWidth", double(0.02));
	}

	void JetClustering::process()
	{
		// clearing old jets
		for(unsigned int n=0;n<_jets->size();n++)
			delete (*_jets)[n];
		_jets->clear();

		Event *event = Event::Instance();

		JetConfig jetCfg;
		jetCfg.nJet = _njets;
		JetFinder* jetFinder = new JetFinder(jetCfg);
		double ycut = _ycut;

		// select vertices
		vector<Vertex *> selectedVertices;
		vector<Track *> residualTracks = event->getTracks();
		for(vector<Vertex *>::const_iterator it = _vertices->begin(); it != _vertices->end();it++){
			Vertex *v = *it;
			double mass = 0.;
			double k0mass = .498;
			if(v->getTracks().size() == 2)
				mass = (*(TLorentzVector *)(v->getTracks()[0]) + *(TLorentzVector *)(v->getTracks()[1])).M();
			if((mass < k0mass - _vsK0MassWidth/2 || mass > k0mass + _vsK0MassWidth/2) && v->getPos().Mag() < _vsMaxDist && v->getPos().Mag() > _vsMinDist){
				selectedVertices.push_back(v);
				for(vector<Track *>::const_iterator itt = v->getTracks().begin(); itt != v->getTracks().end(); itt++){
					vector<Track *>::iterator itt2 = remove_if(residualTracks.begin(), residualTracks.end(), bind2nd(equal_to<Track *>(), *itt));
					residualTracks.erase(itt2, residualTracks.end());
				}
			}
		}
		*_jets = jetFinder->run(residualTracks, event->getNeutrals(), selectedVertices, &ycut, _useMuonID);

		/*
		cout << "JetClustering: _vertices.size() = " << _vertices->size() << endl;
		cout << "JetClustering: selectedVertices.size() = " << selectedVertices.size() << endl;

		for (vector<Jet*>::const_iterator it = _jets->begin(); it != _jets->end(); ++it) {
			Jet* jet = *it;
			cout << "   Jet nvtx = " << jet->getVertices().size() << endl;
		}
		*/
	}

	void JetClustering::end() {
		//cout << "ENDO" << endl;
	}


	/*
	void JetVertexRefiner::init() {
		LcfiplusAlgorithm::init(param);

		string vcolname = param->get("JetVertexRefinerInputJetCollectionName",string("VertexJets"));
		Event::Instance()->Get(vcolname.c_str(), _inputJets);
		string jcolname = param->get("JetVertexRefinerOutputJetCollectionName",string("RefinedJets"));
		Event::Instance()->Register(jcolname.c_str(), _outputJets, EventStore::PERSIST | EventStore::JET_EXTRACT_VERTEX);
	}

	void JetVertexRefiner::process() {
		// clearing old jets
		for(unsigned int n=0;n<_outputJets->size();n++)
			delete (*_outputJets)[n];
		_outputJets->clear();

		Event *event = Event::Instance();

		vector<Track*> residualTracks = event->getTracks();
		for(vector<Jet*>::const_iterator it = _inputJets->begin(); it != _inputJets->end();++it){
			Jet* jet = *it;

			const vector<Track*>& tracks = jet->getTracks();
			const vector<Vertex*>& vertices = jet->getVertices();
			for (vector<Vertex*>::const_iterator it2 = vertices.begin(); it2 != vertices.end(); ++it2) {
			}
		}
	}
	*/
}

