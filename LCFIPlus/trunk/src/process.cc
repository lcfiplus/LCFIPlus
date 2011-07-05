#include <assert.h>
#include <string>

#include "EventStore.h"
#include "LcfiInterface.h"
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

using namespace flavtag;
using namespace flavtag::algoEtc;

namespace flavtag{

	void BuildUpVertex::init(FlavtagParameters *param){
		FlavtagAlgorithm::init(param);

		// register collection
//		string vcolname = (*param)["VertexCollectionName"];
		string vcolname = param->get("VertexCollectionName",string("BuildUpVertex"));
		EventStore::Instance()->Register(vcolname.c_str(), _vertices, EventStore::PERSIST);

		// configuration
		// track cut
		_secVtxCfg = new SecondaryVertexConfig; // todo: memory leak

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
		Event evt;

		// clearing old vertices
		for(unsigned int n=0;n<_vertices->size();n++)
			delete (*_vertices)[n];
		_vertices->clear();

		LcfiInterface interface(&evt);

		// cut bad tracks
		const vector<Track *> &tracks = evt.getTracks();
		vector<Track *> passedTracks;

		for(unsigned int i=0;i<tracks.size();i++){
			if(interface.passesCut(tracks[i], *_secVtxCfg))
				passedTracks.push_back(tracks[i]);
		}

		// build up vertexing
		VertexFinderSuehara::buildUp(passedTracks, *_vertices, _chi2thpri, _chi2thsec, _massth, _posth, _chi2orderinglimit);
		if(_doassoc)
			VertexFinderSuehara::associateIPTracks(*_vertices,_minimumdist, 0, _chi2ratio);
	}

	void JetClustering::init(FlavtagParameters *param){
		FlavtagAlgorithm::init(param);

		string vcolname = param->get("VertexCollectionName",string("BuildUpVertex"));
		EventStore::Instance()->Get(vcolname.c_str(), _vertices);
		string jcolname = param->get("JetCollectionName",string("Jets"));
		EventStore::Instance()->Register(jcolname.c_str(), _jets, EventStore::PERSIST | EventStore::JET_EXTRACT_VERTEX);

		_njets = param->get("NJetsRequested",int(6));
	}

	void JetClustering::process()
	{
		// clearing old jets
		for(unsigned int n=0;n<_jets->size();n++)
			delete (*_jets)[n];
		_jets->clear();

		Event evt;

		JetConfig jetCfg;
		jetCfg.nJet = _njets;
		JetFinder* jetFinder = new JetFinder(jetCfg);
		double ycut = 0;

		// select vertices
		vector<Vertex *> selectedVertices;
		vector<Track *> residualTracks = evt.getTracks();
		for(vector<Vertex *>::const_iterator it = _vertices->begin(); it != _vertices->end();it++){
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
		*_jets = jetFinder->run(residualTracks, evt.getNeutrals(), selectedVertices, &ycut);
	}

	void JetClustering::end()
	{
		cout << "ENDO" << endl;
	}

}

