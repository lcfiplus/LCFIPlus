#include <assert.h>
#include <string>

#include "EventStore.h"
#include "TrackSelector.h"
#include "VertexSelector.h"
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

namespace lcfiplus {

void PrimaryVertexFinder::init(Parameters* param) {
  Algorithm::init(param);

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

  _secVtxCfg->maxD0 = param->get("PrimaryVertexFinder.TrackMaxD0", 20.);
  _secVtxCfg->maxZ0 = param->get("PrimaryVertexFinder.TrackMaxZ0", 20.);
  _secVtxCfg->minVtxPlusFtdHits = param->get("PrimaryVertexFinder.TrackMinVtxFtdHits", 1);
  _secVtxCfg->minTpcHits = param->get("PrimaryVertexFinder.TrackMinTpcHits", 999999);
  _secVtxCfg->minTpcHitsMinPt = param->get("PrimaryVertexFinder.TrackMinTpcHitsMinPt", 999999);
  _secVtxCfg->minFtdHits = param->get("PrimaryVertexFinder.TrackMinFtdHits", 999999);
  _secVtxCfg->minVtxHits = param->get("PrimaryVertexFinder.TrackMinVxdHits", 999999);

  _chi2th = param->get("PrimaryVertexFinder.Chi2Threshold", 25.);

  _beamspotConstraint = param->get("PrimaryVertexFinder.BeamspotConstraint", true);
  _beamspotSmearing = param->get("PrimaryVertexFinder.BeamspotSmearing", true);
}

void PrimaryVertexFinder::process() {
  bool verbose = false;
  Event* event = Event::Instance();

  // clearing old vertices
  if (_vertex->size()>0) {
    delete (*_vertex)[0];
    _vertex->clear();
  }

  // cut bad tracks
  TrackVec& tracks = event->getTracks();
  TrackVec passedTracks = TrackSelector() (tracks, *_secVtxCfg);
  if (verbose)
    cout << "PrimaryVertexFinder / track selection: " << passedTracks.size() << "/" << tracks.size() << " accepted." << endl;

  // primary vertex finder
  Vertex* vtx = findPrimaryVertex(passedTracks,_chi2th,_beamspotConstraint,_beamspotSmearing);
  _vertex->push_back(vtx);
  if (verbose)
    cout << "PrimaryVertexFinder: " << vtx->getTracks().size() << " tracks associated to the primary vertex." << endl;
}

void PrimaryVertexFinder::end() {
  delete _secVtxCfg;
}


void BuildUpVertex::init(Parameters* param) {
  Algorithm::init(param);

  // register collection
  string vcolname = param->get("BuildUpVertexCollectionName",string("BuildUpVertex"));
  Event::Instance()->Register(vcolname.c_str(), _vertices, EventStore::PERSIST);

  _primvtxcolname = param->get("PrimaryVertexCollectionName",string("PrimaryVertex"));
  Event::Instance()->setDefaultPrimaryVertex(_primvtxcolname.c_str());

  _v0vertices = 0;
  string v0vcolname = param->get("BuildUpVertex.V0VertexCollectionName",string(""));
  if (v0vcolname != "") {
    Event::Instance()->Register(v0vcolname.c_str(), _v0vertices, EventStore::PERSIST);
  }

  // configuration
  // track cut
  _secVtxCfg = new TrackSelectorConfig;

  _secVtxCfg->maxD0 = param->get("BuildUpVertex.TrackMaxD0", 10.);
  _secVtxCfg->maxZ0 = param->get("BuildUpVertex.TrackMaxZ0", 20.);
  _secVtxCfg->minPt = param->get("BuildUpVertex.TrackMinPt", 0.1);
  _secVtxCfg->maxInnermostHitRadius = 1e10;

  _secVtxCfg->maxD0Err = param->get("BuildUpVertex.TrackMaxD0Err", .1);
  _secVtxCfg->maxZ0Err = param->get("BuildUpVertex.TrackMaxZ0Err", .1);

  _secVtxCfg->minTpcHits = param->get("BuildUpVertex.TrackMinTpcHits", 999999);
  _secVtxCfg->minTpcHitsMinPt = param->get("BuildUpVertex.TrackMinTpcHitsMinPt", 999999);
  _secVtxCfg->minFtdHits = param->get("BuildUpVertex.TrackMinFtdHits", 999999);
  _secVtxCfg->minVtxHits = param->get("BuildUpVertex.TrackMinVxdHits", 999999);
  _secVtxCfg->minVtxPlusFtdHits = param->get("BuildUpVertex.TrackMinVxdFtdHits", 1);

  // buildup parameters
  _chi2thpri = param->get("BuildUpVertex.PrimaryChi2Threshold", 25.);
  _chi2thsec = param->get("BuildUpVertex.SecondaryChi2Threshold", 9.);
  _massth = param->get("BuildUpVertex.MassThreshold", 10.);
  _posth = param->get("BuildUpVertex.MinDistFromIP", 0.3);
  _chi2orderinglimit = param->get("BuildUpVertex.MaxChi2ForDistOrder", 1.0);

  _doassoc = param->get("BuildUpVertex.AssocIPTracks", 1);
  _minimumdist = param->get("BuildUpVertex.AssocIPTracksMinDist", 0.);
  _chi2ratio = param->get("BuildUpVertex.AssocIPTracksChi2RatioSecToPri", 2.0);
  _v0sel = param->get("BuildUpVertex.UseV0Selection", 1);

  //flag of AVF/chi2
  _avf = param->get("BuildUpVertex.UseAVF", true);
  //AVF parameter
  _temperature = param->get("BuildUpVertex.AVFTemperature", 5.0);
}

void BuildUpVertex::process() {
  Event* event = Event::Instance();

  // clearing old vertices
  for (unsigned int n=0; n<_vertices->size(); n++)
    delete (*_vertices)[n];
  _vertices->clear();

  // cut bad tracks
  TrackVec& tracks = event->getTracks();
  TrackVec passedTracks = TrackSelector() (tracks, *_secVtxCfg);

  //cout << "BuildUpVertex / track selection: " << passedTracks.size() << "/" << tracks.size() << " accepted." << endl;

  Vertex* primvtx;

  try {
    primvtx = new Vertex(*event->getPrimaryVertex(_primvtxcolname.c_str()));
  } catch (lcfiplus::Exception& e) {
    cout << "BuildUpVertex::process(): primary vertex not found - invoking primary vertex finder internally." << endl;
    primvtx = 0;
  }

  VertexFinderSuehara::VertexFinderSueharaConfig cfg;
  cfg.chi2th = _chi2thsec;
  cfg.massth = _massth;
  cfg.v0selVertex.minpos = _posth;
  cfg.chi2orderinglimit = _chi2orderinglimit;
  if (!_v0sel) {
    cfg.v0selTrack.setNoV0Cut();
    cfg.v0selVertex.setNoV0Cut();
  }
  cfg.minimumdistIP = _minimumdist;
  cfg.chi2ratioIP = _chi2ratio;

  //flag AVF/chi2
  cfg.avf=_avf;
  cfg.temperature=_temperature;

  // build up vertexing
  vector<Vertex*> v0tmp;
  VertexFinderSuehara::buildUp(passedTracks, *_vertices, (_v0vertices ? *_v0vertices : v0tmp), _chi2thpri, cfg, primvtx);
  // TODO: deletion of v0tmp

  if (_doassoc){
    //VertexFinderSuehara::associateIPTracks(*_vertices,primvtx, cfg);
    if(!_avf) VertexFinderSuehara::associateIPTracks(*_vertices,primvtx, cfg);
    else VertexFinderSuehara::associateIPTracksAVF(*_vertices,primvtx, cfg);
  }

  /*
  for(unsigned int n=0;n<_vertices->size();n++){
  	(*_vertices)[n]->Print();
  }
  */
}

void BuildUpVertex::end() {
  delete _secVtxCfg;
}

void JetClustering::init(Parameters* param) {
  Algorithm::init(param);

  _vcolname = param->get("JetClustering.InputVertexCollectionName",string("BuildUpVertex"));
  Event::Instance()->Get(_vcolname.c_str(), _vertices);
  if (!_vertices)
    cout << "JetClustering::init: vertex collection not found; will try later." << endl;

  vector<string> jcolnames;
  param->fetchArray("JetClustering.OutputJetCollectionName",jcolnames);
  param->fetchArray("JetClustering.NJetsRequested",_njets);
  param->fetchArray("JetClustering.YCut",_ycut);
  _algo = param->get("JetClustering.JetAlgorithm", string("DurhamVertex")); // DurhamVertex, Durham, KtVertex, Kt
  _useBeamJets = param->get("JetClustering.UseBeamJets", int(0));
  _rParameter = param->get("JetClustering.RParameter", double(1.0));
  _alphaParameter = param->get("JetClustering.AlphaParameter", double(1.0));
  _betaParameter = param->get("JetClustering.BetaParameter", double(1.0));
  _outputVertexStoresVertex = param->get("JetClustering.OutputJetStoresVertex",int(0));

  // checks
  if (jcolnames.size() == 0)
    throw(Exception(
            "JetClustering::init: output jet collection is not specified. please include OutputJetCollectionName parameter."));

  if (_njets.size() == 0 && _ycut.size() == 0)
    throw(Exception(
            "JetClustering::init: please specify at least either NJetsRequested or YCut parameter."));

  if (_njets.size() > 1 && _ycut.size() > 1)
    throw(Exception(
            "JetClustering::init: cannot accept multiple NJetsRequested with multiple YCut."));

  if (_njets.size() > jcolnames.size() || _ycut.size() > jcolnames.size())
    throw(Exception(
            "JetClustering::init: please specify enough number of OutputJetCollectionName to meet number of NJetsRequested or YCut."));

  // Need to add some exceptions for the new variables _useBeamJets and _rParameter

  if (_njets.size() == 0)_njets.push_back(0);
  if (_ycut.size() == 0)_ycut.push_back(0.);

  // sort njetsrequested / ycut
  if (_njets.size() > 1) {
    for (unsigned int n=0; n<_njets.size(); n++) {
      vector<Jet*>* jets;
      Event::Instance()->Register(jcolnames[n].c_str(), jets, EventStore::PERSIST | (_outputVertexStoresVertex ? EventStore::JET_WRITE_VERTEX : 0));
      _jetsmap[(double)_njets[n]] = jets;
    }
    // sort njets
    sort(_njets.begin(), _njets.end(), std::greater<int>());
  } else {
    for (unsigned int n=0; n<_ycut.size(); n++) {
      vector<Jet*>* jets;
      Event::Instance()->Register(jcolnames[n].c_str(), jets, EventStore::PERSIST | (_outputVertexStoresVertex ? EventStore::JET_WRITE_VERTEX : 0));
      _jetsmap[_ycut[n]] = jets;
    }
    // sort njets
    sort(_ycut.begin(), _ycut.end(), std::less<double>());
  }

  _yaddVV = param->get("JetClustering.YAddedForJetVertexVertex", double(100.));
  _yaddVL = param->get("JetClustering.YAddedForJetVertexLepton", double(100.));
  _yaddLL = param->get("JetClustering.YAddedForJetLeptonLepton", double(100.));

  _useMuonID = param->get("JetClustering.UseMuonID", int(1));
  _muonIDExternal = param->get("JetClustering.MuonIDExternal", int(0)); // default is internal, conservatively
  _muonIDMinEnergy = param->get("JetClustering.MuonIDMinimumEnergy", double(5.));
  _muonIDMinD0Sig = param->get("JetClustering.MuonIDMinimumD0Significance", double(5.));
  _muonIDMinZ0Sig = param->get("JetClustering.MuonIDMinimumZ0Significance", double(5.));
  _muonIDMaxDist = param->get("JetClustering.MuonIDMaximum3DImpactParameter", double(5.));
  _muonIDMinProb = param->get("JetClustering.MuonIDMinimumProbability", double(0.5));

  _vsMinDist = param->get("JetClustering.VertexSelectionMinimumDistance", double(0.3));
  _vsMaxDist = param->get("JetClustering.VertexSelectionMaximumDistance", double(30.));
  _vsK0MassWidth = param->get("JetClustering.VertexSelectionK0MassWidth", double(0.02));

  _maxYth = param->get("JetClustering.MaxNumberOfJetsForYThreshold", int(10));
}

void JetClustering::process() {
  // clearing old jets
  map<double, vector<Jet*> * >::iterator it;
  for (it = _jetsmap.begin(); it != _jetsmap.end(); it++) {
    vector<Jet*>* jets = it->second;
    for (unsigned int n=0; n<jets->size(); n++)
      delete (*jets)[n];
    jets->clear();
  }

  if (!_vertices) {
    // retry
    Event::Instance()->Get(_vcolname.c_str(), _vertices);
    if (!_vertices) {
      cout << "JetClustering::Process: Vertex not found, clustering without vertices..." << endl;
    } else
      cout << "JetClustering::Process: Vertex found." << endl;
  }

  Event* event = Event::Instance();

  JetConfig jetCfg;
  jetCfg.nJet = _njets[0];
  jetCfg.Ycut = _ycut[0];
  jetCfg.algo = _algo;
  jetCfg.useBeamJets = _useBeamJets;
  jetCfg.rParameter = _rParameter;
  jetCfg.alphaParameter = _alphaParameter;
  jetCfg.betaParameter = _betaParameter;
  jetCfg.YaddVV = _yaddVV;
  jetCfg.YaddVL = _yaddVL;
  jetCfg.YaddLL = _yaddLL;
  jetCfg.useMuonID = _useMuonID;
  jetCfg.muonIDExternal = _muonIDExternal;
  jetCfg.muonIDMinEnergy = _muonIDMinEnergy;
  jetCfg.muonIDMinD0Sig = _muonIDMinD0Sig;
  jetCfg.muonIDMinZ0Sig = _muonIDMinZ0Sig;
  jetCfg.muonIDMinProb = _muonIDMinProb;
  jetCfg.muonIDMaxDist = _muonIDMaxDist;
  JetFinder* jetFinder = new JetFinder(jetCfg);

  double* ymin = new double[_maxYth];
  memset(ymin, 0, sizeof(double) * _maxYth);

  // select vertices
  vector<const Vertex*> selectedVertices;
  vector<const Track*> residualTracks = event->getTracks();
  if (_vertices) {
    VertexSelectorConfig vscfg;
    vscfg.rejectdist = true;
    vscfg.minpos = _vsMinDist;
    vscfg.maxpos = _vsMaxDist;
    vscfg.rejectk0 = true;
    vscfg.k0width = _vsK0MassWidth/2;

    selectedVertices = VertexSelector()(*_vertices, vscfg, residualTracks,false);
  }

  int nVertexJets;
  vector<Jet*> curjets = jetFinder->prerun(residualTracks, event->getNeutrals(), selectedVertices, &nVertexJets);
  if (_njets.size() > 1) {

    for (unsigned int n=0; n<_njets.size(); n++) {
      //cout << "JetFinder: number of jets = " << _njets[n] << endl;

      jetCfg.nJet = _njets[n];
      jetFinder->Configure(jetCfg);

      if (nVertexJets > _njets[n]) {
        //cout << "JetFinder: number of vertex jets is larger than njet: reprocess from prerun..." << endl;

        // clearing curjets...
        for (unsigned int j=0; j<curjets.size(); j++) {
          delete curjets[j];
        }

        // rerun prerun
        curjets = jetFinder->prerun(residualTracks, event->getNeutrals(), selectedVertices, &nVertexJets);
      }

      vector<Jet*>& jets = *(_jetsmap[_njets[n]]);
      jets = jetFinder->run(curjets, ymin, _maxYth);
      curjets.clear();

      // copy jets to curjets
      for (unsigned int j=0; j<jets.size(); j++) {
        curjets.push_back(new Jet(*jets[j]));
      }

      /*
      if (_outputVertexStoresVertex) {
      	vector<Vertex*> &jetvtx = *(_jetvtxmap[_njets[n]]);
      	jetvtx.clear();
      	for(unsigned int j=0;j<jets.size(); j++){
      		vector<const Vertex*> v = jets[j]->getVertices();
      		for(unsigned int k=0; k<v.size(); ++k) {
      			jetvtx.push_back( const_cast<Vertex*>(v[k]) );
      			//jetvtx.push_back( new Vertex(*v[k]) );
      		}
      	}
      }
      */
    }
  } else {

    for (unsigned int n=0; n<_ycut.size(); n++) {
      //cout << "JetFinder: YCut = " << _ycut[n] << endl;

      jetCfg.Ycut = _ycut[n];
      jetFinder->Configure(jetCfg);

      vector<Jet*>& jets = *(_jetsmap[_ycut[n]]);
      jets = jetFinder->run(curjets, ymin, _maxYth);
      curjets.clear();

      // copy jets to curjets
      for (unsigned int j=0; j<jets.size(); j++) {
        curjets.push_back(new Jet(*jets[j]));
      }

      /*
      if (_outputVertexStoresVertex) {
      	vector<Vertex*> &jetvtx = *(_jetvtxmap[_ycut[n]]);
      	jetvtx.clear();
      	for(unsigned int j=0;j<jets.size(); j++){
      		vector<const Vertex*> v = jets[j]->getVertices();
      		for(unsigned int k=0; k<v.size(); ++k) {
      			jetvtx.push_back( const_cast<Vertex*>(v[k]) );
      			//jetvtx.push_back( new Vertex(*v[k]) );
      		}
      	}
      }
      */
    }
  }

  if (_vertices) {
    //cout << "JetClustering: _vertices.size() = " << _vertices->size() << endl;
    //cout << "JetClustering: selectedVertices.size() = " << selectedVertices.size() << endl;
  }

  // writing ymins to Jet collection

  // jetclustering to end for ymin array
  jetCfg.nJet = 1;
  jetCfg.Ycut = 0;
  jetFinder->Configure(jetCfg);

  vector<Jet*> jets = jetFinder->run(curjets, ymin, _maxYth);
  // delete jets
  for (unsigned int i=0; i<jets.size(); i++)
    delete jets[i];
  jets.clear();

  Parameters param;
  for (int i=0; i<_maxYth; i++) {
    param.add(TString::Format("y%d%d",i,i+1), (double)ymin[i]);

    //cout << "y" << i << i+1 << " = " << ymin[i] << endl;
  }

  map<double, vector<Jet*> * >::iterator itjmap;
  for (itjmap = _jetsmap.begin(); itjmap != _jetsmap.end(); itjmap ++) {
    vector<Jet*>* pjlist = itjmap->second;

    for (unsigned int i=0; i<pjlist->size(); i++) {
      Jet* j = (*pjlist)[i];
      j->addParam("yth", param);
    }
  }

}

void JetClustering::end() {
  //cout << "ENDO" << endl;
}



void JetVertexRefiner::init(Parameters* param) {
  Algorithm::init(param);

  // collections - jets

  _jincolname = param->get("JetVertexRefiner.InputJetCollectionName",string("VertexJets"));
  string jcolname = param->get("JetVertexRefiner.OutputJetCollectionName",string("RefinedJets"));
  Event::Instance()->Register(jcolname.c_str(), _outputJets, EventStore::PERSIST | EventStore::JET_WRITE_VERTEX);

  // collections - vertices
  _vprimcolname = param->get("JetVertexRefiner.PrimaryVertexCollectionName",string("PrimaryVertex"));
  _vincolname = param->get("JetVertexRefiner.InputVertexCollectionName",string("BuildUpVertex"));
  _vv0colname = param->get("JetVertexRefiner.V0VertexCollectionName",string("BuildUpVertex_V0"));
  string voutcolname = param->get("JetVertexRefiner.OutputVertexCollectionName",string("RefinedVertices"));
  Event::Instance()->Register(voutcolname.c_str(), _outvertices, EventStore::PERSIST);

  // parameters - single track finder
  _cfg.minPosSingle = param->get("JetVertexRefiner.MinPosSingle", double(0.3));
  _cfg.maxPosSingle = param->get("JetVertexRefiner.MaxPosSingle", double(30.));
  _cfg.minEnergySingle = param->get("JetVertexRefiner.MinEnergySingle", double(1.));
  _cfg.maxAngleSingle = param->get("JetVertexRefiner.MaxAngleSingle", double(.5));
  _cfg.maxSeparationPerPosSingle = param->get("JetVertexRefiner.MaxSeparationPerPosSingle", double(.1));
  _cfg.mind0SigSingle = param->get("JetVertexRefiner.mind0sigSingle", double(5.));
  _cfg.minz0SigSingle = param->get("JetVertexRefiner.minz0sigSingle", double(5.));
  _oneVtxProbTh = param->get("JetVertexRefiner.OneVertexProbThreshold", double(0.001));
  _maxCharmFlightLengthPerJetEnergy = param->get("JetVertexRefiner.MaxCharmFlightLengthPerJetEnergy", double(0.1));

  _invertices = 0;
  _v0vertices = 0;
  _inputJets = 0;

  //bness setup
  _cfg.useBNess = param->get("JetVertexRefiner.useBNess", bool(0)); 
  _cfg.cutBNess = param->get("JetVertexRefiner.BNessCut", double(-0.80));
  _cfg.cutBNessE1 = param->get("JetVertexRefiner.BNessCutE1", double(-0.15));

  if(_cfg.useBNess){
    string bnessname = param->get("JetVertexRefiner.BNessWeightFileName", string("lcfiweights/TMVAClassification_BDTG_bnesstagger_bjet_noPID.weights.xml"));
    string bnessname1 = param->get("JetVertexRefiner.BNessWeightFileNameE1", string("lcfiweights/TMVAClassification_BDTG_bnesstagger_E1_bjet_noPID.weights.xml"));
    _bnessbookname = param->get("JetVertexRefiner.BNessBookName", string("BDTG_bnesstagger_bjet"));
    _bnessbookname1 = param->get("JetVertexRefiner.BNessBookNameE1", string("BDTG_bnesstagger_E1_bjet"));
    _bness=new TMVA::Reader( "!Color:Silent" );
    
    _bness->AddVariable( "trke/jete", &_var[0]);
    _bness->AddVariable( "trkjttheta", &_var[1]);
    _bness->AddVariable( "trkptrel", &_var[2]);
    _bness->AddVariable( "D0sig", &_var[3]);
    _bness->AddVariable( "Z0sig", &_var[4]);
    _bness->AddVariable( "D0", &_var[5]);
    _bness->AddVariable( "Z0", &_var[6]);
    //_bness->AddVariable( "parttype", &_var[7]);
    _bness->BookMVA( _bnessbookname.c_str(), bnessname.c_str() ); 
    _bness->BookMVA( _bnessbookname1.c_str(), bnessname1.c_str() ); 
    //_bness->BookMVA( "BDTG_bnesstagger_bcjet", "lcfiweights/TMVAClassification_BDTG_bnesstagger_bcjet_noPID.weights.xml" ); 
    //_bness->BookMVA( "BDTG_bnesstagger_E1_bcjet", "lcfiweights/TMVAClassification_BDTG_bnesstagger_E1_bcjet_noPID.weights.xml" ); 
  }
}

void JetVertexRefiner::process() {
  Event* event = Event::Instance();

  if (!_invertices) {
    event->Get(_vincolname.c_str(), _invertices);
    if (!_invertices)throw(Exception("JetVertexRefiner: input vertex collection is invalid."));
  }

  if (!_v0vertices) {
    event->Get(_vv0colname.c_str(), _v0vertices);
    if (!_v0vertices)throw(Exception("JetVertexRefiner: input V0 vertex collection is invalid."));
  }
  if (!_inputJets) {
    event->Get(_jincolname.c_str(), _inputJets);
    if (!_inputJets)throw(Exception("JetVertexRefiner: input jet collection is invalid."));
  }

  //check should use BNess
  bool tmpub = _cfg.useBNess;
  if(_cfg.useBNess){
    _cfg.useBNess = false; 
    if((*_invertices).size()>0){
      for(unsigned int i=0;i<(*_invertices).size();i++){
	string vname = (*_invertices)[i]->getVertexingName();
	//cout << "vertexingname: " << vname.c_str() << endl;
	if(vname.find("AVF")==string::npos) _cfg.useBNess = false;
	else{ 
	  _cfg.useBNess = true; 
	  break;
	}
      }
    }
  }

  /*
  		// clearing old jets : but should be done in EventStore
  		for(unsigned int n=0;n<_outputJets->size();n++)
  			delete (*_outputJets)[n];
  		_outputJets->clear();
  		for(unsigned int n=0;n<_outvertices->size();n++)
  			delete (*_outvertices)[n];
  		_outvertices->clear();
  */

  // copy jet with extracting vertices

  // obtain jetvertices
  const Vertex* ip = event->getPrimaryVertex(_vprimcolname.c_str());

  vector<Vertex*> inVerticesCopy;
  for (unsigned int n=0; n<_invertices->size(); n++) {
    inVerticesCopy.push_back(new Vertex(*((*_invertices)[n])));
  }

  vector<vector<Vertex*> > jetVertices;
  vector<vector<const Track*> > jetResidualTracks;
  algoEtc::connectVerticesToJets(*_inputJets, inVerticesCopy, jetVertices, jetResidualTracks,ip);

  // invoke singletrackfinder
  for (unsigned int j=0; j<_inputJets->size(); j++) {
    // new jet creation
    Jet* nj = new Jet(*(*_inputJets)[j], true);

    if(tmpub){   //flag for BNess tagger
      for(unsigned int ntr =0; ntr<nj->getTracks().size();ntr++){
	double bness=-1.0;   //,cness=-1.0;  for future use
	if(nj->getTracks()[ntr]->E()>=0.0){
	  //evaluate BNess of each track
	  //make sign 
	  double sign=1.0;
	  if(TMath::Cos(nj->Vect().Angle(nj->getTracks()[ntr]->getPos()))<0.0) sign=-1.0;
	  
	  _var[0]=nj->getTracks()[ntr]->E()/nj->E();
	  _var[1]=nj->Vect().Angle(nj->getTracks()[ntr]->Vect());
	  _var[2]=nj->getTracks()[ntr]->P()*TMath::Sin(_var[1]);
	  _var[3]=sign*fabs(nj->getTracks()[ntr]->getD0())/sqrt(nj->getTracks()[ntr]->getCovMatrix()[tpar::d0d0]);
	  _var[4]=sign*fabs(nj->getTracks()[ntr]->getZ0())/sqrt(nj->getTracks()[ntr]->getCovMatrix()[tpar::z0z0]);
	  _var[5]=sign*fabs(nj->getTracks()[ntr]->getD0());
	  _var[6]=sign*fabs(nj->getTracks()[ntr]->getZ0());
	  //_var[7]=tracks[ntr]->getPID();
	  if(nj->getTracks()[ntr]->E()>=1.0) bness=_bness->EvaluateMVA(_bnessbookname.c_str()); 
	  else bness=_bness->EvaluateMVA(_bnessbookname1.c_str()); 

	  //cness for future use
	  //if(nj->getTracks()[ntr]->E()>=1.0) cness=_bness->EvaluateMVA("BDTG_bnesstagger_bcjet"); 
	  //else cness=_bness->EvaluateMVA("BDTG_bnesstagger_E1_bcjet"); 
	}      
	
	const_cast<Track*> (nj->getTracks()[ntr])->setBNess(bness); 
	//for future use
	//const_cast<Track*> (nj->getTracks()[ntr])->setCNess(cness); 
	//cout << "bness: " << nj->getTracks()[ntr]->getBNess() << endl;; 
      }
    }

    vector<Vertex*> singleVtcs = VertexFinderSuehara::makeSingleTrackVertices(constVector(jetVertices[j]), jetResidualTracks[j], *_v0vertices, ip, _cfg);
    if(!_cfg.useBNess) VertexFinderSuehara::recombineVertices(jetVertices[j], singleVtcs);
    else VertexFinderSuehara::recombineVertices(jetVertices[j], singleVtcs, _cfg);

    // final v0 selection
    VertexSelector()(jetVertices[j], _cfg.v0selVertex);

    // selection on one vertex probability
    if (jetVertices[j].size() > 1) {
      vector<const Track*> singletracklist = jetVertices[j][0]->getTracks();
      singletracklist.insert(singletracklist.end(), jetVertices[j][1]->getTracks().begin(), jetVertices[j][1]->getTracks().end());

      Vertex* single = VertexFitterSimple_V()(singletracklist.begin(), singletracklist.end());
      if (single->getProb() > _oneVtxProbTh) {
        delete jetVertices[j][0];
        delete jetVertices[j][1];
        jetVertices[j].clear();
        jetVertices[j].push_back(single);
      } else {
        // caching single vertex probability
        Parameters para;
        para.add("SingleVertexProbability", (double)single->getProb());
        nj->addParam("RefinedVertex", para);

        delete single;
      }
    }else if(jetVertices[j].size() == 1){  
      if(_cfg.useBNess){
	//check bness
	vector<const Track*> tracklist = jetVertices[j][0]->getTracks();
	vector<const Track*> newlist;
	int hbtr=-1;
	double okbness=-1.0;
	
	for(unsigned int ntr=0;ntr<tracklist.size();ntr++){
	  if(tracklist[ntr]->getBNess()>okbness){
	    okbness = tracklist[ntr]->getBNess();
	    hbtr=ntr;
	  }
	  
	  if(tracklist[ntr]->E()>=1.0 && tracklist[ntr]->getBNess()<_cfg.cutBNess) continue; 
	  if(tracklist[ntr]->E()<1.0 && tracklist[ntr]->getBNess()<_cfg.cutBNessE1) continue; 
	  
	  newlist.push_back(tracklist[ntr]);
	}
	
	Vertex* vtx=NULL;
	if(newlist.size()>1){
	  vtx = VertexFitterSimple_V()(newlist.begin(),newlist.end());
	}
	/*else{  //make 2 track vertex using highest score track
	  double okchi2=1.0e+100;
	  for(unsigned int ntr=0;ntr<tracklist.size();ntr++){
	    newlist.clear();
	    if((int)ntr==hbtr) continue;
	    newlist.push_back(tracklist[hbtr]);
	    newlist.push_back(tracklist[ntr]);
	    // if(tracklist.size()>2){
	    //   if(1.0*((*newlist[0])+(*newlist[1])).M() >min(newlist[0]->E(), newlist[1]->E())) continue;
	    // }
	    
	    Vertex *tmpvtx = VertexFitterSimple_V()(newlist.begin(),newlist.end());
	    double tmpmaxchi2 = max(tmpvtx->getChi2Track(newlist[0]), tmpvtx->getChi2Track(newlist[1]));
	    if(okchi2>tmpmaxchi2){
	      okchi2 = tmpmaxchi2;
	      if(vtx!=NULL) delete vtx;
	      vtx = tmpvtx;
	    }else delete tmpvtx;   
	    }
	}*/
	
	if(vtx != NULL){
	  delete jetVertices[j][0];
	  jetVertices[j][0] = vtx;
	}else{
	  delete jetVertices[j][0];
	  jetVertices[j].resize(0);
	}
      }
    }

    // flight length of charm
    if (jetVertices[j].size() > 1) {
      double cflt = (jetVertices[j][1]->getPos() - jetVertices[j][0]->getPos()).Mag();
      if (cflt > _maxCharmFlightLengthPerJetEnergy * (*_inputJets)[j]->E()) {
        // eleminate second vtx
        delete jetVertices[j][1];
        jetVertices[j].resize(1);
      }
    }

    // TODO: split vertices with low probability single vertex

    // set vertex to jets

    if (jetVertices[j].size() > 0) {
      nj->add(jetVertices[j][0]);
      _outvertices->push_back(jetVertices[j][0]);
    }

    if (jetVertices[j].size() > 1) {
      nj->add(jetVertices[j][1]);
      _outvertices->push_back(jetVertices[j][1]);
    }

    _outputJets->push_back(nj);
  }

  //for(unsigned int n=0;n<_outputJets->size();n++){ (*_outputJets)[n]->Print(); }

  _cfg.useBNess = tmpub;
}

void JetVertexRefiner::end() {
  //cout << "ENDO" << endl;
  delete _bness;
}


}

