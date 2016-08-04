#include <map>

#include "VertexMassRecovery.h"

using namespace lcfiplus;

namespace lcfiplus{
  void VertexMassRecovery::init(Parameters *param){
    Algorithm::init(param);
    // collections - jets
    _jincolname = param->get("VertexMassRecovery.InputJetCollectionName",string("RefinedJets"));
    string jcolname = param->get("VertexMassRecovery.OutputJetCollectionName",string("RefinedJets_AttachPi0s"));
    Event::Instance()->Register(jcolname.c_str(), _outputJets, EventStore::PERSIST | EventStore::JET_WRITE_VERTEX);
    
    // collections - vertices
    _vprimcolname = param->get("VertexMassRecovery.PrimaryVertexCollectionName",string("PrimaryVertex"));
    _vincolname = param->get("VertexMassRecovery.InputVertexCollectionName",string("RefinedVertex"));
    string voutcolname = param->get("VertexMassRecovery.OutputVertexCollectionName",string("RefinedVertices_AttachPi0s"));
    Event::Instance()->Register(voutcolname.c_str(), _outvertices, EventStore::PERSIST);
    
    // pi0reco pdf names(to be corrected)
    string pi0pdfname = param->get("VertexMassRecovery.Pi0PDFName",string("pdf/pdf_for_pi0ID_ok.root"));
    string pi0pdfname2d = param->get("VertexMassRecovery.Pi0PDFName2D",string("pdf/pdf_for_pi0ID_2D_ok.root"));

    vector<string> weightfiles;
    vector<string> booknames;
    string tmpstr;
    //weight file names
    tmpstr = param->get("VertexMassRecovery.WeightFileNameForSecondary1",string("TMVAClassification_BDTG_findpi0_1vtx_LCFIPlus.weights.xml"));
    weightfiles.push_back(tmpstr);
    tmpstr = param->get("VertexMassRecovery.WeightFileNameForSecondary2",string("TMVAClassification_BDTG_findpi0_2vtx_second_LCFIPlus.weights.xml"));
    weightfiles.push_back(tmpstr);
    tmpstr = param->get("VertexMassRecovery.WeightFileNameForTertiary",string("TMVAClassification_BDTG_findpi0_2vtx_third_LCFIPlus.weights.xml"));
    weightfiles.push_back(tmpstr);

    //book names
    tmpstr = param->get("VertexMassRecovery.BookNameForSecondary1",string("BDTG_findpi0_1vtx_LCFIPlus"));
    booknames.push_back(tmpstr);
    tmpstr = param->get("VertexMassRecovery.BookNameForSecondary2",string("BDTG_findpi0_2vtx_second_LCFIPlus"));
    booknames.push_back(tmpstr);
    tmpstr = param->get("VertexMassRecovery.BookNameForTertiary",string("BDTG_findpi0_2vtx_third_LCFIPlus"));
    booknames.push_back(tmpstr);

    //operation points
    vector<double> opp;
    double tmpopp;
    tmpopp =  param->get("VertexMassRecovery.CutForSecondary1",double(0.80));
    opp.push_back(tmpopp);
    tmpopp =  param->get("VertexMassRecovery.CutForSecondary2",double(0.60));
    opp.push_back(tmpopp);
    tmpopp =  param->get("VertexMassRecovery.CutForTertiary",double(0.70));
    opp.push_back(tmpopp);

    //pi0 vertex finder
    pi0vtxfinder = new Pi0VertexFinder(weightfiles, booknames, opp, pi0pdfname, pi0pdfname2d);

    _invertices = 0;
    //_v0vertices = 0;
    _inputJets = 0;

  }
  
  void VertexMassRecovery::process(){
    Event *event = Event::Instance();
    
    if(!_invertices){
      event->Get(_vincolname.c_str(), _invertices);
      if(!_invertices)throw(Exception("VertexMassRecovery: input vertex collection is invalid."));
    }
    if(!_inputJets){
      event->Get(_jincolname.c_str(), _inputJets);
      if(!_inputJets)throw(Exception("VertexMassRecovery: input jet collection is invalid."));
    }
    
    // copy jet with extracting vertices
    
    // obtain jetvertices
    const Vertex *ip = event->getPrimaryVertex(_vprimcolname.c_str());
    
    //define pdg table
    map<int, int> pidx;
    pidx.insert( map<int, int>::value_type( 11, 0 ) );
    pidx.insert( map<int, int>::value_type( 13, 1 ) );
    pidx.insert( map<int, int>::value_type( 211, 2 ) );
    pidx.insert( map<int, int>::value_type( 321, 3 ) );
    pidx.insert( map<int, int>::value_type( 2212, 4 ) );

    for(unsigned int k=0;k<_inputJets->size();k++){
      vector<Vertex* > JetVertices;
      
      for(unsigned int j=0; j<(*_inputJets)[k]->getVertices().size();j++){
	JetVertices.push_back(new Vertex(*((*_inputJets)[k]->getVertices()[j])));
      }
      
      //start to attaching pi0s to vertices
      //define variables
      TVector3 vtx, vtxdir;
      TLorentzVector vtxvect,corrvect;
      int nparts[6]={0,0,0,0,0,0};
      //if tertary vertex exists, it is first
      if(JetVertices.size() > 1){
	//-------------- recover third vertex ----------------
	//get vertex position
	vtx = JetVertices[1]->getPos();
	vtxdir = JetVertices[1]->getPos() - JetVertices[0]->getPos();
	//get vertex momentum and particle types on the vertex
	//vertex momentum
	vtxvect = JetVertices[1]->getFourMomentum();
	const vector <const Track*> & vtxtrk = JetVertices[1]->getTracks();
	//particle type
	for(unsigned int i=0;i<6;i++) nparts[i]=0;
	for(unsigned int i=0;i<vtxtrk.size();i++){
	  nparts[0]++;
	  int pid=vtxtrk[i]->getPDG();
	  pid = pidx[pid];

	  if(pid>=0) nparts[pid+1]++;
	  else nparts[3]++; //rejected tracks: move to pion
	}
	
	//pi0 reconstruction
	//first, get neutral pfos in the jet
	const vector <const Neutral*> & jneut = (*_inputJets)[k]->getNeutrals();
	
	//initialize
	pi0vtxfinder->Initialize();
	//get neutrals
	for(unsigned int i=0;i<jneut.size();i++){
	  const Neutral* nt = jneut[i];
	  pi0vtxfinder->Make_GammaVector(nt);
	}
	
	//start to reconstrcut pi0s
	pi0vtxfinder->Reco_Pi0s(vtx, vtxdir);
	
	//attach pi0s: get recovered 4-momentum!
	corrvect=pi0vtxfinder->Corr_VtxVect(2, vtxvect, nparts, vtx, vtxdir);
	//done.
	//set variables to jetvertex
	//4-momentum
	JetVertices[1]->setRecoveredFourMomentum(corrvect);
	//vertex mass
	JetVertices[1]->setRecoveredVertexMass();
	//pi0 4-momentum
	JetVertices[1]->setPi0sFourMomentum(corrvect-vtxvect);
	//num. of pi0s attached
	JetVertices[1]->setNPi0(pi0vtxfinder->Get_nPi0());
	
	//cout << "third vertex done." << endl;
	
	//-------------- recover secondary vertex ----------------
	//get vertex position
	vtx = JetVertices[0]->getPos();
	vtxdir = JetVertices[0]->getPos() - ip->getPos();
	//get vertex momentum and particle types on the vertex
	//vertex momentum
	vtxvect = JetVertices[0]->getFourMomentum();
	const vector <const Track*> & vtxtrk2 = JetVertices[0]->getTracks();
	//particle type
	for(unsigned int i=0;i<6;i++) nparts[i]=0;
	for(unsigned int i=0;i<vtxtrk2.size();i++){
	  nparts[0]++;
	  int pid=vtxtrk2[i]->getPDG();
	  pid = pidx[pid];

	  if(pid>=0) nparts[pid+1]++;
	  else nparts[3]++; //rejected tracks: move to pion
	}
	
	//pi0 reconstruction - gamma vector is already constructed(avoid double use of gammas)   
	//start to reconstrcut pi0s
	pi0vtxfinder->Reco_Pi0s(vtx, vtxdir);
	
	//attach pi0s: get recovered 4-momentum!
	corrvect=pi0vtxfinder->Corr_VtxVect(1, vtxvect, nparts, vtx, vtxdir);
	//done.
	//set variables to jetvertex
	//4-momentum
	JetVertices[0]->setRecoveredFourMomentum(corrvect);
	//vertex mass
	JetVertices[0]->setRecoveredVertexMass();
	//pi0 4-momentum
	JetVertices[0]->setPi0sFourMomentum(corrvect-vtxvect);
	//num. of pi0s attached
	JetVertices[0]->setNPi0(pi0vtxfinder->Get_nPi0());
	
	//output pi0 attached vertices
	Jet *nj= new Jet(*(*_inputJets)[k], true);
	nj->add(JetVertices[0]);
	nj->add(JetVertices[1]);
	_outvertices->push_back(JetVertices[0]);
	_outvertices->push_back(JetVertices[1]);
	_outputJets->push_back(nj);

	//cout << "secondary vertex done." << endl;
      }else if(JetVertices.size()==1){  //just 1 vtx in the jet
	//-------------- recover secondary vertex ----------------
	//get vertex position
	vtx= JetVertices[0]->getPos();
	vtxdir = JetVertices[0]->getPos() - ip->getPos();
	//get vertex momentum and particle types on the vertex
	//vertex momentum
	vtxvect = JetVertices[0]->getFourMomentum();
	const vector <const Track*> & vtxtrk = JetVertices[0]->getTracks();
	//particle type
	for(unsigned int i=0;i<6;i++) nparts[i]=0;
	for(unsigned int i=0;i<vtxtrk.size();i++){
	  nparts[0]++;
	  int pid=vtxtrk[i]->getPDG();
	  pid = pidx[pid];

	  if(pid>=0) nparts[pid+1]++;
	  else nparts[3]++; //rejected tracks: move to pion
	  //cout << "check energy: " << vtxtrk[i]->E() << " " << vtxtrk[i]->M() << endl;
	}
	
	//pi0 reconstruction 
	//first, get neutral pfos in the jet
	const vector <const Neutral*> & jneut = (*_inputJets)[k]->getNeutrals();
	
	//initialize
	pi0vtxfinder->Initialize();
	//get neutrals
	for(unsigned int i=0;i<jneut.size();i++){
	  const Neutral* nt = jneut[i];
	  pi0vtxfinder->Make_GammaVector(nt);
	}
	
	//start to reconstrcut pi0s
	pi0vtxfinder->Reco_Pi0s(vtx, vtxdir);
	
	//attach pi0s: get recovered 4-momentum!
	corrvect=pi0vtxfinder->Corr_VtxVect(0, vtxvect, nparts, vtx, vtxdir);
	//done.
	//set variables to jetvertex
	//4-momentum
	JetVertices[0]->setRecoveredFourMomentum(corrvect);
	//vertex mass
	JetVertices[0]->setRecoveredVertexMass();
	//pi0 4-momentum
	JetVertices[0]->setPi0sFourMomentum(corrvect-vtxvect);
	//num. of pi0s attached
	JetVertices[0]->setNPi0(pi0vtxfinder->Get_nPi0());

	//output pi0 attached vertices
	Jet *nj= new Jet(*(*_inputJets)[k], true);
	nj->add(JetVertices[0]);
	_outvertices->push_back(JetVertices[0]);
	_outputJets->push_back(nj);
      }else if(JetVertices.size()==0){   //do nothing. just add jet
	Jet *nj= new Jet(*(*_inputJets)[k], true);
	_outputJets->push_back(nj);	
      }
    }
    
    //cout << "attaching pi0s done." << endl;
    
    return;
  }
  
  void VertexMassRecovery::end(){
    delete pi0vtxfinder;
  
  }
}
