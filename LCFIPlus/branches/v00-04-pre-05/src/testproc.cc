#include <string>

#include "TFile.h"
#include "TNtuple.h"
#include "TNtupleD.h"

#include "lcfiplus.h"
#include "process.h"
#include "testproc.h"
#include "VertexSelector.h"
#include "algoEtc.h"
#include "VertexFinderSuehara.h"
#include "VertexFitterSimple.h"

using namespace lcfiplus;

namespace lcfiplus{

const Jet * JetMCMatch(JetVec &jets, const MCParticle *mcp, vector<const Track *> &assignedTracks, vector<const Track *> &residualTracks)
{
	const vector<const Track *> *pTracks;
	pTracks = &(Event::Instance()->getTracks());

	vector<const Track *> bTracks;

	vector<int> nTrackInJet;
	nTrackInJet.resize(jets.size());
	vector<int> nVertexTrackInJet;
	nVertexTrackInJet.resize(jets.size());

	int nvtx = 0;
	// get tracks
	for(unsigned int i=0;i<pTracks->size();i++){
		const MCParticle *mcpc = (*pTracks)[i]->getMcp();

		if(mcpc==0)continue;
		if(mcpc->isParent(mcp)){
			bTracks.push_back((*pTracks)[i]);
			for(unsigned int j=0;j<jets.size();j++){
				for(unsigned int k=0;k<jets[j]->getVertices().size();k++){
					const vector<const Track *> &vtr = jets[j]->getVertices()[k]->getTracks();
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

	vector<const Track *> jetTracks = jets[ntijMaxIndex]->getTracks();
	for(unsigned int i=0;i<jets[ntijMaxIndex]->getVertices().size();i++){
		const Vertex *vtx = jets[ntijMaxIndex]->getVertices()[i];
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

	void TestAlgoV0::init(Parameters *param){
		Algorithm::init(param);
		string filename = param->get("FileName",string("testv0.root"));
		_vtxname = param->get("VertexCollectionName",string("BuildUpVertex"));
		_file = new TFile(filename.c_str(),"RECREATE");
		_ntp = new TTree("v0","v0");
		_vertices = 0;

		VtxData &d = _data;
		_ntp->Branch("x",&d.x,"x/D");
		_ntp->Branch("y",&d.y,"y/D");
		_ntp->Branch("z",&d.z,"z/D");
		_ntp->Branch("r",&d.r,"r/D");
		_ntp->Branch("cs",&d.cs,"cs/D");
		_ntp->Branch("phi",&d.phi,"phi/D");
		_ntp->Branch("chrg",&d.chrg,"chrg/D");
		_ntp->Branch("dirdot",&d.dirdot,"dirdot/D");
		_ntp->Branch("dirdot2",&d.dirdot2,"dirdot2/D");
		_ntp->Branch("ntrk",&d.ntrk,"ntrk/I");
		_ntp->Branch("mks",&d.mks,"mks/D");
		_ntp->Branch("ml0",&d.ml0,"ml0/D");
		_ntp->Branch("mconv",&d.mconv,"mconv/D");
		_ntp->Branch("mks2",&d.mks2,"mks2/D");
		_ntp->Branch("ml02",&d.ml02,"ml02/D");
		_ntp->Branch("v0",&d.v0,"v0/I");
		_ntp->Branch("ks",&d.ks,"ks/I");
		_ntp->Branch("l0",&d.l0,"l0/I");
		_ntp->Branch("conv",&d.conv,"conv/I");
		_ntp->Branch("mcpdg1",&d.mcpdg1,"mcpdg1/I");
		_ntp->Branch("mcpdg2",&d.mcpdg2,"mcpdg2/I");
		_ntp->Branch("mcppdg1",&d.mcppdg1,"mcppdg1/I");
		_ntp->Branch("mcppdg2",&d.mcppdg2,"mcppdg2/I");
		_ntp->Branch("mcpp1",&d.mcpp1,"mcpp1/I");
		_ntp->Branch("mcpp2",&d.mcpp2,"mcpp2/I");
	}

	void TestAlgoV0::process(){
		if(!_vertices){ Event::Instance()->Get(_vtxname.c_str(), _vertices); }

		const VertexVec& vtx_list = *_vertices;
		for(unsigned int i=0; i < vtx_list.size(); ++i) {
			const Vertex* vtx = vtx_list[i];
			if (vtx->isPrimary()) continue;

			memset(&_data,0,sizeof(_data));

			_data.x = vtx->getX();
			_data.y = vtx->getY();
			_data.z = vtx->getZ();
			TVector3 pos = vtx->getPos();
			_data.r = pos.Mag();
			_data.cs = pos.CosTheta();
			_data.phi = pos.Phi();

			TVector3 mom;
			TVector3 mom2;
			for (unsigned int j=0; j<vtx->getTracks().size(); ++j) {
				mom += vtx->getTracks()[j]->Vect();
				mom2 += vtx->getTracks()[j]->momentumAtVertex(vtx);
				_data.chrg += vtx->getTracks()[j]->getCharge();
			}

			_data.dirdot = mom.Unit().Dot( pos.Unit() );
			_data.dirdot2 = mom2.Unit().Dot( pos.Unit() );
			_data.ntrk = vtx->getTracks().size();

			// compute ks mass
			if (_data.ntrk == 2 && _data.chrg == 0) {
				const Track* trk1 = vtx->getTracks()[0];
				const Track* trk2 = vtx->getTracks()[1];
				TVector3 mom1 = trk1->Vect();
				TVector3 mom2 = trk2->Vect();
				TVector3 mom1v = trk1->momentumAtVertex(vtx);
				TVector3 mom2v = trk2->momentumAtVertex(vtx);

				TLorentzVector lvec1;
				TLorentzVector lvec2;
				lvec1.SetVectM( mom1, 0.1396 );
				lvec2.SetVectM( mom2, 0.1396 );
				_data.mks = (lvec1+lvec2).M();

				lvec1.SetVectM( mom1v, 0.1396 );
				lvec2.SetVectM( mom2v, 0.1396 );
				_data.mks2 = (lvec1+lvec2).M();

				// compute l0 mass
				TLorentzVector protonForLambda;
				TLorentzVector pionForLambda;
				if (mom1.Mag() > mom2.Mag()) {
					protonForLambda.SetVectM( mom1, 0.9383 );
					pionForLambda.SetVectM( mom2, 0.1396 );
				} else {
					protonForLambda.SetVectM( mom2, 0.9383 );
					pionForLambda.SetVectM( mom1, 0.1396 );
				}
				_data.ml0 = (protonForLambda+pionForLambda).M();

				if (mom1v.Mag() > mom2v.Mag()) {
					protonForLambda.SetVectM( mom1v, 0.9383 );
					pionForLambda.SetVectM( mom2v, 0.1396 );
				} else {
					protonForLambda.SetVectM( mom2v, 0.9383 );
					pionForLambda.SetVectM( mom1v, 0.1396 );
				}
				_data.ml02 = (protonForLambda+pionForLambda).M();

				// compute photon mass
				double ang1 = atan( trk1->getTanLambda() );
				double ang2 = atan( trk2->getTanLambda() );
				_data.mconv = sqrt( mom1.Mag()*mom2.Mag()*(1-cos(ang1-ang2)) );


				const MCParticle* mcp1 = trk1->getMcp();
				const MCParticle* mcp2 = trk2->getMcp();

				if (mcp1 && mcp2) {
					_data.mcpdg1 = mcp1->getPDG();
					_data.mcpdg2 = mcp2->getPDG();
					_data.mcppdg1 = mcp1->getParent()->getPDG();
					_data.mcppdg2 = mcp2->getParent()->getPDG();
				}

				if (mcp1 && mcp2) {
					const MCParticle* parent1 = mcp1->getParent();
					const MCParticle* parent2 = mcp2->getParent();
					const MCParticle* parent = mcp1->getParent();

					_data.mcpp1 = (int)((long long) parent1 );
					_data.mcpp2 = (int)((long long) parent2 );

					if ( abs(mcp1->getPDG())==11 && abs(mcp2->getPDG())==11 && parent->getPDG()==22 ) _data.conv = 1;
					if ( abs(mcp1->getPDG())==211 && abs(mcp2->getPDG())==211 && parent->getPDG()==310 ) _data.ks = 1;
					if ( abs(mcp1->getPDG())==211 && abs(mcp2->getPDG())==2212 && abs(parent->getPDG())==3122 ) _data.l0 = 1;
					if ( abs(mcp2->getPDG())==211 && abs(mcp1->getPDG())==2212 && abs(parent->getPDG())==3122 ) _data.l0 = 1;
				}

				_data.v0 = _data.ks || _data.l0 || _data.conv;
			}
			_ntp->Fill();
		}
	}

	void TestAlgoV0::end() {
		_file->Write();
		_file->Close();
	}

	void ZHHAlgo::init(Parameters *param){
		Algorithm::init(param);

		string filename = param->get("FileName",string("test.root"));
		_jetname = param->get("JetCollectionName6",string("RefinedJets_6"));
		_jetname4 = param->get("JetCollectionName4",string("RefinedJets_4"));

		_file = new TFile(filename.c_str(),"RECREATE");
		_tree = new TTree("tree","tree");

		// mc info
		_tree->Branch("mchdecaypdg",&_d.mchdecaypdg,"mchdecaypdg[2]/I");
		_tree->Branch("mcnb",&_d.mcnb,"mcnb/I");

		// non-jet variables
		_tree->Branch("ycuts",&_d.ycuts,"ycuts[10]/D");

		_tree->Branch("thrust",&_d.thrust,"thrust/D");
		_tree->Branch("thaxis",&_d.thaxis,"thaxis[3]/D");

		// 6-jet variables
		_tree->Branch("bcat",&_d.bcat,"bcat[6]/D");
		_tree->Branch("btag",&_d.btag,"btag[6]/D");
		_tree->Branch("ctag",&_d.ctag,"ctag[6]/D");
		_tree->Branch("ejet",&_d.ejet,"ejet[6]/D");
		_tree->Branch("pxjet",&_d.pxjet,"pxjet[6]/D");
		_tree->Branch("pyjet",&_d.pyjet,"pyjet[6]/D");
		_tree->Branch("pzjet",&_d.pzjet,"pzjet[6]/D");
		_tree->Branch("ntrjet",&_d.ntrjet,"ntrjet[6]/D");

		// combined variables for compatibility
		_tree->Branch("mass",&_d.mass,"mass[15]/D");
		_tree->Branch("ntrjetmin",&_d.ntrjetmin,"ntrjetmin/D");
		_tree->Branch("pmiss",&_d.pmiss,"pmiss[3]/D");
		_tree->Branch("emiss",&_d.emiss,"emiss/D");

		// 4-jet variables
		_tree->Branch("bcat4",&_d.bcat4,"bcat4[6]/D");
		_tree->Branch("btag4",&_d.btag4,"btag4[6]/D");
		_tree->Branch("ctag4",&_d.ctag4,"ctag4[6]/D");
		_tree->Branch("ejet4",&_d.ejet4,"ejet4[6]/D");
		_tree->Branch("pxjet4",&_d.pxjet4,"pxjet4[6]/D");
		_tree->Branch("pyjet4",&_d.pyjet4,"pyjet4[6]/D");
		_tree->Branch("pzjet4",&_d.pzjet4,"pzjet4[6]/D");
		_tree->Branch("ntrjet4",&_d.ntrjet4,"ntrjet4[6]/D");

		_jets = 0;
		_jets4 = 0;
	}



	void ZHHAlgo::process(){
		if(!_jets){
			Event::Instance()->Get(_jetname.c_str(), _jets);
		}
		if(!_jets4){
			Event::Instance()->Get(_jetname4.c_str(), _jets4);
		}

		// check higgs decay & nbs
		const MCParticleVec & mcps = Event::Instance()->getMCParticles();

		_d.mcnb = 0;
		_d.mchdecaypdg[0] = _d.mchdecaypdg[1] = 0;

		int hcount = 0;
		for(unsigned int i=0;i<mcps.size();i++){
			int abspdg = abs(mcps[i]->getPDG());
			int parpdg = 0;
			if(mcps[i]->getParent())parpdg = abs(mcps[i]->getParent()->getPDG());
			if(((abspdg > 500 && abspdg < 600) || (abspdg > 5000 && abspdg < 6000)) && parpdg < 100)
				_d.mcnb ++;

			if(mcps[i]->getPDG() == 25){
				// higgs
				if(mcps[i]->getDaughters().size() != 2){
					cout << "ERR: # of higgs daughters = " << mcps[i]->getDaughters().size() << endl;
					break;
				}
				if(hcount == 2){
					cout << "Too many higgs found!, ignore decay" << endl;
					break;
				}
				_d.mchdecaypdg[hcount++] = abs(mcps[i]->getDaughters()[0]->getPDG());
			}
		}

		// thrust
		vector<TVector3> v;
		const TrackVec &tracks = Event::Instance()->getTracks();
		const NeutralVec &neutrals = Event::Instance()->getNeutrals();

		for(unsigned int n=0;n<tracks.size();n++){
			v.push_back(tracks[n]->Vect());
		}
		for(unsigned int n=0;n<neutrals.size();n++){
			v.push_back(neutrals[n]->Vect());
		}

		TVector3 taxis;
		_d.thrust = algoEtc::calcThrust(v, taxis);
		_d.thaxis[0] = taxis.x();
		_d.thaxis[1] = taxis.y();
		_d.thaxis[2] = taxis.z();

		// sorting btag
		vector<const Jet *> jets;
		jets = *_jets;
		sort(jets.begin(), jets.end(), sortBtag);

		int nmass = 0;
		TLorentzVector totp;
		_d.ntrjetmin = 10000.;

		for(unsigned int nj = 0; nj < 6; nj ++){
			Jet *j = const_cast<Jet *>(jets[nj]); // TODO: bad boy...
			j->recalcFourMomentum();

			_d.btag[nj] = j->getParam("lcfiplus")->get<double>("BTag");
			_d.bcat[nj] = j->getParam("lcfiplus")->get<int>("Category");
			_d.ctag[nj] = j->getParam("lcfiplus")->get<double>("CTag");
			_d.ejet[nj] = j->E();
			_d.pxjet[nj] = j->Px();
			_d.pyjet[nj] = j->Py();
			_d.pzjet[nj] = j->Pz();

			totp += *j;

			unsigned int ntr = j->getAllTracks().size();
			_d.ntrjet[nj] = ntr;
			if(_d.ntrjetmin > ntr)_d.ntrjetmin = ntr;

			// ycut values
			if(nj == 0){
				for(int i=0;i<10;i++){
					TString s;
					s.Form("y%d%d",i,i+1);
					_d.ycuts[i] = j->getParam("yth")->get<double>(s);
				}
			}

			// masses
			if(nj == 5)continue;
			for(unsigned int nj2 = nj + 1; nj2 < 6; nj2 ++){
				const Jet *j2 = jets[nj2];
				TLorentzVector v = *j;
				v += *j2;
				_d.mass[nmass++] = v.M();
			}

		}
		_d.emiss = 500 - totp.E();
		_d.pmiss[0] = -totp.Px();
		_d.pmiss[1] = -totp.Py();
		_d.pmiss[2] = -totp.Pz();

		// sorting btag for
		vector<const Jet *> jets4;
		jets4 = *_jets4;
		sort(jets4.begin(), jets4.end(), sortBtag);

		for(unsigned int nj = 0; nj < 4; nj ++){
			Jet *j = const_cast<Jet *>(jets4[nj]); // TODO: bad boy...
			j->recalcFourMomentum();

			_d.btag4[nj] = j->getParam("lcfiplus")->get<double>("BTag");
			_d.bcat4[nj] = j->getParam("lcfiplus")->get<int>("Category");
			_d.ctag4[nj] = j->getParam("lcfiplus")->get<double>("CTag");
			_d.ejet4[nj] = j->E();
			_d.pxjet4[nj] = j->Px();
			_d.pyjet4[nj] = j->Py();
			_d.pzjet4[nj] = j->Pz();

			unsigned int ntr = j->getAllTracks().size();
			_d.ntrjet4[nj] = ntr;

		}

		_tree->Fill();
	}

	void ZHHAlgo::end()
	{
		_file->Write();
		_file->Close();
	}

	void TestAlgo::init(Parameters *param){
		Algorithm::init(param);

		string filename = param->get("FileName",string("test.root"));
		_privtxname = param->get("PrimaryVertexCollectionName",string("PrimaryVertex"));
		_vtxname = param->get("BuildUpVertexCollectionName",string("BuildUpVertex"));
		_v0vtxname = param->get("V0VertexCollectionName",string("BuildUpVertex_V0"));
		_jetname = param->get("JetCollectionName",string("Durham_2Jets"));
		_vtxsel = param->get("VertexSelection",(int)0);
		_refine = param->get("PerformRefining",(int)0);
		_bbhh = param->get("IsBBHH",int(0));

		_file = new TFile(filename.c_str(),"RECREATE");
		_ntJet = new TNtupleD("ntJet","ntJet","nvtx:1vtxprob:2vtxprob:cflt:ecvtx:ejet:vangle:vmass:esingle");

		_vertices = 0;
		_jets = 0;
	}



	void TestAlgo::process(){
		if(!_vertices){
			Event::Instance()->Get(_vtxname.c_str(), _vertices);
			//cout << "vtx name: " << _vtxname << ", pointer = " << (unsigned int)_vertices << endl;
		}
		if(!_v0vertices){
			Event::Instance()->Get(_v0vtxname.c_str(), _v0vertices);
			//cout << "vtx name: " << _vtxname << ", pointer = " << (unsigned int)_vertices << endl;
		}

		if(!_jets){
			Event::Instance()->Get(_jetname.c_str(), _jets);
			//cout << "jet name: " << _jetname << ", pointer = " << (unsigned int)_jets << endl;
		}

		// check bbbbbb (reject H->WW etc.)
		const MCParticleVec & mcps = Event::Instance()->getMCParticles();
		if(_bbhh){
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
			if(hcount < 2)return;
		}

		// select vertices
		vector<const Track *> residualTracks;
		vector<const Vertex *> selectedVertices;
		const vector<const Vertex *> * pVertices;
		if(_vertices && _vtxsel){
			VertexSelectorConfig vscfg;
			vscfg.rejectdist = true;
			vscfg.minpos = .3;
			vscfg.maxpos = 30.;
			vscfg.rejectk0 = true;
			vscfg.k0width = .01;

			selectedVertices = VertexSelector()(*_vertices, vscfg, residualTracks,false);
			pVertices = &selectedVertices;
		}else{
			pVertices = _vertices;
		}

		cout << "# jet = " << _jets->size() << endl;

		// copy vertices
		vector<Vertex *> vtcs2;
		const Vertex *ip = Event::Instance()->getPrimaryVertex(_privtxname.c_str());
		if(ip == 0)throw(Exception("IP not found!"));

		for(unsigned int v=0;v<pVertices->size(); v++){
			if(!(*pVertices)[v]->isPrimary())
				vtcs2.push_back(new Vertex(*(*pVertices)[v]));
		}		

		cout << "# sec vtx = " << vtcs2.size() << endl;

		vector<vector<Vertex *> > jetVertices;
		vector<vector<const Track *> > jetResidualTracks;
		// jet-vtx association
		lcfiplus::algoEtc::connectVerticesToJets(*_jets, vtcs2, jetVertices, jetResidualTracks,ip);

		VertexFinderSuehara::VertexFinderSueharaConfig cfg;

		for(unsigned int j=0;j<_jets->size();j++){

			// single track probability
			double singleprob = 0;
			double twoprob = 0;
			double cflt = 0;
			double ecvtx = 0;
			double vangle = 0;
			double vmass = 0;
			double esingle = 0;

			if(_refine){
				vector<Vertex *> singleVtcs = VertexFinderSuehara::makeSingleTrackVertices(constVector(jetVertices[j]), jetResidualTracks[j], *_v0vertices, ip, cfg);

				if(jetVertices[j].size() + singleVtcs.size() >= 2){
					cout << "Before recombination:" << endl;
					for(unsigned int k=0;k<jetVertices[j].size();k++)
						jetVertices[j][k]->Print();
					for(unsigned int k=0;k<singleVtcs.size();k++)
						singleVtcs[k]->Print();
				}

				VertexFinderSuehara::recombineVertices(jetVertices[j], singleVtcs);

				// v0 selection again
				VertexSelector()(jetVertices[j], cfg.v0selVertex);

				vector<const Track *> singletracklist;
				if(jetVertices[j].size() > 1){
					twoprob = jetVertices[j][0]->getProb() * jetVertices[j][1]->getProb();
					singletracklist.resize(jetVertices[j][0]->getTracks().size() + jetVertices[j][1]->getTracks().size());
					std::copy(jetVertices[j][0]->getTracks().begin(), jetVertices[j][0]->getTracks().end(), singletracklist.begin());
					std::copy(jetVertices[j][1]->getTracks().begin(), jetVertices[j][1]->getTracks().end(), singletracklist.begin() + jetVertices[j][0]->getTracks().size());

					Vertex *single = VertexFitterSimple_V()(singletracklist.begin(), singletracklist.end());

					singleprob = single->getProb();

					cout << "twoprob: " << jetVertices[j][0]->getProb() << " " << jetVertices[j][1]->getProb() << " oneprob: " << singleprob << endl;

					delete single;

					cflt = (jetVertices[j][1]->getPos() - jetVertices[j][0]->getPos()).Mag();

					// looking for near vertex
					int nnear = (jetVertices[j][0]->getPos().Mag() > jetVertices[j][1]->getPos().Mag() ? 0 : 1);

					for(unsigned int ntr=0; ntr<jetVertices[j][nnear]->getTracks().size(); ntr++){
						const Track *tr = jetVertices[j][nnear]->getTracks()[ntr];
						ecvtx += tr->E();
					}
					cout << "cflt = " << cflt << ", ecvtx = " << ecvtx << ", ejet = " << (*_jets)[j]->E() << endl;

					// single track investigation
					int idx = -1;
					if(jetVertices[j][0]->getTracks().size() == 1) idx = 0;
					if(jetVertices[j][1]->getTracks().size() == 1) idx = 1;

					if(idx >= 0){
						const Track *tr = jetVertices[j][idx]->getTracks()[0];
						TVector3 vpos = jetVertices[j][idx]->getPos();
						vangle = vpos.Angle(tr->Vect());
						vmass = 2 * tr->E() * tr->E() * (1 - cos(vpos.Angle(tr->Vect())));
						esingle = tr->E();
					}

				}
				else if(jetVertices[j].size() == 1){
					singleprob = jetVertices[j][0]->getProb();
				}
			}

			if(jetVertices[j].size() >= 2){
				for(unsigned int k=0;k<jetVertices[j].size();k++)
					jetVertices[j][k]->Print();
			}

			_file->cd();
			_ntJet->Fill((int)jetVertices[j].size(),singleprob, twoprob, cflt, ecvtx, (*_jets)[j]->E(), vangle, vmass, esingle);
		}
	}

#if 0
	void TestAlgo::init(Parameters *param){
		Algorithm::init(param);


		string filename = param->get("FileName",string("test.root"));
		_jetname = param->get("JetCollectionName",string("Durham_6Jets2"));
		_bbhh = param->get("IsBBHH",int(0));

		_file = new TFile(filename.c_str(),"RECREATE");

		_ntJet2 = new TNtuple("ntJet2","ntJet2","nev:njet:nbjetmc:nvtx:nvtxjet:ngoodvtx:fracgoodvtxtrack:ycut:nbjet:fracgoodtrack");
		_nbJet = new TNtuple("nbJet", "number of b tracks in eachjet", "nev:nb1:nb2:nb3:nb4:nb5:nb6:nb11:nb12:nb13:nb14:nb15:nb16");

		_jets = 0;
	}

	void TestAlgo::process(){
		if(!_jets){
			Event::Instance()->Get(_jetname.c_str(), _jets);
			//cout << "jet name: " << _jetname << ", pointer = " << (unsigned int)_jets << endl;
		}
		JetVec &jets = *_jets;
		unsigned int nj = jets.size();

		Event *event = Event::Instance();
		MCParticleVec &mcps = event->getMCParticles();

		// check bbbbbb (reject H->WW etc.)
		if(_bbhh){
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
			if(hcount < 2)return;
		}

		// semistable B
		vector<const MCParticle *> blist;
		for(unsigned int i=0;i<mcps.size();i++){
			if(mcps[i]->isSemiStableB()){
				blist.push_back(mcps[i]);
			}
		}
		cout << "Number of semistable B: " << blist.size() << endl;

//	TNtuple *ntResidual = new TNtuple("ntResidual","ResidualTracks","nev:bid:btracks:mcvx:mcvy:mcvz:d0:d0err:z0:z0err:tre");
		// calculate btracks
/*		int *btracks = new int[blist.size()];
		memset(btracks,0,sizeof(int)*blist.size());
		for(unsigned int i=0;i<tracks.size();i++){
			for(unsigned int k=0;k<blist.size();k++){
				if(tracks[i]->getMcp()->isParent(blist[k]))btracks[k] ++;
			}
		}
*/
		vector<const Track *> assignedTracks;
		vector<const Track *> residualTracks;

		map<const Jet *, int > nbmap;
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
				const Track *tr = jets[nj]->getTracks()[i];
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
					const Track *tr = jets[nj]->getVertices()[nv]->getTracks()[i];
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
		_nbJet->Fill(0,nbjet0[0], nbjet0[1], nbjet0[2], nbjet0[3], nbjet0[4], nbjet0[5],
										 nbjet1[0], nbjet1[1], nbjet1[2], nbjet1[3], nbjet0[4], nbjet0[5]);

		for(unsigned int ib=0;ib<blist.size();ib++){
			vector<const Track *> aTracks, rTracks;
//Jet * JetMCMatch(vector<Jet *> &jets, MCParticle *mcp, vector<Track *> &assignedTracks, vector<Track *> &residualTracks)
			const Jet *jet = JetMCMatch(jets, blist[ib], aTracks, rTracks);
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
		_ntJet2->Fill(0, jets.size(), blist.size(), 0, nvtxjet, 0, 0, 0., nbjet, fracgoodtrack);
	}

#endif

	void TestAlgo::end() {
		_file->Write();
		_file->Close();
	}

}

