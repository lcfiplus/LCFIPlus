// VertexFinderSuehara.cc

#include "lcfiplus.h"
#include "VertexFinderSuehara.h"
#include "VertexFinderTearDown.h"
#include "LcfiInterface.h"

#include "algoEtc.h"
#include "algoSigProb.h"

using namespace lcfiplus;
using namespace lcfiplus::VertexFinderSuehara;

#if 0
vector<lcfiplus::Vertex*> * lcfiplus::VertexFinderSuehara::findCandidates(TrackVec &tracks, double chi2th, double massth, double ipchi2th) {
	vector<lcfiplus::Vertex *> * pvertices;
	pvertices = new vector<lcfiplus::Vertex*>;

	// primary vertex
	Vertex *priVertex = findPrimaryVertex(tracks, 25);

	// copy tracks in the jet into a list
	//const vector<Track *> &v = jets[nj]->getTracks();
	list<const Track *> trackList,trackList2;
	trackList.resize(tracks.size());
	copy(tracks.begin(), tracks.end(), trackList.begin());

	// sort list in IP order
	trackList.sort(SortTracksByIPSig(priVertex));

	LcfiInterface lcfi;

	// create pair-vertex
	cout << "TrackList size: " << trackList.size() << endl;
	list<const Track *>::iterator trkit1, trkit2, trkit3, trkit4;
// 	for(trkit1 = trackList.begin(); trkit1 != trackList.end(); trkit1 ++)
// 	{
// //		double ipchi2 = lcfi.getChi2TrackVtx(priVertex, *trkit1);
// 		double d0sig = trackD0Significance(*trkit1, priVertex);
// 		double z0sig = trackZ0Significance(*trkit1, priVertex);
// 		double ipsig = sqrt(pow(d0sig,2) + pow(z0sig,2));
// 		cout << "Track #" << (*trkit1)->getId() << ", ipsig = " << ipsig;
// 		cout << ", D0 = " << (*trkit1)->getD0() << ", Z0 = " << (*trkit1)->getZ0();
// 		cout << ", D0sig = " << d0sig << ", Z0sig = " << z0sig << endl;
// 	}
// 	trackList.sort(SortTracksByIPSig(priVertex));
// 	for(trkit1 = trackList.begin(); trkit1 != trackList.end(); trkit1 ++)
// 	{
// //		double ipchi2 = lcfi.getChi2TrackVtx(priVertex, *trkit1);
// //		double ipsig = sqrt(pow(trackD0Significance(*trkit1, priVertex),2) + pow(trackZ0Significance(*trkit1,priVertex),2));
// 		(*trkit1)->setFlightLength(0);
// 		double d0sig = trackD0Significance(*trkit1, priVertex);
// 		double flt = (*trkit1)->getFlightLength();
// 		double z0sig = trackZ0Significance(*trkit1, priVertex);
// 		double flt2 = (*trkit1)->getFlightLength();
// 		double ipsig = sqrt(pow(d0sig,2) + pow(z0sig,2));
// 		cout << "Track #" << (*trkit1)->getId() << ", ipsig = " << ipsig;
// 		cout << ", D0 = " << (*trkit1)->getD0() << ", Z0 = " << (*trkit1)->getZ0();
// 		double d0sig2 = trackD0Significance(*trkit1, priVertex);
// 		double z0sig2 = trackZ0Significance(*trkit1, priVertex);
// 		cout << ", D0sig = " << d0sig << "," << d0sig2 << ", Z0sig = " << z0sig << "," << z0sig2;
// 		cout << ", flt = " << flt << ", " << flt2 << endl;
// 	}
// 	for(trkit1 = trackList.begin(); trkit1 != trackList.end(); trkit1 ++)
// 	{
// 		double ipchi2 = lcfi.getChi2TrackVtx(priVertex, *trkit1);
// 		double ipsig = sqrt(pow(trackD0Significance(*trkit1, priVertex),2) + pow(trackZ0Significance(*trkit1,priVertex),2));
// 		cout << "Track #" << (*trkit1)->getId() << ", ipchi2 = " << ipchi2 << ", ipsig = " << ipsig;
// 		cout << ", D0 = " << (*trkit1)->getD0() << ", Z0 = " << (*trkit1)->getZ0();
// 		cout << ", D0sig = " << trackD0Significance(*trkit1, priVertex) << ", Z0sig = " << trackZ0Significance(*trkit1,priVertex) << endl;
// 	}
	for(trkit1 = trackList.begin(); trkit1 != trackList.end(); trkit1 ++)
	{
		double ipchi2 = lcfi.getChi2TrackVtx(priVertex, *trkit1);
		double ipsig = sqrt(pow((*trkit1)->getD0() / sqrt((*trkit1)->getCovMatrix()[tpar::d0d0]),2)
				+ pow((*trkit1)->getZ0() / sqrt((*trkit1)->getCovMatrix()[tpar::z0z0]),2));
//		double ipsig = sqrt(pow(trackD0Significance(*trkit1, priVertex),2) + pow(trackZ0Significance(*trkit1,priVertex),2));
		cout << "Track #" << (*trkit1)->getId() << ", ipchi2 = " << ipchi2 << ", ipsig = " << ipsig;
		cout << ", D0 = " << (*trkit1)->getD0() << ", Z0 = " << (*trkit1)->getZ0();
		cout << ", D0err = " << sqrt((*trkit1)->getCovMatrix()[tpar::d0d0]) << ", Z0err = " << sqrt((*trkit1)->getCovMatrix()[tpar::z0z0]) << endl;

//		if(ipchi2 < ipchi2th)continue;

		trackList2.clear();
		// remove tracks which are already used with trkit1, from trackList2
		for(trkit2 = trkit1, trkit2 ++; trkit2 != trackList.end(); trkit2 ++)
		{
			unsigned int i;
			for(i=0;i<pvertices->size();i++){
				TrackVec & trs = (*pvertices)[i]->getTracks();
				if(find(trs.begin(), trs.end(), *trkit1) != trs.end() && find(trs.begin(), trs.end(), *trkit2) != trs.end())
					break;
			}
	
			if(i == pvertices->size())trackList2.push_back(*trkit2);
		}

		cout << "TrackList2 size: " << trackList2.size() << endl;

		// remove tracks which are already used with trkit1, from trackList2
		for(trkit2 = trackList2.begin(); trkit2 != trackList2.end(); trkit2 ++)
		{
			// mass threshold
			TLorentzVector v1 = **trkit1;
			TLorentzVector v2 = **trkit2;
			double mass = (v1+v2).M();

			if(mass > massth)continue;

			// obtain vertex
			vector<const Track *> vttmp;
			vttmp.push_back(*trkit1);
			vttmp.push_back(*trkit2);

			// make initial point
			HelixClass h1, h2;
// 			HelixClass::Initialize_Canonical(float phi0, float d0, float z0, float omega, 
//                               float tanlambda, float B);
			h1.Initialize_Canonical((*trkit1)->getPhi(), (*trkit1)->getD0(), (*trkit1)->getZ0(), (*trkit1)->getOmega(), (*trkit1)->getTanLambda(), Globals::Instance()->getBField());
			h2.Initialize_Canonical((*trkit2)->getPhi(), (*trkit2)->getD0(), (*trkit2)->getZ0(), (*trkit2)->getOmega(), (*trkit2)->getTanLambda(), Globals::Instance()->getBField());
// float HelixClass::getDistanceToHelix(HelixClass * helix, float * pos, float * mom);
			float hpos[3], hmom[3];
			float dist = h1.getDistanceToHelix(&h2, hpos,hmom);

			float cov[6] = {1,0,1,0,0,1};
			Vertex vtxtmp(1,0.3,hpos[0], hpos[1], hpos[2], cov);
			Vertex *vtx = VertexFitterLCFI_V() (vttmp.begin(), vttmp.end(), &vtxtmp);
//			Vertex *vtx = VertexFitterLCFI_V() (vttmp.begin(), vttmp.end(), NULL);

			// chi2 / primary threshold
// 			if(vtx->getChi2Track(*trkit1) > ipchi2 || vtx->getChi2() > chi2th){
// 				delete vtx;
// 				continue;
// 			}

			cout << "Vertex found: ";
			for(unsigned int i=0;i<vtx->getTracks().size();i++){
				cout << vtx->getTracks()[i]->getId() << ", ";
			}
			cout << "Chi2 = " << vtx->getChi2() << ", Pos = (" << vtx->getX() << ", " << vtx->getY() << ", " << vtx->getZ() << "), ";
			cout << "DIST: " << dist << ", POCA: " << hpos[0] << ", " << hpos[1] << ", " << hpos[2] << ", ";

			const MCParticle *mcp1 = vtx->getTracks()[0]->getMcp();
			const MCParticle *mcp2 = vtx->getTracks()[1]->getMcp();

			const MCParticle *p1 = mcp1->getSemiStableParent();
			const MCParticle *p2 = mcp2->getSemiStableParent();

			bool sameVertex = false;
			list<const MCParticle *>palist;

			const MCParticle *p = p1;
			while(p){
				palist.push_back(p);
				p = p->getSemiStableParent();
			}
			p = p2;
			while(p){
				if(find(palist.begin(), palist.end(), p) != palist.end())sameVertex = true;
				p = p->getSemiStableParent();
			}

			cout << mcp1->getFlavorTagCategory() << mcp2->getFlavorTagCategory() << (p1==p2 ? 1 : 0) << (p1==p2 || sameVertex ? 1 : 0) << ", " ;
			if(p1)
				cout << "TV1: (" << p1->getEndVertex().x() << ", " << p1->getEndVertex().y() << ", " << p1->getEndVertex().z() << ") ";
			if(p2)
				cout << "TV2: (" << p2->getEndVertex().x() << ", " << p2->getEndVertex().y() << ", " << p2->getEndVertex().z() << ") ";

			cout << endl;

/*
			double privChi2 = vtx->getChi2();
			// attach other tracks
			for(trkit3 = trkit2, trkit3 ++; trkit3 != trackList2.end();trkit3 ++){
				TLorentzVector v3 = **trkit3;
				mass = (v1+v3).M();
				if(mass > massth)continue;

				vttmp.push_back(*trkit3);
				Vertex *oldvtx = vtx;
				vtx = VertexFitterLCFI_V() (vttmp.begin(), vttmp.end(), NULL);
				double newchi2 = vtx->getChi2();
				double ndf = vtx->getTracks().size() * 2 - 3;
				if(vtx->getChi2Track(*trkit1) > ipchi2 || newchi2/ndf > chi2th || newchi2 - privChi2 > chi2th){
					vttmp.pop_back();
					delete vtx;
					vtx = oldvtx;
				}
				else{
					delete oldvtx;
					privChi2 = newchi2;

					trkit4 = trkit3--;
					trackList2.erase(trkit4);
					cout << "Vertex updated: number of vertices = " << vtx->getTracks().size() << ": ";
					for(unsigned int i=0;i<vtx->getTracks().size();i++){
						cout << vtx->getTracks()[i]->getId() << ", ";
					}
					cout << "Chi2 = " << vtx->getChi2() << ", Pos = (" << vtx->getX() << ", " << vtx->getY() << ", " << vtx->getZ() << ")" << endl;

				}
			}
*/
/*			cout << "Vertex found: number of vertices = " << vtx->getTracks().size() << ": ";
			for(unsigned int i=0;i<vtx->getTracks().size();i++){
				cout << vtx->getTracks()[i]->getId() << ", ";
			}
			cout << "Chi2 = " << vtx->getChi2() << ", Pos = (" << vtx->getX() << ", " << vtx->getY() << ", " << vtx->getZ() << ")" << endl;*/
			// passed all selections
			pvertices->push_back(vtx);
		}
	}

	return pvertices;
}

lcfiplus::Vertex* VertexFinderSuehara::findOne(list<const Track *> &tracks, double chi2th, double massth, bool removeTracks){
	// copy tracks in the jet into a list
	//const vector<Track *> &v = jets[nj]->getTracks();
	bool verbose = false;

	LcfiInterface lcfi;

	if(verbose)
		cout << "Tracks: " << tracks.size() << endl;
	list<const Track *>::iterator trkit1, trkit2, trkit3, trkit4;

	const Track *tr1, *tr2;
	double minimumchi2 = chi2th;
	Vertex *curvtx = 0;
	TLorentzVector vsum;

	for(trkit1 = tracks.begin(); trkit1 != tracks.end(); trkit1 ++)
	{
		for(trkit2 = trkit1, trkit2 ++; trkit2 != tracks.end(); trkit2 ++)
		{
			// mass threshold
			TLorentzVector v1 = **trkit1;
			TLorentzVector v2 = **trkit2;
			double mass = (v1+v2).M();

			if(mass > massth)continue;

			// obtain vertex
			vector<const Track *> vttmp;
			vttmp.push_back(*trkit1);
			vttmp.push_back(*trkit2);

			// make initial point
			HelixClass h1, h2;
// 			HelixClass::Initialize_Canonical(float phi0, float d0, float z0, float omega, 
//                               float tanlambda, float B);
			h1.Initialize_Canonical((*trkit1)->getPhi(), (*trkit1)->getD0(), (*trkit1)->getZ0(),
					(*trkit1)->getOmega(), (*trkit1)->getTanLambda(), Globals::Instance()->getBField());
			h2.Initialize_Canonical((*trkit2)->getPhi(), (*trkit2)->getD0(), (*trkit2)->getZ0(),
					(*trkit2)->getOmega(), (*trkit2)->getTanLambda(), Globals::Instance()->getBField());
// float HelixClass::getDistanceToHelix(HelixClass * helix, float * pos, float * mom);
			float hpos[3], hmom[3];
			h1.getDistanceToHelix(&h2, hpos,hmom);

			float cov[6] = {1,0,1,0,0,1};
			Vertex vtxtmp(1,0.3,hpos[0], hpos[1], hpos[2], cov);
			Vertex *vtx = VertexFitterLCFI_V() (vttmp.begin(), vttmp.end(), &vtxtmp);

			double chi2 = vtx->getChi2();
			if(minimumchi2 > chi2){
				minimumchi2 = chi2;
				tr1 = *trkit1;
				tr2 = *trkit2;
				delete curvtx;
				curvtx = vtx;
				vsum = v1+v2;
			}
			else{
				delete vtx;
			}
		}
	}

	// pair with < chi2th not found
	if(curvtx == 0)return 0;

	if(verbose)
		cout << "Chi2-minimum: " << curvtx->getChi2() << ", " << tr1->getId() << ", " << tr2->getId() << endl;

	if(removeTracks)tracks.remove(tr1);
	if(removeTracks)tracks.remove(tr2);

	vector<const Track *> vt;
	vt.push_back(tr1);
	vt.push_back(tr2);

	do{
		const Track *tra = 0;
		minimumchi2 = chi2th;
		for(trkit1 = tracks.begin(); trkit1 != tracks.end(); trkit1 ++)
		{
			if(find(vt.begin(), vt.end(), *trkit1) != vt.end())continue;

			TLorentzVector va = **trkit1;
			double mass = (vsum + va).M();
			if(mass > massth)continue;
	
			vt.push_back(*trkit1);
			Vertex *vtx = VertexFitterLCFI_V() (vt.begin(), vt.end(), curvtx);
			int ntr = vtx->getTracks().size();
			double rchi2 = vtx->getChi2() / (ntr*2-3);
			if(minimumchi2 > rchi2){
				minimumchi2 = rchi2;
				tra = *trkit1;
				delete curvtx;
				curvtx = vtx;
				vsum += va;
			}
			else{
				delete vtx;
			}
			vt.pop_back();
		}
		if(tra == 0) return curvtx;

		vt.push_back(tra);
		if(removeTracks)tracks.remove(tra);

		if(verbose)
			cout << "Track added. chi2: " << curvtx->getChi2() << ", " << tra->getId() << endl;
	}while(tracks.size()>0);

	return curvtx;
}

lcfiplus::Vertex* VertexFinderSuehara::findOne2(list<const Track *> &tracks, double chi2th, double massth, bool removeTracks){
	// copy tracks in the jet into a list
	//const vector<Track *> &v = jets[nj]->getTracks();
	bool verbose = false;

	if(verbose)
		cout << "Tracks: " << tracks.size() << endl;
	list<const Track *>::iterator trkit1, trkit2;

	double minchi2_2tr = chi2th;
	double mindist_2tr = 1e+300;
	double mindist_3tr = 1e+300;
	Vertex *curvtx = 0;
	TLorentzVector vsum;

	for(trkit1 = tracks.begin(); trkit1 != tracks.end(); trkit1 ++)
	{
		for(trkit2 = trkit1, trkit2 ++; trkit2 != tracks.end(); trkit2 ++)
		{
			// mass threshold
			TLorentzVector v1 = **trkit1;
			TLorentzVector v2 = **trkit2;
			double mass = (v1+v2).M();

			if(mass > massth)continue;

			// obtain vertex
			vector<const Track *> vttmp;
			vttmp.push_back(*trkit1);
			vttmp.push_back(*trkit2);

			Vertex *vtx = VertexFitterSimple_V() (vttmp.begin(), vttmp.end(), 0);

			double chi2 = max(vtx->getChi2Track(*trkit1), vtx->getChi2Track(*trkit2));
			if(chi2 < chi2th){	// trying 3+ vertex
				Vertex *vtx2 = associateTracks(vtx, tracks, chi2th, massth);
				if(vtx2 != vtx){ // 3+ tracks
					delete vtx;
					vtx = vtx2;
					minchi2_2tr = 0.;
					if(mindist_3tr > vtx->getPos().Mag()){
						mindist_3tr = vtx->getPos().Mag();
						delete curvtx;
						curvtx = vtx;
					}
					else{
						delete vtx;
					}
				}
				else{ // 2 tracks
					if(minchi2_2tr > 0){ // no 3-track vertex
						double dist = vtx->getPos().Mag();
						if((chi2 < 1. && mindist_2tr > dist) || (minchi2_2tr > chi2)){
							minchi2_2tr = max(1., chi2);
							if(chi2<1.)
								mindist_2tr = dist;

							delete curvtx;
							curvtx = vtx;	
						}
						else
							delete vtx;
					}
					else
						delete vtx;
				}
			}
			else{// bad vertex
				delete vtx;
			}
		}
	}

	if(curvtx && removeTracks){
		for(unsigned int i=0;i<curvtx->getTracks().size();i++){
			tracks.remove(curvtx->getTracks()[i]);
		}
	}

	if(verbose && curvtx){
		cout << "Vertex found with " << curvtx->getTracks().size() << " vertex, prob = " << curvtx->getProb() << ", mass = " << curvtx->getVertexMass();
		cout << ", pos = ( " << curvtx->getX() << ", " << curvtx->getY() << ", " << curvtx->getZ() << ")" << endl;
	}

	return curvtx;
}
#endif

// vertex compare functions
bool VertexFinderSuehara::VertexNearer(const Vertex *vtx1, const Vertex *vtx2)
{
	return vtx1->getPos().Mag() < vtx2->getPos().Mag();
}
bool VertexFinderSuehara::VertexProbLarger(const Vertex *vtx1, const Vertex *vtx2)
{
	return vtx1->getProb() > vtx2->getProb();
}

void VertexFinderSuehara::GetVertexList(list<const Track *> &tracks, const Vertex *ip, vector<Vertex *> &vtx, vector<Vertex *> &v0vtx, VertexFinderSueharaConfig &cfg)
{
	// lists of found vertices
	list<Vertex *> tr3list;
	list<Vertex *> tr2list;

	list<const Track *>::iterator trkit1, trkit2;
	list<Vertex *>::iterator vit;

	bool verbose = false;

	list<const Track *> v0tracks;

	// find all vertices
	
	int nv = 0;
	int ntr = tracks.size();
	int ntrmax = (ntr-1) * ntr / 2;

	for(trkit1 = tracks.begin(); trkit1 != tracks.end(); trkit1 ++)
	{
		for(trkit2 = trkit1, trkit2 ++; trkit2 != tracks.end(); trkit2 ++)
		{
			nv ++;
			if(verbose){
				cerr << "Vertex producing: " << nv << "/" << ntrmax << endl;
			}

			// v0tracks rejection
			if(find(v0tracks.begin(), v0tracks.end(), *trkit1)!=v0tracks.end())break;
			if(find(v0tracks.begin(), v0tracks.end(), *trkit2)!=v0tracks.end())continue;

			// mass threshold
			TLorentzVector v1 = **trkit1;
			TLorentzVector v2 = **trkit2;

			double mass = (v1+v2).M();

			// 110216 do not accept tracks with opposite direction
//			if(v1.Vect().Dot(v2.Vect()) < 0.)continue;
			if(mass > min(v1.E(), v2.E()))continue;

			if(mass > cfg.massth)continue;

			// obtain vertex
			vector<const Track *> vttmp;
			vttmp.push_back(*trkit1);
			vttmp.push_back(*trkit2);

			Vertex *vtx = VertexFitterSimple_V() (vttmp.begin(), vttmp.end(), 0);

			double chi2 = max(vtx->getChi2Track(*trkit1), vtx->getChi2Track(*trkit2));
			TVector3 vpos = vtx->getPos();

			if(chi2 < cfg.chi2thV0SelTrack && !VertexSelector().passesCut(vtx, cfg.v0selTrack,ip)){
				//cout << "V0 found at ( " << vtx->getPos().x() << " " << vtx->getPos().y() << " " << vtx->getPos().z() << ") : 2 tracks removed." << endl;
				v0tracks.push_back(*trkit1);
				v0tracks.push_back(*trkit2);
				v0vtx.push_back(vtx);
				break;
			}

			if(!VertexSelector().passesCut(vtx, cfg.v0selVertex,ip)){
				delete vtx;
				continue;
			}

			// direction & chi2 condition
			if(vpos.Dot((v1+v2).Vect()) > 0 && chi2 < cfg.chi2th){	// trying 3+ vertex
				if(verbose)
					cout << "Vertex accepted." << endl;
				Vertex *vtx2 = associateTracks(vtx, constVector(v0vtx), tracks, cfg);
				if(vtx2 != vtx){ // 3+ tracks
					delete vtx;
					tr3list.push_back(vtx2);
				}
				else{ // 2 tracks
					tr2list.push_back(vtx);
				}
			}
			else{// bad vertex
				delete vtx;
			}
		}
	}

	// v0 rejection again
	for(vit = tr3list.begin(); vit != tr3list.end();){
		bool deleted = false;
		for(unsigned int n = 0; n < (*vit)->getTracks().size(); n++){
			const Track *tr = (*vit)->getTracks()[n];
			if(find(v0tracks.begin(), v0tracks.end(), tr) != v0tracks.end()){
				// new vertex with smaller # of tracks
				vector<const Track *> vttmp = (*vit)->getTracks();
				vector<const Track *>::iterator itt = find(vttmp.begin(), vttmp.end(), tr);
				vttmp.erase(itt);

				Vertex *vtx = VertexFitterSimple_V() (vttmp.begin(), vttmp.end(), 0);
				if(vttmp.size()>2)tr3list.push_back(vtx);
				else tr2list.push_back(vtx);

				delete *vit;
				vit = tr3list.erase(vit);
				deleted = true;
				//cout << "Vertex in tr3list removed due to v0 rejection." << endl;
				break;
			}
		}
		if(!deleted)vit++;
	}
	for(vit = tr2list.begin(); vit != tr2list.end();){
		bool deleted = false;
		for(unsigned int n = 0; n < (*vit)->getTracks().size(); n++){
			if(find(v0tracks.begin(), v0tracks.end(), (*vit)->getTracks()[n]) != v0tracks.end()){
				vit = tr2list.erase(vit);
				deleted = true;
				cout << "Vertex in tr2distlist removed due to v0 rejection." << endl;
				break;
			}
		}
		if(!deleted)vit++;
	}

	if(verbose){
		cerr << "Vertex produced. 3tr: " << tr3list.size() << ", 2tr: " << tr2list.size() << endl;
	}

	while(true){
		// sort found vertices
		Vertex *curvtx = 0;
		if(tr3list.size()>0){
//			tr3list.sort(VertexNearer);
			tr3list.sort(VertexProbLarger);
			curvtx = tr3list.front();
			tr3list.pop_front();
		}
		else if(tr2list.size()>0){
			tr2list.sort(VertexProbLarger);
			curvtx = tr2list.front();
			tr2list.pop_front();
		}
		else
			break; // all vertices gone

		if(verbose)
			cerr << "Sort finished." << endl;

		// vertex selected
		vtx.push_back(curvtx);
		if(verbose){
			cout << "Vertex found with " << curvtx->getTracks().size() << " tracks, prob = " << curvtx->getProb() << ", mass = " << curvtx->getVertexMass();
			cout << ", pos = ( " << curvtx->getX() << ", " << curvtx->getY() << ", " << curvtx->getZ() << ")" << endl;
			cout << "err = ( xx: " << curvtx->getCov()[Vertex::xx] << ", yy: " << curvtx->getCov()[Vertex::yy] << ", zz:" << curvtx->getCov()[Vertex::zz];
			cout << ", xy: " << curvtx->getCov()[Vertex::xy] << ", xz: " << curvtx->getCov()[Vertex::xz] << ", yz:" << curvtx->getCov()[Vertex::yz] << ")" << endl;
		}

		// remove vertices/tracks
		list<Vertex *>::iterator itv;

		// tr3list
		if(tr3list.size()>0){
			for(itv = tr3list.begin();itv != tr3list.end();){
				bool deleted = false;
				Vertex *v = (*itv);
				for(unsigned int itr=0;itr<curvtx->getTracks().size();itr++){
					const Track *tr = curvtx->getTracks()[itr];
					if(find(v->getTracks().begin(), v->getTracks().end(), tr) != v->getTracks().end()){
						// new vertex with smaller # of tracks
						vector<const Track *> vttmp = v->getTracks();
						vector<const Track *>::iterator itt = find(vttmp.begin(), vttmp.end(), tr);
						vttmp.erase(itt);

						Vertex *vtx = VertexFitterSimple_V() (vttmp.begin(), vttmp.end(), 0);
						if(vttmp.size()>2)tr3list.push_back(vtx);
						else tr2list.push_back(vtx);

						// vertex removed
						delete v;
						itv = tr3list.erase(itv);
						deleted = true;
						break;
					}
				}
				if(tr3list.size()==0)break;
				if(!deleted)
					itv ++;
			}
		}

		// tr2list
		if(tr2list.size()>0){
			for(itv = tr2list.begin();itv != tr2list.end();){
				bool deleted = false;
				Vertex *v = (*itv);
				for(unsigned int itr=0;itr<curvtx->getTracks().size();itr++){
					const Track *tr = curvtx->getTracks()[itr];
					if(find(v->getTracks().begin(), v->getTracks().end(), tr) != v->getTracks().end()){
						// vertex removed
						delete v;
						itv = tr2list.erase(itv);
						deleted = true;
						break;
					}
				}
				if(tr2list.size()==0)break;
				if(!deleted)
					itv ++;
			}
		}
	}

}


Vertex * VertexFinderSuehara::associateTracks(Vertex *vertex, const VertexVec &v0vtx, list<const Track *> &tracks, VertexFinderSueharaConfig &cfg, list<const Track *> *residualTracks)
{
	vector<const Track *> vt;
	TLorentzVector vsum;

	if(residualTracks){
		residualTracks->resize(tracks.size());
		copy(tracks.begin(), tracks.end(), residualTracks->begin());
	}

	vector<const Track *> v0tracks = algoEtc::extractTracks(v0vtx);

	for(unsigned int i=0;i<vertex->getTracks().size();i++){
		vt.push_back(vertex->getTracks()[i]);
		vsum += *(TLorentzVector *)(vertex->getTracks()[i]);
		if(residualTracks)
			residualTracks->remove(vertex->getTracks()[i]);
	}

	Vertex * curvtx = vertex;

	list<const Track *>::iterator trkit1;
	do{
		const Track *tra = 0;
		double minchi2 = cfg.chi2th;
		for(trkit1 = tracks.begin(); trkit1 != tracks.end(); trkit1 ++)
		{
			if(find(vt.begin(), vt.end(), *trkit1) != vt.end())continue;

			// rejecting v0tracks
			if(find(v0tracks.begin(), v0tracks.end(), *trkit1) != v0tracks.end())continue;

			TLorentzVector va = **trkit1;
			double mass = (vsum + va).M();

			// 110216 do not accept tracks with opposite direction to current tracks
//			if(va.Vect().Dot(vsum.Vect()) < 0.)continue;
			if(mass - vsum.M() > min(vsum.E(), va.E()))continue;

			// 110216 do not accept tracks with opposite direction to vertex
			if(va.Vect().Dot(curvtx->getPos())< 0.)continue;

			if(mass > cfg.massth)continue;
	
			vt.push_back(*trkit1);
			Vertex *vtx = VertexFitterSimple_V() (vt.begin(), vt.end(), curvtx,true);
			int ntr = vtx->getTracks().size();

			double maxchi2 = 0;
			for(int i=0;i<ntr;i++){
				double chi2 = vtx->getChi2Track(vtx->getTracks()[i]);
				if(maxchi2 < chi2)maxchi2 = chi2;
			}
			if(minchi2 > maxchi2){
				minchi2 = maxchi2;
				tra = *trkit1;
				if(curvtx != vertex)
					delete curvtx;
				curvtx = vtx;
			}
			else{
				delete vtx;
			}
			vt.pop_back();
		}
		if(tra == 0) return curvtx;

		vsum += *(TLorentzVector *)tra;

		vt.push_back(tra);
		if(residualTracks)residualTracks->remove(tra);
	}while(tracks.size()>0);

	return curvtx;
}

void VertexFinderSuehara::associateIPTracks(vector<Vertex *> &vertices, Vertex *ip, VertexFinderSueharaConfig &cfg)
{
	bool verbose = false;

	vector<const Track *>::const_iterator it;

	Vertex *vbeam;
	algoEtc::makeBeamVertex(vbeam);

	for(unsigned int i=0;i<vertices.size();i++){
		if(vertices[i]->getPos().Mag()<cfg.minimumdistIP)continue;
		TLorentzVector v;
		for(unsigned int j=0;j<vertices[i]->getTracks().size();j++){
			v += *(vertices[i]->getTracks()[j]);
		}

		vector<const Track *> iptracks;
		for(it = ip->getTracks().begin(); it != ip->getTracks().end();it++){

			TLorentzVector v2 = **it;
			if((v + v2).M() - v.M() > min(v.E(), v2.E())){iptracks.push_back(*it);continue;}
			if(v2.Vect().Dot(vertices[i]->getPos())< 0.){iptracks.push_back(*it);continue;}

			double chi2ip = ip->getChi2Track(*it);
			// vertex fitter
/*			vector<Track *> tracks = vertices[i]->getTracks();
			tracks.push_back(*it);
			Vertex *vtx = VertexFitterSimple_V()(tracks.begin(), tracks.end(), vertices[i],true);
			double chi2new = vtx->getChi2Track(*it);*/

			VertexFitterSimple_V vf;
			double chi2new = vf.getChi2(vertices[i], *it, 0/*chi2mode*/);

			if(chi2ip > chi2new * cfg.chi2ratioIP){
				// invoke vertex fitter
				vector<const Track *> tracks = vertices[i]->getTracks();
				tracks.push_back(*it);
				Vertex *vtx = VertexFitterSimple_V()(tracks.begin(), tracks.end(), vertices[i],true);

				// move track into new vertex
				delete vertices[i];
				vertices[i] = vtx;

				if(verbose)
					cout << "Track # " << (*it)->getId() << " moved to vertex " << i << ", chi2ip = " << chi2ip << ", chi2new = " << chi2new << endl;
			}
			else
				iptracks.push_back(*it);
		}
		if(iptracks.size() < ip->getTracks().size()){
			// tracks removed
			delete ip;
			ip = VertexFitterSimple_V()(iptracks.begin(), iptracks.end(), vbeam);
		}
	}

	delete vbeam;
}

void VertexFinderSuehara::buildUp(TrackVec &tracks, vector<Vertex *> &vtx, vector<Vertex *> &v0vtx, double chi2thpri, VertexFinderSueharaConfig &cfg, Vertex *ip)
{

	// obtain primary vertex
	Vertex *nip = 0;
	if(!ip){
		ip = findPrimaryVertex(tracks, chi2thpri);
//		vtx.push_back(ip);
	}else{
		nip = new Vertex(ip->getChi2(), ip->getProb(), ip->getX(), ip->getY(), ip->getZ(), ip->getCov(), true);
//		vtx.push_back(nip);
	}


	// pickup residuals
	list<const Track *> residualTracks;
	for(unsigned int i=0;i<tracks.size();i++){
		if(find(ip->getTracks().begin(), ip->getTracks().end(), tracks[i]) == ip->getTracks().end()){
			residualTracks.push_back(tracks[i]);
		}else if(nip){
			nip->add(tracks[i], ip->getChi2Track(tracks[i]));
		}
	}

	//cout << "buildUp: secondary tracks " << residualTracks.size() << "/" << tracks.size() << endl;

	// secondary vertex
	GetVertexList(residualTracks, ip, vtx, v0vtx, cfg);
/*
	Vertex *v;
	do{
//	lcfiplus::Vertex* VertexFinderSuehara::findOne2(list<Track *> &tracks, double chi2th, double massth, bool removeTracks){
		v = findOne2(residualTracks, 9, massth, true);
		if(v){
			vtx.push_back(v);
		}
	}while(v);
*/
}

void VertexFinderSuehara::buildUpForJetClustering(TrackVec &tracks, vector<Vertex *> &vtx)
{
	VertexFinderSueharaConfig cfg;
	cfg.chi2th = 25.;

	vector<Vertex *> v0vtx;
	buildUp(tracks, vtx, v0vtx, cfg.chi2th, cfg);
}

vector<Vertex *> VertexFinderSuehara::makeSingleTrackVertices
	(Jet *jet, TrackVec &tracks, VertexVec &v0vtx, const Vertex *ip, VertexFinderSueharaConfig &cfg)
{
	return makeSingleTrackVertices(jet->getVertices(), tracks, v0vtx, ip, cfg);
}

vector<Vertex *> VertexFinderSuehara::makeSingleTrackVertices
	(VertexVec &vtcs, TrackVec &tracks, VertexVec &v0vtx, const Vertex *ip, VertexFinderSueharaConfig &cfg)
{
	vector<Vertex *> singlevtcs;

	vector<const Track *> v0tracks = algoEtc::extractTracks(v0vtx);

	for(unsigned int ntr = 0; ntr < tracks.size(); ntr ++){
		const Track *track = tracks[ntr];
		if(find(v0tracks.begin(), v0tracks.end(), track) != v0tracks.end())continue;
		if(track->E() < cfg.minEnergySingle)continue;

		// d0/z0 cut
		double d0 = track->getD0();
		double d0err = sqrt(track->getCovMatrix()[tpar::d0d0]);
		double z0 = track->getZ0();
		double z0err = sqrt(track->getCovMatrix()[tpar::z0z0]);

		if(fabs(d0/d0err)<cfg.mind0SigSingle && fabs(z0/z0err)<cfg.minz0SigSingle) continue;

		Helix hel(track);

		for(unsigned int nvtcs = 0; nvtcs < vtcs.size(); nvtcs ++){
			const Vertex *vtx = vtcs[nvtcs];

			// rejecting opposite direction
			if(vtx->getPos().Dot(track->Vect()) < 0)continue;

			// angular preselection
			double angle = vtx->getPos().Angle(track->Vect());
			if(angle > cfg.maxAngleSingle)continue;

			// calculate closest point
			VertexLine line(ip->getPos(), vtx->getPos());
			double linedist = 0;
			TVector3 pos = hel.ClosePoint(line, &linedist);
	
			// selection cuts
			if(pos.Mag() < cfg.minPosSingle || pos.Mag() > cfg.maxPosSingle)continue; 
			// rejecting opposite vtx position
			if(pos.Dot(vtx->getPos()) < 0.) continue;

			if(linedist / pos.Mag() > cfg.maxSeparationPerPosSingle)continue;

			// all selection passed: make single track vertex
			double cov[6] = {0.,0.,0.,0.,0.,0.};
			Vertex *newvtx = new Vertex(0,1,pos.x(), pos.y(), pos.z(), cov,false);
			newvtx->add(track);

			singlevtcs.push_back(newvtx);
			break; // end searching for this track
		}
	}

	//cout << "makeSingleTrackVertices: " << singlevtcs.size() << " vertices found." << endl;
	return singlevtcs;
}

void VertexFinderSuehara::optimizeTwoVertices(Vertex *&v1, Vertex *&v2, int nvr)
{
	bool _verbose = false;
	double chi2sum = 0.;

	vector<const Track *> n1tracks = v1->getTracks();
	vector<const Track *> n2tracks = v2->getTracks();

	if(n1tracks.size() > 1) chi2sum += v1->getChi2();
	if(n2tracks.size() > 1) chi2sum += v2->getChi2();

	bool move = true;
	Vertex *oldv1 = v1, *oldv2 = v2;

	// vertex recombination
	for(int n=0;n<nvr && move;n++){
		move = false;
		for(unsigned int ntr = 0; ntr< v1->getTracks().size(); ntr++){
			if(n1tracks.size()<=1)break;

			const Track *curt = v1->getTracks()[ntr];

			n2tracks.push_back(curt);
			Vertex *v2d = VertexFitterSimple_V()(n2tracks.begin(), n2tracks.end());

			if(v2d->getChi2Track(curt) < v1->getChi2Track(curt)){
				if(v2 != oldv2)
					delete v2;
				v2 = v2d;
				n1tracks.erase(find(n1tracks.begin(), n1tracks.end(), curt));
				if(v1 != oldv1)
					delete v1;
				v1 = VertexFitterSimple_V()(n1tracks.begin(), n1tracks.end());
				move = true;
				if (_verbose) cout << "Track " << ntr << " moved from v1 to v2. v1 chi2: " << v1->getChi2Track(curt) << ", v2 chi2: " << v2->getChi2Track(curt) << endl;

				ntr --;
			}else{
				n2tracks.pop_back();
				delete v2d;
			}
		}
		for(unsigned int ntr = 0; ntr< v2->getTracks().size(); ntr++){
			if(n2tracks.size()<=1)break;

			const Track *curt = v2->getTracks()[ntr];
	
			n1tracks.push_back(curt);
			Vertex *v1d = VertexFitterSimple_V()(n1tracks.begin(), n1tracks.end());
	
			if(v1d->getChi2Track(curt) < v2->getChi2Track(curt)){
				if(v1 != oldv1)
					delete v1;
				v1 = v1d;
				n2tracks.erase(find(n2tracks.begin(), n2tracks.end(), curt));
				if(v2 != oldv2)
					delete v2;
				v2 = VertexFitterSimple_V()(n2tracks.begin(), n2tracks.end());
				move = true;
				if (_verbose) cout << "Track " << ntr << " moved from v2 to v1. v1 chi2: " << v1->getChi2Track(curt) << ", v2 chi2: " << v2->getChi2Track(curt) << endl;
			}else{
				n1tracks.pop_back();
				delete v1d;
			}
		}

		if (_verbose) {
			cout << "[VR" << n << "]";
			cout << "Tracks = " << n1tracks.size() << ", " << n2tracks.size();
			cout << ", prob = " << v1->getProb() << ", " << v2->getProb();
			cout << ", chi2 = " << v1->getChi2() << ", " << v2->getChi2() << endl;
		}

		if(!move)break;

		double chi2sumold = chi2sum;

		chi2sum = 0;
		if(n1tracks.size() > 1) chi2sum += v1->getChi2();
		if(n2tracks.size() > 1) chi2sum += v2->getChi2();

		if(chi2sumold < chi2sum){
			if (_verbose) cout << "Chi2 by current VR is bigger than old one; rollback to old one." << endl;
			if(v1 != oldv1){
				delete v1;
				v1 = oldv1;
			}if(v2 != oldv2){
				delete v2;
				v2 = oldv2;
			}
			break;
		}else{
			delete oldv1;
			delete oldv2;
			oldv1 = v1;
			oldv2 = v2;
		}
	}
}


void VertexFinderSuehara::recombineVertices(vector<Vertex *> &vertices, vector<Vertex *> &singleVertices)
{
	bool _verbose = false;
	if(vertices.size() + singleVertices.size() == 0){
		if (_verbose) cout << "VertexFinderSuehara::recombineVertices: no vertices & singleVertices found. no-op." << endl;
		return;
	}
	if(vertices.size() + singleVertices.size() == 1){
		if (_verbose) cout << "VertexFinderSuehara::recombineVertices: just one vertices + singleVertices found. Currently no-op. (should be divided...)" << endl;

		if(vertices.size() == 0)vertices.push_back(singleVertices[0]);
		return;
	}

	// firstly prepare combined vertex collection
	vector<Vertex *> v = vertices;
	v.insert(v.end(), singleVertices.begin(), singleVertices.end());

	if(v.size() > 2){
		if (_verbose) cout << "VertexFinderSuehara::recombineVertices: " << v.size() << " vertices + singleVertices found. Trying to combine..." << endl;

		int nvtx = v.size();

		double chi2min = 1e+300;
		Vertex *v1opt = 0;
		Vertex *v2opt = 0;

		for(int n1 = 0; n1 < nvtx - 1; n1 ++){
			for(int n2 = n1 + 1; n2 < nvtx; n2 ++){
				vector<const Track *> n1tracks = v[n1]->getTracks();
				vector<const Track *> n2tracks = v[n2]->getTracks();
				// n1: first reserved vertex; n2: second reserved vertex

				for(int n = 0; n < nvtx; n++){
					if(n == n1 || n == n2)continue;

					const vector<const Track *> &ntracks = v[n]->getTracks();
					for(unsigned int ntr = 0;ntr < ntracks.size();ntr++){
						n1tracks.push_back(ntracks[ntr]);
						Vertex * n1v = VertexFitterSimple_V()(n1tracks.begin(), n1tracks.end());
						double n1chi2 = n1v->getChi2Track(ntracks[ntr]);
						delete n1v;

						n2tracks.push_back(ntracks[ntr]);
						Vertex * n2v = VertexFitterSimple_V()(n2tracks.begin(), n2tracks.end());
						double n2chi2 = n2v->getChi2Track(ntracks[ntr]);
						delete n2v;

						if(n1chi2 < n2chi2){ // select n1chi2
							n2tracks.pop_back(); // n1tracks remain added
						}else{
							n1tracks.pop_back(); // n2tracks remain added
						}
					}
				}
				// n1tracks / n2tracks now obtained
				Vertex *v1 = VertexFitterSimple_V()(n1tracks.begin(), n1tracks.end());
				Vertex *v2 = VertexFitterSimple_V()(n2tracks.begin(), n2tracks.end());

				// optimizeTwoVertices here is removed: no better effect, comparing to do last

				double prob = 1.;
				if(n1tracks.size() > 1) prob *= v1->getProb();
				if(n2tracks.size() > 1) prob *= v2->getProb();

				if (_verbose) {
					cout << "n1 = " << n1 << ", n2 = " << n2 << ", tracks = " << n1tracks.size() << ", " << n2tracks.size();
					cout << ", prob = " << v1->getProb() << ", " << v2->getProb();
					cout << ", chi2 = " << v1->getChi2() << ", " << v2->getChi2() << endl;
				}

				// save current vertices
				double chi2 = v1->getChi2() + v2->getChi2();
				if(chi2 < chi2min){
					chi2min = chi2;
					if(v1opt)
						delete v1opt;
					if(v2opt)
						delete v2opt;
					v1opt = v1;
					v2opt = v2;
				}
				else{
					delete v1;
					delete v2;
				}
			}
		}

		for(unsigned int n=0;n<vertices.size();n++){
			delete vertices[n];
		}
		vertices.clear();
		for(unsigned int n=0;n<singleVertices.size();n++){
			delete singleVertices[n];
		}
		singleVertices.clear();

		optimizeTwoVertices(v1opt,v2opt,3);

		if(v1opt->getPos().Mag() < v2opt->getPos().Mag()){
			vertices.push_back(v1opt);
			vertices.push_back(v2opt);
		}else{
			vertices.push_back(v2opt);
			vertices.push_back(v1opt);
		}
	}else{ // vertex number = 2
		optimizeTwoVertices(v[0], v[1],3);

		vertices.clear();
		singleVertices.clear();

		if(v[0]->getPos().Mag() < v[1]->getPos().Mag()){
			vertices.push_back(v[0]);
			vertices.push_back(v[1]);
		}else{
			vertices.push_back(v[1]);
			vertices.push_back(v[0]);
		}
	}

}
