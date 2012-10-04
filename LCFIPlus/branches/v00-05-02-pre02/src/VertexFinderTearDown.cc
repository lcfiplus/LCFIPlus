// VertexFinderTearDown.cc

#include "lcfiplus.h"
#include "VertexFinderTearDown.h"

#include "algoEtc.h"

using namespace lcfiplus;
using namespace lcfiplus::algoEtc;

vector<lcfiplus::Vertex*> * lcfiplus::findTearDownVertices(const Event& evt, const Jet& jet) {
	double chi2 = 9.0;
	bool verbose(false);
	vector<lcfiplus::Vertex *> * pvertices;
	pvertices = new vector<lcfiplus::Vertex*>;

	// copy tracks in the jet into a list
	//const vector<Track *> &v = jets[nj]->getTracks();
	const vector<const Track *> &v = jet.getTracks();
	list<const Track *> tracksInJet;
	tracksInJet.resize(v.size());
	copy(v.begin(), v.end(), tracksInJet.begin());

	while(tracksInJet.size() >= 2){
		list<const Track *> residuals;
		lcfiplus::Vertex *secvtx = lcfiplus::VertexFinderTearDown<list>()(tracksInJet,0, chi2, &residuals);
		if(!secvtx)break;

		pvertices->push_back(secvtx);

		if(verbose)
			cout << "    Secondary vertex found! pos = (" << secvtx->getX() << "," << secvtx->getY() << "," << secvtx->getZ() << "), chi2 = "
				<< secvtx->getChi2() << endl;
		for(unsigned int i=0;i<secvtx->getTracks().size();i++){
			const Track * tr = secvtx->getTracks()[i];
			if(verbose) {
				cout << "        Track #" << i << ": p = (" << tr->Px() << "," << tr->Py() << "," << tr->Pz() << "), chi2 = "
					<< secvtx->getChi2Track(tr) << ", mcFlavorType = "
					<< evt.getMCParticle(tr)->getFlavorTagCategory() << "\n";
			}
		}
		tracksInJet = residuals;
	}

	return pvertices;
}

lcfiplus::Vertex * lcfiplus::findPrimaryVertex(TrackVec &tracks, double chi2, int useBeamConstraint)
{
	Vertex *ip;
	makeBeamVertex(ip);

	Vertex *ret;

	if (useBeamConstraint == 1) {
	  ret =  VertexFinderTearDown<vector, VertexFitterSimple>()(tracks, 0, chi2, 0, ip);
	} else {
	  ret =  VertexFinderTearDown<vector, VertexFitterSimple>()(tracks, 0, chi2, 0, 0);
	}

	// Vertex * ret =  VertexFinderTearDown<vector, VertexFitterSimple>()(tracks, 0, chi2, 0, ip);
	//Vertex * ret =  VertexFinderTearDown<vector>()(tracks, 0, chi2, 0, ip);
	if(ret)
		ret->setPrimary(true);

	if (ret == 0) return ip; // FIXME: this is safety procedure in case primary vertex is not found; need to confirm this is good behavior

	delete ip;
	return ret;
}

// lcfiplus::Vertex * lcfiplus::findPrimaryVertex(const vector<Track *> &tracks, const vector<Track *> &beamTracks, double chi2)
// {
// 	// make point constraint with beam tracks for initial condition
// 	Vertex * ipv = VertexFitterSimple_V()(beamTracks.begin(), beamTracks.end());
// 
// 	// fit with ipv, without beam tracks
// 	Vertex * ret =  VertexFinderTearDown<vector, VertexFitterSimple>()(tracks, 0, chi2, 0, ipv);
// 	delete ipv;
// 	return ret;
// }
