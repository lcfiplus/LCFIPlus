// TrackSelector.h

#ifndef TrackSelector_h
#define TrackSelector_h 1

#include "lcfiplus.h"
#include <vector>

using namespace std;
using namespace lcfiplus;

namespace lcfiplus {

	class TrackSelectorConfig{
	public:
		// cuts which are combined using the AND scheme
		double maxD0;
		double maxD0Err;
		double maxZ0;
		double maxZ0Err;
		double minPt;
		double maxInnermostHitRadius;
		// cuts which are combined using the OR scheme, then AND'd with the AND schemes above
		int minTpcHits;
		int minFtdHits;
		int minVtxHitsWithoutTpcFtd;
		int minVtxPlusFtdHits;

		TrackSelectorConfig(){
			maxD0 = 1e+300;
			maxD0Err = 1e+300;
			maxZ0 = 1e+300;
			maxZ0Err = 1e+300;
			minPt = 0.;
			maxInnermostHitRadius = 1e+300;

			minTpcHits = 0;
			minFtdHits = 0;
			minVtxHitsWithoutTpcFtd = 0;
			minVtxPlusFtdHits = 0;
		}
	};

	class TrackSelector {
	public:
		vector<const Track *> operator () (const vector<const Track *> &tracks, TrackSelectorConfig & config){
			vector<const Track *> ret;

			for(unsigned int i=0;i<tracks.size();i++){
				if(passesCut(tracks[i], config))
					ret.push_back(tracks[i]);
			}

			return ret;
		}

		bool passesCut(const Track *trk, const TrackSelectorConfig &cfg){
			// AND cuts

			if (fabs(trk->getD0()) > cfg.maxD0) return false;
			if (trk->getCovMatrix()[tpar::d0d0] > cfg.maxD0Err) return false;

			if (fabs(trk->getZ0()) > cfg.maxZ0) return false;
			if (trk->getCovMatrix()[tpar::z0z0] > cfg.maxZ0Err) return false;

			if (trk->Pt() < cfg.minPt) return false;
			if (trk->getRadiusOfInnermostHit() > cfg.maxInnermostHitRadius) return false;

			// OR cuts
			if (trk->getTpcHits() >= cfg.minTpcHits) return true;
			if (trk->getFtdHits() >= cfg.minFtdHits) return true;
			if (trk->getVtxHits() >= cfg.minVtxHitsWithoutTpcFtd) return true;
			if (trk->getVtxHits() + trk->getFtdHits() >= cfg.minVtxPlusFtdHits) return true;

			return false;
		}

		//c-tor / d-tor
		TrackSelector(){}
		~TrackSelector(){}
	};
}

#endif //TrackSelector_h
