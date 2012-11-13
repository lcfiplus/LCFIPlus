// TrackSelector.h

#ifndef TrackSelector_h
#define TrackSelector_h 1

#include "lcfiplus.h"
#include <vector>
#include <cmath>

using namespace std;
using namespace lcfiplus;

namespace lcfiplus {

	class TrackSelectorConfig{
	public:
		// cuts which are combined using the AND scheme
		double minD0;
		double maxD0;
		double minD0Err;
		double maxD0Err;
		double minD0Sig;
		double maxD0Sig;
		double minZ0;
		double maxZ0;
		double minZ0Err;
		double maxZ0Err;
		double minZ0Sig;
		double maxZ0Sig;
		double minD0Z0Sig;
		double maxD0Z0Sig;
		double minPt;
		double maxInnermostHitRadius;
		// cuts which are combined using the OR scheme, then AND'd with the AND schemes above
		int minTpcHits;
		double minTpcHitsMinPt;
		int minFtdHits;
		int minVtxHits;
		int minVtxPlusFtdHits;

		TrackSelectorConfig(){
			minD0 = 0.;
			maxD0 = 1e+300;
			minD0Err = 0.;
			maxD0Err = 1e+300;
			minD0Sig = 0.;
			maxD0Sig = 1e+300;
			minZ0 = 0.;
			maxZ0 = 1e+300;
			minZ0Err = 0.;
			maxZ0Err = 1e+300;
			minZ0Sig = 0.;
			maxZ0Sig = 1e+300;
			minD0Z0Sig = 0.;
			maxD0Z0Sig = 1e+300;
			minPt = 0.;
			maxInnermostHitRadius = 1e+300;

			minTpcHits = 999999;
			minTpcHitsMinPt = 999999;
			minFtdHits = 999999;
			minVtxHits = 999999;
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

			if (fabs(trk->getD0()) < cfg.minD0) return false;
			if (fabs(trk->getD0()) > cfg.maxD0) return false;
			if (trk->getCovMatrix()[tpar::d0d0] < cfg.minD0Err) return false;
			if (trk->getCovMatrix()[tpar::d0d0] > cfg.maxD0Err) return false;
			double d0sig = fabs(trk->getD0()) / sqrt(trk->getCovMatrix()[tpar::d0d0]);
			if ( d0sig < cfg.minD0Sig) return false;
			if ( d0sig > cfg.maxD0Sig) return false;

			if (fabs(trk->getZ0()) < cfg.minZ0) return false;
			if (fabs(trk->getZ0()) > cfg.maxZ0) return false;
			if (trk->getCovMatrix()[tpar::z0z0] < cfg.minZ0Err) return false;
			if (trk->getCovMatrix()[tpar::z0z0] > cfg.maxZ0Err) return false;
			double z0sig = fabs(trk->getZ0()) / sqrt(trk->getCovMatrix()[tpar::z0z0]);
			if ( z0sig < cfg.minZ0Sig) return false;
			if ( z0sig > cfg.maxZ0Sig) return false;

			if(sqrt(d0sig * d0sig + z0sig * z0sig) < cfg.minD0Z0Sig)return false;
			if(sqrt(d0sig * d0sig + z0sig * z0sig) > cfg.maxD0Z0Sig)return false;

			if (trk->Pt() < cfg.minPt) return false;
			if (trk->getRadiusOfInnermostHit() > cfg.maxInnermostHitRadius) return false;

			// OR cuts
			if (trk->getFtdHits() >= cfg.minFtdHits) return true;
			if (trk->getVtxHits() >= cfg.minVtxHits) return true;
			if (trk->getVtxHits() + trk->getFtdHits() >= cfg.minVtxPlusFtdHits) return true;
			if (trk->getTpcHits() >= cfg.minTpcHits && trk->Pt() > cfg.minTpcHitsMinPt) return true;

			return false;
		}

		//c-tor / d-tor
		TrackSelector(){}
		~TrackSelector(){}
	};
}

#endif //TrackSelector_h
