#include "lcfiplus.h"
#include "JetFinder.h"
#include <cmath>

#include "algoSigProb.h"

namespace lcfiplus{
namespace algoSigProb{

float rpars[7] = {
  1.19312,
  2.17842e-05,
  0.131316,
  6.27001e-06,
  0.630934,
  0.00012063,
  0.0165772 };
float zpars[7] = {
  1.28694,
  2.64263e-05,
  0.121609,
  7.32101e-06,
  0.624213,
  0.000161151,
  0.0138247 };

float trackD0Significance(const Track* trk, const Vertex* pri) {
	// take primary vertex error before minimization 
	// because this is what LCFI does
	trk->setFlightLength(0);
	float x0 = trk->getX();
	float y0 = trk->getY();
	float priErr = ( pri->getCov()[Vertex::xx]*x0*x0 + 2.0*pri->getCov()[Vertex::xy]*x0*y0 + pri->getCov()[Vertex::yy]*y0*y0 ) / (x0*x0+y0*y0);

	/*
	printf("pri ntrk = %d\n",(int)pri->getTracks().size());
	printf("pri: x=%e, y=%e, z=%e\n",
			pri->getX(),
			pri->getY(),
			pri->getZ());

	printf("pri: xx=%e, yy=%e, zz=%e, xy=%e, yz=%e, zx=%e\n",
			pri->getCov()[Vertex::xx],
			pri->getCov()[Vertex::yy],
			pri->getCov()[Vertex::zz],
			pri->getCov()[Vertex::xy],
			pri->getCov()[Vertex::yz],
			pri->getCov()[Vertex::xz]
			);
	 */

	// include error from primary vertex
	TrackPocaXY pocaXY(trk,pri);
	trk->setFlightLength( pocaXY.getFlightLength() );

	float x = trk->getX();
	float y = trk->getY();

	//cout << "x:" << x << ",y:" << y << endl;

	//float d0 = trk->par[Track::d0];
	float d0 = sqrt( pow(x-pri->getX(),2)+pow(y-pri->getY(),2) );
	float d0err = sqrt( trk->getCovMatrix()[tpar::d0d0] + priErr );
	float d0sig = d0/d0err;

	if ( d0 != d0 ) printf("d0 nan\n");
	float d0cov = trk->getCovMatrix()[tpar::d0d0];
	if ( d0cov != d0cov ) printf("d0cov nan\n");
	if ( d0cov < 0 ) printf("d0cov %f\n",d0cov);

	if ( priErr != priErr ) printf("priErr nan\n");
	if ( priErr != priErr ) printf("priErr nan\n");
	if ( d0err != d0err ) printf("d0err nan, d0cov %f, priErr=%f\n", d0cov, priErr);
	/*
	if ( d0sig != d0sig ) {
		printf("d0sig nan\n");
		printf("d0=%f, d0err=%f\n",d0,d0err);
	}
	*/

	/*
	printf("(pri_x,pri_y,pri_z) = (%f,%f,%f)\n",pri->getX(),pri->getY(),pri->getZ());
	printf("(x0,y0) = (%f,%f)\n",x0,y0);
	printf("priErr = %e\n",priErr);
	printf("set flt = %e\n", pocaXY.getFlightLength());
	printf("(x,y) = (%f,%f)\n",x,y);
	printf("d0/d0err/d0sig = %f/%f/%f\n",d0,d0err,d0sig);
	*/

	return d0sig;
}

// z0 significance at the poca taken in the x-y plane
float trackZ0Significance(const Track* trk, const Vertex* pri) {
	trk->setFlightLength(0);
	float priErr = pri->getCov()[Vertex::zz];

	// include error from primary vertex
	TrackPocaXY pocaXY(trk,pri);
	trk->setFlightLength( pocaXY.getFlightLength() );

	float z = trk->getZ();

	//float z0 = trk->par[Track::z0];
	float z0 = fabs( z-pri->getZ() );
	float z0err = sqrt( trk->getCovMatrix()[tpar::z0z0] + priErr );
	float z0sig = z0/z0err;

	return z0sig;
}

float signedD0Significance(const Track* trk, const Jet* jet, const Vertex* pri, bool updateFlt) {
	TVector3 jet2d( jet->Vect().X(), jet->Vect().Y(), 0);
	trk->setFlightLength(0);

	// take primary vertex error before minimization because this is what LCFI does
	float x0 = trk->getX();
	float y0 = trk->getY();
	float priErr = ( pri->getCov()[Vertex::xx]*x0*x0
			+ 2.0*pri->getCov()[Vertex::xy]*x0*y0
			+ pri->getCov()[Vertex::yy]*y0*y0 ) / (x0*x0+y0*y0);

	if (updateFlt) {
		TrackPocaXY pocaXY(trk,pri);
		trk->setFlightLength( pocaXY.getFlightLength() );
	}

	// determine sign of significance relative to jet direction
	TVector3 pca( -trk->getD0()*sin(trk->getPhi()),
			trk->getD0()*cos(trk->getPhi()), trk->getZ0() );
	float signd0(1);
	if (pca.Dot(jet2d) < 0) signd0 = -1;

	float d0 = sqrt( pow(trk->getX()-pri->getX(),2)+pow(trk->getY()-pri->getY(),2) );
	float d0errsq = trk->getCovMatrix()[tpar::d0d0] + priErr;
	float d0sig = sqrt( d0*d0/d0errsq )*signd0;

	return d0sig;
}

float signedZ0Significance(const Track* trk, const Jet* jet, const Vertex* pri, bool updateFlt) {
	TVector3 jetz( 0, 0, jet->Vect().Z() );

	if (updateFlt) {
		TrackPocaXY pocaXY(trk,pri);
		trk->setFlightLength( pocaXY.getFlightLength() );
	}

	// determine sign of significance relative to jet direction
	TVector3 pca( -trk->getD0()*sin(trk->getPhi()),
			trk->getD0()*cos(trk->getPhi()), trk->getZ0() );
	float signz0(1);
	if (pca.Dot(jetz) < 0) signz0 = -1;
	float z0 = trk->getZ0() - pri->getZ();
	float z0errsq = trk->getCovMatrix()[tpar::z0z0] + pri->getCov()[Vertex::zz];
	float z0sig = sqrt( z0*z0/z0errsq )*signz0;
	return z0sig;
}


/////////////////////////////////////////////////
// for joint probability calculation
/////////////////////////////////////////////////
float prob1D(float sig, float maxsig, float* pars) {
	float prob = TMath::Erfc(    sig / ( sqrt( double(2) ) * pars[0] ) )
							-TMath::Erfc( maxsig / ( sqrt( double(2) ) * pars[0] ) );
	prob += pars[1]*( exp(-pars[2] * sig ) - exp(-pars[2] * maxsig ) );
	prob += pars[3]*( exp(-pars[4] * sig ) - exp(-pars[4] * maxsig ) );
	prob += pars[5]*( exp(-pars[6] * sig ) - exp(-pars[6] * maxsig ) );
	return prob;
}

float trackProbD0(const Track* trk, const Vertex* pri) {
	float sig = fabs( trackD0Significance(trk,pri) );
	return prob1D(sig,200,rpars)/prob1D(0,200,rpars);
}

float trackProbZ0(const Track* trk, const Vertex* pri) {
	float sig = fabs( trackZ0Significance(trk,pri) );
	return prob1D(sig,200,zpars)/prob1D(0,200,zpars);
}

float jointProbD0(const Jet* jet, const Vertex* pri) {
	float maxd0sig = 200.;

	float prod(1);
	int ntrk(0);
	float hiprob(0);

	TrackVec & tracks = jet->getTracks();
	for (TrackVecIte it = tracks.begin(); it != tracks.end(); ++it) {
		const Track* trk = *it;
		if (trk->getD0() < 5 && trk->getZ0() < 5
				&& trk->getVtxHits() >= 5) {

			float sig = fabs( trackD0Significance(trk,pri) );
			if (sig<maxd0sig) {
				float prob = prob1D(sig,maxd0sig,rpars)/prob1D(0,maxd0sig,rpars);
				if (prob > hiprob) hiprob = prob;
				prod *= prob;
				++ntrk;
			}
		}
	}

	if (ntrk == 0) {
		return 0;
	}

	if (hiprob == 0) {
		return 0;
	}

	prod *= 1./hiprob;
	--ntrk;

	float jprob(0);
	float factorial(1);

	for (int k=0; k<ntrk; ++k) {
		if (k>0) factorial *= k;
		jprob += pow( -log(prod), k )/factorial;
	}

	jprob *= prod;

	return jprob;
}

float jointProbZ0(const Jet* jet, const Vertex* pri) {
	float maxz0sig = 200.;

	float prod(1);
	int ntrk(0);
	float hiprob(0);

	TrackVec &tracks = jet->getTracks();
	for (TrackVecIte it = tracks.begin(); it != tracks.end(); ++it) {
		const Track* trk = *it;
		if (trk->getD0() < 5 && trk->getZ0() < 5
				&& trk->getVtxHits() >= 5) {
			float sig = fabs( trackZ0Significance(trk,pri) );
			if (sig < maxz0sig) {
				float prob = prob1D(sig,maxz0sig,zpars)/prob1D(0,maxz0sig,zpars);
				if (prob > hiprob) hiprob = prob;
				prod *= prob;
				++ntrk;
			}
		}
	}

	if (ntrk == 0) {
		return 0;
	}

	if (hiprob == 0) {
		return 0;
	}

	prod *= 1./hiprob;
	--ntrk;

	float jprob(0);
	float factorial(1);

	for (int k=0; k<ntrk; ++k) {
		if (k>0) factorial *= k;
		jprob += pow( -log(prod), k )/factorial;
	}

	jprob *= prod;

	return jprob;
}

///////////////////////////////////////////////////

void findMostSignificantTrack(const Jet* jet, const Vertex* pri, float sigVec[6]) {
	float trk1d0sig(-1e3);
	float trk2d0sig(-1e3);
	const Track* trk1(0);
	const Track* trk2(0);
	const vector<const Track*>& tracks = jet->getTracks();

	for (int i=0; i<6; ++i) {
		sigVec[i] = 0;
	}

	if (tracks.size()<2) {
		printf("Number of tracks fewer than two, skipping significance calculation\n");
		return;
	}

	for (TrackVecIte iter = tracks.begin(); iter != tracks.end(); ++iter) {
		const Track* trk = *iter;

		int nHitCut = 5;
		int nVTX = trk->getVtxHits();
		//printf("nVTX: %d\n",nVTX);
		if (nVTX < nHitCut-1) continue;
		double mom = trk->Vect().Mag();
		if (nVTX == nHitCut-1 && mom < 2) continue;
		if (nVTX >= nHitCut && mom < 1) continue;

		float d0sig = signedD0Significance(trk,jet,pri,true);

		if (d0sig > trk1d0sig) {
			trk2d0sig = trk1d0sig;
			trk1d0sig = d0sig;
			trk2 = trk1;
			trk1 = trk;
		} else if (d0sig > trk2d0sig) {
			trk2d0sig = d0sig;
			trk2 = trk;
		}
	}

	if (trk1 == 0 || trk2 == 0) {
		return;
	}

	float trk1z0sig = signedZ0Significance(trk1,jet,pri,false);
	float trk2z0sig = signedZ0Significance(trk2,jet,pri,false);
	float trk1pt = trk1->Vect().Pt();
	float trk2pt = trk2->Vect().Pt();

	sigVec[0] = trk1d0sig;
	sigVec[1] = trk2d0sig;
	sigVec[2] = trk1z0sig;
	sigVec[3] = trk2z0sig;
	sigVec[4] = trk1pt;
	sigVec[5] = trk2pt;
}

}}

