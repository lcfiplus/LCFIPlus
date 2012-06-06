#include "lcfiplus.h"
#include "JetFinder.h"
#include <cmath>

#include "algoSigProb.h"

namespace lcfiplus{
namespace algoSigProb{

double rpars[7] = {
  1.19312,
  2.17842e-05,
  0.131316,
  6.27001e-06,
  0.630934,
  0.00012063,
  0.0165772 };
double zpars[7] = {
  1.28694,
  2.64263e-05,
  0.121609,
  7.32101e-06,
  0.624213,
  0.000161151,
  0.0138247 };

double trackD0Significance(const Track* trk, const Vertex* pri) {
	// take primary vertex error before minimization 
	// because this is what LCFI does
	trk->setFlightLength(0);
	double x0 = trk->getX();
	double y0 = trk->getY();
	//double priErr = ( pri->getCov()[Vertex::xx]*x0*x0 + 2.0*pri->getCov()[Vertex::xy]*x0*y0 + pri->getCov()[Vertex::yy]*y0*y0 ) / (x0*x0+y0*y0);
	double priErr = 0;

	/*
	printf("pri ntrk = %d\n",(int)pri->getAllTracks().size());
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

	double x = trk->getX();
	double y = trk->getY();

	//cout << "x:" << x << ",y:" << y << endl;

	//double d0 = trk->par[Track::d0];
	double d0 = sqrt( pow(x-pri->getX(),2)+pow(y-pri->getY(),2) );
	double d0err = sqrt( trk->getCovMatrix()[tpar::d0d0] + priErr );
	double d0sig = d0/d0err;

	//if ( d0 != d0 ) printf("d0 nan\n");
	double d0cov = trk->getCovMatrix()[tpar::d0d0];
	//if ( d0cov != d0cov ) printf("d0cov nan\n");
	//if ( d0cov < 0 ) printf("d0cov %f\n",d0cov);

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
double trackZ0Significance(const Track* trk, const Vertex* pri) {
	trk->setFlightLength(0);
	//double priErr = pri->getCov()[Vertex::zz];
	double priErr = 0;

	// include error from primary vertex
	TrackPocaXY pocaXY(trk,pri);
	trk->setFlightLength( pocaXY.getFlightLength() );

	double z = trk->getZ();

	//double z0 = trk->par[Track::z0];
	double z0 = fabs( z-pri->getZ() );
	double z0err = sqrt( trk->getCovMatrix()[tpar::z0z0] + priErr );
	double z0sig = z0/z0err;

	return z0sig;
}

double signedD0Significance(const Track* trk, const Jet* jet, const Vertex* pri, bool updateFlt) {
	TVector3 jet2d( jet->Vect().X(), jet->Vect().Y(), 0);
	trk->setFlightLength(0);

	// take primary vertex error before minimization because this is what LCFI does
	double x0 = trk->getX();
	double y0 = trk->getY();
	double priErr = 0.;
	/*
	double priErr = ( pri->getCov()[Vertex::xx]*x0*x0
			+ 2.0*pri->getCov()[Vertex::xy]*x0*y0
			+ pri->getCov()[Vertex::yy]*y0*y0 ) / (x0*x0+y0*y0);
			*/

	if (updateFlt) {
		TrackPocaXY pocaXY(trk,pri);
		trk->setFlightLength( pocaXY.getFlightLength() );
	}

	// determine sign of significance relative to jet direction
	TVector3 pca( -trk->getD0()*sin(trk->getPhi()),
			trk->getD0()*cos(trk->getPhi()), trk->getZ0() );
	double signd0(1);
	if (pca.Dot(jet2d) < 0) signd0 = -1;

	double d0 = sqrt( pow(trk->getX()-pri->getX(),2)+pow(trk->getY()-pri->getY(),2) );
	double d0errsq = trk->getCovMatrix()[tpar::d0d0] + priErr;
	double d0sig = sqrt( d0*d0/d0errsq )*signd0;

	if (fabs(d0sig)<1e-6) {
		cout << "SMALL D0SIG::::: d0=" << d0 << ", d0errsq=" << d0errsq << ", d0sig=" << d0sig << endl;
	}

	return d0sig;
}

double signedD0(const Track* trk, const Jet* jet, const Vertex* pri, bool updateFlt) {
	TVector3 jet2d( jet->Vect().X(), jet->Vect().Y(), 0);
	trk->setFlightLength(0);

	if (updateFlt) {
		TrackPocaXY pocaXY(trk,pri);
		trk->setFlightLength( pocaXY.getFlightLength() );
	}

	// determine sign of significance relative to jet direction
	TVector3 pca( -trk->getD0()*sin(trk->getPhi()),
			trk->getD0()*cos(trk->getPhi()), trk->getZ0() );
	double signd0(1);
	if (pca.Dot(jet2d) < 0) signd0 = -1;

	double d0 = sqrt( pow(trk->getX()-pri->getX(),2)+pow(trk->getY()-pri->getY(),2) );

	return signd0 * d0;
}

double signedZ0Significance(const Track* trk, const Jet* jet, const Vertex* pri, bool updateFlt) {
	TVector3 jetz( 0, 0, jet->Vect().Z() );

	if (updateFlt) {
		TrackPocaXY pocaXY(trk,pri);
		trk->setFlightLength( pocaXY.getFlightLength() );
	}

	// determine sign of significance relative to jet direction
	TVector3 pca( -trk->getD0()*sin(trk->getPhi()),
			trk->getD0()*cos(trk->getPhi()), trk->getZ0() );
	double signz0(1);
	if (pca.Dot(jetz) < 0) signz0 = -1;
	double z0 = trk->getZ0() - pri->getZ();
	//double z0errsq = trk->getCovMatrix()[tpar::z0z0] + pri->getCov()[Vertex::zz];
	double z0errsq = trk->getCovMatrix()[tpar::z0z0];
	double z0sig = sqrt( z0*z0/z0errsq )*signz0;
	return z0sig;
}

double signedZ0(const Track* trk, const Jet* jet, const Vertex* pri, bool updateFlt) {
	TVector3 jetz( 0, 0, jet->Vect().Z() );

	if (updateFlt) {
		TrackPocaXY pocaXY(trk,pri);
		trk->setFlightLength( pocaXY.getFlightLength() );
	}

	// determine sign of significance relative to jet direction
	TVector3 pca( -trk->getD0()*sin(trk->getPhi()),
			trk->getD0()*cos(trk->getPhi()), trk->getZ0() );
	double signz0(1);
	if (pca.Dot(jetz) < 0) signz0 = -1;
	double z0 = trk->getZ0() - pri->getZ();
	return signz0 * z0;
}

/////////////////////////////////////////////////
// for joint probability calculation
/////////////////////////////////////////////////
double prob1D(double sig, double maxsig, double* pars) {
	double prob = TMath::Erfc(    sig / ( sqrt( double(2) ) * pars[0] ) )
							-TMath::Erfc( maxsig / ( sqrt( double(2) ) * pars[0] ) );
	prob += pars[1]*( exp(-pars[2] * sig ) - exp(-pars[2] * maxsig ) );
	prob += pars[3]*( exp(-pars[4] * sig ) - exp(-pars[4] * maxsig ) );
	prob += pars[5]*( exp(-pars[6] * sig ) - exp(-pars[6] * maxsig ) );
	return prob;
}

double trackProbD0(const Track* trk, const Vertex* pri) {
	double sig = fabs( trackD0Significance(trk,pri) );
	return prob1D(sig,200,rpars)/prob1D(0,200,rpars);
}

double trackProbZ0(const Track* trk, const Vertex* pri) {
	double sig = fabs( trackZ0Significance(trk,pri) );
	return prob1D(sig,200,zpars)/prob1D(0,200,zpars);
}

double jointProbD0(const Jet* jet, const Vertex* pri, double maxd0sigcut) {
	double maxd0sig = 200.;

	double prod(1);
	int ntrk(0);
	double hiprob(0);

	TrackVec & tracks = jet->getAllTracks(true);
	for (TrackVecIte it = tracks.begin(); it != tracks.end(); ++it) {
		const Track* trk = *it;
		if (trk->getD0() < 5 && trk->getZ0() < 5
				&& trk->getVtxHits() >= 5) {

			double sig = fabs( trackD0Significance(trk,pri) );
			if (sig>maxd0sigcut)continue;
			if (sig<maxd0sig) {
				double prob = prob1D(sig,maxd0sig,rpars)/prob1D(0,maxd0sig,rpars);
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

	double jprob(0);
	double factorial(1);

	for (int k=0; k<ntrk; ++k) {
		if (k>0) factorial *= k;
		jprob += pow( -log(prod), k )/factorial;
	}

	jprob *= prod;

	return jprob;
}

double jointProbZ0(const Jet* jet, const Vertex* pri, double maxz0sigcut) {
	double maxz0sig = 200.;

	double prod(1);
	int ntrk(0);
	double hiprob(0);

	TrackVec &tracks = jet->getAllTracks(true);
	for (TrackVecIte it = tracks.begin(); it != tracks.end(); ++it) {
		const Track* trk = *it;
		if (trk->getD0() < 5 && trk->getZ0() < 5
				&& trk->getVtxHits() >= 5) {
			double sig = fabs( trackZ0Significance(trk,pri) );
			if (sig>maxz0sigcut)continue;
			if (sig < maxz0sig) {
				double prob = prob1D(sig,maxz0sig,zpars)/prob1D(0,maxz0sig,zpars);
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

	double jprob(0);
	double factorial(1);

	for (int k=0; k<ntrk; ++k) {
		if (k>0) factorial *= k;
		jprob += pow( -log(prod), k )/factorial;
	}

	jprob *= prod;

	return jprob;
}

///////////////////////////////////////////////////

void findMostSignificantTrack(const Jet* jet, const Vertex* pri, double sigVec[6]) {
	double trk1d0sig(-1e3);
	double trk2d0sig(-1e3);
	const Track* trk1(0);
	const Track* trk2(0);
	const vector<const Track*>& tracks = jet->getAllTracks(true);

	for (TrackVecIte iter = tracks.begin(); iter != tracks.end(); ++iter) {
		const Track* trk = *iter;

		int nHitCut = 5;
		int nVTX = trk->getVtxHits();
		int nFTD = trk->getFtdHits();
		int nhits = nVTX+nFTD;

		if (nhits < nHitCut-1) continue;
		double mom = trk->Vect().Mag();
		if (nhits == nHitCut-1 && mom < 2) continue;
		if (nhits >= nHitCut && mom < 1) continue;

		double d0sig = signedD0Significance(trk,jet,pri,true);

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

	if (trk1) {
		double trk1z0sig = signedZ0Significance(trk1,jet,pri,false);
		double trk1pt = trk1->Vect().Pt();
		sigVec[0] = trk1d0sig;
		sigVec[2] = trk1z0sig;
		sigVec[4] = trk1pt;
	} else {
		sigVec[0] = 0;
		sigVec[2] = 0;
		sigVec[4] = 0;
	}

	if (trk2) {
		double trk2z0sig = signedZ0Significance(trk2,jet,pri,false);
		double trk2pt = trk2->Vect().Pt();
		sigVec[1] = trk2d0sig;
		sigVec[3] = trk2z0sig;
		sigVec[5] = trk2pt;
	} else {
		sigVec[1] = 0;
		sigVec[3] = 0;
		sigVec[5] = 0;
	}
}

}}

