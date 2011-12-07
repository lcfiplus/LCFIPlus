// geometry.cc

#include "lcfiplus.h"
#include "geometry.h"

#include "Math/Functor.h"
//#include "Math/GSLMinimizer1D.h"
#include "Minuit2/Minuit2Minimizer.h"

#include "TError.h"

#include "algo.h"

namespace lcfiplus{

	static bool verbose = false;

	// singleton operations /////////////////////////////////////////////////////
	GeometryHandler * GeometryHandler::_theInstance = NULL;
	GeometryHandler * GeometryHandler::Instance()
	{
		if(_theInstance == NULL) _theInstance = new GeometryHandler;
		return _theInstance;
	}
	// ctor/dtor
	GeometryHandler::GeometryHandler(){
		gErrorIgnoreLevel = kWarning; // to silence Minuit warnings
	}
	GeometryHandler::~GeometryHandler()
	{
		// free the singleton
		_theInstance = NULL;
	}

	// initialization ///////////////////////////////////////////////////////////

	Point::Point(Vertex *vtx){
		_pos = SVector3(vtx->getX(), vtx->getY(), vtx->getZ());
		_err(0,0) = vtx->getCov()[Vertex::xx];
		_err(0,1) = vtx->getCov()[Vertex::xy];
		_err(0,2) = vtx->getCov()[Vertex::xz];
		_err(1,1) = vtx->getCov()[Vertex::yy];
		_err(1,2) = vtx->getCov()[Vertex::yz];
		_err(2,2) = vtx->getCov()[Vertex::zz];
		/* isn't it better to save _errInverse here? */
	}

	// likelihood ///////////////////////////////////////////////////////////////

	// obtain likelihood that this point is at p
	double Point::LogLikelihood(TVector3 &p)const{

		// difference between my and given positions
		SVector3 ps(p.X(), p.Y(), p.Z());
		SVector3 dif = ps - _pos;

		int ret;
		// invert matrix;
		SMatrixSym3 inverr = _err.Inverse(ret);

		// variance: tP x E^-1 x P
		double variance = ROOT::Math::Similarity(dif, inverr);

		if (variance < 0) {
			fprintf(stderr,"error: variance is negative (%f)\n",variance);
			return -1e10;
		}

//		cout << "variance: " << sqrt(variance) << endl;

		// variance to likelihood: currently Gaussian
		return -variance;
		//return TMath::Log10(variance);
	}

	void Point::LogLikelihoodDeriv(TVector3 &p,double* output)const{

		// difference between my and given positions
		SVector3 ps(p.X(), p.Y(), p.Z());
		SVector3 dif = ps - _pos;

		SVector3 difx(p.X(),     0,     0);
		SVector3 dify(    0, p.Y(),     0);
		SVector3 difz(    0,     0, p.Z());

		int ret;
		SMatrixSym3 inverr = _err.Inverse(ret);

		output[0] = -( Dot( Transpose(inverr)*difx, dif ) + Dot( Transpose(inverr)*dif, difx ) );
		output[1] = -( Dot( Transpose(inverr)*dify, dif ) + Dot( Transpose(inverr)*dif, dify ) );
		output[2] = -( Dot( Transpose(inverr)*difz, dif ) + Dot( Transpose(inverr)*dif, difz ) );
	}

	double Helix::Variance(TVector3 &p, double t)const
	{
		SVector3 xyz;
		SMatrixSym3 errXyz;

		// 5-3 conversion
		GetPosErr(t, xyz, errXyz);

		SVector3 dif(p.X() - xyz(0), p.Y() - xyz(1), p.Z() - xyz(2));

		//int ret;
		//errXyz.Inverse(ret);
		//if (ret != 0) fprintf(stderr,"error: matrix inversion failed!\n");
		bool success = errXyz.Invert();
		if (success == false) fprintf(stderr,"error: matrix inversion failed!\n");
		double variance = ROOT::Math::Similarity(dif,errXyz);

		if (variance < 0) {
			fprintf(stderr,"helix variance is negative: %f\n",variance);
			fprintf(stderr,"ref-point=(%f,%f,%f)\n",p.X(),p.Y(),p.Z());
			fprintf(stderr,"pos=(%f,%f,%f) at t=%f\n",xyz(0),xyz(1),xyz(2),t);
			fprintf(stderr,"errXyz=\n");
			fprintf(stderr,"[ %f, %f, %f ] \n",errXyz(0,0),errXyz(0,1),errXyz(0,2));
			fprintf(stderr,"[ %f, %f, %f ] \n",errXyz(1,0),errXyz(1,1),errXyz(1,2));
			fprintf(stderr,"[ %f, %f, %f ] \n",errXyz(2,0),errXyz(2,1),errXyz(2,2));
			return 1e10;
		}

		//cout << "VARIANCE: " << t << ", " << variance << endl;

		return variance;
	}

	double Helix::VarianceDeriv(TVector3 &p, double t)const
	{
		SVector3 xyz;
		SMatrixSym3 errXyz;
		SVector3 dxyz;
		SMatrixSym3 derrXyz;

		// 5-3 conversion
		GetPosErr(t, xyz, errXyz);
		GetPosErrDeriv(t, dxyz, derrXyz);

		SVector3 dif(p.X() - xyz(0), p.Y() - xyz(1), p.Z() - xyz(2));
		SVector3 ddif( -dxyz(0), -dxyz(1), -dxyz(2));

		int ret,ret2,ret3;
		double dvariance =
			ROOT::Math::Similarity(dif,derrXyz.Inverse(ret))
			+ Dot( Transpose( errXyz.Inverse(ret2) )*dif, ddif )
			+ Dot( Transpose( errXyz.Inverse(ret3) )*ddif, dif );

//		cout << "DVARIANCE: " << t << ", " << dvariance << endl;

		return dvariance;
	}

	double Helix::VarianceDeriv2(TVector3 &p, double t)const
	{
		SVector3 xyz;
		SMatrixSym3 errXyz;
		SVector3 dxyz;
		SMatrixSym3 derrXyz;
		SVector3 ddxyz;
		SMatrixSym3 dderrXyz;

		// 5-3 conversion
		GetPosErr(t, xyz, errXyz);
		GetPosErrDeriv(t, dxyz, derrXyz);
		GetPosErrDeriv2(t, ddxyz, dderrXyz);

		SVector3 dif(p.X() - xyz(0), p.Y() - xyz(1), p.Z() - xyz(2));
		SVector3 ddif( -dxyz(0), -dxyz(1), -dxyz(2));
		SVector3 dddif( -ddxyz(0), -ddxyz(1), -ddxyz(2));

		int ret,ret2,ret3;
		double ddvariance =
			    Dot( Transpose(  errXyz.Inverse(ret2) )*dddif,   dif )
			+ 2*Dot( Transpose( derrXyz.Inverse(ret3) )* ddif,   dif )
			+ 2*ROOT::Math::Similarity(ddif,  errXyz.Inverse(ret))
			+   ROOT::Math::Similarity( dif,dderrXyz.Inverse(ret))
			+ 2*Dot( Transpose( derrXyz.Inverse(ret3) )*  dif,  ddif )
			+   Dot( Transpose(  errXyz.Inverse(ret2) )*  dif, dddif );

//		cout << "D^2(VARIANCE): " << t << ", " << dvariance << endl;

		return ddvariance;
	}

	// obtain likelihood that this helix passes p
	double Helix::LogLikelihood(TVector3 &p, double &tmin)const{

//		const double pi = TMath::Pi();

//		double d0 = _hel(id0);
		double z0 = _hel(iz0);
//		double phi0 = _hel(iph);
		double r = fabs(1./_hel(iom));
		double tanlambda = _hel(itd);

// 		// initial minimization: using circle geometry
// 		double x0 = -d0 * sin(phi0);
// 		double y0 =  d0 * cos(phi0);
// 
// 		double xcenter = x0 + r * cos(phi0 - pi/2. * _charge);
// 		double ycenter = y0 + r * sin(phi0 - pi/2. * _charge);
// 
// 		double diff2dx = p.X() - xcenter;
// 		double diff2dy = p.Y() - ycenter;
// 
// 		// nearest angle in the track circle
// 		double phi = atan2(diff2dy, diff2dx);
// 
// 		// z minimization
// 		double phic = atan2(y0 - ycenter, x0 - xcenter);
// 		double zleast = z0 + r * (phic - phi) * _charge * tanlambda;
// 		double zperiod = 2. * pi * r * tanlambda;
		double zleast, zperiod;
		FindZCross(p.X(), p.Y(), zleast,zperiod);
		int nz = int((p.Z()-zleast)/zperiod);

		double z1 = zleast + zperiod * nz;
		double z2 = zleast + zperiod * (nz+1);

		// initial parameter t of helix
		double t1 = (z1 - z0) / (r * tanlambda);
		double t2 = (z2 - z0) / (r * tanlambda);

		double v1 = Variance(p, t1);
		double v2 = Variance(p, t2);

//		cout << "phi: " << phi << ", variance 1: " << t1 << ", " << v1 << ", variance 2: " << t2 << ", " << v2 << endl;

		// single initial point
		double t = (v1 < v2 ? t1 : t2);

#if 0
		// Newton-Raphson method to find the zero of the first derivative of the variance
		// (i.e. minimum of the variance)

		double tolerance=1e-5;
		int maxstep=10000;

		for (int i=0; i<maxstep; ++i) {
			double dv = VarianceDeriv(p,t);
			double d2v = VarianceDeriv2(p,t);
			double delta = dv/d2v;
			t -= delta;
			if (fabs(delta) < tolerance) break;

			//printf("(%3d) t=%f, var=%f\n",i,t,Variance(p,t));
		}
		tmin = t;
		return -Variance(p,tmin);
#endif

#if 0
		VarianceFunctor vf(this,p);
		double bdy = fabs(1.0*3.1415926/_hel(iom));
		golden(t-bdy,t,t+bdy, &vf, 1e-4, tmin);
		return -Variance(p,tmin);
#endif

#if 0
		VarianceFunctor vf(this,p);
		double bdy = fabs(1.0*3.1415926/_hel(iom));
		brent(t-bdy,t,t+bdy, &vf, 1e-8, tmin);
		return -Variance(p,tmin);
#endif

#if 0
		VarianceFunctor vf(this,p);
		VarianceDerivFunctor dvf(this,p);
		double bdy = fabs(1.0*3.1415926/_hel(iom));
		dbrent(t-bdy,t,t+bdy, &vf, &dvf, 1e-4, tmin);
		double variance = Variance(p,tmin);
		if (variance < 0) {
			fprintf(stderr,"error: variance is negative. var=%f, p=(%f,%f,%f), t=%f)\n",variance,p[0],p[1],p[2],tmin);
			return -1e10;
		}
		return -variance;
#endif

#if 1
		// one-dimensional minimization
		VarianceFunctor vf(this,p);
		ROOT::Math::Functor f(vf,1);
//		ROOT::Math::GSLMinimizer1D minimizer;
//		minimizer.SetFunction(f, t, t-(t2-t1)/2, t+(t2-t1)/2);
//		bool b = minimizer.Minimize(100, 1e-12,1e-7);
		ROOT::Minuit2::Minuit2Minimizer min(ROOT::Minuit2::kMigrad);
		//ROOT::Minuit2::Minuit2Minimizer min(ROOT::Minuit2::kSimplex);
		min.SetFunction(f);
//		min.SetMaxFunctionCalls(1000);
		min.SetMaxFunctionCalls(10000);
		min.SetMaxIterations(100);
//		min.SetTolerance(1e-5);
		min.SetTolerance(1e-6);
		min.SetPrintLevel(0);
		// TODO: primary vertex finder sometimes fails
		double tllimit = -1. / r;
		double tulimit = 1000. / r;

		// todo: t<0 treatment
		min.SetLimitedVariable(0,"t",t, 1e-6,tllimit, tulimit);

		min.Minimize();

		tmin = min.X()[0];
		if(verbose){
			cout << "Minimizer result: p = ( " << p.x() << ", " << p.y() << ", " << p.z() << "), ";
			cout << "tmin = " << tmin << ", val = " << vf(&tmin) << endl;
		}

		double variance = vf(&tmin);

//		cout << "Minimizer result: " << b << ", p = ( " << p.x() << ", " << p.y() << ", " << p.z() << "), ";
//		cout << "tmin = " << minimizer.XMinimum() << ", val = " << minimizer.FValMinimum() << endl;

//		return -Variance(p, minimizer.XMinimum());

		if (variance < 0) {
			fprintf(stderr,"error: variance is negative. var=%f, p=(%f,%f,%f), t=%f)\n",variance,p[0],p[1],p[2],tmin);
			return -1e10;
		}
		return -variance;


#endif
	}

	// estimate how the NLL changes as p changes (compute partial derivative)
	// by using the tangent line at the current flight length (t)
	void Helix::LogLikelihoodDeriv(TVector3 &p, double* output)const{

		/*
		// difference between my and given positions
		SVector3 ps(p.X(), p.Y(), p.Z());
		SVector3 dif = ps - _pos;

		SVector3 difx(p.X(),     0,     0);
		SVector3 dify(    0, p.Y(),     0);
		SVector3 difz(    0,     0, p.Z());

		int ret;
		SMatrixSym3 inverr = _err.Inverse(ret);

		output[0] = -(
				Dot( Transpose(inverr)*difx, dif )
				+ Dot( Transpose(inverr)*dif, difx ) );
		output[1] = -(
				Dot( Transpose(inverr)*dify, dif )
				+ Dot( Transpose(inverr)*dif, dify ) );
		output[2] = -(
				Dot( Transpose(inverr)*difz, dif )
				+ Dot( Transpose(inverr)*dif, difz ) );
				*/
	}


	Helix::Helix(const Track *trk){
		_hel(id0) = trk->getD0();
		_hel(iz0) = trk->getZ0();
		_hel(iph) = trk->getPhi();
		_hel(iom) = trk->getOmega();
		_hel(itd) = trk->getTanLambda();

		const float *cov = trk->getCovMatrix();
		/*
		for (int i=0; i<15; ++i) {
			printf("err (%d) = %e\n",i, cov[i]);
		}
		*/
		_err(id0,id0) = cov[tpar::d0d0];
		_err(id0,iz0) = cov[tpar::d0z0];
		_err(id0,iph) = cov[tpar::d0ph];
		_err(id0,iom) = cov[tpar::d0om];
		_err(id0,itd) = cov[tpar::d0td];
		_err(iz0,iz0) = cov[tpar::z0z0];
		_err(iz0,iph) = cov[tpar::z0ph];
		_err(iz0,iom) = cov[tpar::z0om];
		_err(iz0,itd) = cov[tpar::z0td];
		_err(iph,iph) = cov[tpar::phph];
		_err(iph,iom) = cov[tpar::phom];
		_err(iph,itd) = cov[tpar::phtd];
		_err(iom,iom) = cov[tpar::omom];
		_err(iom,itd) = cov[tpar::omtd];
		_err(itd,itd) = cov[tpar::tdtd];

		//std::cout << _err << std::endl;

		_charge = (int)trk->getCharge();
	}

	TVector3 Helix::GetPos(double t)const
	{
		const double pi = TMath::Pi();

		double d0 = _hel(id0);
		double z0 = _hel(iz0);
		double phi0 = _hel(iph);
		double r = fabs(1./_hel(iom));
		double tanlambda = _hel(itd);

		double x = -d0 * sin(phi0) + r * cos(phi0 - pi/2 * _charge) + r * cos(-_charge * t + phi0 + pi/2 * _charge);
		double y =  d0 * cos(phi0) + r * sin(phi0 - pi/2 * _charge) + r * sin(-_charge * t + phi0 + pi/2 * _charge);
		double z =  z0 + r * t * tanlambda;

		return TVector3(x,y,z);
	}

	void Helix::GetPosErr(double t, SVector3 &pos, SMatrixSym3 &err)const
	{
		const double pi = TMath::Pi();

		double d0 = _hel(id0);
		double z0 = _hel(iz0);
		double phi0 = _hel(iph);
		double r = fabs(1./_hel(iom));
		double tanlambda = _hel(itd);

		double x = -d0 * sin(phi0) + r * cos(phi0 - pi/2 * _charge) + r * cos(-_charge * t + phi0 + pi/2 * _charge);
		double y =  d0 * cos(phi0) + r * sin(phi0 - pi/2 * _charge) + r * sin(-_charge * t + phi0 + pi/2 * _charge);
		double z =  z0 + r * t * tanlambda;

		pos(0) = x, pos(1) = y, pos(2) = z;

		// 5-to-3 error-matrix conversion
		SMatrix53 trackToXyz;
		trackToXyz(0,0) = -sin(phi0);			// dx/dd0
		trackToXyz(0,1) = cos(phi0);			// dy/dd0
		trackToXyz(0,2) = 0;							// dz/dd0

		trackToXyz(1,0) = 0;							// dx/dz0
		trackToXyz(1,1) = 0;							// dy/dz0
		trackToXyz(1,2) = 1;							// dz/dz0

		trackToXyz(2,0) = -d0 * cos(phi0) - r * sin(phi0 - pi/2 * _charge) - r * sin(-_charge * t + phi0 + pi/2 * _charge);	// dx/dphi0
		trackToXyz(2,1) = -d0 * sin(phi0) + r * cos(phi0 - pi/2 * _charge) + r * cos(-_charge * t + phi0 + pi/2 * _charge); // dy/dphi0
		trackToXyz(2,2) = 0;							// dz/dphi0
		
		trackToXyz(3,0) = -r*r * (cos(phi0 - pi/2 * _charge) + cos(-_charge * t + phi0 + pi/2 * _charge));									// dx/domega
		trackToXyz(3,1) = -r*r * (sin(phi0 - pi/2 * _charge) + sin(-_charge * t + phi0 + pi/2 * _charge));									// dy/domega
		trackToXyz(3,2) = -r*r * t * tanlambda;// dz/domega

		trackToXyz(4,0) = 0;							// dx/dtanlambda
		trackToXyz(4,1) = 0;							// dy/dtanlambda
		trackToXyz(4,2) = r*t;						// dz/dtanlambda

		/*
		fprintf(stderr,"track matrix *********\n");
		fprintf(stderr,"( %.10e, %.10e, %.10e, %.10e, %.10e )\n",_err(0,0),_err(0,1),_err(0,2),_err(0,3),_err(0,4));
		fprintf(stderr,"( %.10e, %.10e, %.10e, %.10e, %.10e )\n",_err(1,0),_err(1,1),_err(1,2),_err(1,3),_err(1,4));
		fprintf(stderr,"( %.10e, %.10e, %.10e, %.10e, %.10e )\n",_err(2,0),_err(2,1),_err(2,2),_err(2,3),_err(2,4));
		fprintf(stderr,"( %.10e, %.10e, %.10e, %.10e, %.10e )\n",_err(3,0),_err(3,1),_err(3,2),_err(3,3),_err(3,4));
		fprintf(stderr,"( %.10e, %.10e, %.10e, %.10e, %.10e )\n",_err(4,0),_err(4,1),_err(4,2),_err(4,3),_err(4,4));
		//*/

		// xyz error matrix
		err = ROOT::Math::SimilarityT(trackToXyz, _err);
	}

	void Helix::GetPosErrDeriv(double t, SVector3 &dpos, SMatrixSym3 &derr)const
	{
		const double pi = TMath::Pi();

		double d0 = _hel(id0);
		double z0 = _hel(iz0);
		double phi0 = _hel(iph);
		double r = fabs(1./_hel(iom));
		double tanlambda = _hel(itd);

		double dx =  _charge * r * sin(-_charge * t + phi0 + pi/2 * _charge);
		double dy = -_charge * r * cos(-_charge * t + phi0 + pi/2 * _charge);
		double dz =  z0 + r * tanlambda;

		dpos(0) = dx, dpos(1) = dy, dpos(2) = dz;

		// 5-to-3 error-matrix conversion
		SMatrix53 trackToXyz;
		trackToXyz(0,0) = -sin(phi0);			// dx/dd0
		trackToXyz(0,1) = cos(phi0);			// dy/dd0
		trackToXyz(0,2) = 0;							// dz/dd0

		trackToXyz(1,0) = 0;							// dx/dz0
		trackToXyz(1,1) = 0;							// dy/dz0
		trackToXyz(1,2) = 1;							// dz/dz0

		trackToXyz(2,0) = -d0 * cos(phi0) - r * sin(phi0 - pi/2 * _charge) - r * sin(-_charge * t + phi0 + pi/2 * _charge);	// dx/dphi0
		trackToXyz(2,1) = -d0 * sin(phi0) + r * cos(phi0 - pi/2 * _charge) + r * cos(-_charge * t + phi0 + pi/2 * _charge); // dy/dphi0
		trackToXyz(2,2) = 0;							// dz/dphi0
		
		trackToXyz(3,0) = -r*r * (cos(phi0 - pi/2 * _charge) + cos(-_charge * t + phi0 + pi/2 * _charge));									// dx/domega
		trackToXyz(3,1) = -r*r * (sin(phi0 - pi/2 * _charge) + sin(-_charge * t + phi0 + pi/2 * _charge));									// dy/domega
		trackToXyz(3,2) = -r*r * t * tanlambda;// dz/domega

		trackToXyz(4,0) = 0;							// dx/dtanlambda
		trackToXyz(4,1) = 0;							// dy/dtanlambda
		trackToXyz(4,2) = r*t;						// dz/dtanlambda

		// derivs
		SMatrix53 dtrackToXyz;
		dtrackToXyz(0,0) = 0;							// d^2x/dd0dt
		dtrackToXyz(0,1) = 0;							// dy/dd0dt
		dtrackToXyz(0,2) = 0;							// dz/dd0dt

		dtrackToXyz(1,0) = 0;							// d^2x/dz0dt
		dtrackToXyz(1,1) = 0;							// d^2y/dz0dt
		dtrackToXyz(1,2) = 0;							// d^2z/dz0dt

		dtrackToXyz(2,0) = _charge * r * cos(-_charge * t + phi0 + pi/2 * _charge);	// d^2x/dphi0dt
		dtrackToXyz(2,1) = _charge * r * sin(-_charge * t + phi0 + pi/2 * _charge);  // d^2y/dphi0dt
		dtrackToXyz(2,2) = 0;							// d^2z/dphi0dt
		
		dtrackToXyz(3,0) = -_charge * r*r * sin(-_charge * t + phi0 + pi/2 * _charge);									// d^2x/domegadt
		dtrackToXyz(3,1) =  _charge * r*r * cos(-_charge * t + phi0 + pi/2 * _charge);									// d^2y/domegadt
		dtrackToXyz(3,2) = -r*r * tanlambda;	// d^2z/domegadt

		dtrackToXyz(4,0) = 0;							// d^2x/dtanlambdadt
		dtrackToXyz(4,1) = 0;							// d^2y/dtanlambdadt
		dtrackToXyz(4,2) = r; 							// d^2z/dtanlambdadt

		SMatrix35  trackToXyzTranspose = ROOT::Math::Transpose(  trackToXyz );
		SMatrix35 dtrackToXyzTranspose = ROOT::Math::Transpose( dtrackToXyz );

		// xyz error matrix
		SMatrix3 derrtmp =
			  trackToXyzTranspose * _err * dtrackToXyz
			+dtrackToXyzTranspose * _err *  trackToXyz;


		if ( derrtmp(0,1)-derrtmp(1,0)>1e-3
				|| derrtmp(0,2)-derrtmp(2,0)>1e-3
				|| derrtmp(1,2)-derrtmp(2,1)>1e-3) {
			printf("error matrix not symmetric ****\n");
			printf("%f : %f : %f\n", derrtmp(0,0), derrtmp(0,1), derrtmp(0,2) );
			printf("%f : %f : %f\n", derrtmp(1,0), derrtmp(1,1), derrtmp(1,2) );
			printf("%f : %f : %f\n", derrtmp(2,0), derrtmp(2,1), derrtmp(2,2) );
		}

		derr(0,0) = derrtmp(0,0);
		derr(1,1) = derrtmp(1,1);
		derr(2,2) = derrtmp(2,2);
		derr(0,1) = derrtmp(0,1);
		derr(0,2) = derrtmp(0,2);
		derr(1,2) = derrtmp(1,2);
	}

	void Helix::GetPosErrDeriv2(double t, SVector3 &ddpos, SMatrixSym3 &dderr)const
	{
		const double pi = TMath::Pi();

		double d0 = _hel(id0);
		//double z0 = _hel(iz0);
		double phi0 = _hel(iph);
		double r = fabs(1./_hel(iom));
		double tanlambda = _hel(itd);

		double ddx = -r * cos(-_charge * t + phi0 + pi/2 * _charge);
		double ddy = -r * sin(-_charge * t + phi0 + pi/2 * _charge);
		double ddz =  0;

		ddpos(0) = ddx, ddpos(1) = ddy, ddpos(2) = ddz;

		// 5-to-3 error-matrix conversion
		SMatrix53 trackToXyz;
		trackToXyz(0,0) = -sin(phi0);			// dx/dd0
		trackToXyz(0,1) = cos(phi0);			// dy/dd0
		trackToXyz(0,2) = 0;							// dz/dd0

		trackToXyz(1,0) = 0;							// dx/dz0
		trackToXyz(1,1) = 0;							// dy/dz0
		trackToXyz(1,2) = 1;							// dz/dz0

		trackToXyz(2,0) = -d0 * cos(phi0) - r * sin(phi0 - pi/2 * _charge) - r * sin(-_charge * t + phi0 + pi/2 * _charge);	// dx/dphi0
		trackToXyz(2,1) = -d0 * sin(phi0) + r * cos(phi0 - pi/2 * _charge) + r * cos(-_charge * t + phi0 + pi/2 * _charge); // dy/dphi0
		trackToXyz(2,2) = 0;							// dz/dphi0
		
		trackToXyz(3,0) = -r*r * (cos(phi0 - pi/2 * _charge) + cos(-_charge * t + phi0 + pi/2 * _charge));									// dx/domega
		trackToXyz(3,1) = -r*r * (sin(phi0 - pi/2 * _charge) + sin(-_charge * t + phi0 + pi/2 * _charge));									// dy/domega
		trackToXyz(3,2) = -r*r * t * tanlambda;// dz/domega

		trackToXyz(4,0) = 0;							// dx/dtanlambda
		trackToXyz(4,1) = 0;							// dy/dtanlambda
		trackToXyz(4,2) = r*t;						// dz/dtanlambda

		// derivs
		SMatrix53 dtrackToXyz;
		dtrackToXyz(0,0) = 0;							// d^2x/dd0dt
		dtrackToXyz(0,1) = 0;							// dy/dd0dt
		dtrackToXyz(0,2) = 0;							// dz/dd0dt

		dtrackToXyz(1,0) = 0;							// d^2x/dz0dt
		dtrackToXyz(1,1) = 0;							// d^2y/dz0dt
		dtrackToXyz(1,2) = 0;							// d^2z/dz0dt

		dtrackToXyz(2,0) = _charge * r * cos(-_charge * t + phi0 + pi/2 * _charge);	// d^2x/dphi0dt
		dtrackToXyz(2,1) = _charge * r * sin(-_charge * t + phi0 + pi/2 * _charge);  // d^2y/dphi0dt
		dtrackToXyz(2,2) = 0;							// d^2z/dphi0dt
		
		dtrackToXyz(3,0) = -_charge * r*r * sin(-_charge * t + phi0 + pi/2 * _charge);									// d^2x/domegadt
		dtrackToXyz(3,1) =  _charge * r*r * cos(-_charge * t + phi0 + pi/2 * _charge);									// d^2y/domegadt
		dtrackToXyz(3,2) = -r*r * tanlambda;	// d^2z/domegadt

		dtrackToXyz(4,0) = 0;							// d^2x/dtanlambdadt
		dtrackToXyz(4,1) = 0;							// d^2y/dtanlambdadt
		dtrackToXyz(4,2) = r; 							// d^2z/dtanlambdadt

		// second derivs
		SMatrix53 ddtrackToXyz;
		ddtrackToXyz(0,0) = 0;							// d^2x/dd0dt^2
		ddtrackToXyz(0,1) = 0;							// dy/dd0dt^2
		ddtrackToXyz(0,2) = 0;							// dz/dd0dt^2

		ddtrackToXyz(1,0) = 0;							// d^2x/dz0dt^2
		ddtrackToXyz(1,1) = 0;							// d^2y/dz0dt^2
		ddtrackToXyz(1,2) = 0;							// d^2z/dz0dt^2

		ddtrackToXyz(2,0) =  r * sin(-_charge * t + phi0 + pi/2 * _charge);	// d^2x/dphi0dt^2
		ddtrackToXyz(2,1) = -r * cos(-_charge * t + phi0 + pi/2 * _charge);  // d^2y/dphi0dt^2
		ddtrackToXyz(2,2) = 0;							// d^2z/dphi0dt^2
		
		ddtrackToXyz(3,0) = r*r * cos(-_charge * t + phi0 + pi/2 * _charge);									// d^2x/domegadt^2
		ddtrackToXyz(3,1) = r*r * sin(-_charge * t + phi0 + pi/2 * _charge);									// d^2y/domegadt^2
		ddtrackToXyz(3,2) = 0;	// d^2z/domegadt^2

		ddtrackToXyz(4,0) = 0;							// d^2x/dtanlambdadt^2
		ddtrackToXyz(4,1) = 0;							// d^2y/dtanlambdadt^2
		ddtrackToXyz(4,2) = 0; 							// d^2z/dtanlambdadt^2

		SMatrix35 trackToXyzTranspose = ROOT::Math::Transpose( trackToXyz );
		SMatrix35 dtrackToXyzTranspose = ROOT::Math::Transpose( dtrackToXyz );
		SMatrix35 ddtrackToXyzTranspose = ROOT::Math::Transpose( ddtrackToXyz );

		// xyz error matrix
		SMatrix3 dderrtmp =
					trackToXyzTranspose*_err*ddtrackToXyz
					+ddtrackToXyzTranspose*_err*trackToXyz
					+ 2*ROOT::Math::SimilarityT(dtrackToXyz, _err);

		assert(dderrtmp(0,1)-dderrtmp(1,0)<1e-3);
		assert(dderrtmp(0,2)-dderrtmp(2,0)<1e-3);
		assert(dderrtmp(1,2)-dderrtmp(2,1)<1e-3);

		dderr(0,0) = dderrtmp(0,0);
		dderr(1,1) = dderrtmp(1,1);
		dderr(2,2) = dderrtmp(2,2);
		dderr(0,1) = dderrtmp(0,1);
		dderr(0,2) = dderrtmp(0,2);
		dderr(1,2) = dderrtmp(1,2);
	}

	void Helix::GetCenter(double &x, double &y)const
	{
		const double pi = TMath::Pi();

		double d0 = _hel(id0);
		double phi0 = _hel(iph);
		double r = fabs(1./_hel(iom));

		x = -d0 * sin(phi0) + r * cos(phi0 - pi/2 * _charge);
		y =  d0 * cos(phi0) + r * sin(phi0 - pi/2 * _charge);
	}

	// find z position in x/y direction
	void Helix::FindZCross(double x, double y, double &zi, double &zp)const
	{
		const double pi = TMath::Pi();

		double d0 = _hel(id0);
		double z0 = _hel(iz0);
		double phi0 = _hel(iph);
		double r = fabs(1./_hel(iom));
		double tanlambda = _hel(itd);

		double x0 = -d0 * sin(phi0);
		double y0 =  d0 * cos(phi0);

		// center of the xy circle
		double xcenter = x0 + r * cos(phi0 - pi/2. * _charge);
		double ycenter = y0 + r * sin(phi0 - pi/2. * _charge);

		// direction in the circle
		double diff2dx = x - xcenter;
		double diff2dy = y - ycenter;
		double phi = atan2(diff2dy, diff2dx);

		double phic = atan2(y0 - ycenter, x0 - xcenter); // initial direction
		zi = z0 + r * (phic - phi) * _charge * tanlambda;
		zp = 2. * pi * r * tanlambda;
	}


	TVector3 Helix::ClosePoint(const Helix &hel)const
	{

		// 2d closest points
		double xc=0,yc=0,xc2=0,yc2=0; // xy of closest points
		int nc;	// # closest points

		double r1 = fabs(1./_hel(iom));
		double r2 = fabs(1./hel._hel(iom));

		//fprintf(stderr,"_hel(iom)=%f,  hel._hel(iom)=%f\n",_hel(iom),hel._hel(iom));
		if(r2>r1) return hel.ClosePoint(*this);
		assert(r1>=r2);

		// r1>r2 only.

		double x1,y1,x2,y2;
		GetCenter(x1,y1);
		hel.GetCenter(x2,y2);

		double dist = sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1));
		double phi = atan2(y2-y1, x2-x1);

		if(dist >= r1+r2){
			nc = 1;
			xc = x1 + (r1 + (dist-r1-r2)/2) * cos(phi);
			yc = y1 + (r1 + (dist-r1-r2)/2) * sin(phi);
		}
		else if(dist > r1-r2){
			double theta = acos((dist*dist + r1*r1 - r2*r2)/(2*dist*r1)); // cosine formula
			nc = 2;
			xc = x1 + r1 * cos(phi+theta);
			yc = y1 + r1 * sin(phi+theta);
			xc2 = x1 + r1 * cos(phi-theta);
			yc2 = y1 + r1 * sin(phi-theta);
		}
		else {
			nc = 1;
			double phi = atan2(y2-y1, x2-x1);
			xc = x1 + (r1 - (r1-r2-dist)/2) * cos(phi);
			yc = y1 + (r1 - (r1-r2-dist)/2) * sin(phi);
		}

		// xy closest points obtained
		if(verbose){
			cout << "xy closest point: (" << xc << ", " << yc << ")";
			if(nc==2)
				cout << ", (" << xc2 << ", " << yc2 << ")";
			cout << endl;
		}

		// determine z
		// nearest angle in the track circle
		double zleast[2], zperiod[2];
		FindZCross(xc, yc, zleast[0], zperiod[0]);
		hel.FindZCross(xc, yc, zleast[1], zperiod[1]);
		double z[2] = {0,0};

		int imin=0,jmin=0;
		double difmin = 1e+300;

		const double zth = 3000.;
		for(int i=0;;i>=0?i++:i--){
			z[0] = zleast[0] + i * zperiod[0];

			if(fabs(z[0])>zth){
				if(i<0)break;
				i = -1;
				z[0] = zleast[0] + i * zperiod[0];
			}

			for(int j=0;;j>=0?j++:j--){
				z[1] = zleast[1] + j * zperiod[1];
				if(fabs(z[1])>zth){
					if(j<0)break;
					j = -1;
					z[1] = zleast[1] + j * zperiod[1];
				}

				if(fabs(z[1] - z[0])<difmin){
					difmin = fabs(z[1] - z[0]);
					imin = i; jmin = j;
				}
			}
		}

		double z11 = zleast[0] + imin * zperiod[0];
		double z12 = zleast[1] + jmin * zperiod[1];
		double z1 = (z11 + z12)/2;

		if(verbose)
			cout << "z11: " << z11 << ", z12: " << z12 << ", z1: " << z1 << endl;

		if(nc == 2){
			FindZCross(xc2, yc2, zleast[0], zperiod[0]);
			hel.FindZCross(xc2, yc2, zleast[1], zperiod[1]);
			z[0]=0;z[1]=0;
			difmin = 1e+300;

			for(int i=0;;i>=0?i++:i--){
				z[0] = zleast[0] + i * zperiod[0];

				if(fabs(z[0])>zth){
					if(i<0)break;
					i = -1;
					z[0] = zleast[0] + i * zperiod[0];
				}

				for(int j=0;;j>=0?j++:j--){
					z[1] = zleast[1] + j * zperiod[1];
					if(fabs(z[1])>zth){
						if(j<0)break;
						j = -1;
						z[1] = zleast[1] + j * zperiod[1];
					}

					if(fabs(z[1] - z[0])<difmin){
						difmin = fabs(z[1] - z[0]);
						imin = i; jmin = j;
					}
				}
			}
			double z21 = zleast[0] + imin * zperiod[0];
			double z22 = zleast[1] + jmin * zperiod[1];
			double z2 = (z21 + z22)/2;

			if(verbose)
				cout << "z21: " << z21 << ", z22: " << z22 << ", z2: " << z2 << endl;

//			if(fabs(z22-z21) < fabs(z12-z11)){
			TVector3 p1(xc, yc, z1);
			double ll1 = LogLikelihood(p1) + hel.LogLikelihood(p1);
			TVector3 p2(xc2, yc2, z2);
			double ll2 = LogLikelihood(p2) + hel.LogLikelihood(p2);
			if(ll2 > ll1){
				xc = xc2, yc = yc2, z1 = z2;
			}
		}

		return TVector3(xc, yc, z1);
	}

	// geometry operation ///////////////////////////////////////////////////////
	double GeometryHandler::PointFit(const vector<PointBase *> &points, const TVector3 &initial, Point * result)
	{
#if 1
		// three-dimensional minimization
		PointFitFunctor pf(points);
		ROOT::Minuit2::Minuit2Minimizer min(ROOT::Minuit2::kMigrad);
		ROOT::Math::Functor f(pf,3);
		min.SetFunction(f);
		min.SetMaxFunctionCalls(10000);
		min.SetMaxIterations(100);
		min.SetTolerance(1e-4);
		min.SetValidError(true);

		//min.SetVariable(0,"x",initial.x(), 1e-6);
		//min.SetVariable(1,"y",initial.y(), 1e-6);
		//min.SetVariable(2,"z",initial.z(), 1e-6);
		min.SetVariable(0,"x",initial.x(), 1e-4);
		min.SetVariable(1,"y",initial.y(), 1e-4);
		min.SetVariable(2,"z",initial.z(), 1e-4);
		min.SetPrintLevel(0);
		bool success = min.Minimize();
		if (!success && verbose) {
			printf("minuit status: %d\n", min.Status());
		}
		const double *xx = min.X();
		double maxll = -pf(xx);
#endif

		/*
		PointFitFunctor pf(points);
		PointFitDerivFunctor dpf(points);
		int iter;
		double fmin;
		double p[3];
		p[0]=initial.X(), p[1]=initial.Y(), p[2]=initial.Z();
		dfpmin(p, 3, 1e-3, &iter, &fmin, &pf, &dpf);
		double maxll = -fmin;
		*/

//		min.PrintResults();
		if(verbose){
			cout << "Mimimizer status: " << min.Status() << ", Edm = " << min.Edm() << ", # tracks = " << points.size() << ", Tolerance = " << min.Tolerance()  << endl;
			cout << "Minimum: (" << xx[0] << ", " << xx[1] << ", " << xx[2] << "), l= " << maxll << endl;
		}

		if(result){
			Point::SVector3 pos(xx[0], xx[1], xx[2]);
			Point::SMatrixSym3 err;
			for(int i=0;i<3;i++)for(int j=i;j<3;j++){
				err(i,j) = min.CovMatrix(i,j);
				if(verbose)
					cout << "CovMatrix(" << i << "," << j << "): " << min.CovMatrix(i,j) << endl;
			}

			if (verbose) {
				cout << scientific;
				cout << "minuit result - covariant matrix:" << endl;
				cout << err(0,0) << "  " << err(0,1) << "  " << err(0,2) << endl;
				cout << err(1,0) << "  " << err(1,1) << "  " << err(1,2) << endl;
				cout << err(2,0) << "  " << err(2,1) << "  " << err(2,2) << endl;
				cout << fixed;
			}
			result->SetPosErr(pos,err);
		}

		return maxll;
	}

	double GeometryHandler::HelixPointFit(const vector<Helix *> &helices, Point * result)
	{
		vector<PointBase *> points;
		points.resize(helices.size());
		copy(helices.begin(), helices.end(), points.begin());

		PointFitFunctor pf(points);
		double xyz[3];

		double maxll = -1e+300;
		TVector3 maxllv;

		for(unsigned int i=0;i<helices.size()-1;i++){
			for(unsigned int j=i+1; j<helices.size();j++){
				TVector3 v = helices[i]->ClosePoint(*helices[j]);
				xyz[0] = v.x(), xyz[1] = v.y(), xyz[2] = v.z();
				double ll = -pf(xyz);
				if(verbose)
					cout << "HelixPointFit: i = " << i << ", j = " << j << ", ll = " << ll << endl;
				if(maxll < ll){
					maxllv = v;
					maxll = ll;
				}
			}
		}

		double pfll = PointFit(points, maxllv, result);
		TVector3 pfpos = result->GetPos();

		if(verbose)
			cout << "HelixPointFit: ll = " << maxll << " -> " << pfll << ", pos = (" << maxllv.x() << " " << maxllv.y() << " " << maxllv.z() << ") -> ("
						<< pfpos.x() << " " << pfpos.y() << " " << pfpos.z() << ")" << endl;
		return pfll;
	}

}
