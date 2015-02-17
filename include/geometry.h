// geometry.h

#ifndef geometry_h
#define geometry_h 1

#include "lcfiplus.h"
#include "TVector3.h"

#include "Math/SVector.h"
#include "Math/SMatrix.h"

namespace lcfiplus{

	class PointBase;
	class Point;
	class Helix;
	class VertexLine;

	class PointBase // pure virtual point-base class
	{
	public:
		virtual double LogLikelihood(const TVector3 &p)const  = 0;
		virtual void LogLikelihoodDeriv(const TVector3 &p,double* output)const  = 0;
		virtual ~PointBase(){}

	protected:
		PointBase(){}
	};

	class Point : public PointBase // real 3D-point with error
	{
	public:
		typedef ROOT::Math::SVector<double, 3> SVector3;
		typedef ROOT::Math::SMatrix<double, 3,3,ROOT::Math::MatRepSym<double,3> > SMatrixSym3;

		double LogLikelihood(const TVector3 &p)const;
		void LogLikelihoodDeriv(const TVector3 &p, double* output)const;

		Point(){}
		Point(const SVector3 &pos, const SMatrixSym3 &err){_pos = pos; _err = err;}
		Point(const Point &ref){_pos = ref._pos; _err = ref._err;}
		Point(const Vertex *vtx);
		~Point(){}

		void SetPosErr(const SVector3 &pos, const SMatrixSym3 &err){_pos = pos; _err = err;}
		double GetErr(int i, int j)const{return _err(i,j);}
		TVector3 GetPos()const{return TVector3(_pos(0), _pos(1), _pos(2));}
	private:
		SVector3 _pos;
		SMatrixSym3 _err;
	};

	class Helix : public PointBase // parametrized point for helix
	{
	public:
   	enum par { id0=0, iz0, iph, iom, itd, parN };

		typedef ROOT::Math::SVector<double, 5> SVector5;
		typedef ROOT::Math::SVector<double, 3> SVector3;
		typedef ROOT::Math::SMatrix<double, 5,5,ROOT::Math::MatRepSym<double,5> > SMatrixSym5;
		typedef ROOT::Math::SMatrix<double, 5,3> SMatrix53; // used for helix-xyz conversion
		typedef ROOT::Math::SMatrix<double, 3,5> SMatrix35; // used for helix-xyz conversion
		typedef ROOT::Math::SMatrix<double, 3,3,ROOT::Math::MatRepSym<double,3> > SMatrixSym3; // used for xyz error
		typedef ROOT::Math::SMatrix<double, 3,3> SMatrix3;

		class HelixLineDistance2Functor{
		public:
			HelixLineDistance2Functor(const Helix *hel, const VertexLine *line): _hel(hel), _line(line){}
			double operator() (const double *t);

		private:
			const Helix * _hel;
			const VertexLine * _line;
		};

		class HelixLineDistance2DerivFunctor{
		public:
			HelixLineDistance2DerivFunctor(const Helix *hel, const VertexLine *line): _hel(hel), _line(line){}
			void operator() (const double *t, double *output);

		private:
			const Helix * _hel;
			const VertexLine * _line;
		};

		class VarianceFunctor{
		public:
			VarianceFunctor(const Helix *hel, const TVector3 &p) : _hel(hel), _p(p){}
			double operator() (const double *t){
				return _hel->Variance(_p,*t);
			}
		private:
			const Helix *_hel;
			TVector3 _p;
		};

		class VarianceDerivFunctor{
		public:
			VarianceDerivFunctor(const Helix *hel, const TVector3 &p) : _hel(hel), _p(p){}
			double operator() (const double *t){
				return _hel->VarianceDeriv(_p,*t);
			}
		private:
			const Helix *_hel;
			TVector3 _p;
		};

		virtual double LogLikelihood(const TVector3 &p)const {// likelihood with t-minimization
			double tmin; return LogLikelihood(p, tmin);
		}
		double LogLikelihood(const TVector3 &p, double &tmin)const;// full version
		void LogLikelihoodDeriv(const TVector3 &p, double* output)const;// compute space partial derivatives
		double Variance(const TVector3 &p, double t)const;		// t-fixed version, internally used
		double VarianceDeriv(const TVector3 &p, double t)const;		// t-fixed version, internally used
		double VarianceDeriv2(const TVector3 &p, double t)const;		// t-fixed version, internally used
		TVector3 GetPos(double t)const;
		TVector3 GetPosDerivT(double t)const;
		void GetPosErr(double t, SVector3 &pos, SMatrixSym3 &err)const;
		void GetPosErr(double t, SVector3 &pos, SMatrixSym3 &err, SMatrix53& trackToXyz)const;
		void GetPosErrDeriv(double t, SVector3 &pos, SMatrixSym3 &err)const;
		void GetPosErrDeriv2(double t, SVector3 &pos, SMatrixSym3 &err)const;

		double LongitudinalDeviation(const Vertex *ip, const Vertex *sec);

		Helix(){}
		Helix(const SVector5 &hel, const SMatrixSym5 &err, int charge){_hel = hel, _err = err;_charge = charge;}
		Helix(const Track *trk);
		Helix(const Helix &ref){_hel = ref._hel; _err = ref._err;_charge = ref._charge;}
		~Helix(){}

		void GetCenter(double &x, double &y)const;
		void FindZCross(double x, double y, double &zi, double &zp)const;
		// obtain closest points in x-y plane and choose nearest z position - not the TRUE closest point
		TVector3 ClosePoint(const Helix &hel)const;
		TVector3 ClosePoint(const VertexLine &line, double *distance = 0)const;

	private:
		SVector5 _hel;
		SMatrixSym5 _err;
		int _charge;
	};

	class VertexLine : public PointBase // line with error for IP-vertex line
	{
		typedef ROOT::Math::SVector<double, 3> SVector3;
		typedef ROOT::Math::SMatrix<double, 3,3,ROOT::Math::MatRepSym<double,3> > SMatrixSym3; // used for xyz error

		friend TVector3 Helix::ClosePoint(const VertexLine &line, double *distance)const;
		friend class Helix::HelixLineDistance2Functor;
		friend class Helix::HelixLineDistance2DerivFunctor;

		class VarianceFunctor{
		public:
			VarianceFunctor(const VertexLine &line, const TVector3 &p) : _line(line), _p(p){}
			double operator() (const double *t){
				return _line.Variance(_p,*t);
			}
		private:
			const VertexLine &_line;
			TVector3 _p;
		};

	public:
		virtual double LogLikelihood(const TVector3 &p)const {// likelihood with t-minimization
			double tmin; return LogLikelihood(p, tmin);
		}
		double LogLikelihood(const TVector3 &p, double &tmin)const;
		void LogLikelihoodDeriv(const TVector3 &p, double* output)const;
		double Variance(const TVector3 &p, double t)const;		// t-fixed version, internally used

		VertexLine(){}
		VertexLine(const Vertex *ip, const Vertex *secvtx){_ip = ip; _vertex = secvtx; _origin = _vertex->getPos(); _unit = (_origin - _ip->getPos()).Unit();}
		VertexLine(const TVector3 &origin, const TVector3 &dir){_origin = origin; _unit = dir.Unit();_ip = 0; _vertex = 0;}
		~VertexLine(){}

		void Set(const TVector3 &origin, const TVector3 &dir){_origin = origin; _unit = dir.Unit();}

	private:
		// (x,y,z) = origin + t * unit
		TVector3 _origin;
		TVector3 _unit;

		const Vertex *_ip;
		const Vertex *_vertex;
/*
		double _dispersionNear;
		double _dispersionFar;
*/
	};

	class GeometryHandler
	{
	public:
		static GeometryHandler * Instance();

		class PointFitFunctor{
		public:
			PointFitFunctor(const vector<PointBase *> &points) : _points(points){}
			double operator() (const double *xx){
				TVector3 p(xx[0], xx[1], xx[2]);

				double ll = 0.;
				for(unsigned int i=0;i<_points.size();i++){
					ll += _points[i]->LogLikelihood(p);
				}

				return -ll;
			}
		private:
			const vector<PointBase *> _points;
		};

		class PointFitDerivFunctor{
		public:
			PointFitDerivFunctor(const vector<PointBase *> &points) : _points(points){}
			void operator() (const double *xx, double* output){
				output[0] = 0;
				output[1] = 0;
				output[2] = 0;

				TVector3 p(xx[0], xx[1], xx[2]);

				double tmp[3];

				for(unsigned int i=0;i<_points.size();i++){
					_points[i]->LogLikelihoodDeriv(p,tmp);
					output[0] += tmp[0];
					output[1] += tmp[1];
					output[2] += tmp[2];
				}
			}
		private:
			const vector<PointBase *> _points;
		};

		// obtain cross section of points with errors
		double PointFit(const vector<PointBase *> &points, const TVector3 &initial, Point * result = 0);
		// initialization + PointFit()
		double HelixPointFit(const vector<Helix *> &helices, Point * result = 0);

	private:
		static GeometryHandler * _theInstance;
		GeometryHandler();
		~GeometryHandler();
	};


}

#endif
