#ifndef lcfiplus_algo_h
#define lcfiplus_algo_h 1

#include <math.h>
#include <stdio.h>
#include <vector>
#include <iostream>

using namespace std;

class TrackPocaXY;

namespace lcfiplus {

	const double R=0.61803399;
	const double C=(1.0-R);
	const double GOLD=1.618034;
	const double GLIMIT=100.0;
	const double TINY=1.e-25;
	const double CGOLD=0.3819660;
	//const double ZEPS = 1e-10;
	const double ZEPS = 1e-20;
	const int ITMAX=200;

#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
#define SHFT2(a,b,c) (a)=(b);(b)=(c);
#define SHFT3(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
#define SIGN(a,b) ((b)>=0.0 ? fabs(a) : -fabs(a))
#define MOV3(a,b,c, d,e,f) (a)=(d);(b)=(e);(c)=(f);

  template<class T>
  double golden(double ax, double bx, double cx, T* obj, double tol, double& xmin) {
    double f1,f2,x0,x1,x2,x3;
    x0=ax;
    x3=cx;
    if (fabs(cx-bx) > fabs(bx-ax)) {
      x1=bx;
      x2=bx+C*(cx-bx);
    } else {
      x2=bx;
      x1=bx-C*(bx-ax);
    }
    f1=(*obj)(&x1);
    f2=(*obj)(&x2);
    while (fabs(x3-x0) > tol*(fabs(x1)+fabs(x2))) {
      if (f2 < f1) {
        SHFT3(x0,x1,x2,R*x1+C*x3);
        SHFT2(f1,f2,(*obj)(&x2));
      } else {
        SHFT3(x3,x2,x1,R*x2+C*x0);
        SHFT2(f2,f1,(*obj)(&x1));
      }
    }
    if (f1 < f2) {
      xmin=x1;
      return f1;
    } else {
      xmin=x2;
      return f2;
    }
  }


  template<class T>
  double brent(double ax, double bx, double cx, T* obj, double tol, double& xmin) {
		int iter;
		double a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
		double e=0.0;

		a=(ax < cx ? ax : cx);
		b=(ax > cx ? ax : cx);
		x=w=v=bx;
		fw=fv=fx=(*obj)(&x);
		for (iter=1;iter<=ITMAX;iter++) {
			xm=0.5*(a+b);
			tol2=2.0*(tol1=tol*fabs(x)+ZEPS);
			if (fabs(x-xm) <= (tol2-0.5*(b-a))) {
				xmin=x;
				return fx;
			}
			if (fabs(e) > tol1) {
				r=(x-w)*(fx-fv);
				q=(x-v)*(fx-fw);
				p=(x-v)*q-(x-w)*r;
				q=2.0*(q-r);
				if (q > 0.0) p = -p;
				q=fabs(q);
				etemp=e;
				e=d;
				if (fabs(p) >= fabs(0.5*q*etemp) || p<= q*(a-x) || p >= q*(b-x))
					d=CGOLD*(e=(x >= xm ? a-x : b-x));
				else {
					d=p/q;
					u=x+d;
					if (u-a < tol2 || b-u < tol2)
						d=SIGN(tol1,xm-x);
				}
			} else {
				d=CGOLD*(e=(x >= xm ? a-x : b-x));
			}
			u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
			fu=(*obj)(&u);
			if (fu <= fx) {
				if (u>= x) a=x; else b=x;
				SHFT(v,w,x,u);
				SHFT(fv,fw,fx,fu);
			} else {
				if (u < x) a=u; else b=u;
				if (fu <= fw || w == x) {
					v=w;
					w=u;
					fv=fw;
					fw=fu;
				} else if (fu <= fv || v == x || v == w) {
					v=u;
					fv=fu;
				}
			}
		}
		fprintf(stderr,"too many iterations in brent()\n");
		xmin=x;
		return fx;
	}

  template<class T, class U>
  double dbrent(double ax, double bx, double cx, T* obj, U* dobj, double tol, double& xmin) {
		int iter,ok1,ok2;
		double a,b,d,d1,d2,du,dv,dw,dx,e=0.0;
		double fu,fv,fw,fx,olde,tol1,tol2,u,u1,u2,v,w,x,xm;

		a=(ax < cx ? ax : cx);
		b=(ax > cx ? ax : cx);
		x=w=v=bx;
		fw=fv=fx=(*obj)(&x);
		dw=dv=dx=(*dobj)(&x);
		for (iter=1;iter<=ITMAX;iter++) {
			xm=0.5*(a+b);
			tol1=tol*fabs(x)+ZEPS;
			tol2=2.0*tol1;
			if (fabs(x-xm) <= (tol2-0.5*(b-a))) {
				xmin=x;
				return fx;
			}
			if (fabs(e) > tol1) {
				d1=2.0*(b-a);
				d2=d1;
				if (dw != dx) d1=(w-x)*dx/(dx-dw);
				if (dv != dx) d2=(v-x)*dx/(dx-dv);
				u1=x+d1;
				u2=x+d2;
				ok1 = (a-u1)*(u1-b) > 0.0 && dx*d1 <= 0.0;
				ok2 = (a-u2)*(u2-b) > 0.0 && dx*d2 <= 0.0;
				olde=e;
				e=d;
				if (ok1 || ok2) {
					if (ok1 && ok2)
						d=(fabs(d1) < fabs(d2) ? d1 : d2);
					else if (ok1)
						d=d1;
					else
						d=d2;
					if (fabs(d) <= fabs(0.5*olde)) {
						u=x+d;
						if (u-a < tol2 || b-u < tol2)
							d=SIGN(tol1,xm-x);
					} else {
						d=0.5*(e=(dx >= 0.0 ? a-x : b-x));
					}
				} else {
					d=0.5*(e=(dx >= 0.0 ? a-x : b-x ));
				}
			} else {
				d=0.5*(e=(dx >= 0.0 ? a-x : b-x ));
			}
			if (fabs(d) >= tol1) {
				u=x+d;
				fu=(*obj)(&u);
			} else {
				u=x+SIGN(tol1,d);
				fu=(*obj)(&u);
				if (fu > fx) {
					xmin = x;
					return fx;
				}
			}
			du = (*dobj)(&u);
			if (fu <= fx) {
				if (u >= x) a=x; else b=x;
				MOV3(v,fv,dv, w,fw,dw);
				MOV3(w,fw,dw, x,fx,dx);
				MOV3(x,fx,dx, u,fu,du);
			} else {
				if (u < x) a=u; else b=u;
				if (fu <= fw || w == x) {
					MOV3(v,fv,dv, w,fw,dw);
					MOV3(w,fw,dw, u,fu,du);
				} else if (fu < fv || v == x || v == w) {
					MOV3(v,fv,dv, u,fu,du);
				}
			}
		}
		fprintf(stderr,"too many iterations in dbrent\n");
		return 0.0;
	}

	template<class T>
	void lnsrch(int n, double xold[], double fold, double g[], double p[], double x[],
			double *f, double stpmax, int *check, T* obj)
	{
		const double ALF=1.0e-4;
		const double TOLX=1.0e-7;

		static double maxarg1,maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ? (maxarg1) : (maxarg2))

		int i;
		double a,alam,alam2,alamin,b,disc,f2,fold2,rhs1,rhs2,slope,sum,temp,
					test,tmplam;

		*check=0;
		for (sum=0.0,i=0;i<n;i++) sum += p[i]*p[i];
		sum=sqrt(sum);
		if (sum > stpmax)
			for (i=0;i<n;i++) p[i] *= stpmax/sum;
		for (slope=0.0,i=0;i<n;i++)
			slope += g[i]*p[i];
		test=0.0;
		for (i=0;i<n;i++) {
			temp=fabs(p[i])/FMAX(fabs(xold[i]),1.0);
			if (temp > test) test=temp;
		}
		alamin=TOLX/test;
		alam=1.0;
		for (;;) {
			for (i=0;i<n;i++) x[i]=xold[i]+alam*p[i];
			*f=(*obj)(x);
			if (alam < alamin) {
				for (i=0;i<n;i++) x[i]=xold[i];
				*check=1;
				return;
			} else if (*f <= fold+ALF*alam*slope) return;
			else {
				if (alam == 1.0)
					tmplam = -slope/(2.0*(*f-fold-slope));
				else {
					rhs1 = *f-fold-alam*slope;
					rhs2=f2-fold2-alam2*slope;
					a=(rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
					b=(-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);
					if (a == 0.0) tmplam = -slope/(2.0*b);
					else {
						disc=b*b-3.0*a*slope;
						if (disc<0.0) fprintf(stderr,"Roundoff problem in lnsrch.\n");
						else tmplam=(-b+sqrt(disc))/(3.0*a);
					}
					if (tmplam>0.5*alam)
						tmplam=0.5*alam;
				}
			}
			alam2=alam;
			f2 = *f;
			fold2=fold;
			alam=FMAX(tmplam,0.1*alam);
		}
	}

	template<class T, class U>
  void dfpmin(double p[], int n, double gtol, int *iter, double *fret, T* obj, U* dobj) {
		bool verbose = false;

		const double EPS=3.0e-8;
		const double TOLX=4*EPS;
		const double STPMAX=100.0;

		static double maxarg1,maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ? (maxarg1) : (maxarg2))

		int check,i,its,j;
		double den,fac,fad,fae,fp,stpmax,sum=0.0,sumdg,sumxi,temp,test;
		double *dg, *g, *hdg, **hessin, *pnew, *xi;

		dg     = new double[n];
		g      = new double[n];
		hdg    = new double[n];
		hessin = new double*[n];
		for (i=0; i<n; ++i)
			hessin[i] = new double[n];
		pnew   = new double[n];
		xi     = new double[n];

		fp=(*obj)(p);

		(*dobj)(p,g);

		for (i=0;i<n;i++) {
			for (j=0;j<n;j++) hessin[i][j]=0.0;
			hessin[i][i]=1.0;
			xi[i] = -g[i];
			sum += p[i]*p[i];
		}
		stpmax=STPMAX*FMAX(sqrt(sum),(double)n);

		for (its=1;its<=ITMAX;its++) {
			if(verbose){
				cout << "dfpmin: iteration " << its << ", p = ";
				for(i=0; i<n;i++)cout << p[i] << " ";
				cout << ", xi = ";
				for(i=0; i<n;i++)cout << xi[i] << " ";
			}

			*iter=its;
			lnsrch(n,p,fp,g,xi,pnew,fret,stpmax,&check,obj);
			fp = *fret;
			for (i=0;i<n;i++) {
				xi[i]=pnew[i]-p[i];
				p[i]=pnew[i];
			}
			test=0.0;
			for (i=0;i<n;i++) {
				temp=fabs(xi[i])/FMAX(fabs(p[i]),1.0);
				if (temp > test) test=temp;
			}
			if(verbose) cout << "testtol = " << test;
			if (test < TOLX) {
				if(verbose) cout << endl;

				delete[] dg;
				delete[] g;
				delete[] hdg;
				for (i=0; i<n; ++i) {
					delete[] hessin[i];
				}
				delete[] hessin;
				delete[] pnew;
				delete[] xi;
				return;
			}
			for (i=0;i<n;i++) dg[i]=g[i];
			(*dobj)(p,g);
			test=0.0;
			den=FMAX(*fret,1.0);
			for (i=0;i<n;i++) {
				temp=fabs(g[i])*FMAX(fabs(p[i]),1.0)/den;
				if (temp > test) test=temp;
			}

			if(verbose) cout << " testgtol = " << test << endl;

			if (test < gtol) {
				delete[] dg;
				delete[] g;
				delete[] hdg;
				for (i=0; i<n; ++i) {
					delete[] hessin[i];
				}
				delete[] hessin;
				delete[] pnew;
				delete[] xi;
				return;
			}
			for (i=0;i<n;i++) dg[i]=g[i]-dg[i];
			for (i=0;i<n;i++) {
				hdg[i]=0.0;
				for (j=0;j<n;j++) hdg[i] += hessin[i][j]*dg[j];
			}
			fac=fae=sumdg=sumxi=0.0;
			for (i=0;i<n;i++) {
				fac += dg[i]*xi[i];
				fae += dg[i]*hdg[i];
				sumdg += (dg[i]*dg[i]);
				sumxi += (xi[i]*xi[i]);
			}
			if (fac > sqrt(EPS*sumdg*sumxi)) {
				fac=1.0/fac;
				fad=1.0/fae;
				for (i=0;i<n;i++) dg[i]=fac*xi[i]-fad*hdg[i];
				for (i=0;i<n;i++) {
					for (j=i;j<n;j++) {
						hessin[i][j] += fac*xi[i]*xi[j]
							-fad*hdg[i]*hdg[j]+fae*dg[i]*dg[j];
						hessin[j][i]=hessin[i][j];
					}
				}
			}
			for (i=0;i<n;i++) {
				xi[i]=0.0;
				for (j=0;j<n;j++) xi[i] -= hessin[i][j]*g[j];
			}
		}
		fprintf(stderr,"too many iterations in dfpmin");
		delete[] dg;
		delete[] g;
		delete[] hdg;
		for (i=0; i<n; ++i) {
			delete[] hessin[i];
		}
		delete[] hessin;
		delete[] pnew;
		delete[] xi;
	}

}

#endif
