#include "algo.h"
#include <math.h>

#include "flavtag.h"

//#define R 0.61803399
//#define C (1.0-R)
namespace flavtag {
//#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
//#define SHFT2(a,b,c) (a)=(b);(b)=(c);
//#define SHFT3(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
//const double GOLD = 1.618034;
//const double GLIMIT = 100.0;
//const double TINY = 1.e-25;
//const double CGOLD = 0.3819660:



  void powell(vector<float> p, int n, float ftol, int *iter, float *fret,
      float (*func)(vector<float>))
  {
    // construct set of unit vectors for initial step directions
    vector< vector<float> > xi;
    for (unsigned int i=0; i<p.size(); ++i) {
      vector<float> x(p.size());
      x[i] = 1;
      xi.push_back(x);
    }

    void linmin(vector<float> p, vector<float> xi, int n, float *fret, float (*func)(vector<float>));
    int ibig;
    float del,fp,fptt,t;
    vector<float> pt,ptt,xit;
    pt.resize(n);
    ptt.resize(n);
    xit.resize(n);

    *fret=(*func)(p);
    for (int j=0;j<n;++j) pt[j]=p[j];
    for (*iter=1;;++(*iter)) {
      fp=(*fret);
      ibig=0;
      del=0.0;
      for (int i=0;i<n;++i) {
        for (int j=0;j<n;++j) xit[j]=xi[j][i];
        fptt=(*fret);
        linmin(p,xit,n,fret,func);
        if (fptt-(*fret) > del) {
          del=fptt-(*fret);
          ibig=i;
        }
      }
      if (2.0*(fp-(*fret)) <= ftol*(fabs(fp)+fabs(*fret))+TINY) {
        return;
      }
      if (*iter == ITMAX) {
        fprintf(stderr,"powell: exceeding maximum iterations\n");
        exit(1);
      }
      for (int j=0;j<n;++j) {
        ptt[j]=2.0*p[j]-pt[j];
        xit[j]=p[j]-pt[j];
        pt[j]=p[j];
      }
      fptt=(*func)(ptt);
      if (fptt < fp) {
        t=2.0*(fp-2.0*(*fret)+fptt)*pow(fp-(*fret)-del,2)-del*pow(fp-fptt,2);
        if (t < 0.0) {
          linmin(p,xit,n,fret,func);
          for (int j=0;j<n;++j) {
            xi[j][ibig]=xi[j][n];
            xi[j][n]=xit[j];
          }
        }
      }
    }
  }

#define TOL 2e-4
  int ncom;
  float (*nrfunc)(vector<float>);
  vector<float>pcom,xicom;

  void linmin(vector<float> p, vector<float> xi, int n, float *fret, float (*func)(vector<float>)) {
    float brent(float ax, float bx, float cx, float (*f)(float), float tol, float *xmin);
    float f1dim(float x);
    void mnbrak(float *ax, float *bx, float *cx, float *fa, float *fb,
        float *fc, float (*func)(float));
    int j;
    float xx,xmin,fx,fb,fa,bx,ax;

    ncom=n;
    pcom.resize(n);
    xicom.resize(n);
    nrfunc=func;
    for (j=0;j<n;++j) {
      pcom[j]=p[j];
      xicom[j]=xi[j];
    }
    ax=0.0;
    xx=1.0;
    mnbrak(&ax,&xx,&bx,&fa,&fx,&fb,f1dim);
    *fret=brent(ax,xx,bx,f1dim,TOL,&xmin);
    for (j=0;j<n;++j) {
      xi[j] *= xmin;
      p[j] += xi[j];
    }
  }

  float f1dim(float x) {
    int j;
    float f;
    vector<float> xt;
    xt.resize(ncom);
    for (j=0;j<ncom;++j)
      xt[j]=pcom[j]+x*xicom[j];
    f=(*nrfunc)(xt);
    return f;
  }

  void mnbrak(float *ax, float *bx, float *cx, float *fa, float *fb, float *fc,
      float (*func)(float)){
    float ulim,u,r,q,fu,dum;

    *fa=(*func)(*ax);
    *fb=(*func)(*bx);
    if (*fb > *fa) {
      SHFT(dum,*ax,*bx,dum);
      SHFT(dum,*fb,*fa,dum);
    }
    *cx=(*bx)+GOLD*(*bx-*ax);
    *fc=(*func)(*cx);
    while (*fb > *fc) {
      r=(*bx-*ax)*(*fb-*fc);
      q=(*bx-*cx)*(*fb-*fa);
      u=(*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/
        (2.0*SIGN(std::max(fabs(q-r),TINY),q-r));
      ulim=(*bx)+GLIMIT*(*cx-*bx);
      if ((*bx-u)*(u-*cx) > 0.0) {
        fu=(*func)(u);
        if (fu < *fc) {
          *ax=(*bx);
          *bx=u;
          *fa=(*fb);
          *fb=fu;
          return;
        } else if (fu > *fb) {
          *cx=u;
          *fc=fu;
          return;
        }
        u=(*cx)+GOLD*(*cx-*bx);
        fu=(*func)(u);
      } else if ((*cx-u)*(u-ulim) > 0.0) {
        fu=(*func)(u);
        if (fu < *fc) {
          SHFT(*bx,*cx,u,*cx+GOLD*(*cx-*bx));
          SHFT(*fb,*fc,fu,(*func)(u));
        }
      } else if ((u-ulim)*(ulim-*cx) >= 0.0) {
        u=ulim;
        fu=(*func)(u);
      } else {
        u=(*cx)+GOLD*(*cx-*bx);
        fu=(*func)(u);
      }
      SHFT(*ax,*bx,*cx,u);
      SHFT(*fa,*fb,*fc,fu);
    }
  }

  float brent(float ax, float bx, float cx, float (*f)(float), float tol, float *xmin) {
    int iter;
    float a,b,d=0.0,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
    float e=0.0;

    a=(ax < cx ? ax : cx);
    b=(ax > cx ? ax : cx);
    x=w=v=bx;
    fw=fv=fx=(*f)(x);
    for (iter=1;iter<=ITMAX;++iter) {
      xm=0.5*(a+b);
      tol2=2.0*(tol1=tol*fabs(x)+ZEPS);
      if (fabs(x-xm) <= (tol2-0.5*(b-a))) {
        *xmin=x;
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
        if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
          d=CGOLD*(e=(x >= xm ? a-x : b-x));
        else {
          d=p/q;
          u=x+d;
          if (u-a < tol2 || b-u < tol2)
            d=SIGN(tol1,xm-x);
        }
      }else {
        d=CGOLD*(e=(x >= xm ? a-x : b-x));
      }
      u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
      fu=(*f)(u);
      if (fu <=fx) {
        if (u >=x) a=x; else b=x;
        SHFT(v,w,x,u);
        SHFT(fv,fw,fx,fu);
      } else {
        if (u< x) a=u; else b=u;
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
    fprintf(stderr,"too many iterations in brent");
    exit(1);
    *xmin=x;
    return fx;
  }
}
