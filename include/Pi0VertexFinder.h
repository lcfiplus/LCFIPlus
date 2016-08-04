#include <vector>
#include <string>

//ROOT includes
#include "TMatrix.h" 
#include "TFile.h"  
#include "TH1.h"  
#include "TH2.h" 
#include "TLorentzVector.h" 
#include "TVector3.h"  
#include "TMVA/Reader.h"  

#include "lcfiplus.h"

class Pi0Finder{
public:
  Pi0Finder(string pi0pdfname, string pi0pdfname2d);
  ~Pi0Finder();

  void Initialize() { prior=0.5; }  //set prior 0.5
  double Get_Prob(double *var);
  double Get_g1Prob(double *var);
  double Get_g2Prob(double *var);
  double Get_PairProb(double *var);

private:
  TFile *fpdf,*fpdf2;
  TH2F *pdfs2D[2][2];
  TH1F *pdfs[2][19];

  double getValue(int valtype, int type, double value);
  double getValue2D(int valtype, int type, double valuex, double valuey);
  double *threshold;
  double prior;
};

class Pi0VertexFinder{
public:
  Pi0VertexFinder();
  Pi0VertexFinder(vector<string> weightfiles, vector<string> booknames, vector<double> opp, string pi0pdfname, string pi0pdfname2d);
  ~Pi0VertexFinder();

  void Initialize();
  void Make_GammaVector(const lcfiplus::Neutral* ntrl);
  void Reco_Pi0s(TVector3 vtx, TVector3 vtxdir);
  TLorentzVector Corr_VtxVect(int vtxtype, TLorentzVector vtx1vect, int *npart, TVector3 vtx, TVector3 vtxdir);
  int Get_nPi0();

private:
  void Make_Pi0Vector(TVector3 vtx, TVector3 vtxdir);
  double Get_GlobalMinimum();

  TMVA::Reader *vtxpi0;
  float var[20];    //for MVA variables

  int ng,counter,npi0;

  TMatrixD dm;
  vector< vector<double> > gv;
  vector< vector<double> > pi0vec;

  Pi0Finder *pi0ID;

  vector<string> _booknames;
  vector<double> _opp;
};
