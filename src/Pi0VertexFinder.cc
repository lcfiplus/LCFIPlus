#include "Pi0VertexFinder.h"

using namespace std;

/*string itos(int i)  
{
  stringstream s;
  s << i;
  return s.str();
  }*/

//Pi0Finder Functions
//constructor
Pi0Finder::Pi0Finder(string pi0pdfname, string pi0pdfname2d){
  threshold= new double[3];

  //open pdf files
  string ffstr1= pi0pdfname;   //"pdf/pdf_for_pi0ID_ok.root";
  fpdf=new TFile(ffstr1.c_str());

  ffstr1= pi0pdfname2d;    //"pdf/pdf_for_pi0ID_2D_ok.root";
  fpdf2=new TFile(ffstr1.c_str());
  
  string hname,hname2;
  for(int i=0;i<2;i++){
    stringstream s;
    s << i+1;

    //g1e
    hname="hg1e" + s.str();
    hname2="hg1e" + s.str() +"_2";
    //cout <<hname2.c_str() << endl;
    pdfs[i][0]=(TH1F*)fpdf->Get(hname.c_str())->Clone(hname2.c_str());
    
    //g2e
    hname="hg2e" + s.str();
    hname2="hg2e" + s.str() +"_2";
    //cout <<hname2.c_str() << endl;
    pdfs[i][1]=(TH1F*)fpdf->Get(hname.c_str())->Clone(hname2.c_str());

    //theta
    hname="htheta" + s.str();
    hname2="htheta" + s.str() +"_2";
    //cout <<hname2.c_str() << endl;
    pdfs[i][2]=(TH1F*)fpdf->Get(hname.c_str())->Clone(hname2.c_str());

    //mass
    hname="hmpi0" + s.str();
    hname2="hmpi0" + s.str() +"_2";
    //cout <<hname2.c_str() << endl;
    pdfs[i][3]=(TH1F*)fpdf->Get(hname.c_str())->Clone(hname2.c_str());

    //g1hadem
    hname="hg1hadem" + s.str();
    hname2="hg1hadem" + s.str() +"_2";
    //cout <<hname2.c_str() << endl;
    pdfs[i][4]=(TH1F*)fpdf->Get(hname.c_str())->Clone(hname2.c_str());

    //g1chi2
    hname="hg1chi2" + s.str();
    hname2="hg1chi2" + s.str() +"_2";
    //cout <<hname2.c_str() << endl;
    pdfs[i][5]=(TH1F*)fpdf->Get(hname.c_str())->Clone(hname2.c_str());
        
    //g1showerMax
    hname="hg1showerMax" + s.str();
    hname2="hg1showerMax" + s.str() +"_2";
    //cout <<hname2.c_str() << endl;
    pdfs[i][6]=(TH1F*)fpdf->Get(hname.c_str())->Clone(hname2.c_str());
        
    //g1absLength
    hname="hg1absLength" + s.str();
    hname2="hg1absLength" + s.str() +"_2";
    //cout <<hname2.c_str() << endl;
    pdfs[i][7]=(TH1F*)fpdf->Get(hname.c_str())->Clone(hname2.c_str());
        
    //g1xl20
    hname="hg1xl20" + s.str();
    hname2="hg1xl20" + s.str() +"_2";
    //cout <<hname2.c_str() << endl;
    pdfs[i][8]=(TH1F*)fpdf->Get(hname.c_str())->Clone(hname2.c_str());
        
    //g2hadem
    hname="hg2hadem" + s.str();
    hname2="hg2hadem" + s.str() +"_2";
    //cout <<hname2.c_str() << endl;
    pdfs[i][9]=(TH1F*)fpdf->Get(hname.c_str())->Clone(hname2.c_str());
        
    //g2chi2
    hname="hg2chi2" + s.str();
    hname2="hg2chi2" + s.str() +"_2";
    //cout <<hname2.c_str() << endl;
    pdfs[i][10]=(TH1F*)fpdf->Get(hname.c_str())->Clone(hname2.c_str());
        
    //g2showerMax
    hname="hg2showerMax" + s.str();
    hname2="hg2showerMax" + s.str() +"_2";
    //cout <<hname2.c_str() << endl;
    pdfs[i][11]=(TH1F*)fpdf->Get(hname.c_str())->Clone(hname2.c_str());
        
    //g2absLength
    hname="hg2absLength" + s.str();
    hname2="hg2absLength" + s.str() +"_2";
    //cout <<hname2.c_str() << endl;
    pdfs[i][12]=(TH1F*)fpdf->Get(hname.c_str())->Clone(hname2.c_str());
    
    //g2xl20
    hname="hg2xl20" + s.str();
    hname2="hg2xl20" + s.str() +"_2";
    //cout <<hname2.c_str() << endl;
    pdfs[i][13]=(TH1F*)fpdf->Get(hname.c_str())->Clone(hname2.c_str());
        
    //g2rate
    hname="hg2rate" + s.str();
    hname2="hg2rate" + s.str() +"_2";
    //cout <<hname2.c_str() << endl;
    pdfs[i][14]=(TH1F*)fpdf->Get(hname.c_str())->Clone(hname2.c_str());
        
    //g2rate2
    hname="hg2rate2" + s.str();
    hname2="hg2rate2" + s.str() +"_2";
    //cout <<hname2.c_str() << endl;
    pdfs[i][15]=(TH1F*)fpdf->Get(hname.c_str())->Clone(hname2.c_str());
        
    //g1rate
    hname="hg1rate" + s.str();
    hname2="hg1rate" + s.str() +"_2";
    //cout <<hname2.c_str() << endl;
    pdfs[i][16]=(TH1F*)fpdf->Get(hname.c_str())->Clone(hname2.c_str());
        
    //g2rate2
    hname="hconeenergy" + s.str();
    hname2="hconeenergy" + s.str() +"_2";
    //cout <<hname2.c_str() << endl;
    pdfs[i][17]=(TH1F*)fpdf->Get(hname.c_str())->Clone(hname2.c_str());
        
    //g1rate
    hname="hngamma" + s.str();
    hname2="hngamma" + s.str() +"_2";
    //cout <<hname2.c_str() << endl;
    pdfs[i][18]=(TH1F*)fpdf->Get(hname.c_str())->Clone(hname2.c_str());
  }

  //get 2D likelihood
  for(int i=0;i<2;i++){
    stringstream s;
    s << i+1;

    //get 2Dhistograms
    hname="hetheta" + s.str();
    hname2="hetheta" + s.str() +"_2";
    //cout <<hname2.c_str() << endl;
    pdfs2D[i][0]=(TH2F*)fpdf2->Get(hname.c_str())->Clone(hname2.c_str());
    
    hname="he1e2" + s.str();
    hname2="he1e2" + s.str() +"_2";
    //cout <<hname2.c_str() << endl;
    pdfs2D[i][1]=(TH2F*)fpdf2->Get(hname.c_str())->Clone(hname2.c_str());
  }
  
  //normalize histograms
  double weight=1.0;
  for(int j=0;j<2;j++){
    for(int i=0;i<19;i++){
      //normalize histograms
      weight=pdfs[j][i]->Integral(0,pdfs[j][i]->GetNbinsX()+1,"");
      if(weight!=0.0) pdfs[j][i]->Scale(1.0/weight);

      if(i<2){
	weight=pdfs2D[j][i]->Integral(0,pdfs2D[j][i]->GetNbinsX()+1,0,pdfs2D[j][i]->GetNbinsY()+1,"");
	if(weight!=0.0) pdfs2D[j][i]->Scale(1.0/weight);
      }
    }
  }

  //set threshold
  threshold[0]=TMath::Exp(-3.0);  //g1 threshold
  threshold[1]=TMath::Exp(-4.5);  //g2 threshold
  threshold[2]=TMath::Exp(-0.3);  //pair threshold   -0.25 for old

  return;
}

//destructor
Pi0Finder::~Pi0Finder(){

  //fpdf2->Close();
  //fpdf->Close();

  delete[] threshold;
  return;
}

double Pi0Finder::Get_PairProb(double *var){
  //bayesian approach to get the probability whether these 2 photons are from pi0 

  //get probability
  //set prior
  double okval[2]={0.0,0.0};
  double total=0.0;
  double posterior=0.0;
  int tmpj=0;

  if(TMath::Cos(var[2])<0.75) return -1.0;
  
  for(int j=0;j<19;j++){   //variables
    //j==0 and j==1 are using 2D likelihood
    tmpj=j;
    if(j<2){
      for(int i=0;i<2;i++){
	//get likelihood
	if(j==0) okval[i]=getValue2D(i,tmpj,var[2],var[0]+var[1]);   //likelihood 
	if(j==1) okval[i]=getValue2D(i,tmpj,var[2],var[1]);   //likelihood 
      }
    }else if(j==3 || j==14 || j==15 || j==16 || j==17 || j==18){
      
      if(var[18]==0 && j>=17) continue;

      for(int i=0;i<2;i++){
	//get likelihood
	okval[i]=getValue(i,tmpj,var[tmpj]);   //likelihood
      }
    }else continue;

    //cal. probability
    total=okval[0]*prior+okval[1]*(1.0-prior);
    if(total==0.0) total=1.0e-100;

    //cal. posterior!
    posterior=TMath::Log(okval[0])+TMath::Log(prior)-TMath::Log(total);

    //cout << "check: " << tmpj << " " << var[0] << " " << var[1] << " " << var[tmpj] << " " << okval[0] << " " << okval[1] << " " << TMath::Exp(posterior) << " " << prior << endl; 

    //bayesian updating
    prior=TMath::Exp(posterior);
  }

  //threshold check
  double ret=TMath::Exp(posterior);
  if(posterior<TMath::Log(threshold[2])) ret=-1.0;

  return ret;
}

double Pi0Finder::Get_g1Prob(double *var){
  //bayesian approach to get the probability whether g1( g1 means higher energy neutral) is gamma
 
  //get probability
  //set prior
  double okval[2]={0.0,0.0};
  double total=0.0;
  double posterior=0.0;
  for(int j=0;j<9;j++){   //variables
    if(!(j==0 || j==4 || j==5 || j==6 || j==7 || j==8)) continue;  
    for(int i=0;i<2;i++){
      //get likelihood
      okval[i]=getValue(i,j,var[j]);   //likelihood using g1 p.d.f
    }

    //cal. probability
    total=okval[0]*prior+okval[1]*(1.0-prior);
    if(total==0.0) total=1.0e-100;

    //cal. posterior!
    posterior=TMath::Log(okval[0])+TMath::Log(prior)-TMath::Log(total);

    //cout << "check: " <<  j << " " << var[j] << " " << okval[0] << " " << okval[1] << " " << TMath::Exp(posterior) << " " << prior << endl; 

    //bayesian updating
    prior=TMath::Exp(posterior);
  }

  //threshold check
  double ret=TMath::Exp(posterior);
  if(posterior<TMath::Log(threshold[0])) ret=-1.0;

  return ret;
}

double Pi0Finder::Get_g2Prob(double *var){
  //baysian approach to get the probability whether g2( g2 means lower energy neutral) is gamma
 
  //get probability
  //set prior
  double okval[2]={0.0,0.0};
  double total=0.0;
  double posterior=0.0;
  for(int j=0;j<9;j++){   //variables
    if(!(j==1 || j==4 || j==5 || j==6 || j==7 || j==8)) continue;  
    for(int i=0;i<2;i++){
      //get likelihood
      if(j==1) okval[i]=getValue(i,j,var[j]);   //likelihood using g2 p.d.f
      if(j!=1) okval[i]=getValue(i,j+5,var[j]);   //likelihood using g2 p.d.f
    }

    //cal. probability
    total=okval[0]*prior+okval[1]*(1.0-prior);
    if(total==0.0) total=1.0e-100;

    //cal. posterior!
    posterior=TMath::Log(okval[0])+TMath::Log(prior)-TMath::Log(total);

    //cout << "check: " << j << " " << okval[0] << " " << okval[1] << " " << TMath::Exp(posterior) << " " << prior << endl; 

    //bayesian updating
    prior=TMath::Exp(posterior);
  }

  //threshold check
  double ret=TMath::Exp(posterior);
  if(posterior<TMath::Log(threshold[1])) ret=-1.0;

  return ret;
}

double Pi0Finder::Get_Prob(double *var){
  double prob1=0.0,prob2=0.0;

  //check g1 is photon
  prior=0.5;  //initialize
  double prob=Get_g1Prob(var);
  if(prob<0.0) return prob;
  prob1=prob;
  
  //check g2 is photon
  prior=0.5;
  double var2[9];
  var2[1]=var[1];
  var2[4]=var[9];
  var2[5]=var[10];
  var2[6]=var[11];
  var2[7]=var[12];
  var2[8]=var[13];
  prob=Get_g2Prob(var2);
  if(prob<0.0) return prob;
  prob2=prob;

  //check photon pair is from pi0
  prior=prob1*prob2;
  prob=Get_PairProb(var);
  return prob;
}

double Pi0Finder::getValue(int valtype, int type, double value){
  /*double minval=pdfs[valtype][type]->GetBinLowEdge(1);
  double interval=pdfs[valtype][type]->GetBinWidth(1);
  int nbins=pdfs[valtype][type]->GetNbinsX();
  
  double val=1.0e-30;
  int bin=0;
  if(type!=4 && type!=9 && type!=14 && type!=15 && type!=16){
    if(value<minval){
      bin=0;
    }else if(value>=minval+nbins*interval){
      bin=(int)nbins+1;
    }else{
      for(int i=0;i<(int)nbins;i++){
	if(value>=minval+i*interval && value<minval+(i+1)*interval){
	  bin=i+1;
	}
      }
    }
  }else{
    if(value<=minval){
      bin=0;
    }else if(value>minval+nbins*interval){
      bin=(int)nbins+1;
    }else{
      for(int i=0;i<(int)nbins;i++){
	if(value>minval+i*interval && value<=minval+(i+1)*interval){
	  bin=i+1;
	}
      }
    }
  }
  val=pdfs[valtype][type]->GetBinContent(bin);*/

  double val=1.0e-30;
  int bin=0;

  bin=pdfs[valtype][type]->GetXaxis()->FindBin(value);
  
  // //get probability
  val=pdfs[valtype][type]->GetBinContent(bin);
  
  return val;
}

double Pi0Finder::getValue2D(int valtype, int type, double valuex, double valuey){
  /*double minvalx=pdfs[valtype][type]->GetXaxis()->GetBinLowEdge(1);
  double intervalx=pdfs[valtype][type]->GetXaxis()->GetBinWidth(1);
  int nbinsx=pdfs[valtype][type]->GetNbinsX();
  double minvaly=pdfs[valtype][type]->GetYaxis()->GetBinLowEdge(1);
  double intervaly=pdfs[valtype][type]->GetYaxis()->GetBinWidth(1);
  int nbinsy=pdfs[valtype][type]->GetNbinsY();

  double val=1.0e-30;
  int binx=0,biny=0;

  //get xbin
  if(valuex<minvalx){
    binx=0;
  }else if(valuex>=minvalx+nbinsx*intervalx){
    binx=(int)nbinsx+1;
  }else{
    for(int i=0;i<(int)nbinsx;i++){
      if(valuex>=minvalx+i*intervalx && valuex<minvalx+(i+1)*intervalx){
	binx=i+1;
	break;
      }
    }
  }

  //get ybin
  if(valuey<minvaly){
    biny=0;
  }else if(valuey>=minvaly+nbinsy*intervaly){
    biny=(int)nbinsy+1;
  }else{
    for(int i=0;i<(int)nbinsy;i++){
      if(valuey>=minvaly+i*intervaly && valuey<minvaly+(i+1)*intervaly){
	biny=i+1;
	break;
      }
    }
  }

  //get value
  val=pdfs2D[valtype][type]->GetBinContent(binx,biny);*/

  double val=1.0e-30;
  int binx=0,biny=0;

  binx=pdfs2D[valtype][type]->GetXaxis()->FindBin(valuex);
  biny=pdfs2D[valtype][type]->GetYaxis()->FindBin(valuey);

  //get probability
  val=pdfs2D[valtype][type]->GetBinContent(binx, biny);
  
  return val;
}

//Pi0VertexFinder Functions
//constructor
/*Pi0VertexFinder::Pi0VertexFinder(){
  ng=0;

  //define MVA
  vtxpi0=new TMVA::Reader( "!Color:!Silent" );
  vtxpi0->AddVariable( "pi0e", &var[0]);
  //vtxpi0->AddVariable( "easym", &var[1]);
  vtxpi0->AddVariable( "pi0vtxpangle", &var[2]);
  vtxpi0->AddVariable( "expangle1", &var[3]);
  vtxpi0->AddVariable( "gammaangle", &var[4]);
  vtxpi0->AddVariable( "pi0dirangle", &var[5]);
  vtxpi0->AddVariable("caldist", &var[6]);
  vtxpi0->AddVariable("eratio", &var[7]);
  vtxpi0->AddVariable("vtxmass", &var[8]);
  vtxpi0->AddVariable("mratio", &var[9]);
  //vtxpi0->AddVariable("ntrk", &var[10]);
  vtxpi0->AddVariable("ne", &var[11]);
  vtxpi0->AddVariable("nmu", &var[12]);
  vtxpi0->AddVariable("npi", &var[13]);
  vtxpi0->AddVariable("nk", &var[14]);
  vtxpi0->AddVariable("np", &var[15]);
  vtxpi0->BookMVA( "BDTG_findpi0_1vtx_LCFIPlus", "lcfiweights/TMVAClassification_BDTG_findpi0_1vtx_LCFIPlus.weights.xml" );
  vtxpi0->BookMVA( "BDTG_findpi0_2vtx_second_LCFIPlus", "lcfiweights/TMVAClassification_BDTG_findpi0_2vtx_second_LCFIPlus.weights.xml" );
  vtxpi0->BookMVA( "BDTG_findpi0_2vtx_third_LCFIPlus", "lcfiweights/TMVAClassification_BDTG_findpi0_2vtx_third_LCFIPlus.weights.xml" );

  return;
  }*/

//constructor
Pi0VertexFinder::Pi0VertexFinder(vector<string> weightfiles, vector<string> booknames, vector<double> opp, string pi0pdfname, string pi0pdfname2d){
  ng=0;

  //define MVA
  vtxpi0=new TMVA::Reader( "!Color:!Silent" );
  vtxpi0->AddVariable( "pi0e", &var[0]);
  //vtxpi0->AddVariable( "easym", &var[1]);
  vtxpi0->AddVariable( "pi0vtxpangle", &var[2]);
  vtxpi0->AddVariable( "expangle1", &var[3]);
  vtxpi0->AddVariable( "gammaangle", &var[4]);
  vtxpi0->AddVariable( "pi0dirangle", &var[5]);
  vtxpi0->AddVariable("caldist", &var[6]);
  vtxpi0->AddVariable("eratio", &var[7]);
  vtxpi0->AddVariable("vtxmass", &var[8]);
  vtxpi0->AddVariable("mratio", &var[9]);
  //vtxpi0->AddVariable("ntrk", &var[10]);
  vtxpi0->AddVariable("ne", &var[11]);
  vtxpi0->AddVariable("nmu", &var[12]);
  vtxpi0->AddVariable("npi", &var[13]);
  vtxpi0->AddVariable("nk", &var[14]);
  vtxpi0->AddVariable("np", &var[15]);
  vtxpi0->BookMVA( booknames[0].c_str(), weightfiles[0].c_str() );
  vtxpi0->BookMVA( booknames[1].c_str(), weightfiles[1].c_str() );
  vtxpi0->BookMVA( booknames[2].c_str(), weightfiles[2].c_str() );

  _booknames = booknames;
  _opp = opp;

  //pi0finder
  pi0ID = new Pi0Finder(pi0pdfname, pi0pdfname2d);

  return;
}

//destructor
Pi0VertexFinder::~Pi0VertexFinder(){
  //cout << "vtxpi0: " << vtxpi0 << endl;
  delete vtxpi0;
  delete pi0ID;

  return;
}

//functions
void Pi0VertexFinder::Initialize(){
  //gamma vector initialize
  vector < vector < double > >().swap(gv);
  //pi0 vector initialize
  vector < vector < double > >().swap(pi0vec);
  //num. of gammas initialize
  ng=0;
  return;
}

void Pi0VertexFinder::Make_GammaVector(const lcfiplus::Neutral* neut){
  double g1prob=0.0,g2prob=0.0;
  double vars[9];

  //initilize
  for(int i=0;i<9;i++) vars[i]=0.0;

  //get deposit energy
  double ecal=0.0,hcal=0.0,mucal=0.0;
  EVENT::ClusterVec cluvec=neut->getClusters();
  EVENT::FloatVec shapes;
  if(cluvec.size()!=0){
    for(unsigned int i=0;i<cluvec.size();i++){
      ecal+=cluvec[i]->getSubdetectorEnergies()[0];
      hcal+=cluvec[i]->getSubdetectorEnergies()[1];
      mucal+=cluvec[i]->getSubdetectorEnergies()[2];
    }

    shapes=cluvec[0]->getShape();
  }

  //get neutral and make gamma vector
  //check whether it is gamma candidates
  vars[0]=neut->E();
  vars[1]=neut->E();
  //vars[2]=TMath::Cos(g1.Vect().Angle(g2.Vect()));
  //vars[3]=pi0.M();
  vars[4]=0.0;
  if(ecal+hcal!=0.0) vars[4]=ecal/(ecal+hcal);
  if(shapes.size()!=0){
    vars[5]=shapes[0];
    vars[6]=shapes[5];
    vars[7]=fabs(shapes[3])/shapes[6];
    vars[8]=shapes[15]/(2*3.50);
  }

  //check if it is gamma
  pi0ID->Initialize();  //set prior 0.5;
  g1prob=pi0ID->Get_g1Prob(vars);
  pi0ID->Initialize();  //set prior 0.5;
  g2prob=pi0ID->Get_g2Prob(vars);

  // cout << "check g1 and g2: " << g1prob << " " << g2prob << endl;

  if(!(g1prob<0.0 && g2prob<0.0)){    //add this neutral to gamma vectors
    vector<double> tmpgv;   //13 elements
    tmpgv.push_back(neut->Px());
    tmpgv.push_back(neut->Py());
    tmpgv.push_back(neut->Pz());
    tmpgv.push_back(neut->E());
    for(int i=4;i<9;i++) tmpgv.push_back(vars[i]);
    if(cluvec.size()!=0){
      tmpgv.push_back(cluvec[0]->getPosition()[0]);
      tmpgv.push_back(cluvec[0]->getPosition()[1]);
      tmpgv.push_back(cluvec[0]->getPosition()[2]);
    }else{ //fill 0.0 3 times
      tmpgv.push_back(0.0);
      tmpgv.push_back(0.0);
      tmpgv.push_back(0.0);
    }
    tmpgv.push_back(0.0);  //for flg element
    gv.push_back(tmpgv);

    ng++;
  }
  
  return;
}

void Pi0VertexFinder::Make_Pi0Vector(TVector3 vtx, TVector3 vtxdir){

  //define matrix
  dm.ResizeTo(ng,ng);

  //cal. posterior for all the combinations
  double pvars[19];
  double cene=0.0;
  int nph=0;
  double pi0value=1.0;

  TLorentzVector tmpg1,tmpg2,tmppi0,tmpph;
  TVector3 tmpdir;

  vector < double > tmpgv1,tmpgv2;   //13 elements
  for(int i=0;i<ng;i++){
    for(int k=0;k<13;k++) tmpgv1.push_back(gv.at(i).at(k));
    for(int j=i+1;j<ng;j++){
      for(int k=0;k<13;k++) tmpgv2.push_back(gv.at(j).at(k));

      //correct dir vector
      tmpdir.SetXYZ(tmpgv1.at(9),tmpgv1.at(10),tmpgv1.at(11));
      tmpdir=tmpdir-vtx;
      tmpg1.SetPxPyPzE(tmpdir.Unit().X()*tmpgv1.at(3)
		       ,tmpdir.Unit().Y()*tmpgv1.at(3)
		       ,tmpdir.Unit().Z()*tmpgv1.at(3)
		       ,tmpgv1.at(3));

      tmpdir.SetXYZ(tmpgv2.at(9),tmpgv2.at(10),tmpgv2.at(11));
      tmpdir=tmpdir-vtx;
      tmpg2.SetPxPyPzE(tmpdir.Unit().X()*tmpgv2.at(3)
		       ,tmpdir.Unit().Y()*tmpgv2.at(3)
		       ,tmpdir.Unit().Z()*tmpgv2.at(3)
		       ,tmpgv2.at(3));

      tmppi0=tmpg1+tmpg2;
      
      pvars[2]=tmpg1.Vect().Angle(tmpg2.Vect());
      pvars[3]=tmppi0.M();
      
      if(tmpg1.E()>=tmpg2.E()){
	pvars[0]=tmpg1.E();
	pvars[1]=tmpg2.E();
	pvars[4]=tmpgv1.at(4);
	pvars[5]=tmpgv1.at(5);
	pvars[6]=tmpgv1.at(6);
	pvars[7]=tmpgv1.at(7);
	pvars[8]=tmpgv1.at(8);
	pvars[9]=tmpgv2.at(4);
	pvars[10]=tmpgv2.at(5);
	pvars[11]=tmpgv2.at(6);
	pvars[12]=tmpgv2.at(7);
	pvars[13]=tmpgv2.at(8);
      }else{
	pvars[0]=tmpg2.E();
	pvars[1]=tmpg1.E();
	pvars[4]=tmpgv2.at(4);
	pvars[5]=tmpgv2.at(5);
	pvars[6]=tmpgv2.at(6);
	pvars[7]=tmpgv2.at(7);
	pvars[8]=tmpgv2.at(8);
	pvars[9]=tmpgv1.at(4);
	pvars[10]=tmpgv1.at(5);
	pvars[11]=tmpgv1.at(6);
	pvars[12]=tmpgv1.at(7);
	pvars[13]=tmpgv1.at(8);
      }
      pvars[14]=(tmppi0.E()*(1-tmppi0.Beta()))/(2*pvars[1]);
      pvars[15]=(2*pvars[1])/(tmppi0.E()*(1+tmppi0.Beta()));
      pvars[16]=(2*pvars[0])/(tmppi0.E()*(1+tmppi0.Beta()));

      //cal. photon energy inside the cone
      cene=0.0;
      nph=0;
      for(int m=0;m<ng;m++){
	if(m==i || m==j) continue;  //avoid adding themselves
	//correct direction
	tmpdir.SetXYZ(gv.at(m).at(9),gv.at(m).at(10),gv.at(m).at(11));
	tmpdir=tmpdir-vtx;
	tmpph.SetPxPyPzE(tmpdir.Unit().X()*gv.at(m).at(3),
			 tmpdir.Unit().Y()*gv.at(m).at(3),
			 tmpdir.Unit().Z()*gv.at(m).at(3),
			 gv.at(m).at(3));

	double tmptheta=0.0;
	if(tmpg1.E()>=tmpg2.E()) tmptheta=tmpg1.Vect().Angle(tmpph.Vect());
	else  tmptheta=tmpg2.Vect().Angle(tmpph.Vect());

	if(pvars[2]>tmptheta){
	  cene=cene+gv.at(m).at(3);
	  nph=nph+1;
	}
      }
      pvars[17]=cene;
      pvars[18]=(double)nph;

      //cal. posterior
      pi0value=pi0ID->Get_Prob(pvars);
      if(pi0value<0.0) pi0value=1.0e-100;
      if((int) gv.at(i).at(12)==1 || (int)gv.at(j).at(12)==1) pi0value=1.0e-100;  //check if it is already used
      dm[i][j]=TMath::Log(pi0value);
      dm[j][i]=TMath::Log(pi0value);

      tmpgv2.clear();
    }
    tmpgv1.clear();
  }

  return;
}

void Pi0VertexFinder::Reco_Pi0s(TVector3 vtx, TVector3 vtxdir){
  Make_Pi0Vector(vtx,vtxdir);

  //get global minimum!
  double ic=Get_GlobalMinimum();

  //dummy
  ic = ic+1;

  return;
}

//perform global minimization
double Pi0VertexFinder::Get_GlobalMinimum(){
  //initialize
  counter=0;

  //coefficient of strength of the penalty term
  double coeff=0.030; //for new?      

  //start to look for the best pairing
  //first. making basic combination

  //this variables are global!
  double BIC=0.0,likeli=0.0;
  int nupg=0, *combi;
  combi=new int [ng];  // this is used for combination check

  //first. making basic combination
  for(int k=0;k<(int)ng/2;k++){
    //cout << "dm: " << dm[k][k+1] << endl;
    if(dm[2*k][2*k+1]>-200.0){
      likeli+=dm[2*k][2*k+1];
      combi[2*k]=2*k+1;
      combi[2*k+1]=2*k;
    }else{
      combi[2*k]=-1;
      combi[2*k+1]=-1;
      nupg+=2;
    }
  }

  //when num. of gammas is odd  
  if(ng%2==1){
    combi[ng-1]=-1;
    nupg++;
  }

  BIC=-2.0*likeli+coeff*nupg*TMath::Log(ng);
  //cout << "before exchange: " << likeli << " " << BIC << endl;

  //to reserve candidates
  double okdiff=0.0,okBIC=0.0;
  int oknupg=0,okcombi=0;
  //for tempolary
  double tmpdiff=0.0,tmpBIC=0.0;
  int tmpnupg=0;

  for(int mm=0;mm<2;mm++){   //looping for converging twice is enough?
    for(int k=0;k<ng;k++){   //for first gamma
      okdiff=0.0;
      okBIC=BIC;
      oknupg=nupg;
      okcombi=-1;
      
      for(int l=0;l<ng;l++){   //for second gamma
	if(k==l) continue;
	tmpdiff=0.0;
	tmpnupg=nupg;

	//first, kill bad one
	if(dm[k][l]<-200.0) continue;

	if(combi[k]==-1 && combi[l]==-1){   //when both of gammas aren't used - try to combine them
	  if(dm[k][l]>-200.0){
	    tmpdiff=(dm[k][l]+dm[l][k])/2.0;
	    tmpnupg=nupg-2;
	  }/*else{
	    tmpdiff=0.0;
	    tmpnupg=nupg;
	    }*/		
	}else if(combi[k]!=-1 && combi[l]==-1){
	  tmpdiff=dm[k][l]-dm[k][combi[k]];
	  tmpnupg=nupg;
	  //kill bad one
	  //if(dm[k][l]<-200.0) tmpdiff=0.0;
	}else if(combi[k]==-1 && combi[l]!=-1){
	  tmpdiff=dm[k][l]-dm[l][combi[l]];
	  tmpnupg=nupg;
	  //kill bad one
	  //if(dm[k][l]<-200.0) tmpdiff=0.0;
	}else if(combi[k]!=-1 && combi[l]!=-1){
	  tmpdiff=dm[k][l]+dm[combi[k]][combi[l]]-dm[k][combi[k]]-dm[l][combi[l]];
	  tmpnupg=nupg;
	  //kill bad one
	  if(dm[combi[k]][combi[l]]<-200.0) tmpdiff=0.0;
	}
	
	tmpBIC=-2.0*(likeli+tmpdiff)+coeff*tmpnupg*TMath::Log(ng);
	
	if(fabs(tmpdiff)<200.0 && fabs(tmpdiff)>0.0 && okBIC>tmpBIC){
	  //candidate for swapping - reserve it!
	  okBIC=tmpBIC;
	  okdiff=tmpdiff;
	  oknupg=tmpnupg;
	  okcombi=l;
	}
      }
      
      //before swapping, check if de-connect is good
      if(combi[k]>=0 && dm[k][combi[k]]>-200.0){
	double deBIC=-2.0*(likeli-dm[k][combi[k]])+coeff*(nupg+2)*TMath::Log(ng);
	
	if(okBIC>deBIC){  //deconnect. update global variables
	  likeli-=dm[k][combi[k]];
	  int tmpc=combi[k];	  
	  combi[k]=-1;
	  combi[tmpc]=-1;
	  nupg=nupg+2;
	  BIC=-2.0*likeli+coeff*nupg*TMath::Log(ng);
	  
	  continue;
	}
      }
       
      //update global variables
      if(fabs(okdiff)<200.0 && fabs(okdiff)>0.0){
	int tmpc=combi[k];
	int tmpc2=combi[okcombi];
	combi[k]=okcombi;
	combi[okcombi]=k;
	if(tmpc!=-1) combi[tmpc]=tmpc2;
	if(tmpc2!=-1) combi[tmpc2]=tmpc;
	likeli=likeli+okdiff;
	nupg=oknupg;
	BIC=-2.0*likeli+coeff*nupg*TMath::Log(ng);
      }else if(combi[k]!=-1 && dm[k][combi[k]]<-200.0){  //de-connect combination
	int tmpc=combi[k];	  
	combi[k]=-1;
	combi[tmpc]=-1;
	likeli=likeli-TMath::Log(1.0e-100);
	nupg=nupg+2;
	BIC=-2.0*likeli+coeff*nupg*TMath::Log(ng);
      }
    }
  }

  //cout << "after exchange: " << likeli << " " << BIC << endl;

  //making pi0 vectors
  vector<double> tmppv;  //16 elements
  for(int k=0;k<ng;k++){
    //cout << "combi: " << k << " " << combi[k] << endl;
    //check if it is used
    Bool_t uf=true;
    for(int l=0;l<k;l++){
      if(k==combi[l]){
	uf=false;
	break;
      }
    }
    
    if(uf==true && combi[k]!=-1){   //fill this gamma pair to pi0 vector
      tmppv.push_back(gv.at(k).at(0));
      tmppv.push_back(gv.at(k).at(1));
      tmppv.push_back(gv.at(k).at(2));
      tmppv.push_back(gv.at(k).at(3));
      tmppv.push_back(gv.at(combi[k]).at(0));
      tmppv.push_back(gv.at(combi[k]).at(1));
      tmppv.push_back(gv.at(combi[k]).at(2));
      tmppv.push_back(gv.at(combi[k]).at(3));
      tmppv.push_back(gv.at(k).at(9));
      tmppv.push_back(gv.at(k).at(10));
      tmppv.push_back(gv.at(k).at(11));
      tmppv.push_back(gv.at(combi[k]).at(9));
      tmppv.push_back(gv.at(combi[k]).at(10));
      tmppv.push_back(gv.at(combi[k]).at(11));
      tmppv.push_back(k);
      tmppv.push_back(combi[k]);

      //fill it!
      pi0vec.push_back(tmppv);
      counter++;
     }
     tmppv.clear();
  }
  
  //arrange pi0s in descending order
  for(int k=0;k<counter;k++){
    for(int l=k+1;l<counter;l++){
      double pie1=pi0vec.at(k).at(3)+pi0vec.at(k).at(7);
      double pie2=pi0vec.at(l).at(3)+pi0vec.at(l).at(7);
      if(pie1<pie2){
	swap(pi0vec[k],pi0vec[l]);
      }
    }
  }

  delete[] combi;
 
  return BIC;
}

TLorentzVector Pi0VertexFinder::Corr_VtxVect(int vtxtype, TLorentzVector vtx1vect, int *npart, TVector3 vtx, TVector3 vtxdir){
  TVector3 tmpdr1,tmpdr2,g1,g2,caldist;  //use for position vectors
  TLorentzVector tmpg1,tmpg2,pi0,hvect,cvect,pi0vect;  //use for gammas and pi0s 4-momentum
  TLorentzVector corrg1,corrg2,corrpi0; //for pi0 energy correction. reserved
  double pi0value=1.0,cutpos=0.80;   //use for MVAoutput and MVAcut position
  npi0=0;

  //cout << "num. of pi0s created: " << pi0vec.size() << endl;

  for(int k=0;k<counter;k++){
    //check if this pi0 is used
    if((int)gv.at((int)pi0vec.at(k).at(14)).at(12)==1 || (int)gv.at((int)pi0vec.at(k).at(15)).at(12)==1) continue;
    
    g1.SetXYZ(pi0vec.at(k).at(8),
	      pi0vec.at(k).at(9),
	      pi0vec.at(k).at(10));
    
    g2.SetXYZ(pi0vec.at(k).at(11),
	      pi0vec.at(k).at(12),
	      pi0vec.at(k).at(13));
    caldist=g1-g2;

    //correct direction of photon
    tmpdr1=g1-vtx;
    tmpdr2=g2-vtx;
    
    //correct gamma momentum
    tmpg1.SetPxPyPzE(tmpdr1.Unit().X()*pi0vec.at(k).at(3),
		     tmpdr1.Unit().Y()*pi0vec.at(k).at(3),
		     tmpdr1.Unit().Z()*pi0vec.at(k).at(3),
		     pi0vec.at(k).at(3));
    
    tmpg2.SetPxPyPzE(tmpdr2.Unit().X()*pi0vec.at(k).at(7),
		     tmpdr2.Unit().Y()*pi0vec.at(k).at(7),
		     tmpdr2.Unit().Z()*pi0vec.at(k).at(7),
		     pi0vec.at(k).at(7));

    pi0=tmpg1+tmpg2;
    hvect=vtx1vect+pi0vect;
    cvect=pi0+hvect;
    
    //get pi0 using MVA
    var[0]=pi0.E();
    var[1]=(tmpg1.E()-tmpg2.E())/(tmpg1.E()+tmpg2.E());
    var[2]=hvect.Vect().Angle(pi0.Vect());
    var[3]=TMath::ACos((2.0*tmpg1.E()*tmpg2.E()-0.13497*0.13497)/(2.0*tmpg1.E()*tmpg2.E()))-tmpdr1.Angle(tmpdr2);
    var[4]=tmpg1.Vect().Angle(tmpg2.Vect());
    var[5]=vtxdir.Angle(pi0.Vect());
    var[6]=caldist.Mag();
    var[7]=hvect.E()/(pi0.E()+hvect.E());
    var[8]=hvect.M();
    var[9]=(cvect.M()-hvect.M())/cvect.M();
    var[10]=npart[0];
    var[11]=npart[1];
    var[12]=npart[2];
    var[13]=npart[3];
    var[14]=npart[4];
    var[15]=npart[5];
    
    //set MVA cut position (so far only >=2 vertex is checked) 
    if(npart[0]>1){
      switch(vtxtype){   //1vtx case
      case 0:
	pi0value=vtxpi0->EvaluateMVA(_booknames[0]);
	cutpos = _opp[0];
	
	//if(npart[0]==2) cutpos=0.71;  //bias! need to erase <-so far, no idea...
	//else if(npart[0]==3) cutpos=0.79;   //0.76;  //bias! need to erase <-so far, no idea...
	//else cutpos=0.83;
	//if(npart[0]==2) cutpos =cutpos*0.88;  //bias! sofar, no idea to escape
	//if(npart[0]==3) cutpos =cutpos*0.93;  //bias! sofar, no idea to escape
	break;

      case 1:  //secondary of 2vtx case
	pi0value=vtxpi0->EvaluateMVA(_booknames[1]);
	cutpos = _opp[1];
	
	//if(npart[0]==2) cutpos=0.61;     //bias! need to erase <-so far, no idea...
	//if(npart[0]==2) cutpos =cutpos*0.97;  //bias! sofar, no idea to escape
	break;
      case 2:  //tertiary of 2vtx case
	pi0value=vtxpi0->EvaluateMVA(_booknames[2]);
	cutpos = _opp[2];
	
	//if(npart[0]==2) cutpos=0.74;  //bias! need to erase <-so far, no idea...
	//else cutpos=0.68;
	//if(npart[0]==2) cutpos =cutpos*1.09;  //bias! sofar, no idea to escape
	break;
      }
      //if(k<=3) cout << "pi0 value: " << pi0value << " " << cutpos << " " << var[2] << endl;
      
      if(pi0value>cutpos && var[2]>0.0){   //need optimization!
	//pi0 energy correction!(is it good?)
	corrg1=(0.134976/pi0.M())*tmpg1;
	corrg2=(0.134976/pi0.M())*tmpg2;
	corrpi0=corrg1+corrg2;
	
	npi0++;
	pi0vect+=pi0;   //corrpi0;
	
	//check uflg
	gv.at((int)pi0vec.at(k).at(14)).at(12)=1;
	gv.at((int)pi0vec.at(k).at(15)).at(12)=1;
      }
    }
  }
  
  //cout << "num. of pi0s attached: " << npi0 << endl;

  //reset pi0vec
  vector < vector < double > >().swap(pi0vec);
  
  return vtx1vect+pi0vect;
}

int Pi0VertexFinder::Get_nPi0(){
  return npi0;
}
