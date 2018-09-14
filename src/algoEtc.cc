#include "lcfiplus.h"
#include "algoEtc.h"

#include "TRandom.h"
#include <math.h>


namespace lcfiplus {
namespace algoEtc {
int min(int a, int b) {
  return a > b ? b : a;
}

void makeBeamTracks(Track*& t1, Track*& t2, bool smear) {
  // beam crossing = 14 mrad
  float beamtd = tan( 0.5*(3.1415926 - 14e-3) );
  //float beamsizeX = 639e-6; // 639 nm converted to mm
  //float beamsizeY = 5.7e-6; // 5.7 nm converted to mm
  // size(d)/size(z) = tan(7e-3) -> size(z) = size(d)/tan(7e-3)
  // float beamsizeZ = beamsizeX / tan(0.5 * 14e-3);

  float beamsizeX = Globals::Instance()->getBeamSizeX();
  float beamsizeZ = Globals::Instance()->getBeamSizeZ();

  float d0rand(0), z0rand(0);

  if (smear) {
    d0rand = gRandom->Gaus(0,beamsizeX);
    z0rand = gRandom->Gaus(0,beamsizeZ);
  }

  t1 = new Track;
  t2 = new Track;

  double par[tpar::parN];
  double cov[tpar::covN];

  t1->setId(1000001);

  par[tpar::d0] = d0rand;
  par[tpar::z0] = z0rand;
  par[tpar::om] = 0.3*Globals::Instance()->getBField()/250./beamtd*1000.;
  par[tpar::ph] = 0;
  par[tpar::td] = beamtd;
  t1->setHelix(par);

  memset(cov, 0, sizeof(cov));
  cov[tpar::d0d0] = pow(beamsizeX,2);
  cov[tpar::z0z0] = pow(beamsizeZ,2);
  cov[tpar::phph] = 1;
  cov[tpar::omom] = 1;
  cov[tpar::tdtd] = 1;
  t1->setCovMatrix(cov);

  d0rand = 0;
  z0rand = 0;
  if (smear) {
    d0rand = gRandom->Gaus(0,beamsizeX);
    z0rand = gRandom->Gaus(0,beamsizeZ);
  }

  t2->setId(1000002);

  par[tpar::d0] = d0rand;
  par[tpar::z0] = z0rand;
  par[tpar::om] = 0.3*Globals::Instance()->getBField()/250./beamtd*1000.;
  par[tpar::ph] = 0;
  par[tpar::td] = -beamtd;
  t2->setHelix(par);

  memset(cov, 0, sizeof(cov));
  cov[tpar::d0d0] = pow(beamsizeX,2);
  cov[tpar::z0z0] = pow(beamsizeZ,2);
  cov[tpar::phph] = 1;
  cov[tpar::omom] = 1;
  cov[tpar::tdtd] = 1;
  t2->setCovMatrix(cov);
}

void makeBeamVertex(Vertex*& vtx, bool smear) {
  /*
  float beamsizeX = 639e-6; // 639 nm converted to mm
  float beamsizeY = 5.7e-6; // 5.7 nm converted to mm
  float beamsizeZ = beamsizeX / tan(0.5 * 14e-3);
  */

  float beamsizeX = Globals::Instance()->getBeamSizeX();
  float beamsizeY = Globals::Instance()->getBeamSizeY();
  float beamsizeZ = Globals::Instance()->getBeamSizeZ();

  double cov[6];
  cov[Vertex::xx] = pow(beamsizeX,2);
  cov[Vertex::xy] = 0;
  cov[Vertex::xz] = 0;
  cov[Vertex::yy] = pow(beamsizeY,2);
  cov[Vertex::yz] = 0;
  cov[Vertex::zz] = pow(beamsizeZ,2);

  if (smear)
    vtx = new Vertex(0,1,gRandom->Gaus(0, beamsizeX), gRandom->Gaus(0,beamsizeY), gRandom->Gaus(0,beamsizeZ),cov, true);
  else
    vtx = new Vertex(0,1,0.,0.,0.,cov, true);
}

void connectVerticesToJets(const JetVec& jets, const vector<Vertex*>& vtcs, vector<vector<Vertex*> >& jetVertices, vector<vector<const Track*> >& jetResidualTracks, const Vertex* ip) {
  unsigned int nj = jets.size();
  jetVertices.resize(nj);
  jetResidualTracks.resize(nj);

  for (unsigned int j=0; j<nj; j++) {
    jetResidualTracks[j] = jets[j]->getTracks();

    TrackVec trackVecV = algoEtc::extractTracks(jets[j]->getVertices());
    jetResidualTracks[j].insert(jetResidualTracks[j].end(), trackVecV.begin(), trackVecV.end());

    // remove IP if supplied
    if (ip) {
      for (unsigned int t=0; t<ip->getTracks().size(); t++) {
        vector<const Track*>::iterator itt = find(jetResidualTracks[j].begin(), jetResidualTracks[j].end(), ip->getTracks()[t]);
        if (itt != jetResidualTracks[j].end())
          jetResidualTracks[j].erase(itt);
      }
    }
  }

  for (unsigned int v=0; v<vtcs.size(); v++) {
    Vertex* curv = vtcs[v];
    if (curv->getTracks().size() == 0) {
      cout << "connectVerticesToJets: vertex-jet association failed! No tracks found in a vertex." << endl;
      continue;
    }

    // selecting associating jets with energy fraction
    vector<double> efrac;
    efrac.resize(nj);

    for (unsigned int vt=0; vt<curv->getTracks().size(); vt++) {
      const Track* curt = curv->getTracks()[vt];
      for (unsigned int j2=0; j2<nj; j2++) {
        vector<const Track*>::iterator itt = find(jetResidualTracks[j2].begin(), jetResidualTracks[j2].end(), curt);
        if (itt != jetResidualTracks[j2].end()) {
          jetResidualTracks[j2].erase(itt);
          efrac[j2] += curt->E();
          break;
        }
      }
    }
    double maxefrac = 0.;
    int j = -1;
    for (unsigned int j2=0; j2<nj; j2++) {
      if (efrac[j2] > maxefrac) {
        maxefrac = efrac[j2];
        j = j2;
      }
    }

    // fix to cope with dropped tracks e.g. in the context of running kt algorithm prior to LCFIPlus
    //if(j == -1)throw(Exception("connectVerticesToJets: vertex-jet association failed!"));
    if (j == -1) {
      cout << "connectVerticesToJets: vertex-jet association failed, vertex dropped" << endl;
      continue;
    }

    // from here: curv is determined to belong to jet j
    //cout << "Jet-vertex pairing: vertex " << v << " assigned to jet " << j << endl;

    jetVertices[j].push_back(curv);

  }
}


vector<const Track*> extractTracks(VertexVec& vtx) {
  vector<const Track*>ret;
  for (unsigned int n=0; n<vtx.size(); n++) {
    ret.insert(ret.end(), vtx[n]->getTracks().begin(), vtx[n]->getTracks().end());
  }

  return ret;
}

double calcThrust( vector<TVector3>& list, TVector3& taxis ) {

  const int nwork=11,iFastMax = 4,iGood=2;
  const float dConv=0.0001; // 0.0001
  int sgn;
  double theta=0,phi=0;
  double thp,thps,tds,tmax;
  // double dOblateness;
  vector<TVector3> TAxes(3),Fast(iFastMax+1),Workv(nwork);
  vector<double> Workf(nwork),dThrust(3);
  TVector3 tdi,tpr,mytest;

  tmax = 0;
  for ( unsigned int i=0; i < list.size(); i++)
    tmax += list[i].Mag();

  // pass = 0: find thrust axis
  // pass = 1: find major axis
  for ( int pass=0; pass <= 1; pass++ ) {
    if ( pass == 1 ) {
      phi   = TAxes[0].Phi();
      theta = TAxes[0].Theta();
      for ( unsigned  int i = 0; i < list.size(); i++) {
        list[i].RotateZ(-phi);
        list[i].RotateY(-theta);
      }
      TAxes[0].SetXYZ(0,0,1);
    } // if pass == 1

    // Find the ifast highest momentum particles and
    // put the highest in Fast[0], next in Fast[1],....Fast[iFast-1].
    // Fast[iFast] is just a workspace.

    for ( unsigned  int i = 0; i < Fast.size(); i++ )
      Fast[i].SetXYZ(0,0,0);

    for ( unsigned int i = 0; i < list.size(); i++ ) {
      for ( int ifast = iFastMax -1; ifast >= 0 ; ifast-- ) {
        if (list[i].Mag2() > Fast[ifast].Mag2() ) {
          Fast[ifast + 1] = Fast[ifast];
          if (ifast == 0) Fast[ifast] = list[i];
        } else {
          Fast[ifast + 1] = list[i];
          break;
        } // if p>p_fast
      } // for ifast
    } // for i

    // Find axis with highest thrust (case 0)/ highest major (case 1).

    for ( unsigned int iw = 0; iw < Workv.size(); iw++ ) {
      Workf[iw] = 0.;
    }
    int p = (int) min( iFastMax,list.size() ) - 1 ;
    int nc = 1 << p;
    for ( int n = 0; n < nc; n++ ) {
      tdi.SetXYZ(0,0,0);
      for (int i = 0; i < min(iFastMax,nc) ; i++) {
        if ( (1 << (i+1)) * ( (n + (1<<i)) / (1<<(i+1)) ) >= n+1) { //i+1
          sgn = -1;
        } else {
          sgn = 1;
        }
        tdi += sgn*Fast[i];
        if (pass==1) tdi.SetZ(0);
      } // for i
      tds = tdi.Mag2();
      for ( int iw = (int) min(n,9); iw >= 0; iw-- ) {
        if (tds > Workf[iw]) {
          Workf[iw+1] = Workf[iw];
          Workv[iw+1] = Workv[iw];
          if (iw == 0) {
            Workv[iw] = tdi;
            Workf[iw] = tds;
          }
        } else { // if tds
          Workv[iw+1] = tdi;
          Workf[iw+1] = tds;
        } // if tds
      } // for iw
    } // for n

    // Iterate direction of axis until stable maximum.

    dThrust[pass] = 0;
    int nagree = 0;
    for ( int iw = 0; iw < min(nc,10) && nagree < iGood; iw++ ) {
      thp = 0;
      thps = -99999.;
      while ( thp > thps + dConv ) {
        thps = thp;
        if ( thp <= 1E-10 ) {
          tdi = Workv[iw];
        } else {
          tdi=tpr;
        }
        tpr.SetXYZ(0,0,0);
        for ( unsigned int i = 0; i < list.size(); i++ ) {
          sgn = (tdi.Dot(list[i]) > 0 ? 1 : -1);
          tpr += sgn*list[i];
          if (pass == 1) {
            tpr.SetZ(0);  // ###
          }
        } // for i
        thp = tpr.Mag()/tmax;
      } // while
      // Save good axis. Try new initial axis until enough
      // tries agree.
      if ( thp < dThrust[pass] - dConv ) continue;
      if ( thp > dThrust[pass] + dConv ) {
        nagree = 0;
        // 	      if (myrnd.flat() > 0.49999)
        // 		{sgn = 1;} else {sgn=-1;}
        sgn = 1;
        TAxes[pass] = sgn*tpr*(1./(tmax*thp));
        dThrust[pass] = thp;
      } // if thp
      nagree++;
    } // for iw (2)
  } // for pass ...

  // Find minor axis and value by orthogonality.
  if (gRandom->Rndm() > 0.49999) {
    sgn = 1;
  } else {
    sgn=-1;
  }
  TAxes[2].SetXYZ( -sgn*TAxes[1].y(), sgn*TAxes[1].x(), 0);
  thp = 0.;
  for ( unsigned int i = 0; i < list.size(); i++ ) {
    thp += fabs(TAxes[2].Dot(list[i]) );
  } // for i
  dThrust[2] = thp/tmax;

  // Rotate back to original coordinate system.
  for ( unsigned int i = 0; i < TAxes.size(); i++) {
    TAxes[i].RotateY(theta);
    TAxes[i].RotateZ(phi);
  }
  // dOblateness = dThrust[1] - dThrust[2];

  //_principleThrustValue = dThrust[0];
  //_majorThrustValue     = dThrust[1];
  //_minorThrustValue     = dThrust[2];
  //_principleThrustAxis  =   TAxes[0];
  //_majorThrustAxis      =   TAxes[1];
  //_minorThrustAxis      =   TAxes[2];

  taxis = TAxes[0];
  return dThrust[0];
}

bool SimpleSecMuonFinder(const Track* tr, double d0sigth, double z0sigth, double maxpos, double mudepmin,
                         double ecaldepmin, double ecaldepmax, double hcaldepmin, double hcaldepmax, double maxclusterpertrackenergy,
                         const Vertex* ip) {

  double z0 = fabs(tr->getZ0()-ip->getZ()); 
  double dist = sqrt(tr->getD0() *tr->getD0() + z0 * z0);
  double sigd0 = fabs(tr->getD0()) / sqrt(tr->getCovMatrix()[tpar::d0d0]);
  double sigz0 = z0 / sqrt(tr->getCovMatrix()[tpar::z0z0]);
  double mudep = tr->getCaloEdep()[tpar::yoke];
  double ecaldep = tr->getCaloEdep()[tpar::ecal];
  double hcaldep = tr->getCaloEdep()[tpar::hcal];
  /*
  		if(tr->getMcp() && (fabs(tr->getMcp()->getPDG()) == 13 || mudep > 0 || (tr->E() > 5 && (ecaldep + hcaldep < tr->E() / 2.)))){
  		  TVector3 dirvec = tr->Vect().Unit();
  			cout << "Simple Sec Muon Finder: MCP = " << tr->getMcp()->getPDG() << " E = " << tr->E() << " sigd0 = " << sigd0 << " sigz0 = " << sigz0;
  			cout << " mudep = " << mudep << " dist = " << dist << " cosq = " << dirvec.z() << " ecaldep/sinq = " << ecaldep / dirvec.Pt() << " ecaldep/cosq = " << ecaldep / dirvec.z() << " hcaldep = " << hcaldep << endl;
  		}
  */
  return (sigd0 > d0sigth || sigz0 > z0sigth) && dist < maxpos && mudep > mudepmin
         && ecaldep > ecaldepmin && ecaldep < ecaldepmax && hcaldep > hcaldepmin && hcaldep < hcaldepmax && (ecaldep + hcaldep)/tr->E() < maxclusterpertrackenergy;
}

bool SimpleSecElectronFinder(const Track* tr, double d0sigth, double z0sigth, double maxpos, double emin,
                             double minfracecal, double minecalpertrackenergy, double maxecalpertrackenergy,
                             const Vertex* ip) {

  double z0 = fabs(tr->getZ0()-ip->getZ()); 
  double dist = sqrt(tr->getD0() *tr->getD0() + z0 * z0);
  double sigd0 = fabs(tr->getD0()) / sqrt(tr->getCovMatrix()[tpar::d0d0]);
  double sigz0 = z0 / sqrt(tr->getCovMatrix()[tpar::z0z0]);
  double ecaldep = tr->getCaloEdep()[tpar::ecal];
  double hcaldep = tr->getCaloEdep()[tpar::hcal];
  double fracecal = ecaldep / (ecaldep + hcaldep);
  double ecalpertrackenergy = ecaldep / tr->E();

  return (sigd0 > d0sigth || sigz0 > z0sigth) && dist < maxpos && tr->E() > emin
         && fracecal > minfracecal && ecalpertrackenergy > minecalpertrackenergy && ecalpertrackenergy < maxecalpertrackenergy;

}


}
}
