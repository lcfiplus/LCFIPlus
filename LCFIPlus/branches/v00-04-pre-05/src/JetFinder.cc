#include "EventStore.h"
#include "JetFinder.h"

#include <assert.h>

#include <string.h>
#include <string>
#include <algorithm>
#include <vector>
#include <math.h>
#include <TMatrixTSym.h>
#include <functional>

#include "TVector3.h"
#include "TLorentzVector.h"

using namespace std;

namespace lcfiplus {

  /* Jet resolution functions */
  double JetFinder::funcDurham(Jet& jet1, Jet& jet2, double Evis2){
    /*
       D(I,J) = 2.*MIN(PL(4,I)*PL(4,I),PL(4,J)*PL(4,J))*MAX(0.,(1.-
       (PL(1,I)*PL(1,J)+PL(2,I)*PL(2,J)+PL(3,I)*PL(3,J))/
       (PL(6,I)*PL(6,J))))
     */

    double e1 = jet1.E();
    double e2 = jet2.E();
    TVector3 mom1 = jet1.Vect();
    TVector3 mom2 = jet2.Vect();
    double val = 2*min(e1*e1,e2*e2)*max(0., 1-(mom1.Dot(mom2))/(mom1.Mag()*mom2.Mag()));
    return val/Evis2;
  }

  double JetFinder::funcDurhamVertex(Jet& jet1, Jet& jet2, double Evis2){
		// do not combine vertex-oriented jets
		if(jet1.getVertices().size() > 0 && jet2.getVertices().size() > 0)
			return 1e+300;

		return funcDurham(jet1, jet2, Evis2);
  }


  double JetFinder::funcDurhamCheat(Jet& jet1, Jet& jet2, double Evis2){
    double e1 = jet1.E();
    double e2 = jet2.E();
    TVector3 mom1 = jet1.Vect();
    TVector3 mom2 = jet2.Vect();
    double val = 2*min(e1*e1,e2*e2)*max(0., 1-(mom1.Dot(mom2))/(mom1.Mag()*mom2.Mag()));
    return val/Evis2;
  }

  double JetFinder::funcJade(Jet& jet1, Jet& jet2, double Evis2){
    /*
       JADE(I,J) = 2.*PL(4,I)*PL(4,J)*MAX(0.,(1.-
       (PL(1,I)*PL(1,J)+PL(2,I)*PL(2,J)+PL(3,I)*PL(3,J))/
       (PL(6,I)*PL(6,J))))
     */

    double e1 = jet1.E();
    double e2 = jet2.E();
    TVector3 mom1 = jet1.Vect();
    TVector3 mom2 = jet2.Vect();
    double val = 2*e1*e2*max(0., 1-(mom1.Dot(mom2))/(mom1.Mag()*mom2.Mag()));
    return val/Evis2;
  }

  double JetFinder::funcJadeE(Jet& jet1, Jet& jet2, double Evis2){
    /*
       E(I,J) = MAX(0.,(PL(4,I)+PL(4,J))**2-(PL(1,I)+PL(1,J))**2-
       (PL(2,I)+PL(2,J))**2-(PL(3,I)+PL(3,J))**2)
     */

    double e1 = jet1.E();
    double e2 = jet2.E();
    TVector3 mom1 = jet1.Vect();
    TVector3 mom2 = jet2.Vect();
    double val = max( 0.,  (e1+e2)*(e1+e2) - (mom1+mom2).Mag2() );
    return val/Evis2;
  }

  JetFinder::JetFinder(const JetConfig& cfg) : _cfg(cfg) {
    _Yfunc = JetFinder::funcDurham;
  }
	void JetFinder::Configure(const JetConfig& cfg){_cfg = cfg;}

  vector<Jet*> JetFinder::run(TrackVec &tracks) {
    vector<Jet*> jets;
    for (TrackVecIte iter = tracks.begin(); iter != tracks.end(); ++iter) {
      const Track* trk = *iter;
      Jet* j = new Jet(trk);
      jets.push_back(j);
    }
    return run(jets);
  }

  vector<Jet*> JetFinder::run(TrackVec &tracks,NeutralVec &neutrals,double *pymin, int ynjetmax) {
    vector<Jet*> jets;
    for (TrackVecIte iter = tracks.begin(); iter != tracks.end(); ++iter) {
      const Track* trk = *iter;
      Jet* j = new Jet(trk);
      jets.push_back(j);
    }
    for (NeutralVecIte iter = neutrals.begin(); iter != neutrals.end(); ++iter) {
      const Neutral* neut = *iter;
      Jet* j = new Jet(neut);
      jets.push_back(j);
    }
    return run(jets,pymin, ynjetmax);
  }

  vector<Jet*> JetFinder::run(TrackVec &tracks, NeutralVec &neutrals, VertexVec &vertices_, double *pymin, bool findmu, int ynjetmax)
	{
		vector<Jet*> prejets = prerun(tracks, neutrals, vertices_, findmu);
		return run(prejets, pymin, ynjetmax);
	}

  vector<Jet*> JetFinder::prerun(TrackVec &tracks, NeutralVec &neutrals, VertexVec &vertices_, bool findmu, int *nVertexJets)
	{

		// angle to combine two vertices
		const double vertexcombineangle = 0.2;
		const double vertexcombineanglelepton = 0.3;
		// angle to associate particles to vertex
		const double vertexassocangle = 0.2;

		_Yfunc = JetFinder::funcDurhamVertex;

		// associate tracks/neutrals to vertex jets
		vector<bool> usedTracks(tracks.size(), false);
		vector<bool> usedNeutrals(neutrals.size(), false);

		// copy vertices
		vector<const Vertex *> vertices = vertices_;

		if(findmu){
			// lepton tag
			for(unsigned int i=0;i<tracks.size();i++){
				const Track *tr = tracks[i];
				double dist = sqrt(tr->getD0() *tr->getD0() + tr->getZ0() * tr->getZ0());
				double sigd0 = fabs(tr->getD0()) / sqrt(tr->getCovMatrix()[tpar::d0d0]);
				double sigz0 = fabs(tr->getZ0()) / sqrt(tr->getCovMatrix()[tpar::z0z0]);
				double mudep = tracks[i]->getCaloEdep()[tpar::yoke];
				double ecaldep = tracks[i]->getCaloEdep()[tpar::ecal];
				double hcaldep = tracks[i]->getCaloEdep()[tpar::hcal];

				/*
				if(mudep>0 && tr->getMcp())
					cout << "Muon selection: " << tr->getMcp()->getPDG() << " " << sigd0 << " " << sigz0 << " " << tr->E() << " " << mudep << " " << ecaldep << " " << hcaldep << endl;
				 */
				// muon selection criteria
//				if(dist > 0.1 && dist < 2.0 && tr->E() > 5. && mudep > 0.05 && ecaldep < 1. && hcaldep > 1.5 && hcaldep < 4.){
				if((sigd0 > 5. || sigz0 > 5.) && dist < 5. && mudep > 0.05 && ecaldep < 1. && hcaldep > 1.5 && hcaldep < 5.){
					// treated as a vertex
					float cov[6] = {0,0,0,0,0,0};
					Vertex *fakevtx = new Vertex(0,1,tr->Px()/tr->E(), tr->Py()/tr->E(), tr->Pz()/tr->E(), cov, false);
					fakevtx->add(tr);

					//cout << "A fake vertex created by mu-candidate: " << tr->Px() << " " << tr->Py() << " " << tr->Pz() << " " << tr->E() << endl;

					vertices.push_back(fakevtx);
					usedTracks[i] = true;
				}
			}
		}

		vector<Jet*> jets;
		vector<bool> usedVertices(vertices.size(),false);

		for(unsigned int i=0;i<vertices.size();i++){
			if(!usedVertices[i]){
				Jet *jet = new Jet(vertices[i]);
				jets.push_back(jet);
				TVector3 pos = vertices[i]->getPos();

				for(unsigned int j=i+1; j<vertices.size();j++){
					if(usedVertices[j])continue; // BUGFIX at 120413

					TVector3 pos2 = vertices[j]->getPos();

					//if(vertices[j]->getTracks().size() == 1){ cout << "Angle (" << j << ", " << i << ") : " << pos.Angle(pos2) << endl; }

					double angle = pos.Angle(pos2);
					if(angle < vertexcombineangle || (vertices[j]->getTracks().size() == 1 && angle < vertexcombineanglelepton)){
						jet->add(vertices[j]);
						usedVertices[j] = true;
						//cout << "JetFinder: Vertex " << i << " and " << j << " are combined." << endl;
					}
				}
			}
		}

		while(jets.size() > (unsigned int)_cfg.nJet){
			//cout << "JetFinder: remained # jets with vertex is " << jets.size() << ", performing further combination..." << endl;

			double minimumangle = 3.15;
			int minimumi = -1;
			int minimumj = -1;

			for(unsigned int i=0;i<jets.size();i++){
				for(unsigned int iv=0;iv<jets[i]->getVertices().size();iv++){
					TVector3 pos = jets[i]->getVertices()[iv]->getPos();

					for(unsigned int j=i+1; j<jets.size();j++){
						for(unsigned int jv=0;jv<jets[j]->getVertices().size();jv++){
							TVector3 pos2 = jets[j]->getVertices()[jv]->getPos();
							// if lepton, criteria are loosen
							double angle = pos.Angle(pos2) - (jets[j]->getVertices()[jv]->getTracks().size() == 1) * (vertexcombineanglelepton - vertexcombineangle);
							if(angle < minimumangle){
								minimumangle = angle;
								minimumi = i;
								minimumj = j;
							}
						}
					}
				}
			}
			for(unsigned int jv=0;jv<jets[minimumj]->getVertices().size();jv++){
				jets[minimumi]->add(jets[minimumj]->getVertices()[jv]);
			}
			delete jets[minimumj];
			jets.erase(jets.begin() + minimumj);
			//cout << "JetFinder: Jet " << minimumi << " and " << minimumj << " are combined." << endl;
		}

		if(nVertexJets)
			*nVertexJets = (int)jets.size();

		for(unsigned int i=0;i<tracks.size();i++){
			if(usedTracks[i])continue;

			Jet *jetToAssoc = 0;
			double minimumangle = vertexassocangle;

			for(unsigned int j=0;j<jets.size();j++){
				Jet *jet = jets[j];

				for(unsigned int k=0;k<jet->getVertices().size();k++){
					const Vertex *vtx = jet->getVertices()[k];
					TVector3 vpos = vtx->getPos();

					double angle = vpos.Angle(tracks[i]->Vect());
					if(angle < minimumangle){
						jetToAssoc = jet;
						minimumangle = angle;
					}
				}
			}

			if(jetToAssoc){
				jetToAssoc->add(tracks[i]);
				usedTracks[i] = true;
			}
		}

		for(unsigned int i=0;i<neutrals.size();i++){
			Jet *jetToAssoc = 0;
			double minimumangle = vertexassocangle;

			for(unsigned int j=0;j<jets.size();j++){
				Jet *jet = jets[j];

				for(unsigned int k=0;k<jet->getVertices().size();k++){
					const Vertex *vtx = jet->getVertices()[k];
					TVector3 vpos = vtx->getPos();

					double angle = vpos.Angle(neutrals[i]->Vect());
					if(angle < minimumangle){
						jetToAssoc = jet;
						minimumangle = angle;
					}
				}
			}

			if(jetToAssoc){
				jetToAssoc->add(neutrals[i]);
				usedNeutrals[i] = true;
			}
		}

#if 0
		const double daumassth = 10.0;
		const double pmassth = 20.0;

		// combine vertices by chain
		// CAUTION: vertices must be sorted by distance to IP
		if(jets.size()>=2){
			for(unsigned int i=0;i<jets.size()-1;i++){
				Vertex *v1 = jets[i]->getVertices()[0];
				for(unsigned int j=i+1;j<jets.size();j++){
					Vertex *v2 = jets[j]->getVertices()[0];
					Vertex *vp, *vd;
					if(v1->getPos().Mag() < v2->getPos().Mag()){
						vp = v1; vd = v2;
					}else{
						vp = v2; vd = v1;
					}
					MCParticle *mcp, *mcd;
					mcp=vp->getMcp();
					mcd=vd->getMcp();
					if(!mcp || !mcd){cout << "Stable MCParticle is not found in vertex!" << endl;continue;}
					cout << "PDG: " << mcp->getPDG() << ", " << mcd->getPDG() << endl;
					if(mcp!=mcd && !mcp->isParent(mcd) && !mcd->isParent(mcp)){					
					     continue;
					}
					   /*
					TVector3 vdau = vd->getPos() - vp->getPos();
					double ptdau, pdau;
					double daumass = vd->getVertexMass(0, &vdau, 1.87, &ptdau, &pdau);

					cout << "Vertex mass of c: parent: " << i << ", daughter: " << j << ": " << daumass << endl;

					if(daumass>daumassth || ptdau/pdau > 0.4)continue;

					double pt, p;
					double pmass = vp->getVertexMass(vd,0,1.87, &pt,&p);
					cout << "Vertex mass of b: " << pmass << ", pt/p: " << pt/p << endl;
					if(pmass>pmassth || pt/p > 0.4) continue;
*/
					cout << "Vertices combined." << endl;
					cout << "(" << vp->getX() << ", " << vp->getY() << ", " << vp->getZ() << "),";
					cout << "(" << vd->getX() << ", " << vd->getY() << ", " << vd->getZ() << ")" << endl;

					if(vp==v1){
						jets[i]->add(vd);
						delete jets[j];
					}else{
						jets[j]->add(vd);
						delete jets[i];
						jets[i] = jets[j];
					}
					jets.erase(jets.begin() + j);
					break;
				}
			}
		}

#endif
/*
		// visible energy
		double evis = 0.;
		for(unsigned int i=0;i<tracks.size();i++)evis += tracks[i]->E();
		for(unsigned int i=0;i<neutrals.size();i++)evis += neutrals[i]->E();
		double evis2 = evis * evis;

		vector<Track *> rtracks;
		vector<Neutral *> rneutrals;
		// pre-combining
		double yth = 1e-4;
		for(unsigned int i=0;i<tracks.size();i++){
			double ymin = yth;
			int jymin = -1;
			for(unsigned int j=0;j<jets.size();j++){
				double y = (*_Yfunc)(*jets[j], Jet(tracks[i]), evis2);
				if(y < ymin){
					ymin = y;
					jymin = j;
				}
			}
			if(jymin>=0)
				jets[jymin]->add(tracks[i]);
			else
				rtracks.push_back(tracks[i]);
		}
		for(unsigned int i=0;i<neutrals.size();i++){
			double ymin = yth;
			int jymin = -1;
			for(unsigned int j=0;j<jets.size();j++){
				double y = (*_Yfunc)(jets[j], Jet(neutrals[i]), evis2);
				if(y < ymin){
					ymin = y;
					jymin = j;
				}
			}
			if(jymin>=0)
				jets[jymin]->add(neutrals[i]);
			else
				rneutrals.push_back(neutrals[i]);
		}
*/

		for (unsigned int i=0;i < tracks.size();i++){
			if(usedTracks[i] == false){
				Jet* j = new Jet(tracks[i]);
				jets.push_back(j);
			}
    }
		for (unsigned int i=0;i < neutrals.size();i++){
			if(usedNeutrals[i] == false){
				Jet* j = new Jet(neutrals[i]);
				jets.push_back(j);
			}
    }

		// TODO: delete fakevtx
		return jets;

	}


  vector<Jet*> JetFinder::run(vector<Jet*> jets, double *pymin, int ynjetmax) {

    int nPartons = jets.size();

    double Evis(0.);

    for (vector<Jet*>::iterator iter = jets.begin(); iter != jets.end(); ++iter) {
      Jet* jet = *iter;
      Evis += jet->E();
    }
    double Evis2 = Evis*Evis;

    int njet = nPartons;

    TMatrixTSym<double> matY(njet);

    for (int i1=0; i1<njet; ++i1) {
      for (int i2=i1+1; i2<njet; ++i2) {
        double val = (*_Yfunc)(*jets[i1],*jets[i2],Evis2);
        matY[i1][i2] = val;
      }
    }

    while ( njet > _cfg.nJet ) {

      double Ymin = 1e8;
      int imin=-1;
      int jmin=-1;

      /* compute upper-right triangle of matY
       * for the values i2 > i1
       */
      for (int i1=0; i1<njet; ++i1) {
        for (int i2=i1+1; i2<njet; ++i2) {
          if (matY[i1][i2] < Ymin) {
            Ymin = matY[i1][i2];
            imin = i1;
            jmin = i2;
          }
        }
      }

      if(imin == -1 || jmin == -1){printf("Further combination cannot be done. njet = %d\n", njet);break;}
//      assert(imin != -1);
//      assert(jmin != -1);

			if(pymin && (ynjetmax >= njet))
				pymin[njet-1] = Ymin;

      // check Ymin against Ycut
      if (_cfg.Ycut > 0 && Ymin > _cfg.Ycut) {
        printf("Ycut reached: %e\n", Ymin);
        break;
      }

      if (imin > jmin) {
        int tmp = jmin;
        jmin = imin;
        imin = tmp;
      }

      jets[imin]->add(*jets[jmin]);
      delete jets[jmin];

      jets[jmin] = jets.back();
      jets.pop_back();

      --njet;

      //printf("combining %d with %d, val=%e\n",imin,jmin,Ymin);

      // recompute matY (only the entries that changed from reordering)
      for (int ii=0; ii<njet; ++ii) {
        if (ii != imin) {
          int i = min(ii,imin);
          int j = max(ii,imin);
          matY[i][j] = (*_Yfunc)(*jets[i],*jets[j],Evis2);
        }
        if (ii != jmin) {
          int i = min(ii,jmin);
          int j = max(ii,jmin);
          matY[i][j] = matY[ii][njet];
        }
      }
    }

    // order jets by largest energy first
    sort(jets.begin(),jets.end(),Jet::sort_by_energy);

    /*

       char buf[100];
    // debug printing
    for (int i=0; i<njet && i<10; ++i) {
    double e = jets[i].getEnergy();
    const double *p = jets[i].getTVector3();
    snprintf(buf,100,"(%d) [%-.2f, %-.2f, %-.2f, %-.2f]\n",i,e,p[0],p[1],p[2]);
    SLM << buf;
    }
     */


    return jets;
  }

	/*
  CheatedJetFinder::CheatedJetFinder(const JetConfig& cfg) : _cfg(cfg) {
  }

  vector<MCParticle*> CheatedJetFinder::originalPartons( MCParticleVec & mcps ) {
    vector<MCParticle*> quarks;
    vector<MCParticle*> quarksAfterGluonEmission;

    for (MCParticleVecIte it = mcps.begin(); it != mcps.end(); ++it) {
      const MCParticle* mcp = *it;
      if (mcp->getPDG() == 92) {
        printf(" color singlet: %d\n", mcp->getPDG());
      }
    }

    return quarksAfterGluonEmission;
  }

  vector<Jet*> CheatedJetFinder::run(Event* event) {
    vector<Jet*> jets;

    //const vector<Track*>& tracks  = event->getTracks();
    //const vector<Neutral*>& neutrals = event->getNeutrals();
    //const vector<MCParticle*>& mcps = event->getMCParticles();

    //vector<const Track&> tracksLeft(tracks);
    //vector<const Neutral&> neutralsLeft(neutrals);
    //vector<const MCParticle&> partons = originalPartons(mcps);

    return jets;
  }
	*/

	Jet* convertJetVertex(const Jet* jet) {
		Jet* newjet = new Jet();

		TrackVec & tracks  = jet->getTracks();
		NeutralVec & neutrals = jet->getNeutrals();
		VertexVec & vertices = jet->getVertices();

		for (TrackVecIte iter = tracks.begin(); iter != tracks.end(); ++iter) {
			newjet->add(*iter);
		}

		for (NeutralVecIte iter = neutrals.begin(); iter != neutrals.end(); ++iter) {
			newjet->add(*iter);
		}

		for (VertexVecIte iter = vertices.begin(); iter != vertices.end(); ++iter) {
			const Vertex* vertex = *iter;
			TrackVec & tracks2  = vertex->getTracks();
			for (TrackVecIte iter = tracks2.begin(); iter != tracks2.end(); ++iter) {
				newjet->add(*iter);
			}
		}

		return newjet;
	}
}
