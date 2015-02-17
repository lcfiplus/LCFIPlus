// VertexFitterLCFI.h

#ifndef VertexFitterSimple_h
#define VertexFitterSimple_h 1

#include "lcfiplus.h"
#include "geometry.h"

using namespace std;
using namespace lcfiplus;

namespace lcfiplus {

	template<class Iterator>
		class VertexFitterSimple {
		public:
			Vertex * operator() (Iterator tracksBegin, Iterator tracksEnd, Vertex *pointConstraint = 0, bool pointInitialOnly = false)
			{
				bool verbose = false;

				GeometryHandler *gh = GeometryHandler::Instance();
				if(pointConstraint){
					Point *ip = new Point(pointConstraint);
					vector<PointBase *> tracks;
					if(!pointInitialOnly)
						tracks.push_back(ip);
					int ntracks = 0;
					for(Iterator it = tracksBegin; it != tracksEnd;it++,ntracks++){
						tracks.push_back(new Helix(*it));
					}
					if(verbose)
						cout << "VertexFitterSimple: number of tracks is " << ntracks << endl;

					TVector3 initial = pointConstraint->getPos();
					Point *result = new Point;
					double chi2 = -gh->PointFit(tracks, initial, result);

					TVector3 vresult = result->GetPos();
					double cov[6];
					cov[Vertex::xx] = result->GetErr(0,0);
					cov[Vertex::xy] = result->GetErr(0,1);
					cov[Vertex::xz] = result->GetErr(0,2);
					cov[Vertex::yy] = result->GetErr(1,1);
					cov[Vertex::yz] = result->GetErr(1,2);
					cov[Vertex::zz] = result->GetErr(2,2);

					if(verbose) {
						cout << "Vertex cov matrix:" << endl;
						cout << scientific << cov[Vertex::xx] << "  ";
						cout << cov[Vertex::yy] << "  ";
						cout << cov[Vertex::zz] << "  ";
						cout << cov[Vertex::xy] << "  ";
						cout << cov[Vertex::yz] << "  ";
						cout << cov[Vertex::xz] << endl << fixed;
					}

					if(verbose){
						cout << "VertexFitterSimple: vertex position is " << endl;
						vresult.Print();
					}

					Vertex *vtx = new Vertex(chi2, TMath::Prob(chi2, ntracks*2-3), vresult.x(), vresult.y(), vresult.z(), cov, false);
					for(Iterator it = tracksBegin; it != tracksEnd; it++,ntracks++){
						Helix hel(*it);
						double ll = hel.LogLikelihood(vresult); // need to incorporate vertex error??

						if(verbose)
							cout << "VertexFitterSimple: track loglikelihood is " << ll << endl;

						vtx->add(*it, -ll);
					}
					for(vector<PointBase *>::iterator it = tracks.begin(); it != tracks.end(); it++){
						delete *it;
					}

					delete result;
					return vtx;
				}

				// without point constraint
				vector<Helix *> tracks;
				int ntracks = 0;
				for(Iterator it = tracksBegin; it != tracksEnd;it++,ntracks++){
					tracks.push_back(new Helix(*it));
				}
				if(verbose)
					cout << "VertexFitterSimple: number of tracks is " << ntracks << endl;

				Point *result = new Point;
				double chi2 = -gh->HelixPointFit(tracks, result);

				TVector3 vresult = result->GetPos();
				double cov[6];
				cov[Vertex::xx] = result->GetErr(0,0);
				cov[Vertex::xy] = result->GetErr(0,1);
				cov[Vertex::xz] = result->GetErr(0,2);
				cov[Vertex::yy] = result->GetErr(1,1);
				cov[Vertex::yz] = result->GetErr(1,2);
				cov[Vertex::zz] = result->GetErr(2,2);

				if(verbose){
					cout << "VertexFitterSimple: vertex position is " << endl;
					vresult.Print();
				}

				Vertex *vtx = new Vertex(chi2, (ntracks > 1 ? TMath::Prob(chi2, ntracks*2-3) : 1), vresult.x(), vresult.y(), vresult.z(), cov, false);
				for(Iterator it = tracksBegin; it != tracksEnd;it++, ntracks++){
					Helix hel(*it);
					double ll = hel.LogLikelihood(vresult); // need to incorporate vertex error??
					if(verbose)
						cout << "VertexFitterSimple: track loglikelihood is " << ll << endl;

					vtx->add(*it, -ll);
				}
				for(vector<Helix *>::iterator it = tracks.begin(); it != tracks.end(); it++){
					delete *it;
				}

				delete result;
				return vtx;
			}

			double getChi2(const Vertex *vtx, const Track *trk, int mode=1){
				// 110510 suehara for IPassoc study
				if(mode == 0){
					// mode 0: no fit at all
					Helix hel(trk);
					TVector3 v = vtx->getPos();
					return -hel.LogLikelihood(v);
				}
				else if(mode == 1){
					// mode 1: vertex treated as errored point

					Point ptVtx(vtx);
					Helix hel(trk);
	
					vector<PointBase *> vpt;
					vpt.push_back(&ptVtx);
					vpt.push_back(&hel);
	
					GeometryHandler *gh = GeometryHandler::Instance();
					return -gh->PointFit(vpt, vtx->getPos(), 0);
				}else{
					// mode 2: vertex treated as track list
					vector<PointBase *> vpt;
					for(unsigned int n=0;n<vtx->getTracks().size();n++){
						vpt.push_back(new Helix(vtx->getTracks()[n]));
					}
					vpt.push_back(new Helix(trk));

					GeometryHandler *gh = GeometryHandler::Instance();
					Point ret;
					double ll = -gh->PointFit(vpt, vtx->getPos(), &ret);

					for(unsigned int n=0;n<vpt.size();n++){
						delete vpt[n];
					}
					if(mode == 2){
						return ll;
					}else{
						Helix hel(trk);
						return -hel.LogLikelihood(ret.GetPos());
					}
				}
			}

		};


		typedef VertexFitterSimple<vector<const Track *>::const_iterator > VertexFitterSimple_V;
		typedef VertexFitterSimple<list<const Track *>::const_iterator > VertexFitterSimple_L;
}

#endif
