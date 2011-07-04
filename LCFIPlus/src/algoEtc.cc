#include "flavtag.h"
#include "algoEtc.h"

#include "TRandom.h"
#include <math.h>

namespace flavtag{
namespace algoEtc{

void makeBeamTracks(Track *&t1, Track *&t2)
{
   // beam crossing = 14 mrad
    float beamtd = tan( 0.5*(3.1415926 - 14e-3) );
    float beamsizeX = 639e-6; // 639 nm converted to mm
    //float beamsizeY = 5.7e-6; // 5.7 nm converted to mm
    // size(d)/size(z) = tan(7e-3) -> size(z) = size(d)/tan(7e-3)
    float beamsizeZ = beamsizeX / tan(0.5 * 14e-3);

    float d0rand(0), z0rand(0);

    //d0rand = (_rand.Rndm()*2-1)*beamsizeX;
    //z0rand = (_rand.Rndm()*2-1)*beamsizeZ;
    d0rand = gRandom->Gaus(0,beamsizeX);
    z0rand = gRandom->Gaus(0,beamsizeZ);

		t1 = new Track;
		t2 = new Track;

		float par[tpar::parN];
		float cov[tpar::covN];

    t1->setId(1000001);
		
    par[tpar::d0] = d0rand;
    par[tpar::z0] = z0rand;
    par[tpar::om] = 0.3*3.5/250./beamtd*1000.;
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

    //d0rand = (_rand.Rndm()*2-1)*beamsizeX;
    //z0rand = (_rand.Rndm()*2-1)*beamsizeZ;
    d0rand = gRandom->Gaus(0,beamsizeX);
    z0rand = gRandom->Gaus(0,beamsizeZ);

		t2->setId(1000002);

    par[tpar::d0] = d0rand;
    par[tpar::z0] = z0rand;
    par[tpar::om] = 0.3*3.5/250./beamtd*1000.;
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

void makeBeamVertex(Vertex *&vtx)
{
	float beamsizeX = 639e-6; // 639 nm converted to mm
	float beamsizeY = 5.7e-6; // 5.7 nm converted to mm
	float beamsizeZ = beamsizeX / tan(0.5 * 14e-3);

	float cov[6];
	cov[Vertex::xx] = pow(beamsizeX,2);
	cov[Vertex::xy] = 0;
	cov[Vertex::xz] = 0;
	cov[Vertex::yy] = pow(beamsizeY,2);
	cov[Vertex::yz] = 0;
	cov[Vertex::zz] = pow(beamsizeZ,2);

	vtx = new Vertex(0,1,gRandom->Gaus(0, beamsizeX), gRandom->Gaus(0,beamsizeY), gRandom->Gaus(0,beamsizeZ),cov);
}

}}
