// testproc.h

#ifndef testproc_h
#define testproc_h 1

#include "lcfiplus.h"
#include "TNtuple.h"
#include "TNtupleD.h"

namespace lcfiplus{

	class ZHHAlgo : public Algorithm
	{
	public:
		struct data{
			int mchdecaypdg[2];
			int mcnb;

			double thrust;
			double thaxis[3];
			double ycuts[10];

			double bcat[6];
			double btag[6];
			double ctag[6];
			double ejet[6];
			double pxjet[6];
			double pyjet[6];
			double pzjet[6];
			double ntrjet[6];

			double mass[15];
			double ntrjetmin;
			double pmiss[3];
			double emiss;

			double bcat4[6];
			double btag4[6];
			double ctag4[6];
			double ejet4[6];
			double pxjet4[6];
			double pyjet4[6];
			double pzjet4[6];
			double ntrjet4[6];

		};

		ZHHAlgo(){}
		virtual ~ZHHAlgo(){}

		static bool sortBtag(const Jet *j1, const Jet *j2){
			double btag1 = j1->getParam("lcfiplus")->get<double>("BTag");
			double btag2 = j2->getParam("lcfiplus")->get<double>("BTag");

			return btag1 > btag2;
		}

		void init(Parameters *param);
		void process();
		void end();

		ClassDef(ZHHAlgo,1);

	private:
		TFile *_file;

		string _jetname;
		string _jetname4;

		JetVec * _jets; //!
		JetVec * _jets4; //!
		TTree *_tree;

		data _d;

	};

	class TestAlgo : public Algorithm
	{
	public:
		TestAlgo(){}
		virtual ~TestAlgo(){}

		void init(Parameters *param);
		void process();
		void end();

		ClassDef(TestAlgo,1);

	private:
		TNtuple *_ntJet2;
		TNtuple *_nbJet;
		TFile *_file;

		string _v0vtxname;
		string _privtxname;
		string _jetname;
		bool _bbhh;

		VertexVec * _vertices; //!
		VertexVec * _v0vertices; //!
		JetVec * _jets; //!

		// for old version

		string _vtxname;
		int _vtxsel;
		int _refine;
		TNtupleD *_ntJet;

	};

	class TestAlgoV0 : public Algorithm {
	public:
		TestAlgoV0(){}
		virtual ~TestAlgoV0(){}
		void init(Parameters *param);
		void process();
		void end();
		ClassDef(TestAlgoV0,1);
	private:
		TTree *_ntp;
		TFile *_file;
		VertexVec * _vertices; //!
		string _vtxname;

		struct VtxData {
			double x;
			double y;
			double z;
			double r;
			double cs;
			double phi;
			double chrg;
			double dirdot;
			double dirdot2;
			int ntrk;
			double mks;
			double ml0;
			double mconv;
			double mks2;
			double ml02;
			int v0;
			int ks;
			int l0;
			int conv;
			int mcpdg1;
			int mcpdg2;
			int mcppdg1;
			int mcppdg2;
			int mcpp1;
			int mcpp2;
		};
		VtxData _data;
	};

}

#endif
