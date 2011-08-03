// Driver.h

#ifndef Suehara_h
#define Suehara_h 1

#include "flavtag.h"
#include "JetFinder.h"

//extern int piuToTree(const char *infile = "input.piu", const char *outfile = "test.root");
extern void treeTest(const char *treefile = "share/test.root");
extern void lcioTest(const char *lciofile = "share/test.slcio");
extern void lcioToTree(const char *infile = "share/test.slcio", const char *outfile = "share/test2.root");
extern void testSuehara();
extern void TearDownVertexing(const char *inputfile = "share/zpole.root", const char *outputfile = "share/zpole_v.root",
	double chi2th_primary = 9.0, double chi2th_secondary = 9.0, bool bonly = false, bool verbose = false);

extern void checkMCTearDown(const char *inputfile = "share/zpole.root", const char *outputfile = "share/zpole_v.root", int nStart=0, int nEnd=100);
extern void testSueharaVertex(const char *inputfile = "share/zpole.root", const char *outputfile = "share/zpole_v.root", int nStart=0, int nEnd=100);
extern void bbhhAnalysis(const char *inputfile, const char *outputfile, int nStart, int nEnd);
extern void bbhhAnalysis2(const char *inputfile, const char *outputfile, int nStart, int nEnd);
extern void bbhhAnalysis3(const char *inputfile, const char *outputfile, int nStart, int nEnd, int vtx=1, int bbhh=1);
extern void pointTest();
extern void helixTest();
extern void helixVarianceTest();
extern void simpleAnalysis(const char *inputfile, const char *outputfile, int nStart, int nEnd);
extern void VertexAnalysis110214(const char *inputfile, const char *outputfile, int nStart, int nEnd);
extern void VertexAnalysis110215(const char *inputfile, const char *outputfile, int nStart, int nEnd, int hh);
extern void TrackDist110218(const char *inputfile, const char *outputfile, int nStart, int nEnd, int hh);
extern vector<Jet *> SueharaJetClustering(Event &evt, int njets);
extern void outVertex(const char *inputfile, const char *outputfile, int nStart, int nEnd, int bbhh, int inclVertex);
extern void lcioToLcio(const char *inputfile, const char *outputfile);

class VertexAnalysis{
	public:
		int event;
		double ipchi2[2];
		double vchi2;
		double trchi2[2];
		int tag[2];
		int success;
		double truevtx[3][2];
		double recovtx[3];

	ClassDefNV(VertexAnalysis,1);
};

class MatchMCRecoVertex{
	public:
		int event;
		MCVertex *mcv;
		vector<Vertex *> recov;
		int recovsize;

	ClassDefNV(MatchMCRecoVertex,2);
};

#endif
