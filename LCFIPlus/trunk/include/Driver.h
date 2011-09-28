// Driver.h

#ifndef Driver_h
#define Driver_h 1

#include "lcfiplus.h"
#include "JetFinder.h"

extern void processEvents(const char* input, const char* output, int nStart, int nEnd);
extern void eventDisplay(const char *infile = "share/test.slcio", int start = 0);

extern vector<int> findMcJetFlavor(vector<Jet*> jets, vector<MCParticle*> mcps);
extern vector<MCVertex*> findMcVertex(vector<MCParticle*> mcps);
extern void matchMcVertex( const Event& evt, vector<MCVertex*>& vtxList, map<MCVertex*,int>& table , bool vertexing = true);
extern void matchMcVertexJet(const Event& evt, const vector<MCVertex*>& vtxList, map<MCVertex*,int>& table, const vector<Jet*>& jets);
extern void matchMcVertexReco(const Event& evt, const vector<MCVertex*>& vtxList, map<MCVertex*,int>& table, Vertex* vertex );

extern vector<Track*> findSingleTracks(const Event& evt, const Jet& jet, const vector<lcfiplus::Vertex*>& vtxList);

#endif
