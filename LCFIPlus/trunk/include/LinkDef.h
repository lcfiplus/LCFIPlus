#ifndef EVENTDISPLAY_LINKDEF_H_
#define EVENTDISPLAY_LINKDEF_H_

#ifdef __CINT__

#pragma link C++ class flavtag::TrackData+;
#pragma link C++ class std::vector<flavtag::TrackData>+;
//#pragma link C++ class flavtag::NeutralData+;
//#pragma link C++ class std::vector<flavtag::NeutralData>+;
//#pragma link C++ class flavtag::MCParticleData+;
//#pragma link C++ class std::vector<flavtag::MCParticleData>+;
#pragma link C++ class flavtag::Track+;
#pragma link C++ class std::vector<flavtag::Track*>+;
#pragma link C++ class flavtag::Neutral+;
#pragma link C++ class std::vector<flavtag::Neutral*>+;
#pragma link C++ class flavtag::MCParticle+;
#pragma link C++ class std::vector<flavtag::MCParticle*>+;
#pragma link C++ class flavtag::Vertex+;
#pragma link C++ class std::vector<flavtag::Vertex*>+;
#pragma link C++ class flavtag::MCVertex+;
#pragma link C++ class std::vector<flavtag::MCVertex*>+;
#pragma link C++ class flavtag::Jet+;
#pragma link C++ class std::vector<flavtag::Jet*>+;

#pragma link C++ class flavtag::LCIOStorer;
#pragma link C++ class flavtag::TreeStorer;
#pragma link C++ class flavtag::EventStore;
#pragma link C++ class flavtag::EventNavigator;

#pragma link C++ class flavtag::FlavtagAlgorithm;
#pragma link C++ class flavtag::BuildUpVertex;
#pragma link C++ class flavtag::JetClustering;
#pragma link C++ class flavtag::MakeNtuple;
#pragma link C++ class flavtag::TrainMVA;

//#pragma link C++ class flavtag::Event;

#pragma link C++ function lcioTest;
#pragma link C++ function lcioToTree;
#pragma link C++ function treeTest;
#pragma link C++ function processEvents;
#pragma link C++ function eventDisplay;
#pragma link C++ function TearDownVertexing;
#pragma link C++ function testSuehara;
#pragma link C++ function checkMCTearDown;
#pragma link C++ function testSueharaVertex;
#pragma link C++ function bbhhAnalysis;
#pragma link C++ function bbhhAnalysis2;
#pragma link C++ function bbhhAnalysis3;
#pragma link C++ function pointTest;
#pragma link C++ function helixTest;
#pragma link C++ function helixVarianceTest;
#pragma link C++ function simpleAnalysis;
#pragma link C++ function VertexAnalysis110214;
#pragma link C++ function VertexAnalysis110215;
#pragma link C++ function TrackDist110218;
#pragma link C++ function outVertex;
#pragma link C++ function lcioToLcio;

#pragma link C++ class VertexAnalysis+;
#pragma link C++ class MatchMCRecoVertex+;
#pragma link C++ class std::vector<MatchMCRecoVertex*>+;

#endif

#endif
