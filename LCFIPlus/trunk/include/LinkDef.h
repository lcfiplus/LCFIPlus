#ifndef EVENTDISPLAY_LINKDEF_H_
#define EVENTDISPLAY_LINKDEF_H_

#ifdef __CINT__

//#pragma link C++ class lcfiplus::TrackData+;
//#pragma link C++ class std::vector<lcfiplus::TrackData>+;
//#pragma link C++ class lcfiplus::NeutralData+;
//#pragma link C++ class std::vector<lcfiplus::NeutralData>+;
//#pragma link C++ class lcfiplus::MCParticleData+;
//#pragma link C++ class std::vector<lcfiplus::MCParticleData>+;
#pragma link C++ class lcfiplus::Track+;
#pragma link C++ class std::vector<lcfiplus::Track*>+;
#pragma link C++ class lcfiplus::Neutral+;
#pragma link C++ class std::vector<lcfiplus::Neutral*>+;
#pragma link C++ class lcfiplus::MCParticle+;
#pragma link C++ class std::vector<lcfiplus::MCParticle*>+;
#pragma link C++ class lcfiplus::Vertex+;
#pragma link C++ class std::vector<lcfiplus::Vertex*>+;
#pragma link C++ class lcfiplus::MCVertex+;
#pragma link C++ class std::vector<lcfiplus::MCVertex*>+;
#pragma link C++ class lcfiplus::Jet+;
#pragma link C++ class std::vector<lcfiplus::Jet*>+;

#pragma link C++ class lcfiplus::LCIOStorer;
#pragma link C++ class lcfiplus::TreeStorer;
#pragma link C++ class lcfiplus::EventStore;
#pragma link C++ class lcfiplus::EventNavigator;

#pragma link C++ class lcfiplus::LcfiplusAlgorithm;
//#pragma link C++ class lcfiplus::FlavtagAlgorithm;
#pragma link C++ class lcfiplus::BuildUpVertex;
#pragma link C++ class lcfiplus::JetClustering;
#pragma link C++ class lcfiplus::PrimaryVertexFinder;
#pragma link C++ class lcfiplus::FlavorTag;
#pragma link C++ class lcfiplus::MakeNtuple;
#pragma link C++ class lcfiplus::TrainMVA;

//#pragma link C++ class lcfiplus::Event;

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
