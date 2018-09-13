#ifndef LINKDEF_H_
#define LINKDEF_H_

#ifdef __CINT__

#pragma link C++ class lcfiplus::Track+;
#pragma link C++ class std::vector<const lcfiplus::Track*>+;
#pragma link C++ class std::vector<lcfiplus::Track*>+;
#pragma link C++ class lcfiplus::Neutral+;
#pragma link C++ class std::vector<const lcfiplus::Neutral*>+;
#pragma link C++ class std::vector<lcfiplus::Neutral*>+;
#pragma link C++ class lcfiplus::MCParticle+;
#pragma link C++ class std::vector<const lcfiplus::MCParticle*>+;
#pragma link C++ class std::vector<lcfiplus::MCParticle*>+;
#pragma link C++ class lcfiplus::MCColorSinglet+;
#pragma link C++ class std::vector<const lcfiplus::MCColorSinglet*>+;
#pragma link C++ class std::vector<lcfiplus::MCColorSinglet*>+;
#pragma link C++ class lcfiplus::Vertex+;
#pragma link C++ class std::vector<const lcfiplus::Vertex*>+;
#pragma link C++ class std::vector<lcfiplus::Vertex*>+;
#pragma link C++ class lcfiplus::MCVertex+;
#pragma link C++ class std::vector<const lcfiplus::MCVertex*>+;
#pragma link C++ class std::vector<lcfiplus::MCVertex*>+;
#pragma link C++ class lcfiplus::Jet+;
#pragma link C++ class std::vector<const lcfiplus::Jet*>+;
#pragma link C++ class std::vector<lcfiplus::Jet*>+;

#pragma link C++ class lcfiplus::LCIOStorer;
#pragma link C++ class lcfiplus::TreeStorer;
#pragma link C++ class lcfiplus::EventStore;
#pragma link C++ class lcfiplus::EventStoreObserver;

#pragma link C++ class lcfiplus::Algorithm;
#pragma link C++ class lcfiplus::BuildUpVertex;
#pragma link C++ class lcfiplus::JetClustering;
#pragma link C++ class lcfiplus::JetVertexRefiner;
#pragma link C++ class lcfiplus::PrimaryVertexFinder;
#pragma link C++ class lcfiplus::FlavorTag;
#pragma link C++ class lcfiplus::MakeNtuple;
#pragma link C++ class lcfiplus::TrainMVA;
#pragma link C++ class lcfiplus::ReadMVA;
#pragma link C++ class lcfiplus::TrackNtuple;
#pragma link C++ class lcfiplus::VertexMassRecovery+;
#pragma link C++ class lcfiplus::VertexNtuple;

#pragma link C++ class lcfiplus::ZHHAlgo;
#pragma link C++ class lcfiplus::TestAlgo;
#pragma link C++ class lcfiplus::FlavtagReader;
#pragma link C++ class lcfiplus::VertexAnalysis;
#pragma link C++ class lcfiplus::TestAlgoV0;

#endif

#endif
