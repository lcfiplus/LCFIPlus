// LCIOStorer.h
#ifndef LCIOStorer_h
#define LCIOStorer_h 1

#include "lcfiplus.h"
#include <vector>
#include <string>

#include "TObject.h"

#include "lcio.h"
#include "IO/LCReader.h"
#include "IO/LCWriter.h"
#include "EVENT/LCEvent.h"
#include <EVENT/ReconstructedParticle.h>
#include "EVENT/MCParticle.h"

using namespace std;

namespace lcfiplus {

class LCIOStorer : public TObject, public EventStoreObserver {
 public:
  LCIOStorer(const char* inputfile = NULL, const char* outputfile = NULL, bool autoread = true, bool autowrite = false, const char* outPrefix = 0);
  virtual ~LCIOStorer();

  // in both func the obtained LCEvent is used for following operations
  bool Next(bool autovertex = false, bool autojet = false); // file mode: move to the next event
  void SetEvent(lcio::LCEvent* event); // non-file mode
  void SetColorSinglets(vector<MCParticle*>& mcps, vector<MCColorSinglet*>& mccs);

  // LCIO -> lcfiplus
  // register basic collections to EventStore
  void InitMCPPFOCollections(const char* pfoColName, const char* mcColName, const char* mcpfoColName);
  void InitPFOCollections(const char* pfoColName);

  // initialize misc collections
  void InitVertexCollection(const char* lcioName, const char* flavtagName, bool readnow = true);
  void InitJetCollection(const char* lcioName, const char* flavtagName, bool readnow = true, bool readvtx = true, const char* vtxname = 0);

  /*
  			// register every LCIO collections: used only in standalone
  			void InitVertexCollectionsAuto(lcio::LCEvent *evt);
  			void InitJetCollectionsAuto(lcio::LCEvent *evt);
  */
  // vertices/jets
  void ReadVertices(const char* vtxname, vector<const Vertex*>* lcficol);
  void ReadJets(const char* jetname, vector<const Jet*>* lcficol, const char* vtxrelname = 0);

  // lcfiplus -> LCIO
  void WriteVertices(const char* vertexName, const char* outName = 0, const char* outRPName = 0);
  void WriteVertices(VertexVec* pvvtx, const char* newName, const char* newRPName = 0);
  void WriteJets(const char* jetName, const char* outName = 0, bool writeVertex = true, const char* vtxName = 0, const char* relName = 0);
// 			void ConvertJetWithFlavor(const char *jetName, const char *flavName);

  // convert PID : lcfiplus -> LCIO
  void WritePID(lcio::LCCollection* lciocol, lcio::ReconstructedParticle* lciojet, const lcfiplus::Jet* lcfijet, const char* paramname);
  void WriteAllPIDs(lcio::LCCollection* lciocol, lcio::ReconstructedParticle* lciojet, const lcfiplus::Jet* lcfijet);

  void WriteEvent(); // write to the outputfile
  void AutoConvert(); // auto convert without WriteEvent()

  // helper functions for std::sort()
  //static bool energy_sort_trk(Track *a, Track *b);
  //static bool energy_sort_mc(MCParticle *a, MCParticle *b);
  static bool energy_sort_pfo(lcio::ReconstructedParticle* a, lcio::ReconstructedParticle* b);

  void setReadSubdetectorEnergies(bool flag) {
    _readSubdetectorEnergies = flag;
  }
  void setTrackHitOrdering(int flag) {
    _trackHitOrdering = flag;
  }
  void setUpdateVertexRPDaughters(bool flag) {
    _updateVertexRPDaughters = flag;
  }
  void setIgnoreLackOfVertexRP(bool flag) {
    _ignoreLackOfVertexRP = flag;
  }

  bool getReadSubdetectorEnergies()const {
    return _readSubdetectorEnergies;
  }
  int getTrackHitOrdering()const {
    return _trackHitOrdering;
  }
  bool getUpdateVertexRPDaughters()const {
    return _updateVertexRPDaughters;
  }
  bool getIgnoreLackOfVertexRP()const {
    return _ignoreLackOfVertexRP;
  }

  // callback function from EventStore
  virtual void GetCallback(const char* name, const char* classname);

 private:
  lcio::LCEvent* _event;
  lcio::LCReader* _reader;
  lcio::LCWriter* _writer;

  // collections to import
  map<string, vector<lcfiplus::MCParticle*> *> _importMCPCols;
  map<string, vector<lcfiplus::MCColorSinglet*> *> _importMCCSCols;
  map<string, pair<vector<lcfiplus::Track*> *, vector<lcfiplus::Neutral*> *> >_importPFOCols;
  map<string, pair<string,string> > _importMCPFOLinkCols; // Link, <MC, PFO>
  map<string, vector<const lcfiplus::Vertex*> *> _importVertexCols;
  map<string, vector<const lcfiplus::Jet*> *> _importJetCols;

  // LCIO - lcfiplus relation
  map<lcfiplus::Track*, lcio::ReconstructedParticle*> _trackLCIORel;
  map<lcfiplus::Neutral*, lcio::ReconstructedParticle*> _neutralLCIORel;
  map<lcfiplus::MCParticle*, lcio::MCParticle*> _mcpLCIORel;
  map<const lcfiplus::Vertex*, lcio::Vertex*> _vtxLCIORel;
  map<const lcfiplus::Jet*, lcio::ReconstructedParticle*> _jetLCIORel;

  // reverse direction
  map<lcio::ReconstructedParticle*, lcfiplus::Track*> _trackLCIORel2;
  map<lcio::ReconstructedParticle*, lcfiplus::Neutral*> _neutralLCIORel2;
  map<lcio::MCParticle*, lcfiplus::MCParticle*> _mcpLCIORel2;
  map<lcio::Vertex*, const lcfiplus::Vertex*> _vtxLCIORel2;
  map<lcio::ReconstructedParticle*, const lcfiplus::Jet*> _jetLCIORel2;

  // autosave for output
  bool _autoread;  // for vertex/jet collections: on in LcfiplusProcessor
  bool _autowrite; // for every persistent collections: off in LcfiplusProcessor
  string _savePrefix;

  bool _readSubdetectorEnergies;
  int _trackHitOrdering;
  bool _updateVertexRPDaughters;
  bool _ignoreLackOfVertexRP;

  ClassDef(LCIOStorer,0);
};

}

#endif
