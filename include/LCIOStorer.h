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

namespace lcfiplus{

	class LCIOStorer : public TObject
	{
		public:
			LCIOStorer(const char *inputfile = NULL, const char *outputfile = NULL, bool autoconvert = false, const char *outPrefix = 0);
			virtual ~LCIOStorer();

			// in both func the obtained LCEvent is used for following operations
			bool Next(); // file mode: move to the next event
			void SetEvent(lcio::LCEvent *event); // non-file mode

			// LCIO -> lcfiplus
			// register basic collections to EventStore
			void InitCollections(const char *pfoColName = "PandoraPFOs", const char *mcColName = "MCParticlesSkimmed",
				const char *mcpfoColName = "RecoMCTruthLink",
				const char *trackName = "Tracks", const char *neutralName = "Neutrals", const char *mcpName = "MCParticles");
			void InitCollectionsWithoutMCP(const char *pfoColName = "PandoraPFOs", const char *trackName = "Tracks", const char *neutralName = "Neutrals");

			// initialize misc collections
			void InitVertexCollection(const char *lcioName, const char *flavtagName);
			void InitJetCollection(const char *lcioName, const char *flavtagName);

			// register every LCIO collections
			void InitVertexCollectionsAuto(lcio::LCEvent *evt);
			void InitJetCollectionsAuto(lcio::LCEvent *evt);

			// lcfiplus -> LCIO
			void ConvertVertex(const char *vertexName, const char *outName = 0, const char *outRPName = 0);
			void ConvertJet(const char *jetName, const char *outName = 0, bool extractVertex = false);
// 			void ConvertJetWithFlavor(const char *jetName, const char *flavName);

			// convert PID : lcfiplus -> LCIO
			void WritePID(lcio::LCCollection *lciocol, lcio::ReconstructedParticle *lciojet, const lcfiplus::Jet *lcfijet, const char *paramname);
			void WriteAllPIDs(lcio::LCCollection *lciocol, lcio::ReconstructedParticle *lciojet, const lcfiplus::Jet *lcfijet);

			void WriteEvent(); // write to the outputfile
			void AutoConvert(); // auto convert without WriteEvent()

			// helper functions for std::sort()
			//static bool energy_sort_trk(Track *a, Track *b);
			//static bool energy_sort_mc(MCParticle *a, MCParticle *b);
			static bool energy_sort_pfo(lcio::ReconstructedParticle *a, lcio::ReconstructedParticle *b);

		private:
			lcio::LCEvent *_event;
			lcio::LCReader *_reader;
			lcio::LCWriter *_writer;

			vector<Track*> * _pTracks;
			vector<Neutral*> * _pNeutrals;
			vector<MCParticle*> * _pMCPs;

			string _pfoColName;
			string _mcColName;
			string _mcpfoColName;

			// vertices / jets
			map<string, vector<lcfiplus::Vertex *> *> _importVertexCols;
			map<string, vector<lcfiplus::Jet *> *> _importJetCols;

			// LCIO - lcfiplus relation
/*			map<lcfiplus::Track *, lcio::ReconstructedParticle *> _trackLCIORel;
			map<lcfiplus::Neutral *, lcio::ReconstructedParticle *> _neutralLCIORel;
			map<lcfiplus::MCParticle *, lcio::MCParticle *> _mcpLCIORel;*/
			vector<lcio::ReconstructedParticle *> _trackLCIORel;
			vector<lcio::ReconstructedParticle *> _neutralLCIORel;
			vector<lcio::MCParticle *> _mcpLCIORel;

			// new classes
			map<const lcfiplus::Vertex *, lcio::Vertex *> _vtxLCIORel;
			map<const lcfiplus::Jet *, lcio::ReconstructedParticle *> _jetLCIORel;

			// autosave for output
			bool _autoconvert;
			string _savePrefix;

			bool _useMcp;

			ClassDef(LCIOStorer,0);
	};

}

#endif
