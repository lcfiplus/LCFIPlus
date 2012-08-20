#include "LcfiplusProcessor.h"
#include "lcfiplus.h"

// ----- include for verbosity dependent logging ---------
#include "marlin/VerbosityLevels.h"
#include "marlin/StringParameters.h"
#define SLM streamlog_out(MESSAGE)

// Marlin stuff
#include <marlin/Global.h>
#include <marlin/Exceptions.h>

// GEAR
#include <gear/BField.h>

#include "TROOT.h"
#include "TApplication.h"

using namespace std;
using namespace lcio ;
using namespace marlin ;
using namespace lcfiplus ;

LcfiplusProcessor aLcfiplusProcessor ;

// static object initialization
LCIOStorer * LcfiplusProcessor::_lcio = 0;

LcfiplusProcessor::LcfiplusProcessor() : Processor("LcfiplusProcessor") {
  
	_inInit = false;
	_printPeriod = 0;

  // modify processor description
	_description = "Lcfiplus general processor" ;

	// MCP on/off
	registerProcessorParameter("UseMCP", "Whether MCParticle collection is imported or not", _useMcp, int(0));

	// input collections 
	registerInputCollection(LCIO::RECONSTRUCTEDPARTICLE, "PFOCollection" , "Particle flow output collection",
		_pfoCollectionName, std::string(""));
	registerInputCollection(LCIO::MCPARTICLE, "MCPCollection" , "MC particle collection",
		_mcpCollectionName, std::string(""));
	registerInputCollection(LCIO::LCRELATION, "MCPFORelation", "Relation between MC and PFO particles",
		_mcpfoRelationName, std::string(""));

	registerProcessorParameter("Algorithms", "LCFIPlus algorithms to run", _algonames, vector<string>());
	registerProcessorParameter("ReadSubdetectorEnergies", "Read subdetector energies (ILD)", _readSubdetectorEnergies, int(1));
	registerProcessorParameter("UpdateVertexRPDaughters", "Writing back obtained vertices to input RP collections (which must be writable)",
		_updateVertexRPDaughters, int(1));
	registerProcessorParameter("IgnoreLackOfVertexRP", "Keep running even if vertex RP collection is not present",
			_ignoreLackOfVertexRP, int(0));
	registerProcessorParameter("PrintEventNumber", "Event number printing period in std output: 0 = no printing", _printPeriod, int(0));

	registerOptionalParameter("MagneticField", "Manually set magnetic field, overriding the value from GEAR [T]", _magneticField, float(0.0));
	registerOptionalParameter("BeamSizeX", "Bunch size in the X direction [mm]", _beamSizeX, float(639e-6));
	registerOptionalParameter("BeamSizeY", "Bunch size in the Y direction [mm]", _beamSizeY, float(5.7e-6));
	registerOptionalParameter("BeamSizeZ", "Bunch size in the Z direction [mm]", _beamSizeZ, float(9.13e-2));

	// ROOT object
/*	int argc = 0;
	if(gROOT->GetApplication() == 0){
		//TApplication *theapp =
		new TApplication("LcfiplusProcessor",&argc,0);
		SLM << "TApplication created." << endl;
	}
*/
}

LcfiplusProcessor::~LcfiplusProcessor()
{
/*
	if(gROOT->GetApplication()){
		cout << "We have application" << endl;
		delete gROOT->GetApplication();
	}
*/
}

void LcfiplusProcessor::init() { 

	try{

		streamlog_out(DEBUG) << "   init called  " 
			<< std::endl ;

		// usually a good idea to
		printParameters() ;

		StringParameters *parameter = parameters();
/*
		// obtain algorithm name
		if(!parameter->isParameterSet("algorithm")){
			SLM << "Lcfiplus algorithm not set. run nothing." << endl;
			return;
		}
		vector<string> algos;
		parameter->getStringVals("algorithm", algos);

*/

		// algorithm check
		if(_algonames.size() == 0){
			SLM << "No algorithms given to run. run nothing." << endl;
			return;
		}

		// set globals
		if (_magneticField != 0)
			Globals::Instance()->setBField(_magneticField);
		else
			Globals::Instance()->setBField( marlin::Global::GEAR->getBField().at(gear::Vector3D(0.,0.,0.)).z() );

		Globals::Instance()->setBeamSizeX(_beamSizeX);
		Globals::Instance()->setBeamSizeY(_beamSizeY);
		Globals::Instance()->setBeamSizeZ(_beamSizeZ);

		// register observer
		Event::Instance()->RegisterObserver(this);

		// conversion StringParameters -> Parameters
		_param = new Parameters;
		
		StringVec keys;
		parameter->getStringKeys(keys);

		for(unsigned int i=0;i<keys.size();i++){
			vector<string> vals;
			parameter->getStringVals(keys[i], vals);
			_param->add(keys[i].c_str(), vals);
		}

		// initialize LCIOStorer
		if(!_lcio){
			_lcio = new LCIOStorer(0,0,true,false,0); // no file
			_lcio->setReadSubdetectorEnergies(_readSubdetectorEnergies);
			_lcio->setUpdateVertexRPDaughters(_updateVertexRPDaughters);
			_lcio->setIgnoreLackOfVertexRP(_ignoreLackOfVertexRP);

			_lcioowner = true;
		}
		else{
			_lcioowner = false;

			if((_lcio->getReadSubdetectorEnergies() != _readSubdetectorEnergies)
			 ||(_lcio->getUpdateVertexRPDaughters() != _updateVertexRPDaughters)
			 ||(_lcio->getIgnoreLackOfVertexRP() != _ignoreLackOfVertexRP)){
					throw(lcfiplus::Exception("Global parameters do not match to previous processors: specify the same for all LcfiplusProcessors."));
			}
		}

		// load basic collection
		if(_useMcp)
			_lcio->InitMCPPFOCollections(_pfoCollectionName.c_str(), _mcpCollectionName.c_str(), _mcpfoRelationName.c_str());
		else
			_lcio->InitPFOCollections(_pfoCollectionName.c_str());

		_inInit = true;

		// make algorithm classes and pass param to init
		for(unsigned int i=0;i<_algonames.size();i++){
			string s = "lcfiplus::";
			s += _algonames[i];
			TClass *cl = TClass::GetClass(s.c_str());
			if(!cl || !cl->InheritsFrom("lcfiplus::Algorithm")){
				SLM << "Algorithm " << _algonames[i] << " is not valid. skip." << endl;
				continue;
			}
			Algorithm *newalgo = (Algorithm *)cl->New();
			if(!newalgo){SLM << "Initialization failed!." << endl; break;}
			_algos.push_back(newalgo);
			SLM << "Algorithm " << _algonames[i] << " is being initialized." << endl;
			newalgo->init(_param);

			SLM << "Algorithm " << _algonames[i] << " successfully initialized." << endl;
		}

		_inInit = false;

		// printing EventStore collections
		Event::Instance()->Print();

		_nRun = 0 ;
		_nEvt = 0 ;

	}catch(lcfiplus::Exception &e){
		streamlog_out(ERROR) << e.what() << endl;
		throw(marlin::StopProcessingException(this));
	}
}

void LcfiplusProcessor::processRunHeader( LCRunHeader* run) { 

  _nRun++ ;
} 

void LcfiplusProcessor::processEvent( LCEvent * evt ) { 

	if(_printPeriod && _nEvt % _printPeriod == 0)
		cout << "processEvent: event # " << _nEvt << endl;

	if (_lcio == false) return;

	try{

/*		// auto register collections
		if(_autoVertex)
			_lcio->InitVertexCollectionsAuto(evt);
		if(_autoJet)
			_lcio->InitJetCollectionsAuto(evt);
*/
		// set LCEvent
		if(_lcioowner)
			_lcio->SetEvent(evt);

		// set deafult mcp/track/netural
		Event::Instance()->setDefaultTracks(_pfoCollectionName.c_str());
		Event::Instance()->setDefaultNeutrals(_pfoCollectionName.c_str());
		Event::Instance()->setDefaultMCParticles(_mcpCollectionName.c_str());

		// process registered algorithms
		for(unsigned int i=0;i<_algos.size();i++){
			_algos[i]->process();
		}

		// convert persistent collections
//		_lcio->AutoConvert();
		for(unsigned int i=0;i<_vertexColNamesToWrite.size();i++){
			bool persist = _vertexColNamesToWriteFlags[i] & EventStore::PERSIST;
			if (persist)
				_lcio->WriteVertices(_vertexColNamesToWrite[i].c_str());
		}
		for(unsigned int i=0;i<_jetColNamesToWrite.size();i++){
			bool persist = _jetColNamesToWriteFlags[i] & EventStore::PERSIST;
			bool writeVertex = _jetColNamesToWriteFlags[i] & EventStore::JET_WRITE_VERTEX;
			if (persist)
				_lcio->WriteJets(_jetColNamesToWrite[i].c_str(),0,writeVertex);
		}

		_nEvt ++ ;
	}catch(lcfiplus::Exception &e){
		streamlog_out(ERROR) << e.what() << endl;
		throw(marlin::StopProcessingException(this));
	}
}


void LcfiplusProcessor::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void LcfiplusProcessor::end(){ 
  
	for(unsigned int i=0;i<_algos.size();i++){
		_algos[i]->end();
	}

	Event::Instance()->ClearObjects();

	std::cout << "LcfiplusProcessor::end()  " << name() 
		<< " processed " << _nEvt << " events in " << _nRun << " runs "
		<< std::endl ;

}

void LcfiplusProcessor::RegisterCallback(const char *name, const char *classname, int flags)
{
	if(!_inInit)return; // no-op other than init

	if(string("vector<lcfiplus::Vertex*>") == classname){
		cout << "Registering output LCIO vertex collection: " << name << endl;
		_vertexColNamesToWrite.push_back(name);
		_vertexColNamesToWriteFlags.push_back(flags);
	}
	else if(string("vector<lcfiplus::Jet*>") == classname){
		cout << "Registering output LCIO jet collection: " << name << endl;
		_jetColNamesToWrite.push_back(name);
		_jetColNamesToWriteFlags.push_back(flags);
	}
}
