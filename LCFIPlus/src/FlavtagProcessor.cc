#include "FlavtagProcessor.h"
#include "EventStore.h"

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"
#include "marlin/StringParameters.h"

// Marlin stuff
#include <marlin/Global.h>

#include "TROOT.h"
#include "TApplication.h"

#define SLM streamlog_out(MESSAGE)

using namespace std;
using namespace lcio ;
using namespace marlin ;
using namespace flavtag ;

FlavtagProcessor aFlavtagProcessor ;

// static object initialization
LCIOStorer * FlavtagProcessor::_lcio = 0;

FlavtagProcessor::FlavtagProcessor() : Processor("FlavtagProcessor") {
  
  // modify processor description
	_description = "Flavtag general processor" ;

	// input collections 
	registerInputCollection(LCIO::RECONSTRUCTEDPARTICLE, "PFOCollection" , "Particle flow output collection",
		_pfoCollectionName, std::string("PandoraPFOs"));
	registerInputCollection(LCIO::MCPARTICLE, "MCPCollection" , "MC particle collection",
		_mcpCollectionName, std::string("MCParticlesSkimmed"));
	registerInputCollection(LCIO::LCRELATION, "MCPFORelation", "Relation between MC and PFO particles",
		_mcpfoRelationName, std::string("RecoMCTruthLink"));


	// ROOT object
	int argc = 0;
	if(gROOT->GetApplication() == 0){
		TApplication *theapp = new TApplication("FlavtagProcessor",&argc,0);
		SLM << "TApplication created." << endl;
	}
}

void FlavtagProcessor::init() { 

  streamlog_out(DEBUG) << "   init called  " 
		       << std::endl ;

  // usually a good idea to
  printParameters() ;

	StringParameters *parameter = parameters();

	// obtain algorithm name
	if(!parameter->isParameterSet("algorithm")){
		SLM << "Flavtag algorithm not set. run nothing." << endl;
		return;
	}
	vector<string> algos;
	parameter->getStringVals("algorithm", algos);

	if(algos.size() == 0){
		SLM << "Flavtag algorithm size is 0. run nothing." << endl;
		return;
	}

	// conversion StringParameters -> FlavtagParameters
	_param = new FlavtagParameters;
	
	StringVec keys;
	parameter->getStringKeys(keys);

	for(unsigned int i=0;i<keys.size();i++){
		vector<string> vals;
		parameter->getStringVals(keys[i], vals);
		_param->add(keys[i].c_str(), vals);
	}

	// initialize LCIOStorer
	if(!_lcio){
		_lcio = new LCIOStorer(0,0,true,0); // no file
		_lcio->InitCollections(_pfoCollectionName.c_str(), _mcpCollectionName.c_str(), _mcpfoRelationName.c_str());
		_lcioowner = true;
	}
	else
		_lcioowner = false;

	// make algorithm classes and pass param to init
	for(unsigned int i=0;i<algos.size();i++){
		string s = "flavtag::";
		s += algos[i];
		TClass *cl = TClass::GetClass(s.c_str());
		if(!cl || !cl->InheritsFrom("flavtag::FlavtagAlgorithm")){
			SLM << "Algorithm " << algos[i] << " is not valid. skip." << endl;
			continue;
		}
		FlavtagAlgorithm *newalgo = (FlavtagAlgorithm *)cl->New();
		if(!newalgo){SLM << "Initialization failed!." << endl; break;}
		_algos.push_back(newalgo);
		SLM << "Algorithm " << algos[i] << " is being initialized." << endl;
		newalgo->init(_param);

		SLM << "Algorithm " << algos[i] << " successfully initialized." << endl;
	}

  _nRun = 0 ;
  _nEvt = 0 ;
  
}

void FlavtagProcessor::processRunHeader( LCRunHeader* run) { 

  _nRun++ ;
} 

void FlavtagProcessor::processEvent( LCEvent * evt ) { 

	// set LCEvent
	if(_lcioowner)
		_lcio->SetEvent(evt);

	// process registered algorithms
	for(unsigned int i=0;i<_algos.size();i++){
		_algos[i]->process();
	}

	// convert persistent collections
	_lcio->AutoConvert();

  _nEvt ++ ;
}


void FlavtagProcessor::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void FlavtagProcessor::end(){ 
  
	for(unsigned int i=0;i<_algos.size();i++){
		_algos[i]->end();
	}


   std::cout << "FlavtagProcessor::end()  " << name() 
 	    << " processed " << _nEvt << " events in " << _nRun << " runs "
 	    << std::endl ;

}

