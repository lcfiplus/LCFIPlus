#include "LcfiplusProcessor.h"
#include "EventStore.h"

// ----- include for verbosity dependent logging ---------
#include "marlin/VerbosityLevels.h"
#include "marlin/StringParameters.h"
#define SLM streamlog_out(MESSAGE)

// Marlin stuff
#include <marlin/Global.h>
#include <marlin/Exceptions.h>

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
  
  // modify processor description
	_description = "Lcfiplus general processor" ;

	// MCP on/off
	registerProcessorParameter("UseMCP", "Whether MCParticle collection is imported or not", _useMcp, int(0));

	// input collections 
	registerInputCollection(LCIO::RECONSTRUCTEDPARTICLE, "PFOCollection" , "Particle flow output collection",
		_pfoCollectionName, std::string("PandoraPFOs"));
	registerInputCollection(LCIO::MCPARTICLE, "MCPCollection" , "MC particle collection",
		_mcpCollectionName, std::string("MCParticlesSkimmed"));
	registerInputCollection(LCIO::LCRELATION, "MCPFORelation", "Relation between MC and PFO particles",
		_mcpfoRelationName, std::string("RecoMCTruthLink"));

	registerProcessorParameter("VertexAutoLoad", "Loading LCIO vertices automatically", _autoVertex, int(1));
	registerProcessorParameter("JetAutoLoad", "Loading LCIO jets automatically", _autoJet, int(1));

	registerProcessorParameter("Algorithms", "LCFIPlus algorithms to run", _algonames, vector<string>());

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
			_lcio = new LCIOStorer(0,0,true,0); // no file
			if(_useMcp)
				_lcio->InitCollections(_pfoCollectionName.c_str(), _mcpCollectionName.c_str(), _mcpfoRelationName.c_str());
			else
				_lcio->InitCollectionsWithoutMCP(_pfoCollectionName.c_str());

			_lcioowner = true;
		}
		else
			_lcioowner = false;

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

		// printing EventStore collections
		Event::Instance()->Print();

		_nRun = 0 ;
		_nEvt = 0 ;

	}catch(lcfiplus::Exception &e){
		e.Print();
		throw(marlin::StopProcessingException(this));
	}
}

void LcfiplusProcessor::processRunHeader( LCRunHeader* run) { 

  _nRun++ ;
} 

void LcfiplusProcessor::processEvent( LCEvent * evt ) { 

	cout << "processEvent: event # " << _nEvt << endl;

	if (_lcio == false) return;

	try{

		// auto register collections
		if(_autoVertex)
			_lcio->InitVertexCollectionsAuto(evt);
		if(_autoJet)
			_lcio->InitJetCollectionsAuto(evt);

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
	}catch(lcfiplus::Exception &e){
		e.Print();
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

