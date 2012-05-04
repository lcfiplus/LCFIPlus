#ifndef LcfiplusProcessor_h
#define LcfiplusProcessor_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>

#include "LCIOStorer.h"
#include "EventStore.h"

using namespace lcio ;
using namespace marlin ;

/**
 *  Marlin processor for LCFIPlus.
 * 
 * @author Tomohiko Tanabe, ICEPP, The University of Tokyo
 * @author Taikan Suehara, ICEPP, The University of Tokyo
 * @version $Id$
 */

class LcfiplusProcessor : public Processor, public lcfiplus::EventStoreObserver
{
  
 public:
  
  virtual Processor*  newProcessor() { return new LcfiplusProcessor ; }
  
  LcfiplusProcessor() ;
  virtual ~LcfiplusProcessor() ;
  
  /** Called at the begin of the job before anything is read.
   */
  virtual void init() ;
  
  /** Called for every run.
   */
  virtual void processRunHeader( LCRunHeader* run ) ;
  
  /** Called for every event - the working horse.
   */
  virtual void processEvent( LCEvent * evt ) ; 
  
  
  virtual void check( LCEvent * evt ) ; 
  
  
  /** Called after data processing for clean up.
   */
  virtual void end() ;
  

	virtual void RegisterCallback(const char *name, const char *classname, int flags);

 private:

	// lciostorer singleton
	static lcfiplus::LCIOStorer *_lcio;
	bool _lcioowner;

	int _useMcp;

  /** Input collection name.
   */
	std::string _pfoCollectionName;
	std::string _mcpCollectionName;
	std::string _mcpfoRelationName;
	std::vector<std::string> _algonames;

	std::vector<lcfiplus::Algorithm *> _algos;
	lcfiplus::Parameters * _param;

  int _nRun ;
  int _nEvt ;
	int _printPeriod;

	int _readSubdetectorEnergies;
	int _updateVertexRPDaughters;
	int _ignoreLackOfVertexRP;

	// collections to register
	std::vector<std::string> _vertexColNamesToWrite;
	std::vector<std::string> _jetColNamesToWrite;
	std::vector<int> _vertexColNamesToWriteFlags;
	std::vector<int> _jetColNamesToWriteFlags;

	bool _inInit;
} ;

#endif



