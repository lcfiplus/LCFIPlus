#ifndef LcfiplusProcessor_h
#define LcfiplusProcessor_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>

#include "LCIOStorer.h"

using namespace lcio ;
using namespace marlin ;

/**  Example processor for marlin.
 * 
 *  If compiled with MARLIN_USE_AIDA 
 *  it creates a histogram (cloud) of the MCParticle energies.
 * 
 *  <h4>Input - Prerequisites</h4>
 *  Needs the collection of MCParticles.
 *
 *  <h4>Output</h4> 
 *  A histogram.
 * 
 * @param CollectionName Name of the MCParticle collection
 * 
 * @author F. Gaede, DESY
 * @version $Id: LcfiplusProcessor.h,v 1.1.1.1 2009/06/04 00:16:27 suehara Exp $ 
 */

class LcfiplusProcessor : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new LcfiplusProcessor ; }
  
  
  LcfiplusProcessor() ;
  
  /** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
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
  
 private:

	// lciostorer singleton
	static lcfiplus::LCIOStorer *_lcio;
	bool _lcioowner;

  /** Input collection name.
   */
	std::string _pfoCollectionName;
	std::string _mcpCollectionName;
	std::string _mcpfoRelationName;

	std::vector<lcfiplus::LcfiplusAlgorithm *> _algos;
	lcfiplus::LcfiplusParameters * _param;

  int _nRun ;
  int _nEvt ;

	int _autoVertex;
	int _autoJet;
} ;

#endif



