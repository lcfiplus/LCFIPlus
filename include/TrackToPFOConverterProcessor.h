#ifndef TrackToPFOConverterProcessor_h
#define TrackToPFOConverterProcessor_h 1

#include "marlin/Processor.h"
#include <string>

using namespace marlin ;

/**
 *  Marlin processor for LCFIPlus.
 *
 * @author Tomohiko Tanabe, ICEPP, The University of Tokyo
 * @author Taikan Suehara, ICEPP, The University of Tokyo
 * @version $Id$
 */

class TrackToPFOConverterProcessor : public Processor {

 public:

  virtual Processor*  newProcessor() {
    return new TrackToPFOConverterProcessor ;
  }

  TrackToPFOConverterProcessor() ;
  virtual ~TrackToPFOConverterProcessor() ;

  /** Called at the begin of the job before anything is read.
   */
  virtual void init() ;

  /** Called for every run.
   */
  virtual void processRunHeader( LCRunHeader* run ) ;

  /** Called for every event - the working horse.
   */
  virtual void processEvent( LCEvent* evt ) ;


  virtual void check( LCEvent* evt ) ;


  /** Called after data processing for clean up.
   */
  virtual void end() ;

 private:

  /** Input collection name.
   */
  std::string _inputTrackCollectionName;
  std::string _outputPFOCollectionName;
} ;

#endif



