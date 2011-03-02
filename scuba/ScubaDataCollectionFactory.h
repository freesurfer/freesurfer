/**
 * @file  ScubaDataCollectionFactory.h
 * @brief Factory for creating DataCollection subclasse based on type
 *
 * Creates Scuba specfic DataCollections based on type strings (MRI
 * and MRIS). Also handles Tcl commands to do the same.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:36 $
 *    $Revision: 1.6 $
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
 *
 * Terms and conditions for use, reproduction, distribution and contribution
 * are found in the 'FreeSurfer Software License Agreement' contained
 * in the file 'LICENSE' found in the FreeSurfer distribution, and here:
 *
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
 *
 * Reporting: freesurfer@nmr.mgh.harvard.edu
 *
 */


#ifndef ScubaDataCollectionFactory_h
#define ScubaDataCollectionFactory_h

#include "string_fixed.h"
#include "DebugReporter.h"
#include "TclCommandManager.h"
#include "DataCollection.h"

class ScubaDataCollectionFactory : public DebugReporter, public TclCommandListener {

  friend class ScubaDataCollectionFactoryTester;

public:
  // Gets the static reference to this class.
  static ScubaDataCollectionFactory& GetFactory();

  virtual TclCommandResult
  DoListenToTclCommand ( char* iCommand, int iArgc, char** iArgv );

  DataCollection& MakeDataCollection ( std::string iType );


protected:
  static bool mbAddedTclCommands;
};

#endif
