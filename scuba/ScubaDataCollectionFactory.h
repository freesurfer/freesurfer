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
 *    $Author: kteich $
 *    $Date: 2007/10/22 04:39:28 $
 *    $Revision: 1.5 $
 *
 * Copyright (C) 2002-2007,
 * The General Hospital Corporation (Boston, MA). 
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
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
