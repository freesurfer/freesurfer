/**
 * @file  ScubaDataCollectionFactory.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 02:09:14 $
 *    $Revision: 1.4 $
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
