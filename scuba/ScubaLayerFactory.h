/**
 * @file  ScubaLayerFactory.h
 * @brief A factory for Scuba specific Layers
 *
 * Makes layers  for specific types.  This is a little  different that
 * our  other factories  because it  creates  layers based  on a  type
 * (2DMRI and 2DMRIS). Also handles Tcl commands for doing the same.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: kteich $
 *    $Date: 2007/10/22 04:39:29 $
 *    $Revision: 1.6 $
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


#ifndef ScubaLayerFactory_h
#define ScubaLayerFactory_h

#include "string_fixed.h"
#include "DebugReporter.h"
#include "TclCommandManager.h"
#include "Layer.h"

class ScubaLayerFactory : public DebugReporter, public TclCommandListener {

  friend class ScubaLayerFactoryTester;

public:
  // Gets the static reference to this class.
  static ScubaLayerFactory& GetFactory();

  virtual TclCommandResult
  DoListenToTclCommand ( char* isCommand, int iArgc, char** iArgv );

  Layer& MakeLayer ( std::string iType );


protected:
  static bool mbAddedTclCommands;
};

#endif
