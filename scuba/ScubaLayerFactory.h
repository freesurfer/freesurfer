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
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:38 $
 *    $Revision: 1.7 $
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
