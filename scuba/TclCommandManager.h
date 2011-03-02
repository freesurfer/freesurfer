/**
 * @file  TclCommandManager.h
 * @brief Allows objects to register and respond to Tcl commands
 *
 * Classes should subclass from TclCommandListener to respond to
 * commands, then use the TclCommandManager to register their
 * commands. Override the DoListenToTclCommand() to respond to the
 * commands.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:39 $
 *    $Revision: 1.19 $
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


//
// TclCommandManager.h
//
// Warning: Do not edit the following four lines.  CVS maintains them.
// Revision Author: $Author: nicks $
// Revision Date  : $Date: 2011/03/02 00:04:39 $
// Revision       : $Revision: 1.19 $

#ifndef TclCommandManager_h
#define TclCommandManager_h

#include <errno.h>
#include <stdlib.h>
#include <list>
#include "string_fixed.h"
#include <iostream>
#include <stdexcept>
#include <sstream>

extern "C" {
#define USE_NON_CONST
#include <tcl.h>
#undef USE_NON_CONST
}

#include "DebugReporter.h"

// This class should be subclassed and the DoListenToTclCommand to
// implement the function to be called from Tcl space.  Values can be
// returned by writing a format and value string. See the test code
// for examples.
class TclCommandListener {
  friend class TclCommandManager;

public:
  enum TclCommandResult { ok, error };

  TclCommandResult
  ListenToTclCommand ( char* isCommand, int iArgc, char** iArgv );
  virtual TclCommandResult
  DoListenToTclCommand ( char* isCommand, int iArgc, char** iArgv ) = 0;

  virtual ~TclCommandListener();

protected:
  // Valid format chars are:
  //   i - int object
  //   f - float object
  //   s - string object
  //   L - begin list
  //   l - end list
  std::string sReturnFormat;
  std::string sReturnValues;
  std::string sResult;
};

class TclCommandManager : public DebugReporter, public TclCommandListener {

  friend class TclCommandManagerTester;

protected:

  // A command being managed. Contains the command name and a list of
  // listeners for that command.
  struct Command {
    std::string msCommand;
    std::string msArgHints;
    std::string msDescription;
    int mcArgs;
    std::list<TclCommandListener*> mlListeners;
  };


public:

  // Gets the static reference to this class.
  static TclCommandManager& GetManager();

  // Adds a TclCommandListener the list of listeners for a
  // command. When the command is called in Tcl, the Listener's
  // ListenToTclCommand() function will be called.
  void AddCommand ( TclCommandListener& iListener, char const* isCommand,
                    int icArgs, char const* isArgHints,
                    char const* isDescription );

  // Start managing commands for the given Tcl interp context.
  void Start ( Tcl_Interp const* iInterp );


  // The manager's Tcl callback function. Shouldn't be publically
  // called, but needs to be static and public as it is a C-style
  // callback.
  static int HandleCommand ( ClientData iClientData, Tcl_Interp* iInterp,
                             int argc, char** argv );

  // Remove all instances of a listener.
  void RemoveListener ( TclCommandListener& iListener );

  bool Started () const {
    return mbStarted;
  }

  // We listen to commands too.
  TclCommandResult
  DoListenToTclCommand ( char* iCommand, int iArgc, char** iArgv );

  std::string PrintAllCommands ();

  std::string SendCommand ( std::string isCommand );

  // Lets the tcl parser handle an event.
  void DoTclEvent ();

  void SetCommandLineParameters ( int iArgc, char** iArgv );

  // Just run a script.
  std::string RunTclScript ( char* ifnScript );

  // For converting arguments. Will throw run_time exception if not
  // it's an unconvertible value. The exception's message will
  // specify what the problem with the value was.
  static int ConvertArgumentToInt ( std::string isArg );
  static float ConvertArgumentToFloat ( std::string isArg );
  static bool ConvertArgumentToBoolean ( std::string isArg );
  static std::string ConvertBooleanToReturnValue ( bool ib );

  Tcl_Interp* GetInterp () {
    return mInterp;
  }

protected:

  TclCommandManager();

  Tcl_Obj* ConvertFStringToTclObj( std::stringstream& isFormat,
                                   std::stringstream& isValues,
                                   Tcl_Interp* iInterp );

  // Create the command in the Tcl context and attach our handler.
  void CreateCommand( Command& iCommand );

  // List of all commands.
  std::list<Command*> mlCommands;

  // Information about the Tcl context and whether or not we're
  // 'listening' to one.
  Tcl_Interp* mInterp;
  bool mbStarted;

  int    mArgc;
  char** mArgv;
};


#endif
