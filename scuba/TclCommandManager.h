//
// TclCommandManager.h
//
// Warning: Do not edit the following four lines.  CVS maintains them.
// Revision Author: $Author: kteich $
// Revision Date  : $Date: 2004/04/03 22:04:02 $
// Revision       : $Revision: 1.8 $

#ifndef TclCommandManager_h
#define TclCommandManager_h

#include <stdlib.h>
#include <list>
#include <tcl.h>
#include <string>
#include <iostream>
#include <stdexcept>
#include <sstream>
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
    ListenToTclCommand ( char* iCommand, int iArgc, char** iArgv );
  virtual TclCommandResult
    DoListenToTclCommand ( char* iCommand, int iArgc, char** iArgv ) = 0;
  
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

  bool Started () const { return mbStarted; }

  // We listen to commands too.
  TclCommandResult
    DoListenToTclCommand ( char* iCommand, int iArgc, char** iArgv );

  std::string PrintAllCommands ();

  std::string SendCommand ( std::string isCommand );

  // Lets the tcl parser handle an event.
  void DoTclEvent ();

  void SetCommandLineParameters ( int iArgc, char** iArgv );

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
