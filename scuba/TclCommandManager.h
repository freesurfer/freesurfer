//
// TclCommandManager.h
//
// Warning: Do not edit the following four lines.  CVS maintains them.
// Revision Author: $Author: kteich $
// Revision Date  : $Date: 2003/10/17 16:52:28 $
// Revision       : $Revision: 1.2 $

#ifndef TclCommandManager_h
#define TclCommandManager_h

#include <stdlib.h>
#include <list>
#include <tcl.h>
#include "DebugReporter.h"

class TclCommandListener {
 public:
  virtual void ListenToTclCommand ( char* iCommand, int iArgc, char** iArgv ) = 0;
};

class TclCommandManager : public DebugReporter {

 protected:

  // A command being managed. Contains the command name and a list of
  // listeners for that command.
  struct Command {
    char msCommand[1024];
    std::list<TclCommandListener*> mlListeners;
  };
  
  
 public:

  // Gets the static reference to this class.
  static TclCommandManager& GetManager();
  
  // Adds a TclCommandListener the list of listeners for a
  // command. When the command is called in Tcl, the Listener's
  // ListenToTclCommand() function will be called.
  void AddCommand ( TclCommandListener& iListener, char const* isCommand );
  
  // Start managing commands for the given Tcl interp context.
  void Start ( Tcl_Interp const* iInterp );
  

  // The manager's Tcl callback function. Shouldn't be publically
  // called, but needs to be static and public as it is a C-style
  // callback.
  static int HandleCommand ( ClientData iClientData, Tcl_Interp* iInterp,
			     int argc, char** argv );

 protected:

  TclCommandManager();

  // Create the command in the Tcl context and attach our handler.
  void CreateCommand( Command& iCommand );

  // List of all commands.
  std::list<Command*> mlCommands;

  // Information about the Tcl context and whether or not we're
  // 'listening' to one.
  Tcl_Interp* mInterp;
  bool mbStarted;
};


#endif
