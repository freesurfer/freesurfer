#include <string>
#include <stdexcept>
#include "TclCommandManager.h"

TclCommandManager::TclCommandManager() : DebugReporter() {

  mbStarted = false;
  mInterp = 0;
}

TclCommandManager& 
TclCommandManager::GetManager() {

  static TclCommandManager sManager;

  return sManager;
}

void
TclCommandManager::AddCommand (  TclCommandListener& iListener,
				 char const* isCommand ) {

  std::list<Command*>::iterator tCommand;
  for( tCommand = mlCommands.begin(); 
       tCommand != mlCommands.end(); ++tCommand ) {

    Command* command = *tCommand;
    if( strcmp( command->msCommand, isCommand ) == 0 ) {
      
      command->mlListeners.push_back( &iListener );
      return;
    }
  }

  struct Command* command = new Command;
  strcpy( command->msCommand, isCommand );
  command->mlListeners.push_back( &iListener );
  mlCommands.push_back( command );

  if( mbStarted ) {
    CreateCommand( *command );
  }
}
  
void
TclCommandManager::Start( Tcl_Interp const* iInterp ) {

  mInterp = (Tcl_Interp*)iInterp;
  mbStarted = true;

  std::list<Command*>::iterator tCommand;
  for( tCommand = mlCommands.begin(); 
       tCommand != mlCommands.end(); ++tCommand ) {

    Command* command = *tCommand;
    CreateCommand( *command );
  }
  
}

void
TclCommandManager::CreateCommand( Command& iCommand ) {

  if( 0 == mInterp ) {
    DebugOutput( << "CreateCommand() called with mInterp NULL" );
    throw std::logic_error( "Tried to CreateCommand without an interpretor" );
  }
  
  Tcl_CreateCommand( mInterp, iCommand.msCommand, 
		     TclCommandManager::HandleCommand,
		     (ClientData)NULL, (Tcl_CmdDeleteProc*)NULL );
}

int
TclCommandManager::HandleCommand ( ClientData iClientData, Tcl_Interp* iInterp,
				   int argc, char** argv ) {

  int rTcl = TCL_OK;
  TclCommandManager& commandMgr = TclCommandManager::GetManager(); 

  std::list<Command*>::iterator tCommand;
  for( tCommand = commandMgr.mlCommands.begin(); 
       tCommand != commandMgr.mlCommands.end(); ++tCommand ) {

    Command* command = *tCommand;
    if( strcmp( command->msCommand, argv[0] ) == 0 ) {

      std::list<TclCommandListener*>::iterator tListener;
      for( tListener = command->mlListeners.begin(); 
	   tListener != command->mlListeners.end(); ++tListener ) {

	TclCommandListener* listener = *tListener;
	listener->ListenToTclCommand( argv[0], argc, argv );
      }
    }
  }

  return rTcl;
}
