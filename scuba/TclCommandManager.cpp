#include <string>
#include <stdexcept>
#include "TclCommandManager.h"

using namespace std;

TclCommandListener::TclCommandResult
TclCommandListener::ListenToTclCommand ( char* iCommand, 
					 int iArgc, char** iArgv ) {
  sReturnFormat.clear();
  sReturnValues.clear();
  sResult.clear();

  return this->DoListenToTclCommand( iCommand, iArgc, iArgv );
}

TclCommandListener::~TclCommandListener () {

  TclCommandManager& commandMgr = TclCommandManager::GetManager();
  commandMgr.RemoveListener( *this );
}

TclCommandManager::TclCommandManager() : DebugReporter() {

  mbStarted = false;
  mInterp = 0;
}

TclCommandManager& 
TclCommandManager::GetManager() {

  static TclCommandManager* sManager = NULL;
  if( NULL == sManager ) {
    sManager = new TclCommandManager();

    sManager->AddCommand( *sManager, "PrintAllCommands", 0, "", 
			 "Print all registered commands." );
  }

  return *sManager;
}

void
TclCommandManager::AddCommand (  TclCommandListener& iListener,
				 char const* isCommand,
				 int  icArgs,
				 char const* isArgHints, 
				 char const* isDescription ) {

  std::list<Command*>::iterator tCommand;
  for( tCommand = mlCommands.begin(); 
       tCommand != mlCommands.end(); ++tCommand ) {

    Command* command = *tCommand;
    if( command->msCommand.compare( string(isCommand) ) == 0 ) {
      
      command->mlListeners.push_back( &iListener );
      return;
    }
  }

  struct Command* command = new Command;
  command->msCommand = isCommand;
  command->mcArgs = icArgs;
  command->msArgHints = isArgHints;
  command->msDescription = isDescription;

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

  char sCommand[1024];
  strcpy( sCommand,  iCommand.msCommand.c_str() );
  Tcl_CreateCommand( mInterp, sCommand,
		     TclCommandManager::HandleCommand,
		     (ClientData)NULL, (Tcl_CmdDeleteProc*)NULL );
}

int
TclCommandManager::HandleCommand ( ClientData iClientData, Tcl_Interp* iInterp,
				   int argc, char** argv ) {

  int rTcl = TCL_OK;
  TclCommandListener::TclCommandResult finalResult = ok;
  TclCommandManager& commandMgr = TclCommandManager::GetManager(); 

  std::list<Command*>::iterator tCommand;
  for( tCommand = commandMgr.mlCommands.begin(); 
       tCommand != commandMgr.mlCommands.end(); ++tCommand ) {

    Command* command = *tCommand;
    if( command->msCommand.compare( string(argv[0]) ) == 0 ) {

      // Check for right number of args. mcArgs should be argc-1
      // because argc also includes the command name.
      if( command->mcArgs != argc - 1 ) {
	stringstream ssError;
	ssError << "wrong # args: should be " << command->msCommand << " "
		<< command->msArgHints;
	char* sResult = strdup( ssError.str().c_str() );
	Tcl_SetResult( iInterp, sResult, TCL_VOLATILE );
	free( sResult );
	finalResult = error;

      } else {
	
	// For each listener to this command...
	std::list<TclCommandListener*>::iterator tListener;
	for( tListener = command->mlListeners.begin(); 
	     tListener != command->mlListeners.end(); ++tListener ) {
	  
	  // Give the listener a chance to listen to this command.
	  TclCommandListener* listener = *tListener;
	  TclCommandListener::TclCommandResult result =
	    listener->ListenToTclCommand( argv[0], argc, argv );
	  
	  // If they gave us a return format string, convert the
	  // format string and values string into a tcl object, which
	  // could be a single value object or a list.
	  if( listener->sReturnFormat.length() != 0 ) {
	    stringstream sFormat( listener->sReturnFormat );
	    stringstream sValues( listener->sReturnValues );
	    Tcl_Obj* rObj = 
	      commandMgr.ConvertFStringToTclObj( sFormat, sValues, iInterp );
	    Tcl_SetObjResult( iInterp, rObj );
	  }

	  // If they gave us a result string, change that into a
	  // string and set the result to that.
	  if( listener->sResult.length() != 0 ) {
	    char* sResult = strdup( listener->sResult.c_str() );
	    Tcl_SetResult( iInterp, sResult, TCL_VOLATILE );
	    free( sResult );
	  }
	  
	  // If the result was not ok, save that in the final
	  // result. This way many objects can respond to the command
	  // but an error will be returned if any of them return an
	  // error.
	  if( result != ok ) {
	    finalResult = result;
	  }
	}
      }
    }
  }
  
  switch( finalResult ) {
  case ok:
    rTcl = TCL_OK;
    break;
  case error:
    rTcl = TCL_ERROR;
    break;
  }
  
  return rTcl;
}


Tcl_Obj*
TclCommandManager::ConvertFStringToTclObj( stringstream& isFormat,
					   stringstream& isValues,
					   Tcl_Interp* iInterp ) {

  Tcl_Obj* rObj = NULL;

  char cFormat;
  isFormat >> cFormat;
  switch( cFormat ) {
    
  case 'L': {
    Tcl_Obj* list = Tcl_NewListObj( 0, NULL );
    
    Tcl_Obj* newObject =
      ConvertFStringToTclObj( isFormat, isValues, iInterp );
    
    while( newObject != NULL ) {
      Tcl_ListObjAppendElement( iInterp, list, newObject );
      
      newObject = 
	ConvertFStringToTclObj( isFormat, isValues, iInterp );
    }
    
    rObj = list;
  }
    break;
    
  case 'l':
    rObj = NULL;
    break;
    
  case 'i': {
    
    // Pull a string off the value stream.
    string sValue;
    isValues >> sValue;
    
    // Try to convert it to an int.
    int value = strtol(sValue.c_str(), (char**)NULL, 10);
    if( ERANGE == errno ) {
      DebugOutput( << "Error converting " << sValue << " to int." );
      throw logic_error( "couldn't covnert string to int" );
    }
    
    // Return it.
    rObj = Tcl_NewIntObj( value );
  }
    break;
    
  case 'f': {
    
    // Pull a string off the value stream.
    string sValue;
    isValues >> sValue;

    // Try to convert it to a float.
    double value = strtod(sValue.c_str(), (char**)NULL);
    if( ERANGE == errno ) {
      DebugOutput( << "Error converting " << sValue << " to double." );
      throw logic_error( "couldn't covnert string to double" );
    }

    // Return it.
    rObj = Tcl_NewDoubleObj( value );
  }
    break;

  case 's': {
    
    // Pull a string off the value stream.
    stringstream sReturn;
    string sValue;
    isValues >> sValue;
    sReturn << sValue;

    // If the first char was a ", keep going until we get another ".
    if( sValue[0] == '\"' ) {
      while( sReturn.str().rfind('\"', sReturn.str().length()) == 0 ) {
	isValues >> sValue;
	sReturn << " " << sValue;
      }
    }

    // Strip the quotes if necessary and return it.
    string sReturn2 = sReturn.str();
    if( sReturn2[0] == '\"' ) {
      sReturn2 = sReturn2.substr( 1, sReturn2.length()-2 );
    }
    rObj = Tcl_NewStringObj( sReturn2.c_str(), sReturn2.length() );
  }
    break;

  default:
    DebugOutput( << "Invalid format char " << cFormat );
  }

  return rObj;
}

void
TclCommandManager::RemoveListener ( TclCommandListener& iListener ) {

  std::list<Command*>::iterator tCommand;
  for( tCommand = mlCommands.begin(); 
       tCommand != mlCommands.end(); ++tCommand ) {

    Command* command = *tCommand;
    command->mlListeners.remove( &iListener );
  }
}

TclCommandListener::TclCommandResult
TclCommandManager::DoListenToTclCommand ( char* isCommand, 
					  int iArgc, char** iasArgv ) {

  // PrintAllCommands
  if( 0 == strcmp( isCommand, "PrintAllCommands" ) ) {
    sResult = PrintAllCommands();
  }

  return ok;
}


string
TclCommandManager::PrintAllCommands () {

  stringstream ssResult;

  std::list<Command*>::iterator tCommand;
  for( tCommand = mlCommands.begin(); 
       tCommand != mlCommands.end(); ++tCommand ) {
    
    Command* command = *tCommand;

    ssResult << command->msCommand << " " << command->msArgHints << endl;
    ssResult << command->msDescription << endl << endl;
  }

  return ssResult.str();
}
