#include <string>
#include <stdexcept>
#include "TclCommandManager.h"

using namespace std;

void
TclCommandListener::ListenToTclCommand ( char* iCommand, 
					 int iArgc, char** iArgv ) {
  sReturnFormat.clear();
  sReturnValues.clear();
  sResult.clear();

  this->DoListenToTclCommand( iCommand, iArgc, iArgv );
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
  }

  return *sManager;
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

	if( listener->sReturnFormat.length() != 0 ) {
	  stringstream sFormat( listener->sReturnFormat );
	  stringstream sValues( listener->sReturnValues );

	  Tcl_Obj* rObj = 
	    commandMgr.ConvertFStringToTclObj( sFormat, sValues, iInterp );
	  Tcl_SetObjResult( iInterp, rObj );
	}

	if( listener->sResult.length() != 0 ) {
	  char* sResult = strdup( listener->sResult.c_str() );
	  Tcl_SetResult( iInterp, sResult, TCL_VOLATILE );
	  free( sResult );
	}
      }
    }
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
