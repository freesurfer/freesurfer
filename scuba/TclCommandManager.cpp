/**
 * @file  TclCommandManager.cpp
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
 *    $Revision: 1.31 $
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


#include <stdlib.h>
#include <errno.h>
#include "string_fixed.h"
#include <stdexcept>
#include "TclCommandManager.h"
#ifdef __cplusplus
extern "C"
{
#endif
#include <string.h> // strcmp
#ifdef __cplusplus
}
#endif

using namespace std;

TclCommandListener::TclCommandResult
TclCommandListener::ListenToTclCommand ( char* isCommand,
    int iArgc, char** iArgv ) {
  sReturnFormat.clear();
  sReturnValues.clear();
  sResult.clear();

  return this->DoListenToTclCommand( isCommand, iArgc, iArgv );
}

TclCommandListener::~TclCommandListener () {

  TclCommandManager& commandMgr = TclCommandManager::GetManager();
  commandMgr.RemoveListener( *this );
}

TclCommandManager::TclCommandManager() : DebugReporter() {

  mbStarted = false;
  mInterp = 0;
  mArgc = 0;
  mArgv = NULL;
}

TclCommandManager&
TclCommandManager::GetManager() {

  static TclCommandManager* sManager = NULL;
  if ( NULL == sManager ) {
    sManager = new TclCommandManager();

    sManager->AddCommand( *sManager, "PrintAllCommands", 0, "",
                          "Print all registered commands." );
    sManager->AddCommand( *sManager, "GetArgc", 0, "",
                          "Returns the argc value from the command line." );
    sManager->AddCommand( *sManager, "GetArgv", 0, "",
                          "Returns the argv list from the command line." );
    sManager->AddCommand( *sManager, "DebugOutput", 1, "message",
                          "Prints a message to the debugging output." );
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
  for ( tCommand = mlCommands.begin();
        tCommand != mlCommands.end(); ++tCommand ) {

    Command* command = *tCommand;
    if ( command->msCommand.compare( string(isCommand) ) == 0 ) {

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

  if ( mbStarted ) {
    CreateCommand( *command );
  }
}

void
TclCommandManager::Start( Tcl_Interp const* iInterp ) {

  mInterp = (Tcl_Interp*)iInterp;
  mbStarted = true;

  std::list<Command*>::iterator tCommand;
  for ( tCommand = mlCommands.begin();
        tCommand != mlCommands.end(); ++tCommand ) {

    Command* command = *tCommand;
    CreateCommand( *command );
  }
}

void
TclCommandManager::CreateCommand( Command& iCommand ) {

  if ( 0 == mInterp ) {
    DebugOutput( << "CreateCommand() called with mInterp NULL" );
    throw std::logic_error( "Tried to CreateCommand without an interpretor" );
  }

  char sCommand[1024];
  strcpy( sCommand,  iCommand.msCommand.c_str() );
  Tcl_CreateCommand( mInterp, sCommand,
                     (Tcl_CmdProc *) TclCommandManager::HandleCommand,
                     (ClientData)NULL, (Tcl_CmdDeleteProc*)NULL );
}

int
TclCommandManager::HandleCommand ( ClientData, Tcl_Interp* iInterp,
                                   int argc, char** argv ) {

  int rTcl = TCL_OK;
  TclCommandListener::TclCommandResult finalResult = ok;
  TclCommandManager& commandMgr = TclCommandManager::GetManager();

  std::list<Command*>::iterator tCommand;
  for ( tCommand = commandMgr.mlCommands.begin();
        tCommand != commandMgr.mlCommands.end(); ++tCommand ) {

    Command* command = *tCommand;
    if ( command->msCommand.compare( string(argv[0]) ) == 0 ) {

      // Check for right number of args. mcArgs should be argc-1
      // because argc also includes the command name.
      if ( command->mcArgs != argc - 1 ) {
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
        for ( tListener = command->mlListeners.begin();
              tListener != command->mlListeners.end(); ++tListener ) {

          try {

            // Give the listener a chance to listen to this command.
            TclCommandListener* listener = *tListener;
            TclCommandListener::TclCommandResult result =
              listener->ListenToTclCommand( argv[0], argc, argv );

            // If they gave us a return format string, convert the
            // format string and values string into a tcl object, which
            // could be a single value object or a list.
            if ( listener->sReturnFormat.length() != 0 ) {
              stringstream sFormat( listener->sReturnFormat );
              stringstream sValues( listener->sReturnValues );
              Tcl_Obj* rObj =
                commandMgr.ConvertFStringToTclObj( sFormat, sValues, iInterp );
              Tcl_SetObjResult( iInterp, rObj );
            }

            // If they gave us a result string, change that into a
            // string and set the result to that.
            if ( listener->sResult.length() != 0 ) {
              char* sResult = strdup( listener->sResult.c_str() );
              Tcl_SetResult( iInterp, sResult, TCL_VOLATILE );
              free( sResult );
            }

            // If the result was not ok, save that in the final
            // result. This way many objects can respond to the command
            // but an error will be returned if any of them return an
            // error.
            if ( result != ok ) {
              finalResult = result;
            }
          } catch ( runtime_error& e ) {
            finalResult = error;
            char* sError = strdup( e.what() );
            Tcl_SetResult( iInterp, sError, TCL_VOLATILE );
            free( sError );
          } catch (...) {
            finalResult = error;
            char sError[1024] = "Unknown error";
            Tcl_SetResult( iInterp, sError, TCL_VOLATILE );
          }
        }
      }
    }
  }

  switch ( finalResult ) {
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
  switch ( cFormat ) {

  case 'L': {
    Tcl_Obj* list = Tcl_NewListObj( 0, NULL );

    Tcl_Obj* newObject =
      ConvertFStringToTclObj( isFormat, isValues, iInterp );

    while ( newObject != NULL ) {
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
    if ( ERANGE == errno ) {
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
    if ( ERANGE == errno ) {
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
    if ( sValue[0] == '\"' ) {
      while ( sReturn.str().rfind('\"', sReturn.str().length()) == 0 ) {
        isValues >> sValue;
        sReturn << " " << sValue;
      }
    }

    // Strip the quotes if necessary and return it.
    string sReturn2 = sReturn.str();
    if ( sReturn2[0] == '\"' ) {
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
  for ( tCommand = mlCommands.begin();
        tCommand != mlCommands.end(); ++tCommand ) {

    Command* command = *tCommand;
    command->mlListeners.remove( &iListener );
  }
}

TclCommandListener::TclCommandResult
TclCommandManager::DoListenToTclCommand ( char* isCommand,
    int, char** iasArgv ) {

  // PrintAllCommands
  if ( 0 == strcmp( isCommand, "PrintAllCommands" ) ) {
    sResult = PrintAllCommands();
  }

  // GetArgc
  if ( 0 == strcmp( isCommand, "GetArgc" ) ) {

    // Return argc-1 because tcl's argc doesn't include the program name.
    stringstream ssReturnValues;
    ssReturnValues << mArgc - 1;
    sReturnValues = ssReturnValues.str();
    sReturnFormat = "i";
  }

  // GetArgv
  if ( 0 == strcmp( isCommand, "GetArgv" ) ) {

    stringstream ssFormat;
    stringstream ssResult;
    ssFormat << "L";

    for ( int nArgc = 1; nArgc < mArgc; nArgc++ ) {

      ssFormat << "s";
      ssResult << "\"" << mArgv[nArgc] << "\" ";
    }
    ssFormat << "l";

    sReturnFormat = ssFormat.str();
    sReturnValues = ssResult.str();
  }

  // DebugOutput message
  if ( 0 == strcmp( isCommand, "DebugOutput" ) ) {

    string sMessage = iasArgv[1];
    DebugOutput( << sMessage );
  }

  return ok;
}


string
TclCommandManager::PrintAllCommands () {

  stringstream ssResult;

  std::list<Command*>::iterator tCommand;
  for ( tCommand = mlCommands.begin();
        tCommand != mlCommands.end(); ++tCommand ) {

    Command* command = *tCommand;

    ssResult << command->msCommand << " " << command->msArgHints << endl;
    ssResult << command->msDescription << endl << endl;
  }

  return ssResult.str();
}

string
TclCommandManager::SendCommand ( string isCommand ) {

  if ( mInterp ) {

    char* sCommand = strdup( isCommand.c_str() );
    int rTcl = Tcl_Eval( mInterp, sCommand );
    const char* sTclResult = Tcl_GetStringResult( mInterp );
    stringstream ssError;
    ssError << "Error on cmd: \"" << sCommand << "\", " << sTclResult;
    free( sCommand );
    if ( TCL_OK != rTcl ) {
      throw runtime_error( ssError.str() );
    }
    return string(sTclResult);
  } else {
    return "";
  }
}

void
TclCommandManager::DoTclEvent () {

  if ( mInterp ) {
    Tcl_DoOneEvent( TCL_ALL_EVENTS | TCL_DONT_WAIT );
  }
}

void
TclCommandManager::SetCommandLineParameters ( int iArgc, char** iArgv ) {

  mArgc = iArgc;
  mArgv = iArgv;
}

string
TclCommandManager::RunTclScript ( char* ifnScript ) {

  if ( mInterp ) {

    int rTcl = Tcl_EvalFile( mInterp, ifnScript );
    const char* sTclResult = Tcl_GetStringResult( mInterp );
    if ( TCL_OK != rTcl ) {
      DebugOutput( << "Error on EvalFile: \"" << ifnScript << "\", "
                   << sTclResult );
    }
    return string(sTclResult);
  } else {
    return "";
  }

}

int
TclCommandManager::ConvertArgumentToInt ( std::string isArg ) {

  int theInteger = strtol(isArg.c_str(), (char**)NULL, 10);
  if ( ERANGE == errno ) {
    string sResult = "non-integer value";
    throw runtime_error( sResult );
  }

  return theInteger;
}

float
TclCommandManager::ConvertArgumentToFloat ( std::string isArg ) {

  float theFloat = strtod(isArg.c_str(), (char**)NULL);
  if ( ERANGE == errno ) {
    string sResult = "non-float value";
    throw runtime_error( sResult );
  }

  return theFloat;
}

bool
TclCommandManager::ConvertArgumentToBoolean ( std::string isArg ) {

  bool theBoolean;
  if ( isArg == "true" ||
       isArg == "1" ) {
    theBoolean = true;
  } else if ( isArg == "false" ||
              isArg == "0" ) {
    theBoolean = false;
  } else {
    string sResult;
    sResult = "should be true, 1, false, or 0";
    throw runtime_error( sResult );
  }

  return theBoolean;
}

std::string
TclCommandManager::ConvertBooleanToReturnValue ( bool ib ) {

  stringstream ssValue;
  ssValue << (int)ib;
  return ssValue.str();
}

