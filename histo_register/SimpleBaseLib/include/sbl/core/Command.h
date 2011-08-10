#ifndef _SBL_COMMAND_H_
#define _SBL_COMMAND_H_
#include <sbl/core/Config.h>
#include <sbl/core/Table.h>
namespace sbl {


/*! \file Command.h
    \brief The Command module manages the set of currently registered command functions.
    The command functions can be executed from scripts, the command-line, and optionally
    a graphical user interface.  Each command function has a string name and a description
    determined from the C++ function name.  All command function take a Config object as their
    only argument.  These Config objects are populated with parameters for the command.
*/


// register commands, etc. defined in this module
void initCommand();


//-------------------------------------------
// COMMAND REGISTRY
//-------------------------------------------


/// register a command name and callback; the command description is determined from the callback function name
#define registerCommand( name, callback ) registerCommandInternal( name, callback, #callback ) 
void registerCommandInternal( const String &name, void (*callback)(Config &conf), const String &functionName );


/// get the description of a command from its (short) name
String commandDescription( const String &name );


/// perform command name tab completion starting with the given prefix
String tabCompleteCommand( const String &prefix );


/// generate description from CamelCase name: "someFuncName" -> "some func name"
String descriptionFromName( const String &name );


//-------------------------------------------
// COMMAND HISTORY
//-------------------------------------------


// get a command from the history, at the specified offset from the current list position
String nextHistoryCommand( int offset );


//-------------------------------------------
// COMMAND EXECUTION
//-------------------------------------------


// if true, use multi-pass config with interactive (GUI) editing
void useInteractiveCommand( bool useInteractive );


/// execute a command with command arguments in a Config object
void execCommand( const String &name, Config &conf );


/// execute a command with arguments included in the command string (space separated)
void execCommand( const String &nameAndArgs, bool addToHistory );


/// execute a command of the form: "cmdname&param1=val1&param2=val2&..."
void execURLCommand( const String &request );


/// run a script containing a sequence of commands
void runScript( const String &fileName );


//-------------------------------------------
// CHECK FOR COMMAND CANCELLING
//-------------------------------------------


/// this callback will be called every time a long-running command calls checkCommandEvents
void setCommandEventCallback( void (*commandEventCallback)() );


/// check for events (e.g. update GUI) and return true if command cancel
bool checkCommandEvents();


/// indicate whether to cancel the currently running command (if any)
void setCancelCommand( bool cancel );


//-------------------------------------------
// CLEANUP MANAGEMENT
//-------------------------------------------


/// register a clean-up function (to be called when program terminates)
void registerCleanUp( void (*callback)() );


/// run all registered clean-up functions (call when program terminates)
void runCleanUp();


} // end namespace sbl
#endif // _SBL_COMMAND_H_

