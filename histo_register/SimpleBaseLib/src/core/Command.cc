// Licensed under MIT license; see license.txt.

#include <sbl/core/Command.h>
#include <sbl/core/Pointer.h>
#include <sbl/core/Display.h>
#include <sbl/core/PathConfig.h> // the following are for command history
#include <sbl/math/MathUtil.h>
#include <sbl/system/FileSystem.h>
#include <sbl/system/TimeUtil.h>
#include <sbl/system/Timer.h>
#ifdef USE_GUI
    #include <sbl/gui/ConfigEditor.h>
#endif
namespace sbl {


//-------------------------------------------
// COMMAND REGISTRY
//-------------------------------------------


// information about a command to be stored in command registry
struct Command {
    String name;
    String functionName;
    String description;
    void (*callback)( Config &conf );
};


// get list of registered commands
Array<Command> &commandList() {
    static aptr<Array<Command> > s_commands;
    if (s_commands.get() == NULL) {
        s_commands.reset( new Array<Command>() ); 
    }
    return *s_commands;
}


/// register a command name and callback; the command description is determined from the callback function name
void registerCommandInternal( const String &name, void (*callback)(Config &conf), const String &functionName ) {

    // create new command object
    Command *cmd = new Command;
    cmd->name = name;
    cmd->functionName = functionName;
    cmd->callback = callback;
    commandList().append( cmd );

    // generate description: "someFuncName" -> "some func name"
    cmd->description = descriptionFromName( functionName );
}


/// get the description of a command from its (short) name
String commandDescription( const String &name ) {
    String result;
    for (int i = 0; i < commandList().count(); i++) {
        const Command &cmd = commandList()[ i ];
        if (cmd.name == name)
            result = cmd.description;
    }
    return result;
}


/// perform command name tab completion starting with the given prefix
String tabCompleteCommand( const String &prefix ) {
    String result = prefix;
    int matchCount = 0;
    bool contSet = false;
    if (prefix.length() == 0)
        return result;

    // loop over commands
    const Array<Command> &cmdList = commandList();
    for (int i = 0; i < cmdList.count(); i++) {
        const Command &cmd = cmdList[ i ];

        // if prefix matches
        if (cmd.name.startsWith( prefix )) {
            matchCount++;

            // if continuation has been set
            if (contSet) {

                // find first non-matching char in continuation
                int noMatchPos = 0;
                while (noMatchPos < result.length() && noMatchPos < cmd.name.length()) {
                    if (result.get( noMatchPos ) != cmd.name.get( noMatchPos ))
                        break;
                    noMatchPos++;
                }

                // keep only up until no-match position
                if (noMatchPos < result.length())
                    result = result.leftOf( noMatchPos );

            // create continuation
            } else {
                result = cmd.name; //cmd.name.leftOf( prefix.length() - 1 );
                contSet = true;
            }
        }
    }

    // display matches
    if (matchCount > 20) {
        disp( 0, "command complete: %d matches", matchCount );
    } else if (matchCount > 1) {
        disp( 0, "command complete:" );
        for (int i = 0; i < cmdList.count(); i++) {
            const Command &cmd = cmdList[ i ];
            if (cmd.name.startsWith( prefix )) {
                disp( 1, "%s: %s", cmd.name.c_str(), cmd.description.c_str() );
            }
        }
    } else if (matchCount == 0) {
        disp( 0, "command complete: no matches" );
    }
    return result;
}


/// generate description from CamelCase name: "someFuncName" -> "some func name"
String descriptionFromName( const String &name ) {
    int pos = 0;
    String description( 0, name.length() + 100 );
    for (int i = 0; i < name.length(); i++) {
        unsigned short c = name.get( i );

        // if lower case, just copy it
        if ('a' <= c && c <= 'z') {
            description.set( pos++, c );

        // if upper case, insert space and convert to upper case
        } else if ('A' <= c && c <= 'Z') {
            description.set( pos++, ' ' );
            description.set( pos++, c - 'A' + 'a' );
        }
    }
    return description;
}


//-------------------------------------------
// DEFAULT COMMAND CONFIG
//-------------------------------------------


// the path where we will store a list of recent commands and (for GUI apps) recent command parameters
String commandHistoryPath() {
    static String s_commandHistoryPath;
    static bool s_commandHistoryPathInit = false;
    if (s_commandHistoryPathInit == false) {
        if (pathConfig().entryExists( "commandHistoryPath" ))
            s_commandHistoryPath = pathConfig().readString( "commandHistoryPath" );
        else 
            s_commandHistoryPath = "history";
        createDir( s_commandHistoryPath );
        s_commandHistoryPathInit = true;
    }
    return s_commandHistoryPath;
}


// load the last parameters for this command
void loadRecentCommandConfig( const String &name, Config &conf ) {
    String fileName = commandHistoryPath() + "/" + name + ".conf";
    if (fileExists( fileName )) 
        conf.load( fileName );
}


// save the parameters for this command
void saveRecentCommandConfig( const String &name, Config &conf ) {
    String fileName = commandHistoryPath() + "/" + name + ".conf";
    conf.save( fileName );
}


//-------------------------------------------
// COMMAND HISTORY
//-------------------------------------------


// access list of recent commands
Array<String> &commandHistory() {
    static aptr< Array<String> > s_commandHistory;
    if (s_commandHistory.get() == NULL)
        s_commandHistory.reset( new Array<String>() );
    return *s_commandHistory;
}


// the current offset within the list of recent commands
static int g_commandHistoryIndex = 0;


// add a command to the list of recent commands
void addToCommandHistory( const String &command ) {
    commandHistory().appendCopy( command );
    if (commandHistory().count() > 100)
        commandHistory().remove( 0 );
    g_commandHistoryIndex = commandHistory().count(); // one past the end
}


// load list of recent commands from file
void loadCommandHistory() {
    String fileName = commandHistoryPath() + "/history.txt";
    commandHistory() = loadStrings( fileName, false );
    g_commandHistoryIndex = commandHistory().count(); // one past the end
}


// save list of recent commands to file
void saveCommandHistory() {
    String fileName = commandHistoryPath() + "/history.txt";
    File file( fileName, FILE_WRITE, FILE_TEXT );
    if (file.openSuccess()) {
        for (int i = 0; i < commandHistory().count(); i++) {
            file.writeRawString( commandHistory()[ i ] );
            file.writeRawString( "\n" );
        }
    }
}


// load command history and register function to save history on exit
void initCommandHistory() {
    loadCommandHistory();
    registerCleanUp( saveCommandHistory );    
}


// get a command from the history, at the specified offset from the current list position
String nextHistoryCommand( int offset ) {
    String command;
    int count = commandHistory().count();
    if (count) {
        g_commandHistoryIndex += offset;
        g_commandHistoryIndex = bound( g_commandHistoryIndex, 0, count ); // note: allow to be one past the end
        if (g_commandHistoryIndex < count)
            command = commandHistory()[ g_commandHistoryIndex ];
    }
    return command;
}


//-------------------------------------------
// COMMAND EXECUTION
//-------------------------------------------


// if true, use multi-pass config with interactive (GUI) editing
bool g_useInteractive = false;


// if true, use multi-pass config with interactive (GUI) editing
void useInteractiveCommand( bool useInteractive ) {
    g_useInteractive = useInteractive;
}


/// run a command in the user-interface with interactive parameter editing
void interactiveExecCommand( const Command &cmd, Config &conf ) {
#ifdef USE_GUI

    // run a first-pass of the command to see what parameters it wants
    conf.enableDualPass();
    conf.setInitialPass( true );
    cmd.callback( conf );
    conf.setInitialPass( false );

    // if the command didn't check it for the initial pass, then we've already finished running the command
    if (conf.checkedInitialPass() == false)
        return;

    // if all the arguments were specified in the initial config, just run the command now
    if (conf.missingValue() == false) {
        cmd.callback( conf );
        return;
    }

    // at this point, we have some unspecified parameters that the command would like, 
    // so we will display a dialog asking for those parameters; 
    // first we will load the parameters used last time (if any)
    Config recentConf;
    loadRecentCommandConfig( cmd.name, recentConf );
    
    // check whether all the entries have the same names
    if (recentConf.entryCount() && recentConf.entryCount() == conf.entryCount()) {
        bool allMatch = true;
        for (int i = 0; i < conf.entryCount(); i++) {
            if (conf.entry( i ).name != recentConf.entry( i ).name) 
                allMatch = false;
        }

        // if the recent config has the same set/ordering of paramters, copy them over the command-supplied defaults
        if (allMatch) 
            conf.updateFrom( recentConf );
    }

    // now, show the config editor so that the user can change the parameters
    bool ok = editConfig( conf, String( "Command: " ) + cmd.description );
    if (ok) {

        // save the user's parameter values to use as defaults next time
        if (conf.entryCount()) {
            saveRecentCommandConfig( cmd.name, conf );
            String args;
            for (int i = 0; i < conf.entryCount(); i++) {
                args += conf.entry( i ).value;
                if (i + 1 < conf.entryCount())
                    args += " ";
            }
            disp( 0, "%s %s", cmd.name.c_str(), args.c_str() );
        }

        // run the command
        cmd.callback( conf );
    }
#endif // USE_GUI
}


/// execute a command with command arguments in a Config object
void execCommand( const String &name, Config &conf ) {
    setCancelCommand( false );
    bool found = false;
    for (int i = 0; i < commandList().count(); i++) {
        const Command &cmd = commandList()[ i ];
        if (cmd.name == name || cmd.functionName == name) {

            // compute string of command args (for display purposes)
            String args = join( conf.commandArgs(), " " );
            Array<String> namedArgs;
            for (int j = 0; j < conf.entryCount(); j++) 
                namedArgs.appendCopy( conf.entry( j ).name + "=" + conf.entry( j ).value );
            if (namedArgs.count()) {
                if (args.length())
                    args += " ";
                args += join( namedArgs, " " );
            }

            // prepare to execute command
            found = true;
            Timer timer( true );

            // execute the command
            if (g_useInteractive) {
                interactiveExecCommand( cmd, conf );
            } else {
                disp( 0, "%s %s", name.c_str(), args.c_str() );
                cmd.callback( conf );
            }

            // display command timing
            timer.stop();
            bool cancelled = checkCommandEvents();
            if (cancelled)
                status( "command cancelled; duration: %5.3fs", timer.timeSum() );
            else
                status( "command done; duration: %5.3fs", timer.timeSum() );
            printf( "\n" ); // if not in gui, add new line after command status
        }
    }
    if (found == false) {
        String nameAsScript = name + ".txt";
        Array<String> scriptNames = dirFileList( "scripts", "", "txt" );
        for (int i = 0; i < scriptNames.count(); i++) {
            if (scriptNames[ i ] == nameAsScript) {
                runScript( String( "scripts/" ) + nameAsScript );
                found = true;
                break;
            }
        }
    }
    if (found == false) {
        disp( 0, "command not found: %s", name.c_str() );
    }
}


/// execute one or more commands (semicolon-separated)
void execCommand( const String &commandText, bool addToHistory ) {

    // add to command history
    if (addToHistory && commandText.strip() != "x")
        addToCommandHistory( commandText );

    // exec each command (separated by semicolons)
    Array<String> cmdList = commandText.split( ";" );
    for (int i = 0; i < cmdList.count(); i++) {
        const String &nameAndArgs = cmdList[ i ].strip();
        if (nameAndArgs.length()) {

            // get command args
            Config conf;
            Array<String> split = nameAndArgs.split( " " );
            Array<String> commandArgs;
            for (int i = 1; i < split.count(); i++) // copy all but first
                commandArgs.append( new String( split[ i ] ));
            conf.setCommandArgs( commandArgs );
            conf.allowNewEntries( true );
            conf.allowDefault( true );

            // execute command
            execCommand( split[ 0 ], conf );
        }
    }
}


// execute a command of the form: "cmdname&param1=val1&param2=val2&..."
void execURLCommand( const String &request ) {
    Config conf;
    Array<String> split = request.split( "&" );
    for (int i = 1; i < split.count(); i++) {
        Array<String> subSplit = split[ i ].split( "=" );
        if (subSplit.count() == 2) {
            conf.writeString( subSplit[ 0 ].c_str(), subSplit[ 1 ] );
            disp( 1, "%s = %s", subSplit[ 0 ].c_str(), subSplit[ 1 ].c_str() );
        }
    }
    conf.allowDefault( true );

    // execute the command (command name is first entry in split)
    execCommand( split[ 0 ].c_str(), conf );
}


/// run a script containing a sequence of commands
void runScript( const String &fileName ) {
    Array<String> lines = loadStrings( fileName, false );
    for (int i = 0; i < lines.count(); i++) {
        const String &line = lines[ i ].strip();
        if (line == "stop")
            break;
        if (line.length()) 
            execCommand( line, false );
    }
}


//-------------------------------------------
// CHECK FOR COMMAND CANCELLING
//-------------------------------------------


// true if current command should halt itself
bool g_cancelCommand = false;


// this callback will be called every time a long-running command calls checkCommandEvents
void (*g_commandEventCallback)() = NULL;


/// this callback will be called every time a long-running command calls checkCommandEvents
void setCommandEventCallback( void (*commandEventCallback)() ) {
    g_commandEventCallback = commandEventCallback; 
}


/// check for events (e.g. update GUI) and return true if command cancel
bool checkCommandEvents() {
    if (g_commandEventCallback)
        g_commandEventCallback();
    return g_cancelCommand;
}


/// indicate whether to cancel the currently running command (if any)
void setCancelCommand( bool cancel ) {
    g_cancelCommand = cancel;
}


//-------------------------------------------
// CLEANUP MANAGEMENT
//-------------------------------------------


// a clean-up function
class CleanUpCallback {
public:
    void (*callback)();
};


// a list of module clean-up functions
Array<CleanUpCallback> &cleanUpCallbacks() {
    static aptr<Array<CleanUpCallback> > s_cleanUpCallbacks;
    if (s_cleanUpCallbacks.get() == NULL) {
        s_cleanUpCallbacks.reset( new Array<CleanUpCallback>() ); 
    }
    return *s_cleanUpCallbacks;
}


// run all registered clean-up functions (call when program terminates)
void runCleanUp() {
    for (int i = 0; i < cleanUpCallbacks().count(); i++) 
        cleanUpCallbacks()[ i ].callback();
}


// register a clean-up function (to be called when program terminates)
void registerCleanUp( void (*callback)() ) {
    CleanUpCallback *cleanUp = new CleanUpCallback;
    cleanUp->callback = callback;
    cleanUpCallbacks().append( cleanUp );
}


//-------------------------------------------
// COMMAND-RELATED COMMANDS
//-------------------------------------------


// run a basic test script (list of commands)
void runScript( Config &conf ) {

    // get command parameters
    String fileName = addDataPath( conf.readString( "fileName", dataPath() + "script.txt" ) );
    if (conf.initialPass())
        return;

    // check whether the file exists
    if (fileExists( fileName ) == false) {
        warning( "script not found: %s", fileName.c_str() );

    // if found, run the script
    } else {
        runScript( fileName );
    }
}


// show all registered commands
void showCommandList( Config &conf ) {
    for (int i = 0; i < commandList().count(); i++) {
        const Command &cmd = commandList()[ i ];
        disp( 1, "%s: %s", cmd.name.c_str(), cmd.description.c_str() );
    }
}


// start an interactive command prompt
// fix(later): add tab-completion
void startCommandPrompt( Config &conf ) {
    const int bufLen = 1000;
    char buf[ bufLen ];
    int bufPos = 0;
    printf( ">> " );
    while (1) {
        int c = getchar();
        if (c >= ' ' && c <= '~' && bufPos < bufLen) 
            buf[ bufPos++ ] = c;
        if (c == 10 || c == 13) {
            buf[ bufPos ] = 0;
            if (String( buf ) == "x")
                break;
            execCommand( buf, true );
            bufPos = 0;
            printf( ">> " );
        }
    }
}


// cancel currently running command
void cancelCurrentCommands( Config &conf ) {
    setCancelCommand( true );
}


// set a parameter in the main config file
// fix(clean): move to MainConfig.cc
void setMainConfigParameter( Config &conf ) {
    const String &paramName = conf.readString( "paramName" );
    const String &paramValue = conf.readString( "paramValue" );
    if (conf.initialPass())
        return;
    pathConfig().writeString( paramName, paramValue );
}


//-------------------------------------------
// INIT / CLEAN-UP
//-------------------------------------------


// register commands, etc. defined in this module
void initCommand() {
    initCommandHistory();
    registerCommand( "run", runScript );
    registerCommand( "cmdlist", showCommandList );
    registerCommand( "cmdprompt", startCommandPrompt );
    registerCommand( "cancel", cancelCurrentCommands );    
    registerCommand( "set", setMainConfigParameter );
}


} // end namespace sbl

