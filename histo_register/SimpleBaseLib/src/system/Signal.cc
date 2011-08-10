// Licensed under MIT license; see license.txt.

#include <sbl/system/Signal.h>
#include <sbl/core/Command.h>
#ifdef WIN32
    #include <windows.h> 
#else
    #include <signal.h>
#endif
namespace sbl {


// control-C handler
void handleCtrlC( int signal ) {
    static int s_count = 0;
    s_count++;
    if (s_count >= 3) {
        warning( "three Ctrl-C events; shutting down" );
        exit( 1 );
    } 
    disp( 1, "cancelling..." );
    setCancelCommand( true );
}


// windows control handler
#ifdef WIN32
BOOL ctrlHandler( DWORD fdwCtrlType ) {
    if (fdwCtrlType == CTRL_C_EVENT) {
        handleCtrlC( 0 );
        return TRUE;
    } 
    return FALSE;
} 
#endif


//-------------------------------------------
// INIT / CLEAN-UP
//-------------------------------------------


// register commands, etc. defined in this module
void initSignal() {
#ifdef WIN32
    if (SetConsoleCtrlHandler( (PHANDLER_ROUTINE) ctrlHandler, TRUE ) == false) {
        warning( "error setting control handler" );
    }
#else
    signal( SIGINT, handleCtrlC );
#endif
}


} // end namespace sbl

