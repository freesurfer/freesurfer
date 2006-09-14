#include <stdlib.h>
#include "string_fixed.h"
#include <iostream>
#include <stdexcept>

extern "C" {
#define USE_NON_CONST
#include <tcl.h>
#undef USE_NON_CONST
}

#include "TclCommandManager.h"
#include "Scuba-impl.h"

#define Assert(x,s)   \
  if(!(x)) { \
  stringstream ss; \
  ss << "Line " << __LINE__ << ": " << s; \
  cerr << ss.str().c_str() << endl; \
  throw runtime_error( ss.str() ); \
  }

char* Progname = "test_TclCommandManager";

using namespace std;

class TestListenerCounter : public TclCommandListener {

protected:
  int mID;
  
  struct ListenCount {
    char msCommand[1024];
    int mcHeard;
  };

  list<ListenCount*> mlListenCounts;

public:
  TestListenerCounter( int const iID ) : mID( iID ) { }

  TclCommandResult
  DoListenToTclCommand ( char* isCommand, int iArgc, char** iArgv ) {
    list<ListenCount*>::iterator tListenCount;
    for( tListenCount = mlListenCounts.begin();
	 tListenCount != mlListenCounts.end(); ++tListenCount ) {
      ListenCount* listenCount = *tListenCount;
      if( 0 == strcmp( isCommand, listenCount->msCommand ) ) {
	listenCount->mcHeard++;
      }
    }
    return ok;
  }

  void CountCommand ( char* isCommand ) {
    ListenCount* listenCount = new ListenCount;
    strcpy( listenCount->msCommand, isCommand );
    listenCount->mcHeard = 0;
    mlListenCounts.push_back( listenCount );
  }

  int GetCountForCommand ( char const* isCommand ) {
    list<ListenCount*>::iterator tListenCount;
    for( tListenCount = mlListenCounts.begin();
	 tListenCount != mlListenCounts.end(); ++tListenCount ) {
      ListenCount* listenCount = *tListenCount;
      if( 0 == strcmp( isCommand, listenCount->msCommand ) ) {
	return listenCount->mcHeard;
      }
    }
    return 0;
  }
};

class TestListenerReturner : public TclCommandListener {

public:
  TclCommandResult
  DoListenToTclCommand( char *isCommand, int iArgc, char** iArgv ) {

    // These test the functionality of returning Tcl objects to the
    // Tcl context. Used in conjunction with test_TclCommandManager.tcl.
    if( 0 == strcmp( isCommand, "ReturnSingleInt" ) ) {
      sReturnFormat = "i";
      sReturnValues = "5";
    } else if( 0 == strcmp( isCommand, "ReturnSingleFloat" ) ) {
      sReturnFormat = "f";
      sReturnValues = "5.5";
    } else if( 0 == strcmp( isCommand, "ReturnSingleString" ) ) {
      sReturnFormat = "s";
      sReturnValues = "hello";
    } else if( 0 == strcmp( isCommand, "ReturnSingleTwoWordString" ) ) {
      sReturnFormat = "s";
      sReturnValues = "\"hello world\"";
    } else if( 0 == strcmp( isCommand, "ReturnSingleLongString" ) ) {
      sReturnFormat = "s";
      sReturnValues = "\"to neither love nor reverence will thou be tied\"";
    } else if( 0 == strcmp( isCommand, "ReturnSingleList" ) ) {
      sReturnFormat = "Lifsl";
      sReturnValues = "5 5.5 hello";
    } else if( 0 == strcmp( isCommand, "ReturnNestedList" ) ) {
      sReturnFormat = "LifsLifsll";
      sReturnValues = "5 5.5 hello 6 6.6 \"hello world\"";
    } else if( 0 == strcmp( isCommand, "ReturnNestedList2" ) ) {
      sReturnFormat = "LLsslLssll";
      sReturnValues = "\"Label 1\" \"Value 1\" \"Label 2\" \"Value 2\"";
    } else if( 0 == strcmp( isCommand, "ReturnMessage" ) ) {
      sResult = "This is a result string.";
    } else if( 0 == strcmp( isCommand, "ReturnError" ) ) {
      sResult = "This is an error string.";
      return error;
    } else if( 0 == strcmp( isCommand, "TestThrow" ) ) {
      throw runtime_error( "throw" );
    }
    return ok;
  }
};

class TclCommandManagerTester {
public:
  void Test( Tcl_Interp* iInterp ) {
    
    int rTcl = TCL_OK;
    int const kzListeners = 10;
    int const kzCommands = 10;
    int const kNumberOfCallsToMake = 10;
    
    TestListenerCounter* aListener[kzListeners];
    for( int nListener = 0; nListener < kzListeners; nListener++ ) {
      aListener[nListener] = new TestListenerCounter( nListener );
    }
    
    char asCommandNames[1024][kzCommands]; 
    for( int nCommand = 0; nCommand < kzCommands; nCommand++ ) {
      sprintf( asCommandNames[nCommand], "command%d", nCommand );
    }
    
    TclCommandManager& commandMgr = TclCommandManager::GetManager();
    commandMgr.SetOutputStreamToCerr();
    
    
    // Add half of the commands for each even listener before Start()
    // is called. Don't add them to the odd listeners.
    for( int nListener = 0; nListener < kzListeners; nListener++ ) {
      if( nListener % 2 == 0 ) {
	for( int nCommand = 0; nCommand < kzCommands/2; nCommand++ ) {
	  aListener[nListener]->CountCommand( asCommandNames[nCommand] );
	  commandMgr.AddCommand( *aListener[nListener], 
				 asCommandNames[nCommand], 3, 
				 "arg1 arg2 arg3", "" );
	}
      }
    }
    

    // Call start.
    commandMgr.Start( iInterp );
    
    
    // Call those commands many times each.
    for( int nCommand = 0; nCommand < kzCommands/2; nCommand++ ) {
      char sCommandAndArgs[1024];
      sprintf( sCommandAndArgs, "%s arg1 arg2 arg3", 
	       asCommandNames[nCommand] );
      for( int nCall = 0; nCall < kNumberOfCallsToMake; nCall++ ) {
	rTcl = Tcl_Eval( iInterp, sCommandAndArgs );
	Assert( TCL_OK == rTcl, "Tcl_Eval returned not TCL_OK" );
      }
    }
    
    
    // Check the counts in all listeners. All evens should have the
    // full number of calls, all odds should have 0 counts.
    for( int nListener = 0; nListener < kzListeners; nListener++ ) {
      for( int nCommand = 0; nCommand < kzCommands/2; nCommand++ ) {
	int cCalls = 
	  aListener[nListener]->GetCountForCommand( asCommandNames[nCommand]);
	if( nListener % 2 == 0 ) {
	  if( kNumberOfCallsToMake != cCalls ) {
	    throw logic_error( "Count mismatch" );
	  }
	} else {
	  if( 0 != cCalls ) {
	    throw logic_error( "Count mismatch" );
	  }
	}
      }
    }
    
    
    // Add some commands for each even listener now that Start() has
    // called. Don't add them to the odd listeners.
    for( int nListener = 0; nListener < kzListeners; nListener++ ) {
      if( nListener % 2 == 0 ) {
	for( int nCommand = kzCommands/2; nCommand < kzCommands; nCommand++ ) {
	  aListener[nListener]->CountCommand( asCommandNames[nCommand] );
	  commandMgr.AddCommand( *aListener[nListener], 
				 asCommandNames[nCommand], 
				 3, "arg1 arg2 arg3", "" );
	}
      }
    }
    
    
    // Call those commands many times each.
    for( int nCommand = kzCommands/2; nCommand < kzCommands; nCommand++ ) {
      char sCommandAndArgs[1024];
      sprintf( sCommandAndArgs, "%s arg1 arg2 arg3", 
	       asCommandNames[nCommand] );
      for( int nCall = 0; nCall < kNumberOfCallsToMake; nCall++ ) {
	rTcl = Tcl_Eval( iInterp, sCommandAndArgs );
	Assert( TCL_OK == rTcl, "Tcl_Eval returned not TCL_OK" );
      }
    }
    
    
    // Check the counts in all listeners. All evens should have 100
    // counts, all odds should have 0 counts.
    for( int nListener = 0; nListener < kzListeners; nListener++ ) {
      for( int nCommand = 0; nCommand < kzCommands; nCommand++ ) {
	int cCalls = 
	  aListener[nListener]->GetCountForCommand( asCommandNames[nCommand]);
	if( nListener % 2 == 0 ) {
	  if( kNumberOfCallsToMake != cCalls ) {
	    throw logic_error( "Count mismatch" );
	  }
	} else {
	  if( 0 != cCalls ) {
	    throw logic_error( "Count mismatch" );
	  }
	}
      }
    }


    // Set our command line stuff. We'll test it in the script.
    int argc = 5;
    char **argv = (char**) calloc( argc, sizeof(char*) ); 
    argv[0] = (char*) calloc( 256, sizeof(char) );
    strcpy( argv[0], "test_TclCommandManager" );
    argv[1] = (char*) calloc( 256, sizeof(char) );
    strcpy( argv[1], "param 1" );
    argv[2] = (char*) calloc( 256, sizeof(char) );
    strcpy( argv[2], "param 2" );
    argv[3] = (char*) calloc( 256, sizeof(char) );
    strcpy( argv[3], "param 3" );
    argv[4] = (char*) calloc( 256, sizeof(char) );
    strcpy( argv[4], "param 4" );
    commandMgr.SetCommandLineParameters( argc, argv );

    
    // Add the return testing commands.
    TestListenerReturner* listener = new TestListenerReturner();
    commandMgr.AddCommand( *listener, "ReturnSingleInt", 0, "", "" );
    commandMgr.AddCommand( *listener, "ReturnSingleFloat", 0, "", "" );
    commandMgr.AddCommand( *listener, "ReturnSingleString", 0, "", "" );
    commandMgr.AddCommand( *listener, "ReturnSingleTwoWordString", 0, "", "" );
    commandMgr.AddCommand( *listener, "ReturnSingleLongString", 0, "", "" );
    commandMgr.AddCommand( *listener, "ReturnSingleList", 0, "", "" );
    commandMgr.AddCommand( *listener, "ReturnNestedList", 0, "", "" );
    commandMgr.AddCommand( *listener, "ReturnNestedList2", 0, "", "" );
    commandMgr.AddCommand( *listener, "ReturnMessage", 0, "", "" );  
    commandMgr.AddCommand( *listener, "ReturnError", 0, "", "" );
    commandMgr.AddCommand( *listener, "TestThrow", 0, "", "" );
  
    // Run the script that will test tcl return stuff.
    rTcl = Tcl_EvalFile( iInterp, "test_TclCommandManager.tcl" );
    const char* sTclResult = Tcl_GetStringResult( iInterp );
    if( 0 != strcmp( sTclResult, "" ) ) {
      cerr << "result not null" << endl;
      throw logic_error( iInterp->result );
    }
    if( rTcl != TCL_OK ) {
      cerr << "rTcl not OK" << endl;
      throw logic_error( iInterp->result );
    }
    
    
    // Delete all the listeners and then make sure the the manager's
    // command lists are empty.
    for( int nListener = 0; nListener < kzListeners; nListener++ ) {
      delete aListener[nListener];
    }
    delete listener;

    // Don't check the commands that are listened to by the
    // TclCommandManager, as these are 'global'.
    std::list<TclCommandManager::Command*>::iterator tCommand;
    for( tCommand = commandMgr.mlCommands.begin(); 
	 tCommand != commandMgr.mlCommands.end(); ++tCommand ) {
      TclCommandManager::Command* command = *tCommand;
      if( 0 != command->mlListeners.size() &&
	  command->msCommand != "PrintAllCommands" && 
	  command->msCommand != "GetArgc" && 
	  command->msCommand != "GetArgv" && 
	  command->msCommand != "DebugOutput" ) {
	cerr << "Not all listeners removed for command " 
	     << command->msCommand << " size = " 
	     << command->mlListeners.size() << endl;
      }
    }
    
  }
};

int main ( int iArgc, char** iArgv ) {

  cerr << "Beginning test" << endl;
 
  try {

    Tcl_Interp* interp = Tcl_CreateInterp();
    Assert( interp, "Tcl_CreateInterp returned null" );
    
    int rTcl = Tcl_Init( interp );
    Assert( TCL_OK == rTcl, "Tcl_Init returned not TCL_OK" );

    for( int nTrial = 0; nTrial < 1; nTrial++ ) {
      
      TclCommandManagerTester tester0;
      tester0.Test( interp );
    }
  }

  catch ( char* msg ) {
    cerr << msg << " failed" << endl;
    exit( 1 );
  }
  catch( logic_error& e ) {
    cerr << "failed: " << e.what() << endl;
    exit( 1 );
  }
  catch(...) {
    cerr << "failed." << endl;
    exit( 1 );
  }

  cerr << "Passed." << endl;

  return 0;
}
