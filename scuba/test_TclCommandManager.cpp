#include <stdlib.h>
#include <string>
#include <iostream>
#include <tcl.h>
#include <stdexcept>
#include "TclCommandManager.h"

#define Assert(x,s)   if(!(x)) { throw logic_error( s ); }

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

  void DoListenToTclCommand ( char* isCommand, int iArgc, char** iArgv ) {
    list<ListenCount*>::iterator tListenCount;
    for( tListenCount = mlListenCounts.begin();
	 tListenCount != mlListenCounts.end(); ++tListenCount ) {
      ListenCount* listenCount = *tListenCount;
      if( 0 == strcmp( isCommand, listenCount->msCommand ) ) {
	listenCount->mcHeard++;
      }
    }
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
  void DoListenToTclCommand( char *isCommand, int iArgc, char** iArgv ) {

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
    } else if( 0 == strcmp( isCommand, "ReturnMessage" ) ) {
      sResult = "This is a result string.";
    }
  }
};

int main ( int iArgc, char** iArgv ) {

  int rTcl = TCL_OK;
  int const kzListeners = 100;
  int const kzCommands = 100;
  int const kNumberOfCallsToMake = 100;

  TestListenerCounter* aListener[kzListeners];
  for( int nListener = 0; nListener < kzListeners; nListener++ ) {
    aListener[nListener] = new TestListenerCounter( nListener );
  }

  char asCommandNames[1024][kzCommands]; 
  for( int nCommand = 0; nCommand < kzCommands; nCommand++ ) {
    sprintf( asCommandNames[nCommand], "command%d", nCommand );
  }

  try { 

    cerr << "TclCommandListener::GetManager()" << endl;
    TclCommandManager& commandMgr = TclCommandManager::GetManager();
    commandMgr.SetOutputStreamToCerr();
    
    cerr << "Tcl_CreateInterp()" << endl;
    Tcl_Interp* interp = Tcl_CreateInterp();
    Assert( interp, "Tcl_CreateInterp returned null" );
    
    cerr << "Tcl_Init()" << endl;
    rTcl = Tcl_Init( interp );
    Assert( TCL_OK == rTcl, "Tcl_Init returned not TCL_OK" );


    // Add half of the commands for each even listener before Start()
    // is called. Don't add them to the odd listeners.
    cerr << "Adding " << kzCommands/2 << " commands to even listeners "
      "before Start() is called" << endl;
    for( int nListener = 0; nListener < kzListeners; nListener++ ) {
      if( nListener % 2 == 0 ) {
	for( int nCommand = 0; nCommand < kzCommands/2; nCommand++ ) {
	  aListener[nListener]->CountCommand( asCommandNames[nCommand] );
	  commandMgr.AddCommand( *aListener[nListener], 
				 asCommandNames[nCommand] );
	}
      }
    }


    // Call start.
    cerr << "TclCommandManager::Start()" << endl;
    commandMgr.Start( interp );


    // Call those commands many times each.
    cerr << "Calling " << kzCommands/2 << " commands " 
	 << kNumberOfCallsToMake << " times" << endl;
    for( int nCommand = 0; nCommand < kzCommands/2; nCommand++ ) {
      char sCommandAndArgs[1024];
      sprintf( sCommandAndArgs, "%s arg1 arg2 arg3", 
	       asCommandNames[nCommand] );
      for( int nCall = 0; nCall < kNumberOfCallsToMake; nCall++ ) {
	rTcl = Tcl_Eval( interp, sCommandAndArgs );
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
    cerr << "Adding " << kzCommands/2 << " commands to even listeners "
      "after Start() is called" << endl;
    for( int nListener = 0; nListener < kzListeners; nListener++ ) {
      if( nListener % 2 == 0 ) {
	for( int nCommand = kzCommands/2; nCommand < kzCommands; nCommand++ ) {
	  aListener[nListener]->CountCommand( asCommandNames[nCommand] );
	  commandMgr.AddCommand( *aListener[nListener], 
				 asCommandNames[nCommand] );
	}
      }
    }


    // Call those commands many times each.
    cerr << "Calling " << kzCommands/2 << " commands " 
	 << kNumberOfCallsToMake << " times" << endl;
    for( int nCommand = kzCommands/2; nCommand < kzCommands; nCommand++ ) {
      char sCommandAndArgs[1024];
      sprintf( sCommandAndArgs, "%s arg1 arg2 arg3", 
	       asCommandNames[nCommand] );
      for( int nCall = 0; nCall < kNumberOfCallsToMake; nCall++ ) {
	rTcl = Tcl_Eval( interp, sCommandAndArgs );
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

    
    // Add the return testing commands.
    TestListenerReturner listener;
    commandMgr.AddCommand( listener, "ReturnSingleInt" );
    commandMgr.AddCommand( listener, "ReturnSingleFloat" );
    commandMgr.AddCommand( listener, "ReturnSingleString" );
    commandMgr.AddCommand( listener, "ReturnSingleTwoWordString" );
    commandMgr.AddCommand( listener, "ReturnSingleLongString" );
    commandMgr.AddCommand( listener, "ReturnSingleList" );
    commandMgr.AddCommand( listener, "ReturnNestedList" );
    commandMgr.AddCommand( listener, "ReturnMessage" );

    // Run the script that will test tcl return stuff.
    cerr << "Running tcl script..." << endl;
    rTcl = Tcl_EvalFile( interp, "test_TclCommandManager.tcl" );
    char* sTclResult = Tcl_GetStringResult( interp );
    if( 0 != strcmp( sTclResult, "" ) ) {
      cerr << "result not null" << endl;
      throw logic_error( interp->result );
    }
    if( rTcl != TCL_OK ) {
      cerr << "rTcl not OK" << endl;
      throw logic_error( interp->result );
    }
    
  }
  catch ( char* msg ) {
    cerr << msg << " failed" << endl;
    if( TCL_OK != rTcl ) {
      cerr << "rTcl is " << rTcl << endl;
    }
    exit( 1 );
  }
  catch( exception e ) {
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
