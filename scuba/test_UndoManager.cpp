#include <map>

#include "UndoManager.h"
extern "C" {
#include "macros.h"
}
#include "Scuba-impl.h"

char* Progname = "test_UndoManager";

using namespace std;

#define Assert(x,s)   \
  if(!(x)) { \
  stringstream ss; \
  ss << "Line " << __LINE__ << ": " << s; \
  cerr << ss.str().c_str() << endl; \
  throw runtime_error( ss.str() ); \
  }

#define AssertTclOK(x) \
    if( TCL_OK != (x) ) { \
      ssError << "Tcl_Eval returned not TCL_OK: " << endl  \
	     << "Command: " << sCommand << endl \
	     << "Result: " << iInterp->result; \
      throw runtime_error( ssError.str() ); \
    } \


class UndoManagerTester {
public:
  void Test( Tcl_Interp* iInterp );
};


class TestUndoAction : public UndoAction {
public:
  TestUndoAction ( map<int,int>& iData,
		  int inIndex, int iUndoValue, int iRedoValue ) :
  mData ( iData ) {
    mData = iData;
    mnIndex = inIndex;
    mUndoValue = iUndoValue;
    mRedoValue = iRedoValue;
  }

  virtual void Undo () { 
    mData[mnIndex] = mUndoValue; 
  }
  virtual void Redo () {
    mData[mnIndex] = mRedoValue; 
  }

protected:
  map<int,int>& mData;
  int mnIndex;
  int mUndoValue;
  int mRedoValue;

};

void 
UndoManagerTester::Test ( Tcl_Interp* iInterp ) {

  stringstream ssError;

  try {

    UndoManager& undoList = UndoManager::GetManager();

    map<int,int> data;

    int const kOrigValue = 1;
    int const kNewValue = 2;
    int const kNewValue2 = 3;

    for( int n = 0; n < 10; n ++ ) {
      data[n] = kOrigValue;
    }

    undoList.BeginAction( "Change values" );

    for( int n = 0; n < 10; n ++ ) {
      if( n%2 ) {
	data[n] = kNewValue;
	TestUndoAction* action  = 
	  new TestUndoAction ( data, n, kOrigValue, kNewValue );
	undoList.AddAction( action );
      }
    }

    undoList.EndAction();

    for( int n = 0; n < 10; n ++ ) {
      if( n%2 ) {
	Assert((kNewValue == data[n]), "pre-undo wasn't correct" );
      } else {
	Assert((kOrigValue == data[n]), "pre-undo wasn't correct" );
      }
    }

    undoList.Undo();

    for( int n = 0; n < 10; n ++ ) {
      Assert((kOrigValue == data[n]), "post-undo wasn't correct" );
    }

    undoList.Redo();

    for( int n = 0; n < 10; n ++ ) {
      if( n%2 ) {
	Assert((kNewValue == data[n]), "post-redo wasn't correct" );
      } else {
	Assert((kOrigValue == data[n]), "post-redo wasn't correct" );
      }
    }

    string sTitle = undoList.GetTitle();
    Assert( ("Undo Change values" == sTitle), 
	    "GetTitle failed for undo" );

    undoList.Undo();

    sTitle = undoList.GetTitle();;
    Assert( ("Redo Change values" == sTitle), 
	    "GetTitle failed for redo" );


    // Try the tcl commands.
    char sCommand[1024];
    int rTcl;
    
    undoList.Redo();

    sprintf( sCommand, "GetUndoTitle" );
    rTcl = Tcl_Eval( iInterp, sCommand );
    AssertTclOK( rTcl );
    char* sTclResult = Tcl_GetStringResult( iInterp );
    sTitle = sTclResult;
    ssError << "Tcl GetUndoTitle failed for redo, was " << sTitle;
    Assert( ("Undo Change values" == sTitle), ssError.str() );

    undoList.Undo();
    
    rTcl = Tcl_Eval( iInterp, sCommand );
    AssertTclOK( rTcl );
    sTclResult = Tcl_GetStringResult( iInterp );
    sTitle = sTclResult;
    ssError << "Tcl GetUndoTitle failed for undo, was " << sTitle;
    Assert( ("Redo Change values" == sTitle), ssError.str() );


    // Try a second list but don't add any undo actions to it. Then
    // undo, and make sure the list is the same.
    undoList.BeginAction( "Dummy" );

    for( int n = 0; n < 10; n ++ ) {
      if( n%2 ) {
	data[n] = kNewValue2;
      } else {
	data[n] = kNewValue;
      }
    }

    undoList.EndAction();

    undoList.Undo();

    for( int n = 0; n < 10; n ++ ) {
      if( n%2 ) {
	Assert((kNewValue2 == data[n]), "post dummy undo wasn't correct" );
      } else {
	Assert((kNewValue == data[n]), "post dummy undo wasn't correct" );
      }
    }

    undoList.Redo();

    for( int n = 0; n < 10; n ++ ) {
      if( n%2 ) {
	Assert((kNewValue2 == data[n]), "post dummy redo wasn't correct" );
      } else {
	Assert((kNewValue == data[n]), "post dummy redo wasn't correct" );
      }
    }



    // Do the same undo test we did before but through tcl.
    undoList.BeginAction( "Change values" );

    for( int n = 0; n < 10; n ++ ) {
      if( n%2 ) {
	data[n] = kNewValue;
	TestUndoAction* action  = 
	  new TestUndoAction ( data, n, kOrigValue, kNewValue );
	undoList.AddAction( action );
      } else {
	data[n] = kOrigValue;
      }
    }

    undoList.EndAction();

    rTcl = Tcl_Eval( iInterp, "UndoOrRedo" );
    AssertTclOK( rTcl );

    for( int n = 0; n < 10; n ++ ) {
      Assert((kOrigValue == data[n]), "post-tcl undo wasn't correct" );
    }

    rTcl = Tcl_Eval( iInterp, "UndoOrRedo" );
    AssertTclOK( rTcl );

    for( int n = 0; n < 10; n ++ ) {
      if( n%2 ) {
	Assert((kNewValue == data[n]), "post-tcl redo wasn't correct" );
      } else {
	Assert((kOrigValue == data[n]), "post-tcl redo wasn't correct" );
      }
    }

  }
  catch( runtime_error e ) {
    cerr << "failed with exception: " << e.what() << endl;
    exit( 1 );
  }
  catch(...) {
    cerr << "failed" << endl;
    exit( 1 );
  }
}


int main ( int argc, char** argv ) {

  cerr << "Beginning test" << endl;

  try {

    Tcl_Interp* interp = Tcl_CreateInterp();
    Assert( interp, "Tcl_CreateInterp returned null" );
  
    int rTcl = Tcl_Init( interp );
    Assert( TCL_OK == rTcl, "Tcl_Init returned not TCL_OK" );
    
    TclCommandManager& commandMgr = TclCommandManager::GetManager();
    commandMgr.SetOutputStreamToCerr();
    commandMgr.Start( interp );


    UndoManagerTester tester0;
    tester0.Test( interp );

 
  }
  catch( runtime_error e ) {
    cerr << "failed with exception: " << e.what() << endl;
    exit( 1 );
  }
  catch(...) {
    cerr << "failed" << endl;
    exit( 1 );
  }

  cerr << "Success" << endl;

  exit( 0 );
}

