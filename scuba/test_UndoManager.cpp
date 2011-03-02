/**
 * @file  test_UndoManager.cpp
 * @brief test UndoManager class
 *
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:39 $
 *    $Revision: 1.9 $
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


#include <map>

#include "UndoManager.h"
extern "C" {
#include "macros.h"
}
#include "Scuba-impl.h"

const char* Progname = "test_UndoManager";

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

    for ( int n = 0; n < 10; n ++ ) {
      data[n] = kOrigValue;
    }

    int id = undoList.BeginAction( "Change values" );

    for ( int n = 0; n < 10; n ++ ) {
      if ( n%2 ) {
        data[n] = kNewValue;
        TestUndoAction* action  =
          new TestUndoAction ( data, n, kOrigValue, kNewValue );
        undoList.AddAction( id, action );
      }
    }

    undoList.EndAction( id );

    for ( int n = 0; n < 10; n ++ ) {
      if ( n%2 ) {
        Assert((kNewValue == data[n]), "pre-undo wasn't correct" );
      } else {
        Assert((kOrigValue == data[n]), "pre-undo wasn't correct" );
      }
    }

    {
      stringstream ssErr;
      ssErr << "Pre undo Size of undo list was "
      << undoList.mUndoActions.size();
      Assert( (undoList.mUndoActions.size() == 1), ssErr.str() );
    }

    {
      stringstream ssErr;
      ssErr << "Pre undo Size of redo list was "
      << undoList.mRedoActions.size();
      Assert( (undoList.mRedoActions.size() == 0), ssErr.str() );
    }

    undoList.Undo();

    {
      stringstream ssErr;
      ssErr << "Post undo Size of undo list was "
      << undoList.mUndoActions.size();
      Assert( (undoList.mUndoActions.size() == 0), ssErr.str() );
    }

    {
      stringstream ssErr;
      ssErr << "Post undo Size of redo list was "
      << undoList.mRedoActions.size();
      Assert( (undoList.mRedoActions.size() == 1), ssErr.str() );
    }

    for ( int n = 0; n < 10; n ++ ) {
      Assert((kOrigValue == data[n]), "post-undo wasn't correct" );
    }

    undoList.Redo();

    {
      stringstream ssErr;
      ssErr << "Post redo Size of undo list was "
      << undoList.mUndoActions.size();
      Assert( (undoList.mUndoActions.size() == 1), ssErr.str() );
    }

    {
      stringstream ssErr;
      ssErr << "Post redo Size of redo list was "
      << undoList.mRedoActions.size();
      Assert( (undoList.mRedoActions.size() == 0), ssErr.str() );
    }

    for ( int n = 0; n < 10; n ++ ) {
      if ( n%2 ) {
        Assert((kNewValue == data[n]), "post-redo wasn't correct" );
      } else {
        Assert((kOrigValue == data[n]), "post-redo wasn't correct" );
      }
    }

    string sTitle = undoList.GetUndoTitle();
    Assert( ("Undo Change values" == sTitle),
            "GetUndoTitle failed for undo" );

    undoList.Undo();

    sTitle = undoList.GetRedoTitle();
    ;
    Assert( ("Redo Change values" == sTitle),
            "GetRedoTitle failed for redo" );


    // Try the tcl commands.
    char sCommand[1024];
    int rTcl;

    undoList.Redo();

    sprintf( sCommand, "GetUndoTitle" );
    rTcl = Tcl_Eval( iInterp, sCommand );
    AssertTclOK( rTcl );
    const char* sTclResult = Tcl_GetStringResult( iInterp );
    sTitle = sTclResult;
    {
      stringstream ssErr;
      ssErr << "Tcl GetUndoTitle failed for redo, was " << sTitle << endl;
      Assert( ("Undo Change values" == sTitle), ssErr.str() );
    }

    undoList.Undo();

    sprintf( sCommand, "GetRedoTitle" );
    rTcl = Tcl_Eval( iInterp, sCommand );
    AssertTclOK( rTcl );
    sTclResult = Tcl_GetStringResult( iInterp );
    sTitle = sTclResult;
    {
      stringstream ssErr;
      ssErr << "Tcl GetUndoTitle failed for undo, was " << sTitle << endl;
      Assert( ("Redo Change values" == sTitle), ssErr.str() );
    }


    // Try a second list but don't add any undo actions to it. Then
    // undo, and make sure the list is the same.
    id = undoList.BeginAction( "Dummy" );

    for ( int n = 0; n < 10; n ++ ) {
      if ( n%2 ) {
        data[n] = kNewValue2;
      } else {
        data[n] = kNewValue;
      }
    }

    undoList.EndAction( id );

    undoList.Undo();

    for ( int n = 0; n < 10; n ++ ) {
      if ( n%2 ) {
        Assert((kNewValue2 == data[n]), "post dummy undo wasn't correct" );
      } else {
        Assert((kNewValue == data[n]), "post dummy undo wasn't correct" );
      }
    }

    undoList.Redo();

    for ( int n = 0; n < 10; n ++ ) {
      if ( n%2 ) {
        Assert((kNewValue2 == data[n]), "post dummy redo wasn't correct" );
      } else {
        Assert((kNewValue == data[n]), "post dummy redo wasn't correct" );
      }
    }



    data.clear();
    // Try multiple concurrent IDs. Make 10 ten lists, one for each
    // change, and build them all at the same time, then undo 10
    // times.
    int ids[10];
    for( int n = 0; n < 10; n++ ) {
      data[n] = kNewValue;
      ids[n] = undoList.BeginAction( "Change" );
      TestUndoAction* action  =
	new TestUndoAction ( data, n, kOrigValue, kNewValue );
      undoList.AddAction( ids[n], action );
    }
    for( int n = 0; n < 10; n++ ) {
      undoList.EndAction( ids[n] );
    }
    for( int n = 0; n < 10; n++ ) {
      undoList.Undo();
    }
    for ( int n = 0; n < 10; n ++ ) {
      Assert((kOrigValue == data[n]), "post-undo after multiple undo lists wasn't correct" );
    }


    // Do the same undo test we did before but through tcl.
    id = undoList.BeginAction( "Change values" );

    for ( int n = 0; n < 10; n ++ ) {
      if ( n%2 ) {
        data[n] = kNewValue;
        TestUndoAction* action  =
          new TestUndoAction ( data, n, kOrigValue, kNewValue );
        undoList.AddAction( id, action );
      } else {
        data[n] = kOrigValue;
      }
    }

    undoList.EndAction( id );

    rTcl = Tcl_Eval( iInterp, "Undo" );
    AssertTclOK( rTcl );

    for ( int n = 0; n < 10; n ++ ) {
      Assert((kOrigValue == data[n]), "post-tcl undo wasn't correct" );
    }

    rTcl = Tcl_Eval( iInterp, "Redo" );
    AssertTclOK( rTcl );

    for ( int n = 0; n < 10; n ++ ) {
      if ( n%2 ) {
        Assert((kNewValue == data[n]), "post-tcl redo wasn't correct" );
      } else {
        Assert((kOrigValue == data[n]), "post-tcl redo wasn't correct" );
      }
    }


    // Clear list and check sizes.
    undoList.Clear();

    {
      stringstream ssErr;
      ssErr << "After clearing, undo size was "
      << undoList.mUndoActions.size();
      Assert( (undoList.mUndoActions.size() == 0), ssErr.str() );
    }

    {
      stringstream ssErr;
      ssErr << "After clearing, redo size was "
      << undoList.mRedoActions.size();
      Assert( (undoList.mRedoActions.size() == 0), ssErr.str() );
    }


    // Test multiple undoes. Create an undo action for each value from
    // 1 to 10, then undo them one by one, each time checking the
    // values in the array. Then redo them, checking again.
    for ( unsigned int n = 0; n < 10; n++ ) {
      data[n] = kNewValue;
      stringstream ssTitle;
      ssTitle  << "Change value " << n;
      id = undoList.BeginAction( ssTitle.str() );
      TestUndoAction* action =
        new TestUndoAction( data, n, kOrigValue, kNewValue );
      undoList.AddAction( id, action );
      undoList.EndAction( id );

      {
        stringstream ssErr;
        ssErr << "Pre multi undo size of undo list with n " << n << " was "
        << undoList.mUndoActions.size();
        Assert( (undoList.mUndoActions.size() == n+1), ssErr.str() );
      }

      {
        stringstream ssErr;
        ssErr << "Pre multi redo size of undo list with n " << n << " was "
        << undoList.mRedoActions.size();
        Assert( (undoList.mRedoActions.size() == 0), ssErr.str() );
      }
    }

    for ( int n = 9; n >= 5; n-- ) {

      undoList.Undo();

      {
        stringstream ssErr;
        ssErr << "Data not correct n " << n << ", array: ";
        for ( int m = 0; m < 10; m++ ) {
          ssErr << data[m] << " ";
        }

        for ( int m = 0; m < n; m++ ) {
          Assert( (data[m] == kNewValue), ssErr.str() );
        }
        for ( int m = n; m < 10; m++ ) {
          Assert( (data[m] == kOrigValue), ssErr.str() );
        }
      }

      {
        stringstream ssErr;
        ssErr << "multi undo size of undo list with n " << n << " was "
        << undoList.mUndoActions.size();
        Assert( (undoList.mUndoActions.size() == (unsigned int)n),
                ssErr.str() );
      }

      {
        stringstream ssErr;
        ssErr << "multi redo size of undo list with n " << n << " was "
        << undoList.mRedoActions.size();
        Assert( (undoList.mRedoActions.size() == (unsigned int)(10-n)),
                ssErr.str() );
      }
    }

    for ( int n = 5; n < 10; n++ ) {

      undoList.Redo();


      // Now 0 -> n are new, and n+1 -> 9 are orig
      {
        stringstream ssErr;
        ssErr << "Data not correct n " << n << ", array: ";
        for ( int m = 0; m < 10; m++ ) {
          ssErr << data[m] << " ";
        }

        for ( int m = 0; m <= n; m++ ) {
          Assert( (data[m] == kNewValue), ssErr.str() );
        }
        for ( int m = n+1; m < 10; m++ ) {
          Assert( (data[m] == kOrigValue), ssErr.str() );
        }
      }

      {
        stringstream ssErr;
        ssErr << "multi undo size of undo list with n " << n << " was "
        << undoList.mUndoActions.size();
        Assert( (undoList.mUndoActions.size() == (unsigned int)n+1),
                ssErr.str() );
      }

      {
        stringstream ssErr;
        ssErr << "multi redo size of undo list with n " << n << " was "
        << undoList.mRedoActions.size();
        Assert( (undoList.mRedoActions.size() == (unsigned int)(9-n)),
                ssErr.str() );
      }
    }


    // Test the limit. Get the max, make more than that, and make sure
    // the count is still good. (In the future it would be good to
    // test if the right one is being deleted.)
    undoList.Clear();
    for ( int n = 0; n < undoList.mcMaxActions + 2; n++ ) {
      id = undoList.BeginAction( "test" );
      undoList.EndAction( id );

      Assert( (undoList.mUndoActions.size() + undoList.mRedoActions.size() <=
               (unsigned int)undoList.mcMaxActions),
              "Went above max number of actions." );
    }
    for ( int n = 0; n < 10; n++ ) {
      undoList.Undo();
      Assert( (undoList.mUndoActions.size() + undoList.mRedoActions.size() <=
               (unsigned int)undoList.mcMaxActions),
              "Went above max number of actions." );
    }
    for ( int n = 0; n < 10; n++ ) {
      id = undoList.BeginAction( "test" );
      undoList.EndAction( id );

      Assert( (undoList.mUndoActions.size() + undoList.mRedoActions.size() <=
               (unsigned int)undoList.mcMaxActions),
              "Went above max number of actions." );
    }


  } catch ( runtime_error& e ) {
    cerr << "failed with exception: " << e.what() << endl;
    exit( 1 );
  } catch (...) {
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


  } catch ( runtime_error& e ) {
    cerr << "failed with exception: " << e.what() << endl;
    exit( 1 );
  } catch (...) {
    cerr << "failed" << endl;
    exit( 1 );
  }

  cerr << "Success" << endl;

  exit( 0 );
}

