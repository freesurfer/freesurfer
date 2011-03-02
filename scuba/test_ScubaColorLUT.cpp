/**
 * @file  test_ScubaColorLUT.cpp
 * @brief test ScubaColorLUT class
 *
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:39 $
 *    $Revision: 1.7 $
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


#include "ScubaColorLUT.h"
extern "C" {
#include "macros.h"
}
#include "Scuba-impl.h"

const char* Progname = "test_ScubaColorLUT";

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


class ScubaColorLUTTester {
public:
void Test( Tcl_Interp* iInterp );
};

void
ScubaColorLUTTester::Test ( Tcl_Interp* iInterp ) {

  stringstream ssError;

  try {

    ScubaColorLUT lut;

    lut.UseFile( "testLUT.lut" );

    int color[3];
    string sLabel;
    lut.GetColorAtIndex( 0, color );
    sLabel = lut.GetLabelAtIndex( 0 );
    ssError << "Failed for entry 0: color " << color[0] << ", "
    << color[1] << ", " << color[2] << ", label " << sLabel;
    Assert( (0 == color[0] && 0 == color[1] && 0 == color[2] &&
             sLabel == "Unknown"), ssError.str() );

    lut.GetColorAtIndex( 1, color );
    sLabel = lut.GetLabelAtIndex( 1 );
    Assert( (255 == color[0] && 255 == color[1] && 255 == color[2] &&
             sLabel == "Entry1"),
            "Failed for entry 1" );

    lut.GetColorAtIndex( 2, color );
    sLabel = lut.GetLabelAtIndex( 2 );
    Assert( (1 == color[0] && 2 == color[1] && 3 == color[2] &&
             sLabel == "Entry2"),
            "Failed for entry 2" );

    lut.GetColorAtIndex( 4, color );
    sLabel = lut.GetLabelAtIndex( 4 );
    Assert( (1 == color[0] && 2 == color[1] && 3 == color[2] &&
             sLabel == "Entry4"),
            "Failed for entry 4" );

    // Try the tcl commands.
    char sCommand[1024];
    int rTcl;

    sprintf( sCommand, "MakeNewColorLUT" );
    rTcl = Tcl_Eval( iInterp, sCommand );
    AssertTclOK( rTcl );
    const char* sTclResult = Tcl_GetStringResult( iInterp );
    int lutID = strtol(sTclResult, (char**)NULL, 10);
    Assert( (ERANGE != errno), "MakeNewColorLUT did not return valid ID" );

    ScubaColorLUT& lut1 = ScubaColorLUT::FindByID( lutID );

    sprintf( sCommand, "SetColorLUTFileName %d testLUT.lut", lut1.GetID() );
    rTcl = Tcl_Eval( iInterp, sCommand );
    AssertTclOK( rTcl );

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


    ScubaColorLUTTester tester0;
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

