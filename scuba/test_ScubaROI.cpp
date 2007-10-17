/**
 * @file  test_ScubaROI.cpp
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: kteich $
 *    $Date: 2007/10/17 23:59:49 $
 *    $Revision: 1.6 $
 *
 * Copyright (C) 2002-2007,
 * The General Hospital Corporation (Boston, MA). 
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
 *
 */


#include "ScubaROI.h"
extern "C" {
#include "macros.h"
}
#include "Scuba-impl.h"

char* Progname = "test_ScubaROI";

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


class ScubaROITester {
public:
void Test( Tcl_Interp* iInterp );
};

void
ScubaROITester::Test ( Tcl_Interp* iInterp ) {

  stringstream ssError;

  try {

    ScubaROI roi;

    string sLabel, sTestLabel;
    sLabel = "punk ROI";
    roi.SetLabel( sLabel );
    Assert( (sLabel == roi.msLabel), "SetLabel didn't work" );

    sTestLabel = roi.GetLabel();
    Assert( (sTestLabel == sLabel), "GetLabel didn't work" );

    roi.SetType( ScubaROI::Structure );
    Assert( (ScubaROI::Structure == roi.GetType()),
            "Get/SetType structure didn't work" );

    roi.SetType( ScubaROI::Free );
    Assert( (ScubaROI::Free == roi.GetType()),
            "Get/SetType free didn't work" );

    roi.SetStructure( 5 );
    Assert( (5 == roi.GetStructure()),
            "Get/SetStructure didn't work" );

    roi.SetStructure( -1 );
    Assert( (-1 == roi.GetStructure()),
            "Get/SetStructure didn't work" );

    int color[3];
    color[0] = 1;
    color[1] = 2;
    color[2] = 3;
    roi.SetFreeColor( color );
    Assert( (1 == roi.mFreeColor[0] && 2 == roi.mFreeColor[1] &&
             3 == roi.mFreeColor[2]), "SetFreeColor failed" );

    // Try the tcl commands.
    char sCommand[1024];
    int rTcl;

    int roiID = roi.GetID();

    sprintf( sCommand, "SetROILabel %d \"i hate ROIs\"", roiID );
    rTcl = Tcl_Eval( iInterp, sCommand );
    AssertTclOK( rTcl );
    Assert( ("i hate ROIs" == roi.GetLabel()), "Tcl SetROILabel failed" );

    sprintf( sCommand, "GetROILabel %d", roiID );
    rTcl = Tcl_Eval( iInterp, sCommand );
    AssertTclOK( rTcl );
    const char* sTclResult = Tcl_GetStringResult( iInterp );
    sTestLabel = sTclResult;
    Assert( ("i hate ROIs" == sTestLabel), "Tcl GetROILabel failed" );


    sprintf( sCommand, "SetROIType %d structure", roiID );
    rTcl = Tcl_Eval( iInterp, sCommand );
    AssertTclOK( rTcl );
    Assert( (ScubaROI::Structure == roi.mType),
            "Tcl SetROIType failed" );
    sprintf( sCommand, "SetROIType %d free", roiID );
    rTcl = Tcl_Eval( iInterp, sCommand );
    AssertTclOK( rTcl );
    Assert( (ScubaROI::Free == roi.mType),
            "Tcl SetROIType failed" );

    sprintf( sCommand, "GetROIType %d", roiID );
    rTcl = Tcl_Eval( iInterp, sCommand );
    AssertTclOK( rTcl );
    sTclResult = Tcl_GetStringResult( iInterp );
    string sType = sTclResult;
    ssError << "Tcl GetROIType failed: sType=" << sType << ", should be free";
    Assert( ("free" == sType), ssError.str() );

    roi.SetType( ScubaROI::Structure );
    sprintf( sCommand, "GetROIType %d", roiID );
    rTcl = Tcl_Eval( iInterp, sCommand );
    AssertTclOK( rTcl );
    sTclResult = Tcl_GetStringResult( iInterp );
    sType = sTclResult;
    ssError << "Tcl GetROIType failed: sType=" << sType << ", should be structure";
    Assert( ("structure" == sType), "Tcl GetROIType failed" );


    sprintf( sCommand, "SetROIStructure %d 5", roiID );
    rTcl = Tcl_Eval( iInterp, sCommand );
    AssertTclOK( rTcl );
    Assert( (5 == roi.mStructure), "Tcl SetROIStructure failed" );

    roi.SetStructure( 6 );
    sprintf( sCommand, "GetROIStructure %d", roiID );
    rTcl = Tcl_Eval( iInterp, sCommand );
    AssertTclOK( rTcl );
    sTclResult = Tcl_GetStringResult( iInterp );
    int structure = strtol(sTclResult, (char**)NULL, 10);
    Assert( (ERANGE != errno), "MakeNewColorLUT did not return valid ID" );
    Assert( (6 == structure), "Tcl GetROIStructure failed" );

    sprintf( sCommand, "SetROIColor %d 7 8 9", roiID );
    rTcl = Tcl_Eval( iInterp, sCommand );
    AssertTclOK( rTcl );
    Assert( (7 == roi.mFreeColor[0] && 8 == roi.mFreeColor[1] && 9 == roi.mFreeColor[2]),
            "Tcl SetROIColor failed" );

    color[0] = 10;
    color[1] = 11;
    color[2] = 12;
    roi.SetFreeColor( color );
    sprintf( sCommand, "GetROIColor %d", roiID );
    rTcl = Tcl_Eval( iInterp, sCommand );
    AssertTclOK( rTcl );
    sTclResult = Tcl_GetStringResult( iInterp );
    string sColor = sTclResult;
    Assert( ("10 11 12" == sColor), "Tcl GetROIColor failed" );

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


    ScubaROITester tester0;
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

