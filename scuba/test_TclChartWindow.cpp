/**
 * @file  test_TclChartWindow.cpp
 * @brief test TclChartWindow class
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


#include <stdlib.h>
#include "string_fixed.h"
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <map>
#include "TclScubaKeyCombo.h"
extern "C" {
#include "mri.h"
}
#include "Scuba-impl.h"
#include "tk.h"
#include "tix.h"
#include "blt.h"
#include "TclChartWindow.h"

#define Assert(x,s)   \
  if(!(x)) { \
  stringstream ss; \
  ss << "Line " << __LINE__ << ": " << s; \
  cerr << ss.str().c_str() << endl; \
  throw runtime_error( ss.str() ); \
  }


#define AssertTclOK(x) \
    if( TCL_OK != (x) ) { \
      stringstream ssError; \
      ssError << "Tcl_Eval returned not TCL_OK: " << endl  \
      << "Command: " << sCommand << endl \
      << "Result: " << iInterp->result; \
      throw runtime_error( ssError.str() ); \
    } \


using namespace std;

const char* Progname = "test_TclChartWindow";

// Tester class.
class TclChartWindowTester : public TclCommandListener {
public:

  // Main tester entry function.
  void Test( Tcl_Interp* iInterp );

  // Make a window with a pair of Pass/Fail buttons on it. The buttons
  // will call TestResult with the test ID and the result.
  void MakeTestButtonWindow ( int iTestID );

  // Generate new random data for the chart.
  void NewData ();
  void NewDataMulti ();

  // Tcl command callback function.
  TclCommandResult DoListenToTclCommand ( char* isCommand,
                                          int iArgc, char** iArgv );

protected:
  ChartWindow* mChart;
  map<int,bool> mTestPass;
  bool mNewPass;
};

void
TclChartWindowTester::MakeTestButtonWindow ( int iTestID ) {

  TclCommandManager& commandMgr = TclCommandManager::GetManager();

  stringstream ssCommand;
  ssCommand << "toplevel .w";
  commandMgr.SendCommand( ssCommand.str() );

  ssCommand.str("");
  ssCommand << "button .w.bwPass -text PASS -background green"
  << " -activebackground green"
  << " -command {TestResult " << iTestID
  << " true ; wm withdraw .w}";
  commandMgr.SendCommand( ssCommand.str() );

  ssCommand.str("");
  ssCommand << "button .w.bwFail -text FAIL -background red"
  << " -activebackground red"
  << " -command {TestResult " << iTestID
  << " false ; wm withdraw .w}";
  commandMgr.SendCommand( ssCommand.str() );

  ssCommand.str("");
  ssCommand << "pack .w.bwPass .w.bwFail -fill both -expand yes";
  commandMgr.SendCommand( ssCommand.str() );
}

void
TclChartWindowTester::NewData () {

  mChart->ClearData();

  list<ChartWindow::PointData> lData;

  ChartWindow::PointData data;
  int cData = 100;
  for ( int nData = 0; nData < cData; nData++ ) {

    data.mX = nData;
    data.mY = random() % 100;

    stringstream ssLabel;
    ssLabel << "Item " << nData;
    data.msLabel = ssLabel.str();

    lData.push_back( data );
  }

  mChart->SetPointData( lData );

  list<ChartWindow::MarkerData> lMarkers;
  ChartWindow::MarkerData marker;
  int cMarkers = 5;
  for ( int nMarker = 0; nMarker < cMarkers; nMarker++ ) {

    marker.mValue = random() % 100;

    marker.mColorRGBi[0] = random() % 256;
    marker.mColorRGBi[1] = random() % 256;
    marker.mColorRGBi[2] = random() % 256;

    stringstream ssLabel;
    ssLabel << "Marker " << nMarker;
    marker.msLabel = ssLabel.str();

    lMarkers.push_back( marker );
  }

  mChart->SetXAxisMarkers( lMarkers );

}

void
TclChartWindowTester::NewDataMulti () {

  mChart->ClearData();

  list<ChartWindow::PointData> lData;

  ChartWindow::PointData data;
  int cGroups = 10;
  for ( int nGroup = 0; nGroup < cGroups; nGroup++ ) {

    lData.clear();

    int cData = 20;
    for ( int nData = 0; nData < cData; nData++ ) {

      data.mX = nData;
      data.mY = random() % 100;

      stringstream ssLabel;
      ssLabel << "Item " << nData;
      data.msLabel = ssLabel.str();

      lData.push_back( data );
    }

    mChart->SetPointData( nGroup, lData );

    stringstream ssLabel;
    ssLabel << "Group " << nGroup;
    mChart->SetGroupLabel( nGroup, ssLabel.str() );
    mChart->SetGroupConnected( nGroup, (nGroup % 2) );
  }

}

void
TclChartWindowTester::Test( Tcl_Interp* iInterp ) {

  mChart = ChartWindow::NewChartWindow();

  TclCommandManager& commandMgr = TclCommandManager::GetManager();
  commandMgr.AddCommand( *this, "TestResult", 2, "testID pass", "" );
  commandMgr.AddCommand( *this, "NewData", 0, "", "" );
  commandMgr.AddCommand( *this, "NewDataMulti", 0, "", "" );
  commandMgr.AddCommand( *this, "ToggleShowLegend", 0, "", "" );

  stringstream ssCommand;
  ssCommand << "frame .f";
  commandMgr.SendCommand( ssCommand.str() );

  ssCommand.str("");
  ssCommand << "button .f.bwNewData -text \"New 1 Data\""
  << " -command NewData";
  commandMgr.SendCommand( ssCommand.str() );

  ssCommand.str("");
  ssCommand << "button .f.bwNewDataMulti -text \"New multi Data\""
  << " -command NewDataMulti";
  commandMgr.SendCommand( ssCommand.str() );

  ssCommand.str("");
  ssCommand << "button .f.bwToggleLegend -text \"Toggle Legend\""
  << " -command {ToggleShowLegend}";
  commandMgr.SendCommand( ssCommand.str() );

  ssCommand.str("");
  ssCommand << "pack .f.bwNewData .f.bwNewDataMulti .f.bwToggleLegend -fill both -expand yes";
  commandMgr.SendCommand( ssCommand.str() );

  ssCommand.str("");
  ssCommand << "pack .f -fill both -expand yes";
  commandMgr.SendCommand( ssCommand.str() );

  mChart->SetTitle( "Window title" );
  mChart->SetXAxisLabel( "x Axis" );
  mChart->SetYAxisLabel( "y Axis" );
  mChart->SetInfo( "Info string" );

  NewData();
  mChart->Draw();

  MakeTestButtonWindow( 1 );
  mNewPass = false;
  mTestPass[1] = false;
  while ( 1 && !mNewPass ) {
    Tcl_DoOneEvent( TCL_ALL_EVENTS );
  }

  if ( !mTestPass[1] ) {
    throw runtime_error( "Didn't pass test 1" );
  }
}

TclCommandListener::TclCommandResult
TclChartWindowTester::DoListenToTclCommand ( char* isCommand,
    int iArgc, char** iasArgv ) {

  if ( 0 == strcmp( isCommand, "TestResult" ) ) {

    int testID;
    bool pass;
    testID = TclCommandManager::ConvertArgumentToInt( iasArgv[1] );
    pass = TclCommandManager::ConvertArgumentToBoolean( iasArgv[2] );

    mTestPass[testID] = pass;
    mNewPass = true;

    return ok;

  } else if ( 0 == strcmp( isCommand, "NewData" ) ) {

    static int cNew = 0;
    stringstream ssLabel;
    ssLabel << "Info String " << cNew++;

    NewData();
    mChart->SetInfo( ssLabel.str() );
    mChart->Draw();

    return ok;

  } else if ( 0 == strcmp( isCommand, "NewDataMulti" ) ) {

    static int cNew = 0;
    stringstream ssLabel;
    ssLabel << "Info String " << cNew++;

    NewDataMulti();
    mChart->SetInfo( ssLabel.str() );
    mChart->Draw();

    return ok;

  } else if ( 0 == strcmp( isCommand, "ToggleShowLegend" ) ) {

    mChart->SetShowLegend( !mChart->GetShowLegend() );
    mChart->Draw();

    return ok;
  }



  return ok;
}


int main ( int argc, char** argv ) {

  cerr << "Beginning test" << endl;

  try {

    Tcl_Interp* interp = Tcl_CreateInterp();
    Assert( interp, "Tcl_CreateInterp returned null" );

    int rTcl = Tcl_Init( interp );
    Assert( TCL_OK == rTcl, "Tcl_Init returned not TCL_OK" );

    int rTk = Tcl_Init( interp );
    Assert( TCL_OK == rTk, "Tk_Init returned not TCL_OK" );

    int rTix = Tix_Init( interp );
    Assert( TCL_OK == rTix, "Tix_Init returned not TCL_OK" );

    ChartWindow::SetFactory( new TclChartWindowFactory );

    TclCommandManager& commandMgr = TclCommandManager::GetManager();
    commandMgr.SetOutputStreamToCerr();
    commandMgr.Start( interp );

    commandMgr.SendCommand( "package require Tix" );
    commandMgr.SendCommand( "package require BLT" );

    TclChartWindowTester tester0;
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
