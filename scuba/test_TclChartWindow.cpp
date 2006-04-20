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

char* Progname = "test_TclChartWindow";

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
  for( int nData = 0; nData < cData; nData++ ) {

    data.mX = nData;
    data.mY = random() % 100;

    stringstream ssLabel;
    ssLabel << "Item " << nData;
    data.msLabel = ssLabel.str();

    lData.push_back( data );
  }

  mChart->SetPointData( lData );
  
}

void
TclChartWindowTester::Test( Tcl_Interp* iInterp ) {

  mChart = ChartWindow::NewChartWindow();

  TclCommandManager& commandMgr = TclCommandManager::GetManager();
  commandMgr.AddCommand( *this, "TestResult", 2, "testID pass", "" );
  commandMgr.AddCommand( *this, "NewData", 0, "", "" );
  commandMgr.AddCommand( *this, "ToggleShowLegend", 0, "", "" );

  stringstream ssCommand;
  ssCommand << "frame .f";
  commandMgr.SendCommand( ssCommand.str() );

  ssCommand.str("");
  ssCommand << "button .f.bwNewData -text \"New Data\""
	    << " -command NewData";
  commandMgr.SendCommand( ssCommand.str() );

  ssCommand.str("");
  ssCommand << "button .f.bwToggleLegend -text \"Toggle Legend\""
	    << " -command {ToggleShowLegend}";
  commandMgr.SendCommand( ssCommand.str() );

  ssCommand.str("");
  ssCommand << "pack .f.bwNewData .f.bwToggleLegend -fill both -expand yes";
  commandMgr.SendCommand( ssCommand.str() );

  ssCommand.str("");
  ssCommand << "pack .f -fill both -expand yes";
  commandMgr.SendCommand( ssCommand.str() );

  ssCommand.str("");

  NewData();
  mChart->Draw();

  MakeTestButtonWindow( 1 );
  mNewPass = false;
  mTestPass[1] = false;
  while( 1 && !mNewPass ) {
    Tcl_DoOneEvent( TCL_ALL_EVENTS );
  }

  if( !mTestPass[1] ) {
    throw runtime_error( "Didn't pass test 1" );
  }
}

TclCommandListener::TclCommandResult
TclChartWindowTester::DoListenToTclCommand ( char* isCommand, 
					      int iArgc, char** iasArgv ) { 

  if( 0 == strcmp( isCommand, "TestResult" ) ) {

    int testID;
    bool pass;
    testID = TclCommandManager::ConvertArgumentToInt( iasArgv[1] );
    pass = TclCommandManager::ConvertArgumentToBoolean( iasArgv[2] );
    
    mTestPass[testID] = pass;
    mNewPass = true;

    return ok;

  } else if( 0 == strcmp( isCommand, "NewData" ) ) {
    
    NewData();
    mChart->Draw();
    
    return ok;

  } else if( 0 == strcmp( isCommand, "ToggleShowLegend" ) ) {
    
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

 
  }
  catch( runtime_error& e ) {
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
