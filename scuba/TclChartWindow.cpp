/**
 * @file  TclChartWindow.cpp
 * @brief Implements a Tcl based ChartWindow
 *
 * A Tcl implementation of the ChartWindow interface. Uses BLT code in
 * TclChartWindow.tcl to implement the window itself. The Draw()
 * function uses the current ChartWindow's setting to set up the BLT
 * window.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:39 $
 *    $Revision: 1.14 $
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


#include "string_fixed.h"
#include <sstream>
#include "TclChartWindow.h"
#include "TclCommandManager.h"
#ifdef __cplusplus
extern "C"
{
#endif
#include <string.h> // strcmp
#ifdef __cplusplus
}
#endif

using namespace std;

TclChartWindowStaticListener TclChartWindow::mStaticListener;
bool TclChartWindow::sbInitedTclFile = false;

TclChartWindow::TclChartWindow () :
    ChartWindow() {
  msTitle = "Chart";

  TclCommandManager& commandMgr = TclCommandManager::GetManager();

  // If we're not already inited...
  if ( !sbInitedTclFile ) {
    try {
      // Try using a utility function in scuba.tcl to load the
      // file. If this doesn't work...
      commandMgr.SendCommand( "LoadScubaSupportFile TclChartWindow.tcl" );
    } catch (...) {
      // Just use the source command. If this doesn't work, we can't
      // find the TclChartWindow.tcl and we have a genuine error, so
      // let it be thrown.
      commandMgr.SendCommand( "source TclChartWindow.tcl" );
    }

    // This inits the chart stuff.
    commandMgr.SendCommand( "Chart_Init" );

    sbInitedTclFile = true;
  }

  // Add our commands.
  commandMgr.AddCommand( *this, "GenerateChartReport", 6, "chartID "
                         "includeGroupLabel includePointLabel includeX "
                         "includeY fileName",
                         "Generates a text file with the requested "
                         "information from the specified chart." );

  // Our chart ID is the same as our IDTracker ID.
  int chartID = GetID();

  // Create a window.
  stringstream ssCommand;
  ssCommand << "Chart_NewWindow " << chartID;
  commandMgr.SendCommand( ssCommand.str() );


}

TclChartWindow::~TclChartWindow () {}

void
TclChartWindow::Close () {

  TclCommandManager& manager = TclCommandManager::GetManager();

  int chartID = GetID();

  stringstream ssCommand;
  ssCommand << "Chart_CloseWindow " << chartID;
  manager.SendCommand( ssCommand.str() );
}

void
TclChartWindow::Draw () {

  TclCommandManager& manager = TclCommandManager::GetManager();

  // Our chart ID is the same as our IDTracker ID.
  int chartID = GetID();

  stringstream ssCommand;

  // Call the Chart_ commands from TclChartWindow.tcl to set up our
  // chart.
  ssCommand.str("");
  ssCommand << "Chart_SetWindowTitle " << chartID
  << " \"" << msTitle << "\"";
  manager.SendCommand( ssCommand.str() );

  ssCommand.str("");
  ssCommand << "Chart_SetXAxisLabel " << chartID
  << " \"" << msXLabel << "\"";
  manager.SendCommand( ssCommand.str() );

  ssCommand.str("");
  ssCommand << "Chart_SetYAxisLabel " << chartID
  << " \"" << msYLabel << "\"";
  manager.SendCommand( ssCommand.str() );

  ssCommand.str("");
  ssCommand << "Chart_SetInfo " << chartID
  << " \"" << msInfo << "\"";
  manager.SendCommand( ssCommand.str() );

  ssCommand.str("");
  ssCommand << "Chart_SetShowLegend " << chartID << " "
  << (mbShowLegend ? "true" : "false");
  manager.SendCommand( ssCommand.str() );

  // Clear the existing data.
  ssCommand.str("");
  ssCommand << "Chart_ClearData " << chartID;
  manager.SendCommand( ssCommand.str() );

  // Make some tcl commands with all our data groups and send it.
  map<int,list<PointData> >::iterator tGroup;
  for ( tGroup = mPointData.begin(); tGroup != mPointData.end(); ++tGroup ) {

    int nGroup = tGroup->first;
    list<PointData>& lPoints = tGroup->second;

    // Send the tcl commands to set the group data.
    ssCommand.str("");
    ssCommand << "Chart_SetGroupConnected " << chartID << " "
    << nGroup << " " << mGroupData[nGroup].mbConnected;
    manager.SendCommand( ssCommand.str() );

    if ( mGroupData[nGroup].msLabel != "" ) {
      ssCommand.str("");
      ssCommand << "Chart_SetGroupLabel " << chartID << " "
      << nGroup << " \"" << mGroupData[nGroup].msLabel << "\"";
      manager.SendCommand( ssCommand.str() );
    }

    ssCommand.str("");
    ssCommand << "Chart_SetGroupColor " << chartID << " "
    << nGroup << " "
    << mGroupData[nGroup].mColorRGBi[0] << " "
    << mGroupData[nGroup].mColorRGBi[1] << " "
    << mGroupData[nGroup].mColorRGBi[2];
    manager.SendCommand( ssCommand.str() );

    // Go through the point list and build the points command.
    list<PointData>::iterator tPoint;
    ssCommand.str("");
    ssCommand << "Chart_SetPointData " << chartID << " "
    << nGroup << " [list ";
    for ( tPoint = lPoints.begin(); tPoint != lPoints.end(); ++tPoint ) {
      PointData& point = *tPoint;
      ssCommand << "[list x " << point.mX << " y " << point.mY
      << " label \"" << point.msLabel << "\"] ";
    }
    ssCommand << "]";
    manager.SendCommand( ssCommand.str() );
  }


  if ( mXAxisMarkerData.size() > 0 ) {

    ssCommand.str("");
    ssCommand << "Chart_SetXAxisMarkers " << chartID << " "
    << " [list ";

    list<MarkerData>::iterator tMarker;
    for ( tMarker = mXAxisMarkerData.begin(); tMarker != mXAxisMarkerData.end();
          ++tMarker ) {

      MarkerData& marker = *tMarker;
      /*
      std::string msLabel;
      float mValue;
      int mColorRGBi[3];
      */
      ssCommand << "[list value " << marker.mValue << " label \""
      << marker.msLabel << "\" red " << marker.mColorRGBi[0]
      << " green " << marker.mColorRGBi[1]
      << " blue " << marker.mColorRGBi[2] << "] ";
    }

    ssCommand << "]";
    manager.SendCommand( ssCommand.str() );
  }


  // Make sure the window is showing.
  ssCommand.str("");
  ssCommand << "Chart_ShowWindow " << chartID;
  manager.SendCommand( ssCommand.str() );

}

TclCommandListener::TclCommandResult
TclChartWindow::DoListenToTclCommand ( char* isCommand,
                                       int, char** iasArgv ) {

  // GenerateChartReport <chartID> <includeGroupLabel> <includePointLabel>
  // <includeX> <includeY> <fileName>
  if ( 0 == strcmp( isCommand, "GenerateChartReport" ) ) {
    int chartID;
    try {
      chartID = TclCommandManager::ConvertArgumentToInt( iasArgv[1] );
    } catch ( runtime_error& e ) {
      sResult = string("bad collectionID: ") + e.what();
      return error;
    }

    if ( mID == chartID ) {

      bool bIncludeGroupLabel = false;
      try {
        bIncludeGroupLabel =
          TclCommandManager::ConvertArgumentToBoolean( iasArgv[2] );
      } catch ( runtime_error& e ) {
        sResult = "bad includeGroupLabel \"" + string(iasArgv[2]) + "\"," + e.what();
        return error;
      }

      bool bIncludePointLabel = false;
      try {
        bIncludePointLabel =
          TclCommandManager::ConvertArgumentToBoolean( iasArgv[3] );
      } catch ( runtime_error& e ) {
        sResult = "bad includePointLabel \"" + string(iasArgv[3]) + "\"," + e.what();
        return error;
      }

      bool bIncludeX = false;
      try {
        bIncludeX =
          TclCommandManager::ConvertArgumentToBoolean( iasArgv[4] );
      } catch ( runtime_error& e ) {
        sResult = "bad includeX \"" + string(iasArgv[4]) + "\"," + e.what();
        return error;
      }

      bool bIncludeY = false;
      try {
        bIncludeY =
          TclCommandManager::ConvertArgumentToBoolean( iasArgv[5] );
      } catch ( runtime_error& e ) {
        sResult = "bad includeY \"" + string(iasArgv[5]) + "\"," + e.what();
        return error;
      }

      try {
        GenerateReport( string(iasArgv[6]), bIncludeGroupLabel,
                        bIncludePointLabel, bIncludeX, bIncludeY );
      } catch ( runtime_error& e ) {
        sResult = string("Error generating report: ") + e.what();
        return error;
      }

      return ok;
    }
  }

  return ok;
}

TclChartWindowStaticListener::TclChartWindowStaticListener () {

  TclCommandManager& commandMgr = TclCommandManager::GetManager();
  commandMgr.AddCommand( *this, "DeleteTclChartWindow", 1, "chartID",
                         "Deletes a TclChartWindow." );

}

TclCommandListener::TclCommandResult
TclChartWindowStaticListener::DoListenToTclCommand ( char* isCommand,
    int, char** iasArgv ) {

  // DeleteTclChartWindow <chartID>
  if ( 0 == strcmp( isCommand, "DeleteTclChartWindow" ) ) {

    // Get our chart ID.
    int chartID;
    try {
      chartID = TclCommandManager::ConvertArgumentToInt( iasArgv[1] );
    } catch ( runtime_error& e ) {
      sResult = string("bad chartID: ") + e.what();
      return error;
    }

    // Try to find the chart window.
    TclChartWindow* chart;
    try {
      chart = (TclChartWindow*) &TclChartWindow::FindByID( chartID );
    } catch (...) {
      throw runtime_error( "Couldn't find the chart window." );
    }

    // Delete the chart window object.
    if ( NULL != chart )
      delete chart;
  }

  return ok;
}
