#include "string_fixed.h"
#include <sstream>
#include "TclChartWindow.h"
#include "TclCommandManager.h"

using namespace std;

bool TclChartWindow::sbInitedTclFile = false;

TclChartWindow::TclChartWindow () :
  ChartWindow() {
  msTitle = "Chart";

  TclCommandManager& manager = TclCommandManager::GetManager();
  manager.SendCommand( "Chart_Init" );

  if( !sbInitedTclFile ) {
    
    stringstream ssCommand;

    ssCommand.str("");
    ssCommand << "LoadScubaSupportFile TclChartWindow.tcl";
    manager.SendCommand( ssCommand.str() );

    sbInitedTclFile = true;
  }

}

TclChartWindow::~TclChartWindow () {

}

void
TclChartWindow::Draw () {

  TclCommandManager& manager = TclCommandManager::GetManager();

  int chartID = 1;

  stringstream ssCommand;

  ssCommand.str("");
  ssCommand << "Chart_BuildWindow " << chartID;
  manager.SendCommand( ssCommand.str() );
  
  ssCommand.str("");
  ssCommand << "Chart_SetShowLegend " << chartID << " "
	    << (mbShowLegend ? "true" : "false");
  manager.SendCommand( ssCommand.str() );

  list<PointData>::iterator tPoint;
  ssCommand.str("");
  ssCommand << "Chart_SetData " << chartID << " [list ";
  for( tPoint = mPointData.begin(); tPoint != mPointData.end(); ++tPoint ) {
    PointData& point = *tPoint;
    ssCommand << "[list x " << point.mX << " y " << point.mY
	      << " label \"" << point.msLabel << "\"] ";
  }
  ssCommand << "]";
  manager.SendCommand( ssCommand.str() );

}
