#include <stdexcept>
#include "ChartWindow.h"
#include <fstream>

using namespace std;

ChartWindowFactory* ChartWindow::sFactory = NULL;

void
ChartWindow::SetFactory ( ChartWindowFactory* iFactory ) {
  
  sFactory = iFactory;
}

ChartWindow*
ChartWindow::NewChartWindow () {

  if( NULL != sFactory ) {
    return sFactory->NewChartWindow();
  } else {
    throw runtime_error( "No factory defined" );
  }

  return NULL;
}

ChartWindow::ChartWindow ():
  mbShowLegend( false )
{

}

ChartWindow::~ChartWindow () {

}

void
ChartWindow::ClearData () {

  mPointData.clear();
}

void
ChartWindow::SetPointData ( list<PointData>& iaData ) {

  mPointData = iaData;
}

void
ChartWindow::AddPointData ( PointData& iData ) {

  mPointData.push_back( iData );
}


void
ChartWindow::SetTitle ( string isTitle ) {

  msTitle = isTitle;
}

void
ChartWindow::SetXAxisLabel ( string isLabel ) {

  msXLabel = isLabel;
}

void
ChartWindow::SetYAxisLabel ( string isLabel ) {

  msYLabel = isLabel;
}

void
ChartWindow::SetInfo ( string isInfo ) {

  msInfo = isInfo;
}

void
ChartWindow::GenerateReport ( string ifnReport,
			       bool ibIncludeLabelColumn,
			       bool ibIncludeXColumn,
			       bool ibIncludeYColumn ) {

  try {
    // Check file name first.
    ofstream fReport( ifnReport.c_str(), ios::out );

    list<PointData>::iterator tPoint;
    for( tPoint = mPointData.begin(); tPoint != mPointData.end(); ++tPoint ) {

      PointData& point = *tPoint;

      // Output the data they want.
      if( ibIncludeLabelColumn ) {
	fReport << point.msLabel << "\t";
      }
      if( ibIncludeXColumn ) {
	fReport << point.mX << "\t";
      }
      if( ibIncludeYColumn ) {
	fReport << point.mY;
      }
      fReport << endl;
    }
  }
  catch( exception& e ) {
    throw runtime_error( "Error writing " + ifnReport + ": " + e.what() );
  }
  catch( ... ) {
    throw runtime_error( "Error writing " + ifnReport );
  }
}

