/**
 * @file  ChartWindow.h
 * @brief Abstract chart window class
 *
 * This is a virtual class meant to be reimplemented in a specific
 * widget framework. See TclChartWindow and QtChartWindow. Clients can
 * use this class to display data in a standard chart, plot, or graph.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:36 $
 *    $Revision: 1.10 $
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


#ifndef ChartWindow_h
#define ChartWindow_h

#include "string_fixed.h"
#include <list>
#include "DebugReporter.h"
#include "IDTracker.h"
#include "Broadcaster.h"

class ChartWindowFactory;

class ChartWindow : public DebugReporter,
      public IDTracker<ChartWindow>,
      public Broadcaster  // chartDeleted
{

public:

  virtual ~ChartWindow ();

  // Sets the factory to the correct subclass factory. Use the
  // NewChartWindow() to make sure the correct subclass is being made.
  static void SetFactory ( ChartWindowFactory* iFactory );
  static ChartWindow* NewChartWindow ();

  class PointData {
  public:
    float mX, mY;
    std::string msLabel;
  };

  // How to display a group.
  class GroupData {
  public:
    GroupData() : msLabel(""), mbConnected(false) {
      mColorRGBi[0] = mColorRGBi[1] = mColorRGBi[2] = 0;
    }
    std::string msLabel;
    bool mbConnected;
    int mColorRGBi[3];
  };

  // Display data for a marker.
  class MarkerData {
  public:
    MarkerData() : msLabel(""), mValue(0) {
      mColorRGBi[0] = mColorRGBi[1] = mColorRGBi[2] = 0;
    }
    std::string msLabel;
    float mValue;
    int mColorRGBi[3];
  };

  // Clear the existing data.
  void ClearData ();

  // Provide a list of values to use in the chart. If none is provided
  // for an axis, the default will be used.
  void SetPointData ( std::list<PointData>& iaData );
  void AddPointData ( PointData& iData );

  // Provide a list of values if you want multiple groups. Each can
  // have different color or drawing styles.
  void SetPointData ( int inGroup, std::list<PointData>& iaData );
  void AddPointData ( int inGroup, PointData& iData );

  // Provide drawing style for a group.
  void SetGroupLabel     ( int inGroup, std::string isLabel );
  void SetGroupConnected ( int inGroup, bool ibConnected );
  void SetGroupColor     ( int inGroup, int iColorRGBi[3] );
  void GetGroupColor     ( int inGroup, int oColorRGBi[3] );

  // Set the labels in the graph.
  void SetTitle ( std::string isTitle );
  void SetXAxisLabel ( std::string isLabel );
  void SetYAxisLabel ( std::string isLabel );

  // Vertical markers.
  void SetXAxisMarkers ( std::list<MarkerData>& iaData );
  void AddXAxisMarker  ( MarkerData& iData );

  // An info label under the graph.
  void SetInfo ( std::string isInfo );

  // Configuration options for the graph.
  void SetShowLegend ( bool ibShowLegend ) {
    mbShowLegend = ibShowLegend;
  }
  bool GetShowLegend () {
    return mbShowLegend;
  }

  // Close the chart window.
  virtual void Close () = 0;

  // Draw the graph.
  virtual void Draw () = 0;

  // Output a text file. Specify the columns to generate and other
  // settings.
  void GenerateReport ( std::string ifnReport,
                        bool ibIncludeGroupColumn,
                        bool ibIncludeLabelColumn,
                        bool ibIncludeXColumn,
                        bool ibIncludeYColumn );

protected:

  // Init group data if needed.
  void InitGroupDataIfNotSet ( int inGroup );
  void SetGroupData          ( int inGroup, GroupData& iData );

  // Creates a new chart window. Protected so they have to use the
  // static NewChartWindow function to create a member of the proper
  // subclass.
  ChartWindow ();

  std::map<int,std::list<PointData> > mPointData;
  std::map<int,GroupData> mGroupData;
  std::list<MarkerData> mXAxisMarkerData;
  std::string msTitle;
  std::string msXLabel;
  std::string msYLabel;
  std::string msInfo;
  bool mbShowLegend;

private:

  static ChartWindowFactory* sFactory;
};

// A factory class for the ChartWindow so it can create ChartWindow
// subclasses. ChartWindow subclasses should also have their own
// subclass of the ChartWindowFactory and pass it to the ChartWindow's
// SetFactory().
class ChartWindowFactory {
public:
  ChartWindowFactory() {};
  virtual ~ChartWindowFactory() {};

  virtual ChartWindow* NewChartWindow() {
    throw std::runtime_error( "Default chart window constructor "
                              "is being thrown." );
  }
};


#endif

