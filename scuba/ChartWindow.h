#ifndef ChartWindow_h
#define ChartWindow_h

#include "string_fixed.h"
#include <list>
#include "DebugReporter.h"
#include "IDTracker.h"
#include "Broadcaster.h"

// This is a virtual class meant to be reimplemented in a specific
// widget framework. See TclChartWindow and QtChartWindow. Clients
// can use this class to display data in a standard chart, plot, or
// graph.

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

  // Clear the existing data.
  void ClearData ();

  // Provide a list of values to use in the chart. If none is provided
  // for an axis, the default will be used.
  void SetPointData ( std::list<PointData>& iaData );
  void AddPointData ( PointData& iData );

  // Set the labels in the graph.
  void SetTitle ( std::string isTitle );
  void SetXAxisLabel ( std::string isLabel );
  void SetYAxisLabel ( std::string isLabel );

  // An info label under the graph.
  void SetInfo ( std::string isInfo );

  // Configuration options for the graph.
  void SetShowLegend ( bool ibShowLegend ) { mbShowLegend = ibShowLegend; }
  bool GetShowLegend () { return mbShowLegend; }

  // Close the chart window.
  virtual void Close () = 0;

  // Draw the graph.
  virtual void Draw () = 0;

  // Output a text file. Specify the columns to generate and other
  // settings.
  void GenerateReport ( std::string ifnReport,
			bool ibIncludeLabelColumn,
			bool ibIncludeXColumn,
			bool ibIncludeYColumn );

 protected:

  // Creates a new chart window. Protected so they have to use the
  // static NewChartWindow function to create a member of the proper
  // subclass.
  ChartWindow ();

  std::list<PointData> mPointData;
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

