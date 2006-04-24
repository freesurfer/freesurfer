#ifndef TclChartWindow_h
#define TclChartWindow_h

#include "ChartWindow.h"
#include "string_fixed.h"
#include <list>

// A Tcl implementation of the ChartWindow interface. Uses BLT code in
// TclChartWindow.tcl to implement the window itself. The Draw()
// function uses the current ChartWindow's setting to set up the BLT
// window.

class TclChartWindow : public ChartWindow {

  friend class TclChartWindowFactory;

 public:

  // Sends tcl commands to the TclChartWindow.tcl code.
  void Draw ();

 protected:

  // If not already done so, this will send a Tcl command to load our
  // TclChartWindow.tcl file and call Chart_Init.
  TclChartWindow ();
  virtual ~TclChartWindow ();

  static bool sbInitedTclFile;
};


// This factory will create TclChartWindow objects. Pass it to
// ChartWindow::SetFactory in the program startup.
class TclChartWindowFactory : public ChartWindowFactory {
 public:
  TclChartWindowFactory() {};
  virtual ~TclChartWindowFactory() {};

  virtual ChartWindow* NewChartWindow() { 
    return new TclChartWindow();
  }
};

#endif
