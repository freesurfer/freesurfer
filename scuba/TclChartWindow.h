#ifndef TclChartWindow_h
#define TclChartWindow_h

#include "string_fixed.h"
#include <list>
#include "ChartWindow.h"
#include "TclCommandManager.h"

// A Tcl implementation of the ChartWindow interface. Uses BLT code in
// TclChartWindow.tcl to implement the window itself. The Draw()
// function uses the current ChartWindow's setting to set up the BLT
// window.

class TclChartWindowStaticListener;

class TclChartWindow : public ChartWindow {

  friend class TclChartWindowFactory;

 public:

  virtual ~TclChartWindow ();

  // Sends tcl commands to the TclChartWindow.tcl code.
  void Draw ();
  
  // Sends tcl command to close the window.
  void Close ();

 protected:

  // If not already done so, this will send a Tcl command to load our
  // TclChartWindow.tcl file and call Chart_Init.
  TclChartWindow ();

  static bool sbInitedTclFile;
  static TclChartWindowStaticListener mStaticListener;
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

// This static listener listens for window delete callbacks sent by
// the TclChartWindow.tcl file and will get and delete the proper
// TclChartWindow object.
class TclChartWindowStaticListener : public TclCommandListener {

public:
  TclChartWindowStaticListener ();
  ~TclChartWindowStaticListener () {};

    virtual TclCommandResult
      DoListenToTclCommand ( char* isCommand, int iArgc, char** iArgv );
};




#endif
