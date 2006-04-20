#ifndef TclChartWindow_h
#define TclChartWindow_h

#include "ChartWindow.h"
#include "string_fixed.h"
#include <list>

class TclChartWindow : public ChartWindow {

  friend class TclChartWindowFactory;

 public:

  void Draw ();

 protected:

  TclChartWindow ();
  virtual ~TclChartWindow ();
};

class TclChartWindowFactory : public ChartWindowFactory {
  virtual ~TclChartWindowFactory() {};
  virtual ChartWindow* NewChartWindow() { 
    return new TclChartWindow();
  }
};

#endif
