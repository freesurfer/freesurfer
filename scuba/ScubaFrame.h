#ifndef ScubaFrame_h
#define ScubaFrame_h

#include <string>
#include "ToglManager.h"
#include "TclCommandManager.h"
#include "IDTracker.h"
#include "View.h"


class ScubaFrame : public ToglFrame, public TclCommandListener {

  friend class ScubaFrameTester;

public:
  ScubaFrame( ToglFrame::ID iID );
  virtual ~ScubaFrame();
  
  // View configurations. The cxx numbers spec the number of columns
  // in each row. So c22 is 2 columns in row 0, and 2 in row 1, while
  // c13 is 1 view in row 0 and 3 in row 1.
  enum ViewConfiguration { c11, c22, c44, c13 };
  void SetViewConfiguration( ViewConfiguration iConfig );

  virtual void DoListenToTclCommand( char* iCommand, int iArgc, char** iArgv );

  // Sets the factory to use for creating new frames.
  static void SetViewFactory( ViewFactory* const iFactory ) { 
    mFactory = iFactory; 
  }

protected:

  // Sets the sizes for all of our views according to our current
  // configuration and view size.
  void SizeViewsToConfiguration ();

  // Impelmentations of the ToglFrame callbacks.
  virtual void DoDraw();
  virtual void DoReshape();
  virtual void DoTimer();
  virtual void DoMouseMoved( int inX, int inY, InputState& iState );
  virtual void DoMouseUp( int inX, int inY, InputState& iState );
  virtual void DoMouseDown( int inX, int inY, InputState& iState );
  virtual void DoKeyDown( int inX, int inY, InputState& iState );
  virtual void DoKeyUp( int inX, int inY, InputState& iState );

  // Given a window location, returns a pointer to a view. Or could
  // throw an exception.
  View* FindViewAtWindowLoc( int iWindowX, int iWindowY,
			     int* onCol, int* onRow );

  // Access functions to get/set views at a col/row position. This can
  // be used even if the current view configuration doesn't have the
  // given col/row in range. However if GetView is called and the
  // col/row view hasn't been set yet, it will throw an exception.
  View* GetViewAtColRow( int iCol, int iRow );
  void SetViewAtColRow( int iCol, int iRow, View* const iView );

  // View configuration and the number of rows we have, and for each
  // row, the numbe rof columns in that row.
  ViewConfiguration mViewConfiguration;
  int mcRows;
  std::map<int,int> mcCols;

  // The selected view. Always in ViewConfiguration col/row bounds.
  int mnSelectedViewCol;
  int mnSelectedViewRow;

  // Keys.
  std::string msCycleKey;

  // The map of view pointers. The first is a map for columns and the
  // second is a map of those for rows.
  typedef std::map<int,View*> ViewColMap;
  typedef std::map<int, std::map<int,View*> > ViewRowMap;
  ViewRowMap mViews;

  // Uses this factory to make views.
  static ViewFactory* mFactory;
};

// The factory passed to ToglManager so that this type of Frame is
// created.
class ScubaFrameFactory : public ToglFrameFactory {
public:
  virtual ToglFrame* NewToglFrame( ToglFrame::ID iID ) { 
    ScubaFrame* frame = new ScubaFrame( iID );
    frame->SetOutputStreamToCerr();
    return frame;
  }
};

#endif
