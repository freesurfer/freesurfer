#ifndef ScubaFrame_h
#define ScubaFrame_h

#include "string_fixed.h"
#include "ToglManager.h"
#include "TclCommandManager.h"
#include "IDTracker.h"
#include "View.h"
#include "ScubaToolState.h"

class ScubaFrame : public ToglFrame, public TclCommandListener {

  friend class ScubaFrameTester;

public:
  ScubaFrame( ID iID );
  virtual ~ScubaFrame();
  
  // View configurations. The cxx numbers spec the number of columns
  // in each row. So c22 is 2 columns in row 0, and 2 in row 1, while
  // c13 is 1 view in row 0 and 3 in row 1.
  enum ViewConfiguration { c1, c22, c44, c13 };
  void SetViewConfiguration( ViewConfiguration iConfig );

  virtual TclCommandResult
    DoListenToTclCommand( char* iCommand, int iArgc, char** iArgv );

  // Sets the factory to use for creating new frames.
  static void SetViewFactory( ViewFactory* const iFactory ) { 
    mFactory = iFactory; 
  }

protected:

  // Adjusts window coords for a view.
  void TranslateWindowToView ( int iWindow[2], int inCol, int inRow,
			       int oView[2] );
  
  // Sets the sizes for all of our views according to our current
  // configuration and view size.
  void SizeViewsToConfiguration ();

  // Impelmentations of the ToglFrame callbacks.
  virtual void DoDraw();
  virtual void DoReshape();
  virtual void DoTimer();
  virtual void DoMouseMoved( int iWindow[2], InputState& iInput );
  virtual void DoMouseUp( int iWindow[2], InputState& iInput );
  virtual void DoMouseDown( int iWindow[2], InputState& iInput );
  virtual void DoKeyDown( int iWindow[2], InputState& iInput );
  virtual void DoKeyUp( int iWindow[2], InputState& iInput );

  // Given a window location, returns a pointer to a view. Or could
  // throw an exception.
  View* FindViewAtWindowLoc( int iWindow[2], int* onCol, int* onRow );

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

  // The map of view pointers. The first is a map for columns and the
  // second is a map of those for rows.
  typedef std::map<int,View*> ViewColMap;
  typedef std::map<int, std::map<int,View*> > ViewRowMap;
  ViewRowMap mViews;

  // Uses this factory to make views.
  static ViewFactory* mFactory;

  // Our tool.
  ScubaToolState mTool;

};

// The factory passed to ToglManager so that this type of Frame is
// created.
class ScubaFrameFactory : public ToglFrameFactory {
public:
  virtual WindowFrame* NewWindowFrame( WindowFrame::ID iID ) { 
    ScubaFrame* frame = new ScubaFrame( iID );
    frame->SetViewConfiguration( ScubaFrame::c1 );
    frame->SetOutputStreamToCerr();
    return frame;
  }
};

#endif
