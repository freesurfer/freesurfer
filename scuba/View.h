#ifndef View_h
#define View_h

#include <string>
#include "IDTracker.h"
#include "DebugReporter.h"
#include "InputState.h"
#include "TclCommandManager.h"

// This is a base class used by ScubaFrame as an interface to views,
// or the internal components of the Frame. Each Frame can have
// multiple non-overlapping Views, and each View can have multiple
// layers.

class View : public DebugReporter, public IDTracker<View>, public TclCommandListener {
  
  friend class ViewTester;

public:
  View();
  virtual ~View();

  void Draw();
  void Reshape( int iWidth, int iHeight );
  void Timer();
  void MouseMoved( int inX, int inY, InputState& iState );
  void MouseUp( int inX, int inY, InputState& iState );
  void MouseDown( int inX, int inY, InputState& iState );
  void KeyDown( int inX, int inY, InputState& iState );
  void KeyUp( int inX, int inY, InputState& iState );

  void SetWidth( int iWidth ) { mWidth = iWidth; }
  void SetHeight( int iHeight ) { mHeight = iHeight; }

  void SetLabel( std::string isLabel ) { msLabel = isLabel; }
  std::string GetLabel() { return msLabel; }

  virtual void DoListenToTclCommand ( char* iCommand, int iArgc, char** iArgv );

 protected:
  virtual void DoDraw();
  virtual void DoReshape( int iWidth, int iHeight );
  virtual void DoTimer();
  virtual void DoMouseMoved( int inX, int inY, InputState& iState );
  virtual void DoMouseUp( int inX, int inY, InputState& iState );
  virtual void DoMouseDown( int inX, int inY, InputState& iState );
  virtual void DoKeyDown( int inX, int inY, InputState& iState );
  virtual void DoKeyUp( int inX, int inY, InputState& iState );

  int mWidth;
  int mHeight;
  
  std::string msLabel;
};


class ViewFactory {
 public:
  virtual View* NewView() { 
    return new View();
  }
};



#endif
