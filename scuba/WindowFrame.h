#ifndef WindowManager_h
#define WindowManager_h

#include <map>
#include "DebugReporter.h"
#include "InputState.h"
#include "IDTracker.h"

class WindowFrame : public DebugReporter,
		    public IDTracker<WindowFrame> {

 public:
  typedef int ID;

  // Constructor takes an ID that is ultimately set by the Tcl script.
  WindowFrame( ID iID );
  virtual ~WindowFrame();

  // Accessor for the ID
  ID GetID() const { return mID; }

  // These callabacks are called by the WindowManager. Subclasses
  // can't/shouldn't override these as these setup the environment
  // before the overridable versions are called.
  virtual void Draw();
  virtual void Reshape( int iWidth, int iHeight );
  virtual void Timer();
  virtual void MouseMoved( int iWindow[2], InputState& iInput );
  virtual void MouseUp( int iWindow[2], InputState& iInput );
  virtual void MouseDown( int iWindow[2], InputState& iInput );
  virtual void KeyDown( int iWindow[2], InputState& iInput );
  virtual void KeyUp( int iWindow[2], InputState& iInput );
  
  // These manage flags that the WindowManager will check to see if the
  // frame wants a redisplay. The frame should call RequestRedisplay()
  // to request one. Then WindowManager will ask the frame if it wants
  // one by calling WantRedisplay() on it, and notify that the
  // redisplay has been posted with RedisplayPosted().
  void RequestRedisplay() { mbPostRedisplay = true; }
  bool WantRedisplay() const { return mbPostRedisplay; }
  void RedisplayPosted() { mbPostRedisplay = false; }

  int GetWidth () const { return mWidth; }
  int GetHeight () const { return mHeight; }

 protected:

  // These are the overridable functions that subclass frames can use
  // to implement specific behavior.
  virtual void DoDraw();
  virtual void DoReshape();
  virtual void DoTimer();
  virtual void DoMouseMoved( int iWindow[2], InputState& iInput  );
  virtual void DoMouseUp( int iWindow[2], InputState& iInput );
  virtual void DoMouseDown( int iWindow[2], InputState& iInput );
  virtual void DoKeyDown( int iWindow[2], InputState& iInput );
  virtual void DoKeyUp( int iWindow[2], InputState& iInput );

  // These are set by the Reshape() function. These should not be
  // manually set by subclasses.
  int mHeight;
  int mWidth;

  // Redisplay requested flag.
  bool mbPostRedisplay;

  // Keep track of last mouse moved coords so we interpolate and send
  // moved commands for every point in between.
  int mLastMoved[2];
};



// A factory class for the WindowManager so it can create WindowFrame
// subclasses. WindowFrame subclasses should also have their own
// subclass of the WindowFrameFactory and pass it to the WindowManager's
// SetFrameFactory().
class WindowFrameFactory {
 public:
  virtual ~WindowFrameFactory() {};
  virtual WindowFrame* NewWindowFrame( WindowFrame::ID iID ) { 
    return new WindowFrame( iID );
  }
};

#endif
