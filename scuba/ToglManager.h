#ifndef ToglManager_h
#define ToglManager_h

#include <map>
extern "C" {
#include "togl.h"
}
#include "DebugReporter.h"
#include "InputState.h"

class ToglFrame : public DebugReporter {

 public:
  typedef int ID;

  // Constructor takes an ID that is ultimately set by the Tcl script.
  ToglFrame( ID iID );
  virtual ~ToglFrame();

  // Accessor for the ID
  ID GetID() const { return mID; }

  // These callabacks are called by the ToglManager. Subclasses
  // can't/shouldn't override these as these setup the environment
  // before the overridable versions are called.
  void Draw();
  void Reshape( int iWidth, int iHeight );
  void Timer();
  void MouseMoved( int inX, int inY, InputState& iState );
  void MouseUp( int inX, int inY, InputState& iState );
  void MouseDown( int inX, int inY, InputState& iState );
  void KeyDown( int inX, int inY, InputState& iState );
  void KeyUp( int inX, int inY, InputState& iState );
  
  // These manage flags that the ToglManager will check to see if the
  // frame wants a redisplay. The frame should call RequestRedisplay()
  // to request one. Then ToglManager will ask the frame if it wants
  // one by calling WantRedisplay() on it, and notify that the
  // redisplay has been posted with RedisplayPosted().
  void RequestRedisplay() { mbPostRedisplay = true; }
  bool WantRedisplay() const { return mbPostRedisplay; }
  void RedisplayPosted() { mbPostRedisplay = false; }

 protected:

  // These are the overridable functions that subclass frames can use
  // to implement specific behavior.
  virtual void DoDraw();
  virtual void DoReshape();
  virtual void DoTimer();
  virtual void DoMouseMoved( int inX, int inY, InputState& iState  );
  virtual void DoMouseUp( int inX, int inY, InputState& iState );
  virtual void DoMouseDown( int inX, int inY, InputState& iState );
  virtual void DoKeyDown( int inX, int inY, InputState& iState );
  virtual void DoKeyUp( int inX, int inY, InputState& iState );

  // These are set by the Reshape() function. These should not be
  // manually set by subclasses.
  int mHeight;
  int mWidth;

  // Frame ID.
  ID mID;

  // Redisplay requested flag.
  bool mbPostRedisplay;
};



// A factory class for the ToglManager so it can create ToglFrame
// subclasses. ToglFrame subclasses should also have their own
// subclass of the ToglFrameFactory and pass it to the ToglManager's
// SetFrameFactory().
class ToglFrameFactory {
 public:
  virtual ToglFrame* NewToglFrame( ToglFrame::ID iID ) { 
    return new ToglFrame( iID );
  }
};



class ToglManager {

 public:
  // These are the Togl callbacks that we will register to get event
  // notifications from the Togl system. ToglManager determines to
  // which window they apply and calls ToglFrame functions in a
  // Togl-free way.
  static void DrawCallback ( struct Togl* iTogl );
  static void CreateCallback ( struct Togl* iTogl );
  static void DestroyCallback ( struct Togl* iTogl );
  static void ReshapeCallback ( struct Togl* iTogl );
  static void TimerCallback ( struct Togl* iTogl );
  
  // These are the Togl frame-specific callbacks that are attached to
  // the Togl Tcl/Tk object with the Tk bind command. ToglManager
  // determines to which window they apply, parses the arguments, and
  // calls ToglFrame functions in a Togl-free way.
  static int MouseMotionCallback ( struct Togl* iTogl, int, char* iArgv[] );
  static int MouseDownCallback ( struct Togl* iTogl, int, char* iArgv[] );
  static int MouseUpCallback ( struct Togl* iTogl, int iArgc, char* iArgv[] );
  static int KeyDownCallback ( struct Togl* iTogl, int iArgc, char* iArgv[] );
  static int KeyUpCallback ( struct Togl* iTogl, int iArgc, char* iArgv[] );
  static int ExitCallback ( struct Togl* iTogl, int iArgc, char* iArgv[] );

  // The main entry point should call this function to register all
  // the callbacks.
  static void InitializeTogl ( Tcl_Interp *iInterp );

  // Gets the static instance of the manager.
  static ToglManager& GetManager();

  // Sets the factory to use for creating new frames.
  static void SetFrameFactory( ToglFrameFactory* iFactory ) { 
    mFactory = iFactory; 
  }

 protected:
  // Maps Window IDs to frames.
  static std::map<ToglFrame::ID,ToglFrame*> mFrames;

  static ToglFrameFactory* mFactory;

  // Our modifers for shift, contrl, and alt keys.
  static InputState mState;

};


#endif
