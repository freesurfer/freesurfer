#ifndef ToglManager_h
#define ToglManager_h

extern "C" {
#include "togl.h"
}
#include "DebugReporter.h"
#include <map>


class ToglFrame : public DebugReporter {

 public:
  typedef int ID;

  ToglFrame( ID iID );
  virtual ~ToglFrame();

  void Draw();
  void Reshape( int iWidth, int iHeight );
  void Timer();
  void MouseMoved( int inX, int inY, int iButton, int iModifiers );
  void MouseUp( int inX, int inY, int iButton, int iModifers );
  void MouseDown( int inX, int inY, int iButton, int iModifers );
  void KeyDown( int inX, int inY, std::string isKey, int iModifers );
  void KeyUp( int inX, int inY, std::string isKey, int iModifers );
  
  ID GetID() const { return mID; }
  
 protected:
  virtual void DoDraw();
  virtual void DoReshape();
  virtual void DoTimer();
  virtual void DoMouseMoved( int inX, int inY, int iButton, int iModifiers );
  virtual void DoMouseUp( int inX, int inY, int iButton, int iModifers );
  virtual void DoMouseDown( int inX, int inY, int iButton, int iModifers );
  virtual void DoKeyDown( int inX, int inY, std::string isKey, int iModifers );
  virtual void DoKeyUp( int inX, int inY, std::string isKey, int iModifers );

  int mHeight;
  int mWidth;

  ID mID;
};

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
  // which window they apple and calls ToglFrame functions in a
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
};


#endif
