#ifndef ToglManager_h
#define ToglManager_h

#include <map>
extern "C" {
#include "togl.h"
}
#include "string_fixed.h"
#include "WindowFrame.h"
#include "DebugReporter.h"
#include "InputState.h"

class ToglFrame : public WindowFrame {

 public:
  ToglFrame( ID iID ) : WindowFrame( iID ) {}
  virtual ~ToglFrame() {}

  virtual void Reshape( int iWidth, int iHeight );
};



// A factory class for the ToglManager so it can create ToglFrame
// subclasses. ToglFrame subclasses should also have their own
// subclass of the ToglFrameFactory and pass it to the ToglManager's
// SetFrameFactory().
class ToglFrameFactory : public WindowFrameFactory {
 public:
  virtual WindowFrame* NewWindowFrame( WindowFrame::ID iID ) { 
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
  static void SetFrameFactory( WindowFrameFactory* iFactory ) { 
    mFactory = iFactory; 
  }

 protected:

  static inline int YFlip ( ToglFrame* iFrame, int iY ) {
    return (iFrame->GetHeight() - iY); }

  // Maps Window IDs to frames.
  static std::map<ToglFrame::ID,ToglFrame*> mFrames;

  static WindowFrameFactory* mFactory;

  // Our modifers for shift, contrl, and alt keys.
  static InputState mState;

};


#endif
