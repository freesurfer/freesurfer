#ifndef View_h
#define View_h

#include "string_fixed.h"
#include "IDTracker.h"
#include "DebugReporter.h"
#include "InputState.h"
#include "TclCommandManager.h"
#include "ScubaToolState.h"
#include "Listener.h"
#include "Broadcaster.h"

// This is a base class used by ScubaFrame as an interface to views,
// or the internal components of the Frame. Each Frame can have
// multiple non-overlapping Views, and each View can have multiple
// layers.

class View : public DebugReporter, 
	     public IDTracker<View>,
	     public TclCommandListener,
	     public Broadcaster,
	     public Listener
{
  
  friend class ViewTester;

public:
  View();
  virtual ~View();

  void Draw();
  void Reshape( int iWidth, int iHeight );
  void Timer();
  void MouseMoved( int iWindow[2], InputState& iInput, ScubaToolState& iTool );
  void MouseUp( int iWindow[2], InputState& iInput, ScubaToolState& iTool );
  void MouseDown( int iWinodw[2], InputState& iInput, ScubaToolState& iTool );
  void KeyDown( int iWindow[2], InputState& iInput, ScubaToolState& iTool );
  void KeyUp( int iWindow[2], InputState& iInput, ScubaToolState& iTool );

  void SetWidth( int iWidth ) { mWidth = iWidth; }
  void SetHeight( int iHeight ) { mHeight = iHeight; }

  void SetLabel( std::string isLabel ) { msLabel = isLabel; }
  std::string GetLabel() { return msLabel; }

  virtual TclCommandResult
    DoListenToTclCommand ( char* isCommand, int iArgc, char** iArgv );

  // Handle broadcast messages.
  virtual void
    DoListenToMessage ( std::string isMessage, void* iData );

  // Redisplay posters.
  void RequestRedisplay() { mbPostRedisplay = true; }
  bool WantRedisplay() const { return mbPostRedisplay; }
  void RedisplayPosted() { mbPostRedisplay = false; }

  // Set the frame visiblity flag.
  void SetVisibleInFrame ( bool ibVisible ) { mbVisibleInFrame = ibVisible; }
  bool IsVisibleInFrame () { return mbVisibleInFrame; }

 protected:
  virtual void DoDraw();
  virtual void DoReshape( int iWidth, int iHeight );
  virtual void DoTimer();
  virtual void DoMouseMoved( int iWindow[2], 
			     InputState& iInput, ScubaToolState& iTool );
  virtual void DoMouseUp( int iWindow[2],
			  InputState& iInput, ScubaToolState& iTool );
  virtual void DoMouseDown( int iWindow[2],
			    InputState& iInput, ScubaToolState& iTool );
  virtual void DoKeyDown( int iWindow[2],
			  InputState& iInput, ScubaToolState& iTool );
  virtual void DoKeyUp( int iWindow[2],
			InputState& iInput, ScubaToolState& iTool );

  int mWidth;
  int mHeight;
  
  std::string msLabel;

  // Redisplay requested flag.
  bool mbPostRedisplay;

  // A flag set by the Frame to tell the view whether it's visible or
  // not.
  bool mbVisibleInFrame;
};


class ViewFactory {
 public:
  virtual ~ViewFactory() {};
  virtual View* NewView() { 
    return new View();
  }
};



#endif
