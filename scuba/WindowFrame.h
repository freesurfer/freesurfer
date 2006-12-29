/**
 * @file  WindowFrame.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 02:09:15 $
 *    $Revision: 1.6 $
 *
 * Copyright (C) 2002-2007,
 * The General Hospital Corporation (Boston, MA). 
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
 *
 */


#ifndef WindowFrame_h
#define WindowFrame_h

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
  ID GetID() const {
    return mID;
  }

  // These callabacks are called by the WindowFrame. Subclasses
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

  // These manage flags that the WindowFrame will check to see if the
  // frame wants a redisplay. The frame should call RequestRedisplay()
  // to request one. Then WindowFrame will ask the frame if it wants
  // one by calling WantRedisplay() on it, and notify that the
  // redisplay has been posted with RedisplayPosted().
  void RequestRedisplay() {
    mbPostRedisplay = true;
  }
  bool WantRedisplay() const {
    return mbPostRedisplay;
  }
  void RedisplayPosted() {
    mbPostRedisplay = false;
  }

  int GetWidth () const {
    return mWidth;
  }
  int GetHeight () const {
    return mHeight;
  }

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



// A factory class for the WindowFrame so it can create WindowFrame
// subclasses. WindowFrame subclasses should also have their own
// subclass of the WindowFrameFactory and pass it to the WindowFrame's
// SetFrameFactory().
class WindowFrameFactory {
public:
  virtual ~WindowFrameFactory() {};
  virtual WindowFrame* NewWindowFrame( WindowFrame::ID iID ) {
    return new WindowFrame( iID );
  }
};

#endif
