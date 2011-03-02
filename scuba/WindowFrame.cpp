/**
 * @file  WindowFrame.cpp
 * @brief Top level window class
 *
 * These are created by ToglManager when the togl object in Tcl is
 * created. The Window receives events and manages a redisplay flag
 * that ToglManager uses to know when to post draw events.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:39 $
 *    $Revision: 1.19 $
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
 *
 * Terms and conditions for use, reproduction, distribution and contribution
 * are found in the 'FreeSurfer Software License Agreement' contained
 * in the file 'LICENSE' found in the FreeSurfer distribution, and here:
 *
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
 *
 * Reporting: freesurfer@nmr.mgh.harvard.edu
 *
 */


#include "string_fixed.h"
#include <stdexcept>
#include <math.h>
#include "WindowFrame.h"
#include "Point2.h"

using namespace std;

WindowFrame::WindowFrame() {
  mLastMoved[0] = 0;
  mLastMoved[1] = 1;
  mWidth = 0;
  mHeight = 0;
}

WindowFrame::~WindowFrame() {}

void
WindowFrame::Draw() {

  this->DoDraw();
}

void
WindowFrame::Reshape( int iWidth, int iHeight ) {

  mWidth = iWidth;
  mHeight = iHeight;

  this->DoReshape();
}

void
WindowFrame::Timer() {

  this->DoTimer();
}

void
WindowFrame::MouseMoved( int iWindow[2], InputState& iInput ) {

  if ( iWindow[0] >= 0 && iWindow[0] <= mWidth-1 &&
       iWindow[1] >= 0 && iWindow[1] <= mHeight-1 ) {

    if ( iInput.IsButtonDragEvent() ) {

      // Calculate the delta. Make sure there is one; if not, no mouse
      // moved event.
      float delta[2];
      delta[0] = iWindow[0] - mLastMoved[0];
      delta[1] = iWindow[1] - mLastMoved[1];

      if ( delta[0] != 0 || delta[1] != 0 ) {

        // Find the greater one (absolute value). Divide each delta by
        // this value to get the step; one will be 1.0, and the other
        // will be < 1.0.
        float greater = fabsf(delta[0]) > fabsf(delta[1]) ?
                        fabsf(delta[0]) : fabsf(delta[1]);

        delta[0] /= greater;
        delta[1] /= greater;

        // We step window coords in floats, but transform to ints. Start
        // at the last moved place.
        float windowF[2];
        int   windowI[2];
        windowF[0] = mLastMoved[0];
        windowF[1] = mLastMoved[1];

        // While we're not at the current location...  NOTE - this
        // should work, but jc can't handle simple float comparison,
        // so we have to do wacky 0 comparison instead.
        while ( !(  fabsf((float)iWindow[0] - windowF[0]) < 1.0 &&
                    fabsf((float)iWindow[1] - windowF[1]) < 1.0) ) {

          // Get an integer value and send it to the frame.
          windowI[0] = (int) rint( windowF[0] );
          windowI[1] = (int) rint( windowF[1] );

          this->DoMouseMoved( windowI, iInput );

          // Increment the float window coords.
          windowF[0] += delta[0];
          windowF[1] += delta[1];

          if ( windowF[0] < -10 || windowF[1] < -10 )
            exit( 1 );
        }
      }
    } else {

      this->DoMouseMoved( iWindow, iInput );
    }

    // Save this position.
    mLastMoved[0] = iWindow[0];
    mLastMoved[1] = iWindow[1];
  }

}

void
WindowFrame::MouseUp( int iWindow[2], InputState& iInput ) {

  if ( iWindow[0] > 0 && iWindow[0] < mWidth-1 &&
       iWindow[1] > 0 && iWindow[1] < mHeight-1 ) {
    this->DoMouseUp( iWindow, iInput );
  }
}

void
WindowFrame::MouseDown( int iWindow[2], InputState& iInput ) {

  if ( iWindow[0] > 0 && iWindow[0] < mWidth-1 &&
       iWindow[1] > 0 && iWindow[1] < mHeight-1 ) {
    this->DoMouseDown( iWindow, iInput );
  }
}

void
WindowFrame::KeyDown( int iWindow[2], InputState& iInput ) {

  if ( iWindow[0] > 0 && iWindow[0] < mWidth-1 &&
       iWindow[1] > 0 && iWindow[1] < mHeight-1 ) {
    this->DoKeyDown( iWindow, iInput );
  }
}

void
WindowFrame::KeyUp( int iWindow[2], InputState& iInput ) {

  if ( iWindow[0] > 0 && iWindow[0] < mWidth-1 &&
       iWindow[1] > 0 && iWindow[1] < mHeight-1 ) {
    this->DoKeyUp( iWindow, iInput );
  }
}

void
WindowFrame::DoDraw() {

  DebugOutput( << "WindowFrame " << mID << ": DoDraw()" );
}

void
WindowFrame::DoReshape() {

  DebugOutput( << "WindowFrame " << mID << ": DoReshape()" );
}

void
WindowFrame::DoTimer() {

  DebugOutput( << "WindowFrame " << mID << ": DoTimer()" );
}

void
WindowFrame::DoMouseMoved( int[2], InputState& ) {

  cerr << "WindowFrame::DoMouseMoved " << mID << endl;
  DebugOutput( << "WindowFrame " << mID << ": DoMouseMoved()" );
}

void
WindowFrame::DoMouseUp( int[2], InputState& ) {

  DebugOutput( << "WindowFrame " << mID << ": DoMouseUp()" );
}

void
WindowFrame::DoMouseDown( int[2], InputState& ) {

  DebugOutput( << "WindowFrame " << mID << ": DoMouseDown()" );
}

void
WindowFrame::DoKeyDown( int[2], InputState& ) {

  DebugOutput( << "WindowFrame " << mID << ": DoKeyDown()" );
}

void
WindowFrame::DoKeyUp( int[2], InputState& ) {

  DebugOutput( << "WindowFrame " << mID << ": DoKeyUp()" );
}
