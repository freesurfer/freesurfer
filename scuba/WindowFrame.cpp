#include "string_fixed.h"
#include <stdexcept>
#include "WindowFrame.h"

using namespace std;



WindowFrame::WindowFrame( ID iID ) {
}

WindowFrame::~WindowFrame() {

}

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

  if( iWindow[0] > 0 && iWindow[0] < mWidth-1 &&
      iWindow[1] > 0 && iWindow[1] < mHeight-1 ) {
    this->DoMouseMoved( iWindow, iInput );
  }
}

void
WindowFrame::MouseUp( int iWindow[2], InputState& iInput ) {

  if( iWindow[0] > 0 && iWindow[0] < mWidth-1 &&
      iWindow[1] > 0 && iWindow[1] < mHeight-1 ) {
    this->DoMouseUp( iWindow, iInput );
  }
}

void
WindowFrame::MouseDown( int iWindow[2], InputState& iInput ) {

  if( iWindow[0] > 0 && iWindow[0] < mWidth-1 &&
      iWindow[1] > 0 && iWindow[1] < mHeight-1 ) {
    this->DoMouseDown( iWindow, iInput );
  }
}

void
WindowFrame::KeyDown( int iWindow[2], InputState& iInput ) {

  if( iWindow[0] > 0 && iWindow[0] < mWidth-1 &&
      iWindow[1] > 0 && iWindow[1] < mHeight-1 ) {
    this->DoKeyDown( iWindow, iInput );
  }
}

void
WindowFrame::KeyUp( int iWindow[2], InputState& iInput ) {

  if( iWindow[0] > 0 && iWindow[0] < mWidth-1 &&
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
WindowFrame::DoMouseMoved( int iWindow[2], InputState& iInput ) {

  cerr << "WindowFrame::DoMouseMoved " << mID << endl;
  DebugOutput( << "WindowFrame " << mID << ": DoMouseMoved()" );
}

void
WindowFrame::DoMouseUp( int iWindow[2], InputState& iInput ) {

  DebugOutput( << "WindowFrame " << mID << ": DoMouseUp()" );
}

void
WindowFrame::DoMouseDown( int iWindow[2], InputState& iInput ) {

  DebugOutput( << "WindowFrame " << mID << ": DoMouseDown()" );
}

void
WindowFrame::DoKeyDown( int iWindow[2], InputState& iInput ) {

  DebugOutput( << "WindowFrame " << mID << ": DoKeyDown()" );
}

void
WindowFrame::DoKeyUp( int iWindow[2], InputState& iInput ) {

  DebugOutput( << "WindowFrame " << mID << ": DoKeyUp()" );
}
