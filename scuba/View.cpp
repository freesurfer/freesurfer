#include "View.h"

using namespace std;


template IDTracker<View>;
int IDTracker<View>::mNextID = 0;
map<int,View*> IDTracker<View>::mIDMap;


View::View() {
  mWidth = 0;
  mHeight = 0;
  mbPostRedisplay = false;
}

View::~View() {

}

void
View::Draw () {
  this->DoDraw();
}

void
View::Reshape ( int iWidth, int iHeight ) {
  mWidth = iWidth;
  mHeight = iHeight;
  this->DoReshape( iWidth, iHeight );
}

void 
View::Timer() {

  this->DoTimer();
}

void
View::MouseMoved( int inX, int inY, InputState& iState ) {

  this->DoMouseMoved( inX, inY, iState );
}

void
View::MouseUp( int inX, int inY, InputState& iState ) {

  this->DoMouseUp( inX, inY, iState );
}

void
View::MouseDown( int inX, int inY, InputState& iState ) {

  this->DoMouseDown( inX, inY, iState );
}

void
View::KeyDown( int inX, int inY, InputState& iState ) {

  this->DoKeyDown( inX, inY, iState );
}

void
View::KeyUp( int inX, int inY, InputState& iState ) {

  this->DoKeyUp( inX, inY, iState );
}

void 
View::DoListenToTclCommand ( char* iCommand, int iArgc, char** iArgv ) {

}

void
View::DoDraw() {

  DebugOutput( << "View " << msLabel << ": DoDraw()" );
}

void
View::DoReshape( int iWidth, int iHeight ) {

  DebugOutput( << "View " << msLabel << ": DoReshape()" );
}

void
View::DoTimer() {

  DebugOutput( << "View " << msLabel << ": DoTimer()" );
}

void
View::DoMouseMoved( int inX, int inY, InputState& iState ) {

  DebugOutput( << "View " << msLabel << ": DoMouseMoved()" );
}

void
View::DoMouseUp( int inX, int inY, InputState& iState ) {

  DebugOutput( << "View " << msLabel << ": DoMouseUp()" );
}

void
View::DoMouseDown( int inX, int inY, InputState& iState ) {

  DebugOutput( << "View " << msLabel << ": DoMouseDown()" );
}

void
View::DoKeyDown( int inX, int inY, InputState& iState ) {

  DebugOutput( << "View " << msLabel << ": DoKeyDown()" );
}

void
View::DoKeyUp( int inX, int inY, InputState& iState ) {

  DebugOutput( << "View " << msLabel << ": DoKeyUp()" );
}

