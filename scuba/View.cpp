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
View::MouseMoved( int iWindow[2], InputState& iInput, ScubaToolState& iTool ) {

  this->DoMouseMoved( iWindow, iInput, iTool );
}

void
View::MouseUp( int iWindow[2], InputState& iInput, ScubaToolState& iTool ) {

  this->DoMouseUp( iWindow, iInput, iTool );
}

void
View::MouseDown( int iWindow[2], InputState& iInput, ScubaToolState& iTool ) {

  this->DoMouseDown( iWindow, iInput, iTool );
}

void
View::KeyDown( int iWindow[2], InputState& iInput, ScubaToolState& iTool ) {

  this->DoKeyDown( iWindow, iInput, iTool );
}

void
View::KeyUp( int iWindow[2], InputState& iInput, ScubaToolState& iTool ) {

  this->DoKeyUp( iWindow, iInput, iTool );
}

TclCommandListener::TclCommandResult 
View::DoListenToTclCommand ( char* iCommand, int iArgc, char** iArgv ) {
  return ok;
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
View::DoMouseMoved( int iWindow[2], 
		    InputState& iInput, ScubaToolState& iTool ) {

  DebugOutput( << "View " << msLabel << ": DoMouseMoved()" );
}

void
View::DoMouseUp( int iWindow[2], 
		 InputState& iInput, ScubaToolState& iTool ) {

  DebugOutput( << "View " << msLabel << ": DoMouseUp()" );
}

void
View::DoMouseDown( int iWindow[2], 
		   InputState& iInput, ScubaToolState& iTool ) {

  DebugOutput( << "View " << msLabel << ": DoMouseDown()" );
}

void
View::DoKeyDown( int iWindow[2], 
		 InputState& iInput, ScubaToolState& iTool ) {

  DebugOutput( << "View " << msLabel << ": DoKeyDown()" );
}

void
View::DoKeyUp( int iWindow[2], 
	       InputState& iInput, ScubaToolState& iTool ) {

  DebugOutput( << "View " << msLabel << ": DoKeyUp()" );
}

