#include "View.h"

using namespace std;


template IDTracker<View>;
int IDTracker<View>::mNextID = 0;
map<int,View*> IDTracker<View>::mIDMap;


View::View() {
  mWidth = 0;
  mHeight = 0;
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
View::MouseMoved( int inX, int inY, int iButton, int iModifiers ) {

  this->DoMouseMoved( inX, inY, iButton, iModifiers );
}

void
View::MouseUp( int inX, int inY, int iButton, int iModifers ) {

  this->DoMouseUp( inX, inY, iButton, iModifers );
}

void
View::MouseDown( int inX, int inY, int iButton, int iModifers ) {

  this->DoMouseDown( inX, inY, iButton, iModifers );
}

void
View::KeyDown( int inX, int inY, string isKey, int iModifers ) {

  this->DoKeyDown( inX, inY, isKey, iModifers );
}

void
View::KeyUp( int inX, int inY, string isKey, int iModifers ) {

  this->DoKeyUp( inX, inY, isKey, iModifers );
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
View::DoMouseMoved( int inX, int inY, int iButton, int iModifiers ) {

  DebugOutput( << "View " << msLabel << ": DoMouseMoved()" );
}

void
View::DoMouseUp( int inX, int inY, int iButton, int iModifers ) {

  DebugOutput( << "View " << msLabel << ": DoMouseUp()" );
}

void
View::DoMouseDown( int inX, int inY, int iButton, int iModifers ) {

  DebugOutput( << "View " << msLabel << ": DoMouseDown()" );
}

void
View::DoKeyDown( int inX, int inY, string isKey, int iModifers ) {

  DebugOutput( << "View " << msLabel << ": DoKeyDown()" );
}

void
View::DoKeyUp( int inX, int inY, string isKey, int iModifers ) {

  DebugOutput( << "View " << msLabel << ": DoKeyUp()" );
}
