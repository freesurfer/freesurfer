#include "glut.h"
#include "View.h"

using namespace std;

View::View( string isID ) : FrameView( isID ) {
}

void
View::Draw() {
 
 
  glRasterPos2i( mWidth / 2, mHeight / 2 - msID.length()/2);
		     
  glColor3f( 1, 0, 1 );
  for( int nChar = 0; nChar < msID.length(); nChar++ ) {
    glutBitmapCharacter( GLUT_BITMAP_HELVETICA_18, msID[nChar] );
  }
}

void
View::Reshape() {
}

void
View::MouseMoved( int inX, int inY, int iButton, int iModifiers ) {
}

void
View::MouseUp( int inX, int inY, int iButton, int iModifers ) {

}

void
View::MouseDown( int inX, int inY, int iButton, int iModifers ) {
  cerr << msID << ": click " << endl;
}

void
View::KeyDown( int inX, int inY, std::string isKey, int iModifers ) {

}

void
View::KeyUp( int inX, int inY, std::string isKey, int iModifers ) {

}

