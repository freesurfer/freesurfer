#include "glut.h"
#include "View.h"

using namespace std;

View::View( string isID ) {
  msID = isID;
  mWidth = 0;
  mHeight = 0;
}

void
View::Draw() {
 
 
  glRasterPos2i( mWidth / 2, mHeight / 2 - msID.length()/2);
		     
  glColor3f( 1, 0, 1 );
  for( int nChar = 0; nChar < msID.length(); nChar++ ) {
    glutBitmapCharacter( GLUT_BITMAP_HELVETICA_18, msID[nChar] );
  }
}

