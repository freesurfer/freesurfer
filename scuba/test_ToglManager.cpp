#include <stdexcept>
#include "ToglManager.h"

extern "C" {
#include "glut.h"
}

using namespace std;


class TestFrame : public ToglFrame {
public:
  TestFrame( ToglFrame::ID iID );
  virtual ~TestFrame();
protected:
  virtual void DoDraw();
  virtual void DoReshape();
  virtual void DoTimer();
  virtual void DoMouseMoved( int inX, int inY, int iButton, int iModifiers );
  virtual void DoMouseUp( int inX, int inY, int iButton, int iModifers );
  virtual void DoMouseDown( int inX, int inY, int iButton, int iModifers );
  virtual void DoKeyDown( int inX, int inY, string isKey, int iModifers );
  virtual void DoKeyUp( int inX, int inY, string isKey, int iModifers );

  bool bTimerCalled;

  static int const kcBytesPerPixel;
  GLubyte* mFrameBuffer;
};

int const TestFrame::kcBytesPerPixel = 4;

TestFrame::TestFrame( ToglFrame::ID iID ) : ToglFrame( iID ) {
  DebugOutput( << "Created TestFrame " << iID );
  SetOutputStreamToCerr();
  bTimerCalled = false;
  mFrameBuffer = NULL;
}

TestFrame::~TestFrame() {
  DebugOutput( << "Destroyed TestFrame " << mID );
}

void
TestFrame::DoDraw() {

  char sID[10];
  sprintf( sID, "%d", mID );
  glRasterPos2i( mWidth/2, mHeight/2 );
  glClearColor( 0, 0, 0, 1 );
  glClear( GL_COLOR_BUFFER_BIT );
  for( int nChar = 0; nChar < strlen(sID); nChar++ ) {
    glutBitmapCharacter( GLUT_BITMAP_HELVETICA_18, sID[nChar] );
  }
}

void
TestFrame::DoReshape() {

  if( NULL != mFrameBuffer ) {
    free( mFrameBuffer );
  }

  mFrameBuffer = (GLubyte*) malloc( mWidth * mHeight * kcBytesPerPixel );
  if( NULL == mFrameBuffer ) {
    DebugOutput( << "Allocation of buffer with width " << mWidth
		 << " and height " << mHeight << " failed" );
    throw logic_error( "Couldn't allocate buffer" );
  }

  for( int nPixel = 0; 
       nPixel < mWidth * mHeight * kcBytesPerPixel; 
       nPixel += kcBytesPerPixel ) {
    mFrameBuffer[nPixel] = (GLubyte) 0;
    mFrameBuffer[nPixel+1] = (GLubyte) 0;
    mFrameBuffer[nPixel+2] = (GLubyte) 0;
    mFrameBuffer[nPixel+3] = (GLubyte) 1;
  }
}

void
TestFrame::DoTimer() {
  bTimerCalled = true;
}

void
TestFrame::DoMouseMoved( int inX, int inY, int iButton, int iModifiers ) {

  DebugOutput( << "TestFrame " << mID << ": DoMouseMoved " 
	       << inX << ", " << inY << ", " << iButton << ", " << iModifiers);
}

void
TestFrame::DoMouseUp( int inX, int inY, int iButton, int iModifiers ) {

  DebugOutput( << "TestFrame " << mID << ": DoMouseUp "
	       << inX << ", " << inY << ", " << iButton << ", " << iModifiers);
}

void
TestFrame::DoMouseDown( int inX, int inY, int iButton, int iModifiers ) {

  DebugOutput( << "TestFrame " << mID << ": DoMouseDown "
	       << inX << ", " << inY << ", " << iButton << ", " << iModifiers);
}

void
TestFrame::DoKeyDown( int inX, int inY, string isKey, int iModifiers ) {

  DebugOutput( << "TestFrame " << mID << ": DoKeyDown "
	       << inX << ", " << inY << ", " << isKey << ", " << iModifiers);
}

void
TestFrame::DoKeyUp( int inX, int inY, string isKey, int iModifiers ) {

  DebugOutput( << "TestFrame " << mID << ": DoKeyUp "
	       << inX << ", " << inY << ", " << isKey << ", " << iModifiers);
}


class TestFrameFactory : public ToglFrameFactory {
public:
  virtual ToglFrame* NewToglFrame( ToglFrame::ID iID ) { 
    return new TestFrame( iID );
  }
};


extern "C" {
int Test_toglmanager_Init ( Tcl_Interp* iInterp ) {

  ToglManager& toglMgr = ToglManager::GetManager();

  try {
    toglMgr.InitializeTogl( iInterp );
    toglMgr.SetFrameFactory( new TestFrameFactory );
  }
  catch( ... ) {
    return TCL_ERROR;
  }

  return TCL_OK;
}
}
