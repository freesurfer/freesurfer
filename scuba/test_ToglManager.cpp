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
  virtual void DoMouseMoved( int inX, int inY, InputState& iModifier );
  virtual void DoMouseUp( int inX, int inY, InputState& iState );
  virtual void DoMouseDown( int inX, int inY, InputState& iState );
  virtual void DoKeyDown( int inX, int inY, InputState& iState );
  virtual void DoKeyUp( int inX, int inY, InputState& iState );

  bool bTimerCalled;
};


TestFrame::TestFrame( ToglFrame::ID iID ) : ToglFrame( iID ) {
  DebugOutput( << "Created TestFrame " << iID );
  SetOutputStreamToCerr();
  bTimerCalled = false;
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
}

void
TestFrame::DoTimer() {
  bTimerCalled = true;
}

void
TestFrame::DoMouseMoved( int inX, int inY, InputState& iState ) {

  DebugOutput( << "TestFrame " << mID << ": DoMouseMoved " 
	       << inX << ", " << inY << ", " << iState );
}

void
TestFrame::DoMouseUp( int inX, int inY, InputState& iState ) {

  DebugOutput( << "TestFrame " << mID << ": DoMouseUp "
	       << inX << ", " << inY << ", " << iState );
}

void
TestFrame::DoMouseDown( int inX, int inY, InputState& iState ) {

  DebugOutput( << "TestFrame " << mID << ": DoMouseDown "
	       << inX << ", " << inY << ", " << iState );
}

void
TestFrame::DoKeyDown( int inX, int inY, InputState& iState ) {

  DebugOutput( << "TestFrame " << mID << ": DoKeyDown "
	       << inX << ", " << inY << ", " << iState );
}

void
TestFrame::DoKeyUp( int inX, int inY, InputState& iState ) {

  DebugOutput( << "TestFrame " << mID << ": DoKeyUp "
	       << inX << ", " << inY << ", " << iState );
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


class InputStateTester {
public:
  void Test();
};

void 
InputStateTester::Test() {
  
  try {

    InputState m;
    if( m.IsShiftKeyDown() || m.IsAltKeyDown() || m.IsControlKeyDown() ) {
      cerr << "InputState not init'd properly" << endl;
      throw new logic_error( "InputState not init'd properly" );
    }
    
    m.mbShiftKey = true;
    if( !m.IsShiftKeyDown() || m.IsAltKeyDown() || m.IsControlKeyDown() ) {
      cerr << "InputState shift not set properly" << endl;
      throw new logic_error( "InputState shift not set properly" );
    }

    m.mbAltKey = true;
    if( !m.IsShiftKeyDown() || !m.IsAltKeyDown() || m.IsControlKeyDown() ) {
      cerr << "InputState alt not set properly" << endl;
      throw new logic_error( "InputState alt not set properly" );
    }

    m.mbControlKey = true;
    if( !m.IsShiftKeyDown() || !m.IsAltKeyDown() || !m.IsControlKeyDown() ) {
      cerr << "InputState control not set properly" << endl;
      throw new logic_error( "InputState control not set properly" );
    }

  }
  catch( exception e ) {
    cerr << "InputStateTester failed with exception: " << e.what() << endl;
    exit( 1 );
  }
  catch(...) {
    cerr << "InputStateTester failed" << endl;
    exit( 1 );
  }
  
}

int main( int argc, char** argv ) {

  cerr << "Beginning test" << endl;
 
  try {
    for( int nTrial = 0; nTrial < 50; nTrial++ ) {

      InputStateTester mt0;
      mt0.Test();

      InputStateTester mt1;
      mt1.Test();

    }
  }
  catch( exception e ) {
    cerr << "failed with exception: " << e.what() << endl;
    exit( 1 );
  }
  catch(...) {
    cerr << "failed" << endl;
    exit( 1 );
  }

  cerr << "Success" << endl;

  exit( 0 );
}
