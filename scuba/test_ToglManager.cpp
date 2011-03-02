/**
 * @file  test_ToglManager.cpp
 * @brief test routines
 *
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:39 $
 *    $Revision: 1.14 $
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


#include <stdexcept>
#include "ToglManager.h"
#include "Scuba-impl.h"

extern "C" {
#include "glut.h"
}

using namespace std;

const char* Progname = "test_ToglManager";

class TestFrame : public WindowFrame {
public:
  TestFrame( int iID );
  virtual ~TestFrame();
protected:
  virtual void DoDraw();
  virtual void DoReshape();
  virtual void DoTimer();
  virtual void DoMouseMoved( int iWindow[2], InputState& iModifier );
  virtual void DoMouseUp( int iWindow[2], InputState& iInput );
  virtual void DoMouseDown( int iWindow[2], InputState& iInput );
  virtual void DoKeyDown( int iWindow[2], InputState& iInput );
  virtual void DoKeyUp( int iWindow[2], InputState& iInput );

  bool bTimerCalled;
};


TestFrame::TestFrame( int iID ) : WindowFrame( ) {
  SetOutputStreamToCerr();
  DebugOutput( << "Created TestFrame " << iID );
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
  for ( int nChar = 0; nChar < (int)strlen(sID); nChar++ ) {
    glutBitmapCharacter( GLUT_BITMAP_HELVETICA_18, sID[nChar] );
  }
}

void
TestFrame::DoReshape() {}

void
TestFrame::DoTimer() {
  bTimerCalled = true;
}

void
TestFrame::DoMouseMoved( int iWindow[2], InputState& iInput ) {

  DebugOutput( << "TestFrame " << mID << ": DoMouseMoved "
               << iWindow[0] << ", " << iWindow[1] << ", " << iInput );
}

void
TestFrame::DoMouseUp( int iWindow[2], InputState& iInput ) {

  DebugOutput( << "TestFrame " << mID << ": DoMouseUp "
               << iWindow[0] << ", " << iWindow[1] << ", " << iInput );
}

void
TestFrame::DoMouseDown( int iWindow[2], InputState& iInput ) {

  DebugOutput( << "TestFrame " << mID << ": DoMouseDown "
               << iWindow[0] << ", " << iWindow[1] << ", " << iInput );
}

void
TestFrame::DoKeyDown( int iWindow[2], InputState& iInput ) {

  DebugOutput( << "TestFrame " << mID << ": DoKeyDown "
               << iWindow[0] << ", " << iWindow[1] << ", " << iInput );
}

void
TestFrame::DoKeyUp( int iWindow[2], InputState& iInput ) {

  DebugOutput( << "TestFrame " << mID << ": DoKeyUp "
               << iWindow[0] << ", " << iWindow[1] << ", " << iInput );
}


class TestFrameFactory : public WindowFrameFactory {
public:
  virtual WindowFrame* NewWindowFrame( int iID ) {
    return new TestFrame( iID );
  }
};


#if BUILD_TCL_TEST
extern "C" {
  int Test_toglmanager_Init ( Tcl_Interp* iInterp ) {

    ToglManager& toglMgr = ToglManager::GetManager();

    try {
      toglMgr.InitializeTogl( iInterp );
      toglMgr.SetFrameFactory( new TestFrameFactory );
    } catch ( ... ) {
      return TCL_ERROR;
    }

    return TCL_OK;
  }
}
#endif

class InputStateTester {
public:
  void Test();
};

void
InputStateTester::Test() {

  try {

    InputState m;
    if ( m.IsShiftKeyDown() || m.IsAltKeyDown() || m.IsControlKeyDown() ) {
      cerr << "InputState not init'd properly" << endl;
      throw logic_error( "InputState not init'd properly" );
    }

    m.mbShiftKey = true;
    if ( !m.IsShiftKeyDown() || m.IsAltKeyDown() || m.IsControlKeyDown() ) {
      cerr << "InputState shift not set properly" << endl;
      throw logic_error( "InputState shift not set properly" );
    }

    m.mbAltKey = true;
    if ( !m.IsShiftKeyDown() || !m.IsAltKeyDown() || m.IsControlKeyDown() ) {
      cerr << "InputState alt not set properly" << endl;
      throw logic_error( "InputState alt not set properly" );
    }

    m.mbControlKey = true;
    if ( !m.IsShiftKeyDown() || !m.IsAltKeyDown() || !m.IsControlKeyDown() ) {
      cerr << "InputState control not set properly" << endl;
      throw logic_error( "InputState control not set properly" );
    }

  } catch ( exception& e ) {
    cerr << "InputStateTester failed with exception: " << e.what() << endl;
    exit( 1 );
  } catch (...) {
    cerr << "InputStateTester failed" << endl;
    exit( 1 );
  }

}

int main( int argc, char** argv ) {

  cerr << "Beginning test" << endl;

  try {
    for ( int nTrial = 0; nTrial < 50; nTrial++ ) {

      InputStateTester mt0;
      mt0.Test();

      InputStateTester mt1;
      mt1.Test();

    }
  } catch ( exception& e ) {
    cerr << "failed with exception: " << e.what() << endl;
    exit( 1 );
  } catch (...) {
    cerr << "failed" << endl;
    exit( 1 );
  }

  cerr << "Success" << endl;

  exit( 0 );
}
