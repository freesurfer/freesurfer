/**
 * @file  test_ScubaFrame.cpp
 * @brief test ScubaFrame class
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
#include <sstream>
#include "ScubaFrame.h"
extern "C" {
#include "glut.h"
}
#include "Scuba-impl.h"

const char* Progname = "test_ScubaFrame";

using namespace std;

#define Assert(x,s)   \
  if(!(x)) { \
  stringstream sError; \
  sError << "Line " << __LINE__ << ": " << s; \
  cerr << sError.str().c_str() << endl; \
  throw logic_error( sError.str() ); \
  }

#define AssertTclOK(x) \
    if( TCL_OK != (x) ) { \
      sError << "Tcl_Eval returned not TCL_OK: " << endl  \
      << "Command: " << sCommand << endl \
      << "Result: " << iInterp->result; \
      cerr << sError.str().c_str() << endl; \
      throw logic_error( sError.str() ); \
    } \


// Custom View implementation. ----------------------------------------

class TestView : public View {
public:
  TestView();
protected:
  virtual void DoDraw();
  virtual void DoReshape( int iWidth, int iHeight );
  virtual void DoTimer();
  virtual void DoMouseMoved( int inX, int inY, InputState& iInput );
  virtual void DoMouseUp( int inX, int inY, InputState& iInput );
  virtual void DoMouseDown( int inX, int inY, InputState& iInput );
  virtual void DoKeyDown( int inX, int inY, InputState& iInput );
  virtual void DoKeyUp( int inX, int inY, InputState& iInput );
};

TestView::TestView() {}

void
TestView::DoDraw() {

  stringstream ssLabel;
  ssLabel << mID << ": " << msLabel;
  string sLabel = ssLabel.str();

  glColor3i( 0, 0, 0 );
  glBegin( GL_POLYGON );
  glVertex2d( 0, 0 );
  glVertex2d( 0, mWidth );
  glVertex2d( mHeight, mWidth );
  glVertex2d( mHeight, 0 );
  glVertex2d( 0, 0 );
  glEnd();

  glColor3f( 1, 1, 1 );

  glRasterPos2i( mWidth / 2, mHeight / 2 - sLabel.length()/2);

  glColor3f( 1, 0, 1 );
  for ( int nChar = 0; nChar < (int)sLabel.length(); nChar++ ) {
    glutBitmapCharacter( GLUT_BITMAP_HELVETICA_18, sLabel[nChar] );
  }

  glBegin( GL_LINE_STRIP );
  glVertex2d( 4, 4 );
  glVertex2d( mWidth-4, 4 );
  glVertex2d( mWidth-4, mHeight-4 );
  glVertex2d( 4, mHeight-4 );
  glVertex2d( 4, 4 );
  glEnd ();
}

void
TestView::DoReshape( int iWidth, int iHeight ) {}

void
TestView::DoTimer() {}

void
TestView::DoMouseMoved( int inX, int inY, InputState& iInput ) {}

void
TestView::DoMouseUp( int inX, int inY, InputState& iInput ) {}

void
TestView::DoMouseDown( int inX, int inY, InputState& iInput ) {
  cerr << msLabel << ": click " << endl;
}

void
TestView::DoKeyDown( int inX, int inY, InputState& iInput ) {}

void
TestView::DoKeyUp( int inX, int inY, InputState& iInput ) {}

class TestViewFactory : public ViewFactory {
public:
  virtual View* NewView() {
    return new TestView();
  }
};

// ----------------------------------------------------------------------



// Togl tester ---------------------------------------------------------

#if BUILD_TCL_TEST
extern "C" {
  int Test_scubaframe_Init ( Tcl_Interp* iInterp ) {

    ToglManager& toglMgr = ToglManager::GetManager();

    try {
      toglMgr.InitializeTogl( iInterp );
      toglMgr.SetFrameFactory( new ScubaFrameFactory );
      ScubaFrame::SetViewFactory( new TestViewFactory );

      TclCommandManager& commandMgr = TclCommandManager::GetManager();
      commandMgr.SetOutputStreamToCerr();
      commandMgr.Start( iInterp );
    } catch ( ... ) {
      return TCL_ERROR;
    }

    return TCL_OK;
  }
}
#endif


// Non-togl tester --------------------------------------------------------


class ScubaFrameTester {
public:
  void Test( Tcl_Interp* iInterp );
};

void
ScubaFrameTester::Test( Tcl_Interp* iInterp ) {

  stringstream sError;

  try {

    sError.flush();
    sError << "Setting view factory" << endl;
    ScubaFrame::SetViewFactory( new TestViewFactory );

    sError.flush();
    sError << "Creating frame 0" << endl;
    ScubaFrame* frame = new ScubaFrame( );
    Assert( (NULL != frame), "frame was not created" );

    sError.flush();
    sError << "Setting view config to c22" << endl;
    frame->SetViewConfiguration( ScubaFrame::c22 );

    View* view = NULL;
    for ( int nRow = 0; nRow < 2; nRow++ ) {
      for ( int nCol = 0; nCol < 2; nCol++ ) {
        sError.flush();
        sError << "View " << nCol << ", " << nRow << " was NULL.";
        view = frame->GetViewAtColRow( nCol, nRow );
      }
    }

    // Check the tcl functions that return number of cols and rows.
    char sCommand[1024];
    int rTcl;
    sprintf( sCommand, "GetNumberOfRowsInFrame %d", frame->GetID() );
    rTcl = Tcl_Eval( iInterp, sCommand );
    AssertTclOK( rTcl );
    const char* sTclResult = Tcl_GetStringResult( iInterp );
    string scRows( sTclResult );
    Assert( (scRows == "2"), "tcl function returned wrong number of rows" );
    for ( int nRow = 0; nRow < 2; nRow++ ) {
      sprintf( sCommand, "GetNumberOfColsAtRowInFrame %d %d",
               frame->GetID(), nRow );
      rTcl = Tcl_Eval( iInterp, sCommand );
      AssertTclOK( rTcl );
      const char* sTclResult = Tcl_GetStringResult( iInterp );
      string scCols( sTclResult );
      Assert( (scCols == "2"), "tcl function returned wrong number of cols" );
    }

    sError.flush();
    sError << "deleting frame" << endl;
    delete frame;
  } catch (...) {
    throw logic_error(sError.str());
  }
}


int main( int argc, char** argv ) {

  cerr << "Beginning test" << endl;

  try {
    Tcl_Interp* interp = Tcl_CreateInterp();
    Assert( interp, "Tcl_CreateInterp returned null" );

    int rTcl = Tcl_Init( interp );
    Assert( TCL_OK == rTcl, "Tcl_Init returned not TCL_OK" );

    TclCommandManager& commandMgr = TclCommandManager::GetManager();
    commandMgr.SetOutputStreamToCerr();
    commandMgr.Start( interp );

    for ( int nTrial = 0; nTrial < 50; nTrial++ ) {
      ScubaFrameTester tester0;
      tester0.Test( interp );

      ScubaFrameTester tester1;
      tester1.Test( interp );

      ScubaFrameTester tester2;
      tester2.Test( interp );
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
