#include <stdexcept>
#include "ScubaFrame.h"
extern "C" {
#include "glut.h"
}

using namespace std;

#define Assert(x,s)   \
  if(!(x)) { \
  stringstream sError; \
  sError << "Line " << __LINE__ << ": " << s; \
  cerr << sError.str().c_str() << endl; \
  throw logic_error( sError.str() ); \
  }


// Custom View implementation. ----------------------------------------

class TestView : public View {
public:
  TestView();
protected:
  virtual void DoDraw();
  virtual void DoReshape( int iWidth, int iHeight );
  virtual void DoTimer();
  virtual void DoMouseMoved( int inX, int inY, int iButton, int iModifiers );
  virtual void DoMouseUp( int inX, int inY, int iButton, int iModifers );
  virtual void DoMouseDown( int inX, int inY, int iButton, int iModifers );
  virtual void DoKeyDown( int inX, int inY, std::string isKey, int iModifers );
  virtual void DoKeyUp( int inX, int inY, std::string isKey, int iModifers );
};  

TestView::TestView() {
}

void
TestView::DoDraw() {
 
  if( msLabel != "" ) {
    glRasterPos2i( mWidth / 2, mHeight / 2 - msLabel.length()/2);
    
    glColor3f( 1, 0, 1 );
    for( int nChar = 0; nChar < msLabel.length(); nChar++ ) {
      glutBitmapCharacter( GLUT_BITMAP_HELVETICA_18, msLabel[nChar] );
    }
  }
}

void
TestView::DoReshape( int iWidth, int iHeight ) {
}

void
TestView::DoTimer() {

}

void
TestView::DoMouseMoved( int inX, int inY, int iButton, int iModifiers ) {
}

void
TestView::DoMouseUp( int inX, int inY, int iButton, int iModifers ) {

}

void
TestView::DoMouseDown( int inX, int inY, int iButton, int iModifers ) {
  cerr << msLabel << ": click " << endl;
}

void
TestView::DoKeyDown( int inX, int inY, std::string isKey, int iModifers ) {

}

void
TestView::DoKeyUp( int inX, int inY, std::string isKey, int iModifers ) {

}

class TestViewFactory : public ViewFactory {
public:
  virtual View* NewView() { 
    return new TestView();
  }
};

// ----------------------------------------------------------------------



// Togl tester ---------------------------------------------------------

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
  }
  catch( ... ) {
    return TCL_ERROR;
  }

  return TCL_OK;
}
}



// Non-togl tester --------------------------------------------------------


class ScubaFrameTester {
public:
  void Test();
};

void 
ScubaFrameTester::Test() {

  stringstream sError;

  try {
    
    sError.flush();
    sError << "Setting view factory" << endl;
    ScubaFrame::SetViewFactory( new TestViewFactory );
    
    sError.flush();
    sError << "Creating frame 0" << endl;
    ScubaFrame* frame = new ScubaFrame( 0 );
    Assert( (NULL != frame), "frame was not created" );
    
    sError.flush();
    sError << "Setting view config to c22" << endl;
    frame->SetViewConfiguration( ScubaFrame::c22 );
    
    View* view = NULL;
    for( int nRow = 0; nRow < 2; nRow++ ) {
      for( int nCol = 0; nCol < 2; nCol++ ) {
	sError.flush();
	sError << "View " << nCol << ", " << nRow << " was NULL.";
	view = frame->GetViewAtColRow( nCol, nRow );
      }
    }
    
    sError.flush();
    sError << "deleting frame" << endl;
    delete frame;
  }
  catch(...) {
    throw logic_error(sError.str());
  }
}


int main( int argc, char** argv ) {

  cerr << "Beginning test" << endl;
 
  try {
    for( int nTrial = 0; nTrial < 50; nTrial++ ) {
      ScubaFrameTester tester0;
      tester0.Test();

      ScubaFrameTester tester1;
      tester1.Test();

      ScubaFrameTester tester2;
      tester2.Test();
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
