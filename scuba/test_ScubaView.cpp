#include <stdexcept>
#include <sstream>
extern "C" {
#include "glut.h"
#include "tcl.h"
}
#include "ScubaView.h"
#include "ScubaFrame.h"
#include "ToglManager.h"

using namespace std;

#define Assert(x,s)   \
  if(!(x)) { \
  stringstream sError; \
  sError << "Line " << __LINE__ << ": " << s; \
  cerr << sError.str().c_str() << endl; \
  throw logic_error( sError.str() ); \
  }




// Togl tester ---------------------------------------------------------

extern "C" {
int Test_scubaview_Init ( Tcl_Interp* iInterp ) {

  ToglManager& toglMgr = ToglManager::GetManager();

  try {
    toglMgr.InitializeTogl( iInterp );
    toglMgr.SetFrameFactory( new ScubaFrameFactory );
    ScubaFrame::SetViewFactory( new ScubaViewFactory );

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



class TestLayer : public Layer {

public:
  TestLayer();
  virtual void DrawIntoBuffer( GLbyte* iBuffer, ViewState& iViewState );
  virtual void GetInfoAtRAS ( float inX, float inY, float inZ,
			    std::map<std::string,std::string>& iLabelValues );
  bool WasDrawn() const { return mbWasDrawn; }
protected:
  bool mbWasDrawn;
};

TestLayer::TestLayer() {
  mbWasDrawn = false;
}

void
TestLayer::DrawIntoBuffer( GLbyte* iBuffer, ViewState& iViewState ) {
  mbWasDrawn = true;
}

void 
TestLayer::GetInfoAtRAS ( float inX, float inY, float inZ,
			  std::map<std::string,std::string>& iLabelValues ) {

  string sLabel;
  if( msLabel != "" ) {
    sLabel = msLabel;
  } else {
    stringstream ssLabel;
    ssLabel << mID;
    sLabel = ssLabel.str();
  }

  iLabelValues[sLabel] = "Hello world";
}




class ScubaViewTester {
public:
  void Test( Tcl_Interp* iInterp );
  void TestCoords ( ScubaView& iView, 
		    int iWidth, int iHeight,
		    float iZoomLevel,
		    float iCenterRASX, float iCenterRASY, float iCenterRASZ,
		    ViewState::Plane iInPlane,
		    int iXWindow, int iYWindow,
		    float iRASX, float iRASY, float iRASZ );
};

void 
ScubaViewTester::TestCoords( ScubaView& iView, 
			     int iWidth, int iHeight,
			     float iZoomLevel,
			     float iCenterRASX, float iCenterRASY,
			     float iCenterRASZ,
			     ViewState::Plane iInPlane,
			     int iXWindow, int iYWindow,
			     float iRASX, float iRASY, float iRASZ ) {

    iView.Reshape( iWidth, iHeight );
    iView.Set2DZoomLevel( iZoomLevel );
    float center[3];
    center[0] = iCenterRASX; 
    center[1] = iCenterRASY; 
    center[2] = iCenterRASZ;
    iView.Set2DRASCenter( center );
    iView.Set2DInPlane( iInPlane );
    float rasX, rasY, rasZ;
    iView.TranslateWindowToRAS( iXWindow, iYWindow, rasX, rasY, rasZ );
    if( rasX != iRASX || rasY != iRASY || rasZ != iRASZ ) {
      cerr << "translate didn't work: " << iView.mViewState << ", ras: "
	   << rasX << ", " << rasY << ", " << rasZ << endl;
      throw 0;
    }
}

void 
ScubaViewTester::Test( Tcl_Interp* iInterp ) {

  stringstream sError;

  try {

    ScubaView view;
    
    // Set our view state stuff and check it.
    float center[3] = { 5.0, 5.1, 5.2 };
    view.Set2DRASCenter( center );
    Assert( (view.mViewState.mCenterRAS[0] == center[0] &&
	     view.mViewState.mCenterRAS[1] == center[1] &&
	     view.mViewState.mCenterRAS[2] == center[2]), 
	    "Set2DRASCenter failed" );

    float zoomLevel = 5.0;
    view.Set2DZoomLevel( zoomLevel );
    Assert( (view.mViewState.mZoomLevel == 5.0),
	    "SetZoomLevel failed" );

    ViewState::Plane plane = ViewState::X;
    view.Set2DInPlane( plane );
    Assert( (view.mViewState.mInPlane == ViewState::X),
	    "Set2DInPlane failed" );

    // Add a bunch of layers to different levels and make sure they
    // match our own expectations.
    int const kcLayers = 10;
    TestLayer layer[kcLayers];
    int nLevel = 0;
    map<int,int> levelLayerID;
    for( int nLayer = 0; nLayer < kcLayers; nLayer++ ) {
      int layerID = layer[nLayer].GetID();
      levelLayerID[nLevel] = layerID;
      view.AddLayer( layerID, nLevel );
      nLevel++;
    }

    bool bFailed = false;
    map<int,int>::iterator tLevelLayerID;
    for( tLevelLayerID = view.mLevelLayerIDMap.begin(); 
	 tLevelLayerID != view.mLevelLayerIDMap.end(); ++tLevelLayerID ) {
      int level = (*tLevelLayerID).first;
      int layerID = (*tLevelLayerID).second;
      if( levelLayerID[level] != layerID ) {
	cerr << "ID at level " << level << " incorrect" << endl << ", was "
	     << layerID << ", should be " << levelLayerID[level] << endl;
	bFailed = true;
      }
    }
    Assert( (!bFailed), "Adding IDs failed" );


    // Make sure reshaping works and creates our buffer. Make a few to
    // make sure they dispose of the old one properly. Make one with a
    // negative size to check the error handling.
    view.Reshape( 200, 200 );
    Assert( (NULL != view.mBuffer), "buffer not allocated on reshape" );
    view.Reshape( 200, 200 );
    Assert( (NULL != view.mBuffer), "buffer not allocated on second reshape" );
    try {
      view.Reshape( -1, 0 );
      sError << "error not thrown when reshaping with -1, 0" << endl;
      throw  0;
    }
    catch(...) {}
    try {
      view.Reshape( -1, -1 );
      sError << "error not thrown when reshaping with -1, -1" << endl;
      throw  0;
    }
    catch(...) {}
    view.Reshape( 200, 200 );


    // Draw the view, all our layers should be drawn.
    bFailed = false;
    view.DoDraw();
    for( int nLayer = 0; nLayer < kcLayers; nLayer++ ) {
      if( !layer[nLayer].WasDrawn() ) {
	cerr << "layer ID " << layer[nLayer].GetID() << " not drawn" << endl;
	bFailed = true;
      }
    }
    Assert( (!bFailed), "Drawing failed" );


    // Set labels in all the layers. Call GetInfoAtRAS and check the
    // labelvalue map to make sure they got copied properly.
    map<string,string> labelValueMap;
    bFailed = false;
    for( int nLayer = 0; nLayer < kcLayers; nLayer++ ) {
      stringstream ssLabel;
      ssLabel << nLayer;
      layer[nLayer].SetLabel( ssLabel.str() );
      layer[nLayer].GetInfoAtRAS( 0, 0, 0, labelValueMap );
      if( labelValueMap[ssLabel.str()] != "Hello world" ) {
	cerr << "layer ID " << layer[nLayer].GetID() << " not copied" << endl;
	bFailed = true;
      }
    }
    Assert( (!bFailed), "GetInfoAtRAS failed" );



    // Test the coordinate stuff. At zoom level 1, with size 200x200,
    // inplane z, RASCenter 100,100,100, window coords h,v should
    // equal RAS coords x,y,100 where x==h and y==v.
    TestCoords( view,
		200, 200,   1,   100, 100, 100,   ViewState::Z,
		0, 0,       0, 0, 100 );
    TestCoords( view,
		200, 200,   1,   100, 100, 100,   ViewState::Z,
		199, 199,       199, 199, 100 );

    TestCoords( view,
		200, 200,   2,   100, 100, 100,   ViewState::Z,
		0, 0,       50, 50, 100 );
    TestCoords( view,
		200, 200,   2,   100, 100, 100,   ViewState::Z,
		199, 199,       149.5, 149.5, 100 );

    TestCoords( view,
		100, 100,   1,   0, 0, 0,   ViewState::Z,
		0, 0,       -50, -50, 0 );
    TestCoords( view,
		100, 100,   1,   0, 0, 0,   ViewState::Z,
		99, 99,       49, 49, 0 );

    TestCoords( view,
		200, 200,   1,   100, 100, 100,   ViewState::X,
		0, 0,       100, 0, 0 );
    TestCoords( view,
		200, 200,   1,   100, 100, 100,   ViewState::Y,
		0, 0,       0, 100, 0 );
    TestCoords( view,
		200, 200,   1,   100, 100, 100,   ViewState::Z,
		0, 0,       0, 0, 100 );


    // Test the tcl commands.
    char sCommand[1024];
    int rTcl;

    sprintf( sCommand, "SetViewRASCenter %d 1.2 3.4 5.6", view.GetID() );
    rTcl = Tcl_Eval( iInterp, sCommand );
    if( TCL_OK != rTcl ) {
      sError << "Tcl_Eval returned not TCL_OK: " << endl 
	     << "Command: " << sCommand << endl
	     << "Result: " << iInterp->result;
      throw 0;
    }
    Assert( (view.mViewState.mCenterRAS[0] == (float)1.2 &&
	     view.mViewState.mCenterRAS[1] == (float)3.4 &&
	     view.mViewState.mCenterRAS[2] == (float)5.6), 
	    "SetViewRASCenter didn't set properly" );

    sprintf( sCommand, "SetViewZoomLevel %d 3.4", view.GetID() );
    rTcl = Tcl_Eval( iInterp, sCommand );
    if( TCL_OK != rTcl ) {
      sError << "Tcl_Eval returned not TCL_OK: " << endl 
	     << "Command: " << sCommand << endl
	     << "Result: " << iInterp->result;
      throw 0;
    }
    Assert( (view.mViewState.mZoomLevel == (float)3.4),
	    "SetViewZoomLevel didn't set properly" );

    sprintf( sCommand, "SetViewInPlane %d X", view.GetID() );
    rTcl = Tcl_Eval( iInterp, sCommand );
    if( TCL_OK != rTcl ) {
      sError << "Tcl_Eval returned not TCL_OK: " << endl 
	     << "Command: " << sCommand << endl
	     << "Result: " << iInterp->result;
      throw 0;
    }
    Assert( (view.mViewState.mInPlane == ViewState::X),
	    "SetViewInPlane didn't set properly" );

  }
  catch(...) {
    throw logic_error(sError.str());
  }
}

// -----------------------------------------------------------------------


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

    for( int nTrial = 0; nTrial < 50; nTrial++ ) {
      ScubaViewTester tester0;
      tester0.Test( interp );
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
