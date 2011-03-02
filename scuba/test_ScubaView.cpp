/**
 * @file  test_ScubaView.cpp
 * @brief test ScubaView class
 *
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:39 $
 *    $Revision: 1.16 $
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
extern "C" {
#include "glut.h"
#define USE_NON_CONST
#include "tcl.h"
#undef USE_NON_CONST
}
#include "ScubaView.h"
#include "ScubaFrame.h"
#include "ToglManager.h"
#include "Scuba-impl.h"

const char* Progname = "test_ScubaView";

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

#define fequal(x,y) (fabs((x)-(y))<0.0000001)

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
    } catch ( ... ) {
      return TCL_ERROR;
    }

    return TCL_OK;
  }
}



// Non-togl tester --------------------------------------------------------



class TestLayer : public Layer {

public:
  TestLayer();
  virtual void DrawIntoBuffer( GLubyte* iBuffer, int iWidth, int iHeight,
                               ViewState& iViewState,
                               ScubaWindowToRASTranslator& iTranslator );
  virtual void DrawIntoGL    ( ViewState& iViewState,
                               ScubaWindowToRASTranslator& iTranslator );
  virtual void GetInfoAtRAS ( float inX, float inY, float inZ,
                              std::map<std::string,std::string>& iLabelValues );
  virtual string GetTypeDescription() {
    return "TestLayer";
  }
  bool WasDrawn() const {
    return (mbBufferWasDrawn && mbGLWasDrawn);
  }
  int GetWidth() const {
    return mWidth;
  }
  int GetHeight() const {
    return mHeight;
  }
protected:
  bool mbBufferWasDrawn;
  bool mbGLWasDrawn;
};

TestLayer::TestLayer() {
  mbBufferWasDrawn = false;
  mbGLWasDrawn = false;
}

void
TestLayer::DrawIntoBuffer( GLubyte* iBuffer, int iWidth, int iHeight,
                           ViewState& iViewState,
                           ScubaWindowToRASTranslator& iTranslator ) {
  mbBufferWasDrawn = true;
}

void
TestLayer::DrawIntoGL( ViewState&, ScubaWindowToRASTranslator& ) {
  mbGLWasDrawn = true;
}

void
TestLayer::GetInfoAtRAS ( float inX, float inY, float inZ,
                          std::map<std::string,std::string>& iLabelValues ) {

  string sLabel;
  if ( msLabel != "" ) {
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
  float RAS[3];
  int window[2];
  window[0] = iXWindow;
  window[1] = iYWindow;

  iView.TranslateWindowToRAS( window, RAS );
#if 0
  if ( !fequal(RAS[0],iRASX) ||
       !fequal(RAS[1],iRASY) ||
       !fequal(RAS[2],iRASZ) ) {
    cerr << "translate didn't work: " << endl
    << iView.mViewState << endl
    << "   window: " << window[0] << ", " << window[1] << endl
    << "  got RAS: "<< RAS[0]<< ", " << RAS[1] << ", " << RAS[2] << endl
    << "should be: "<< iRASX<< ", " << iRASY << ", " << iRASZ << endl;

    throw 0;
  }
#endif

  int windowTest[2];
  iView.TranslateRASToWindow( RAS, windowTest );
#if 0
  if ( windowTest[0] != window[0] || windowTest[1] != windowTest[1] ) {
    cerr << "translate didn't work: " << endl
    << iView.mViewState << endl
    << "   window: " << window[0] << ", " << window[1] << endl
    << "  got RAS: "<< RAS[0]<< ", " << RAS[1] << ", " << RAS[2] << endl
    << "  back to: " << windowTest[0] << ", " << windowTest[1] << endl
    << "should be: " << window[0] << ", " << window[1] << endl;

    throw 0;
  }
#endif
}

void
ScubaViewTester::Test( Tcl_Interp* iInterp ) {

  stringstream sError;

  try {

    ScubaView view;
    view.SetFlipLeftRightYZ( false );


    // Make sure they both have the default transform.
    int viewTransformID = view.GetWorldToViewTransform();
    Assert( (viewTransformID == 0), "view didn't get default transform" );

    ScubaView view2;
    int view2TransformID = view2.GetWorldToViewTransform();
    Assert( (view2TransformID == 0), "view2 didn't get default transform" );

    ScubaView view3;
    int view3TransformID = view3.GetWorldToViewTransform();
    Assert( (view3TransformID == 0), "view3 didn't get default transform" );

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


    // Set the view size. When we add layers we'll make sure they got
    // the right size.
    view.Reshape( 123, 456 );

    // Add a bunch of layers to different levels and make sure they
    // match our own expectations.
    int const kcLayers = 10;
    TestLayer aLayer[kcLayers];
    int nLevel = 0;
    map<int,int> levelLayerID;
    for ( int nLayer = 0; nLayer < kcLayers; nLayer++ ) {
      int layerID = aLayer[nLayer].GetID();
      levelLayerID[nLevel] = layerID;
      view.SetLayerAtLevel( layerID, nLevel );
      nLevel++;
    }

    bool bFailed = false;
    map<int,int>::iterator tLevelLayerID;
    for ( tLevelLayerID = view.mLevelLayerIDMap.begin();
          tLevelLayerID != view.mLevelLayerIDMap.end(); ++tLevelLayerID ) {
      int level = (*tLevelLayerID).first;
      int layerID = (*tLevelLayerID).second;
      if ( levelLayerID[level] != layerID ) {
        cerr << "ID at level " << level << " incorrect" << endl << ", was "
        << layerID << ", should be " << levelLayerID[level] << endl;
        bFailed = true;
      }
    }
    Assert( (!bFailed), "Adding IDs failed" );

    // Check layer sizes.
    for ( int nLayer = 0; nLayer < kcLayers; nLayer++ ) {
      Assert( (aLayer[nLayer].GetWidth() == 123), "Width not set correctly" );
      Assert( (aLayer[nLayer].GetHeight() == 456), "Height not set correctly" );
    }

    // Check removing layers.
    int levelToRemove = 0;
    view.RemoveLayerAtLevel( levelToRemove );
    for ( tLevelLayerID = view.mLevelLayerIDMap.begin();
          tLevelLayerID != view.mLevelLayerIDMap.end(); ++tLevelLayerID ) {
      int level = (*tLevelLayerID).first;
      int layerID = (*tLevelLayerID).second;
      if ( levelToRemove == level ||
           layerID == levelLayerID[levelToRemove] ) {
        cerr << "ID at level " << level << " was not removed" << endl;
        bFailed = true;
      }
    }
    Assert( (!bFailed), "Removing level failed" );

    // Add it back.
    view.SetLayerAtLevel( levelLayerID[levelToRemove], levelToRemove );


    // Make sure reshaping works and creates our buffer. Make a few to
    // make sure they dispose of the old one properly. Make one with a
    // negative size to check the error handling. Check the
    // width/height on layers each time.
    view.Reshape( 200, 200 );
    Assert( (NULL != view.mBuffer), "buffer not allocated on reshape" );
    Assert( (200 == aLayer[0].mWidth &&
             200 == aLayer[0].mHeight), "layer height/width not correct" );
    view.Reshape( 200, 200 );
    Assert( (NULL != view.mBuffer), "buffer not allocated on second reshape" );
    view.Reshape( 100, 50 );
    Assert( (NULL != view.mBuffer), "buffer not allocated on reshape" );
    Assert( (100 == aLayer[0].mWidth &&
             50 == aLayer[0].mHeight), "layer width/height not correct" );
    view.DisableOutput();
    try {
      view.Reshape( -1, 0 );
      sError << "error not thrown when reshaping with -1, 0" << endl;
      throw  0;
    } catch (...) {}
    try {
      view.Reshape( -1, -1 );
      sError << "error not thrown when reshaping with -1, -1" << endl;
      throw  0;
    } catch (...) {}
    view.SetOutputStreamToCerr();
    view.Reshape( 200, 200 );

    // Draw the view, all our layers should be drawn.
    bFailed = false;
    view.DoDraw();
    for ( int nLayer = 0; nLayer < kcLayers; nLayer++ ) {
      if ( !aLayer[nLayer].WasDrawn() ) {
        cerr << "layer ID " << aLayer[nLayer].GetID() << " not drawn" << endl;
        bFailed = true;
      }
    }
    Assert( (!bFailed), "Drawing failed" );


    // Set labels in all the layers. Call GetInfoAtRAS and check the
    // labelvalue map to make sure they got copied properly.
    map<string,string> labelValueMap;
    bFailed = false;
    for ( int nLayer = 0; nLayer < kcLayers; nLayer++ ) {
      stringstream ssLabel;
      ssLabel << nLayer;
      aLayer[nLayer].SetLabel( ssLabel.str() );
      aLayer[nLayer].GetInfoAtRAS( 0, 0, 0, labelValueMap );
      if ( labelValueMap[ssLabel.str()] != "Hello world" ) {
        cerr << "layer ID " << aLayer[nLayer].GetID() << " not copied" << endl;
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


    // Test the ScubaView tcl commands.
    char sCommand[1024];
    int rTcl;

    sprintf( sCommand, "SetViewRASCenter %d 1.2 3.4 5.6", view.GetID() );
    rTcl = Tcl_Eval( iInterp, sCommand );
    AssertTclOK( rTcl );
    Assert( (view.mViewState.mCenterRAS[0] == (float)1.2 &&
             view.mViewState.mCenterRAS[1] == (float)3.4 &&
             view.mViewState.mCenterRAS[2] == (float)5.6),
            "SetViewRASCenter didn't set properly" );

    sprintf( sCommand, "SetViewZoomLevel %d 3.4", view.GetID() );
    rTcl = Tcl_Eval( iInterp, sCommand );
    AssertTclOK( rTcl );
    Assert( (view.mViewState.mZoomLevel == (float)3.4),
            "SetViewZoomLevel didn't set properly" );

    sprintf( sCommand, "SetViewInPlane %d X", view.GetID() );
    rTcl = Tcl_Eval( iInterp, sCommand );
    AssertTclOK( rTcl );
    Assert( (view.mViewState.mInPlane == ViewState::X),
            "SetViewInPlane didn't set properly" );


    int level = 30;
    sprintf( sCommand, "SetLayerInViewAtLevel %d %d %d",
             view.GetID(), aLayer[0].GetID(), level );
    rTcl = Tcl_Eval( iInterp, sCommand );
    AssertTclOK( rTcl );
    Assert( (view.mLevelLayerIDMap[level] == aLayer[0].GetID()),
            "SetLayerInViewAtLevel didn't work" );

    sprintf( sCommand, "RemoveLayerFromViewAtLevel %d %d",
             view.GetID(), level );
    rTcl = Tcl_Eval( iInterp, sCommand );
    AssertTclOK( rTcl );
    map<int,int>::iterator tLevelLayer = view.mLevelLayerIDMap.find( level );
    Assert( (tLevelLayer == view.mLevelLayerIDMap.end()),
            "RemoveLayerFromViewAtLevel didn't work" );


    // Test the Layer tcl commands.
    sprintf( sCommand, "GetLayerIDList" );
    rTcl = Tcl_Eval( iInterp, sCommand );
    AssertTclOK( rTcl );
    const char* sTclResult = Tcl_GetStringResult( iInterp );
    string sIDList( sTclResult );
    bFailed = false;
    for ( int nLayer = 0; nLayer < kcLayers; nLayer++ ) {
      int id = aLayer[nLayer].GetID();
      stringstream ssID;
      ssID << id;
      string::size_type position = sIDList.find( ssID.str(), 0 );
      if ( position == string::npos ) {
        cerr << "ID " << id << " not found" << endl;
        bFailed = true;
      }
    }
    stringstream ssIDList( sIDList );
    while ( !ssIDList.eof() ) {
      int id;
      ssIDList >> id;
      bool bIDFound = false;
      for ( int nLayer = 0; nLayer < kcLayers; nLayer++ ) {
        if ( id == aLayer[nLayer].GetID() ) {
          bIDFound = true;
          break;
        }
      }
      if ( !bIDFound ) {
        cerr << "incorrect ID returned: " << id << endl;
        bFailed = true;
      }
    }
    Assert( (!bFailed), "GetLayerIDList failed" );

    Layer& layer = aLayer[0];
    int id = layer.GetID();
    sprintf( sCommand, "GetLayerType %d", id );
    rTcl = Tcl_Eval( iInterp, sCommand );
    AssertTclOK( rTcl );
    sTclResult = Tcl_GetStringResult( iInterp );
    string sLayerType( sTclResult );
    Assert( (sLayerType == "TestLayer"), "GetLayerType didn't work" );

    layer.SetOpacity( 1.5 );
    sprintf( sCommand, "GetLayerOpacity %d", id );
    rTcl = Tcl_Eval( iInterp, sCommand );
    AssertTclOK( rTcl );
    sTclResult = Tcl_GetStringResult( iInterp );
    string sOpacity( sTclResult );
    Assert( (sOpacity == "1.5"), "GetLayerOpacity didn't work" );

    sprintf( sCommand, "SetLayerOpacity %d 5.5", id );
    rTcl = Tcl_Eval( iInterp, sCommand );
    AssertTclOK( rTcl );
    Assert( (layer.GetOpacity() == (float)5.5),
            "SetLayerOpacity didn't work" );

    sprintf( sCommand, "SetLayerLabel %d testLabel", id );
    rTcl = Tcl_Eval( iInterp, sCommand );
    AssertTclOK( rTcl );
    Assert( (layer.GetLabel() == "testLabel"),
            "SetLayerLabel didn't work" );

    sprintf( sCommand, "GetLayerLabel %d", id );
    rTcl = Tcl_Eval( iInterp, sCommand );
    AssertTclOK( rTcl );
    sTclResult = Tcl_GetStringResult( iInterp );
    string sLabel( sTclResult );
    Assert( (sLabel == layer.GetLabel()), "GetLayerLabel didn't work" );

    // Check removing all the layers.
    view.RemoveAllLayers();
    Assert( (view.mLevelLayerIDMap.size() == 0 ), "RemoveAllLayers failed" );

  } catch (...) {
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

    ScubaTransform t;

    // Test constructor/destructor segfaults.
    {
      ScubaView v;
    }

    for ( int nTrial = 0; nTrial < 50; nTrial++ ) {
      ScubaViewTester tester0;
      tester0.Test( interp );
    }

    Tcl_DeleteInterp( interp );
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
