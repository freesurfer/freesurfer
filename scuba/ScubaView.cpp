#include <stdexcept>
#include "glut.h"
#include "ScubaView.h"
#include "PreferencesManager.h"

using namespace std;

int const ScubaView::kBytesPerPixel = 4;

ScubaView::ScubaView() {
  mBuffer = NULL;
  mbPostRedisplay = false;
  mbRebuildOverlayDrawList = true;

  TclCommandManager& commandMgr = TclCommandManager::GetManager();
  commandMgr.AddCommand( *this, "SetViewInPlane" );
  commandMgr.AddCommand( *this, "SetViewZoomLevel" );
  commandMgr.AddCommand( *this, "SetViewRASCenter" );
  commandMgr.AddCommand( *this, "AddLayerToView" );
  commandMgr.AddCommand( *this, "RemoveAllLayersFromView" );
  commandMgr.AddCommand( *this, "RemoveLayerFromViewAtLevel" );

  PreferencesManager& prefsMgr = PreferencesManager::GetManager();
  prefsMgr.UseFile( ".scuba" );

  PreferencesManager::StringPrefValue moveViewLeft( "Left" );
  prefsMgr.RegisterValue( "key-MoveViewLeft", 
			  "Key to move the view to the left.",
			  moveViewLeft );
  msMoveViewLeft = prefsMgr.GetValue( "key-MoveViewLeft" );

  PreferencesManager::StringPrefValue moveViewRight( "Right" );
  prefsMgr.RegisterValue( "key-MoveViewRight", 
			  "Key to move the view to the right.",
			  moveViewRight );
  msMoveViewRight = prefsMgr.GetValue( "key-MoveViewRight" );

  PreferencesManager::StringPrefValue moveViewUp( "Up" );
  prefsMgr.RegisterValue( "key-MoveViewUp", 
			  "Key to move the view up.",
			  moveViewUp );
  msMoveViewUp = prefsMgr.GetValue( "key-MoveViewUp" );

  PreferencesManager::StringPrefValue moveViewDown( "Down" );
  prefsMgr.RegisterValue( "key-MoveViewDown", 
			  "Key to move the view down.",
			  moveViewDown );
  msMoveViewDown = prefsMgr.GetValue( "key-MoveViewDown" );

  PreferencesManager::StringPrefValue moveViewIn( "period" );
  prefsMgr.RegisterValue( "key-MoveViewIn", 
			  "Key to move the view in in plane.",
			  moveViewIn );
  msMoveViewIn = prefsMgr.GetValue( "key-MoveViewIn" );

  PreferencesManager::StringPrefValue moveViewOut( "comma" );
  prefsMgr.RegisterValue( "key-MoveViewOut", 
			  "Key to move the view out in plane.",
			  moveViewOut );
  msMoveViewOut = prefsMgr.GetValue( "key-MoveViewOut" );

  PreferencesManager::StringPrefValue zoomViewIn( "equal" );
  prefsMgr.RegisterValue( "key-ZoomViewIn", 
			  "Key to zoom the view in in plane.",
			  zoomViewIn );
  msZoomViewIn = prefsMgr.GetValue( "key-ZoomViewIn" );

  PreferencesManager::StringPrefValue zoomViewOut( "minus" );
  prefsMgr.RegisterValue( "key-ZoomViewOut", 
			  "Key to zoom the view out in plane.",
			  zoomViewOut );
  msZoomViewOut = prefsMgr.GetValue( "key-ZoomViewOut" );

  PreferencesManager::StringPrefValue inPlaneX( "x" );
  prefsMgr.RegisterValue( "key-InPlaneX", 
			  "Key to change in plane to X in the view.",
			  inPlaneX );
  msInPlaneXKey = prefsMgr.GetValue( "key-InPlaneX" );

  PreferencesManager::StringPrefValue inPlaneY( "y" );
  prefsMgr.RegisterValue( "key-InPlaneY", 
			  "Key to change in plane to Y in the view.",
			  inPlaneY );
  msInPlaneYKey = prefsMgr.GetValue( "key-InPlaneY" );

  PreferencesManager::StringPrefValue inPlaneZ( "z" );
  prefsMgr.RegisterValue( "key-InPlaneZ", 
			  "Key to change in plane to Z in the view.",
			  inPlaneZ );
  msInPlaneZKey = prefsMgr.GetValue( "key-InPlaneZ" );

  
}

ScubaView::~ScubaView() {
}

void
ScubaView::Set2DRASCenter ( float iRASCenter[] ) {

  mViewState.mCenterRAS[0] = iRASCenter[0];
  mViewState.mCenterRAS[1] = iRASCenter[1];
  mViewState.mCenterRAS[2] = iRASCenter[2];
}

void
ScubaView::Set2DZoomLevel ( float iZoomLevel ) {

  mViewState.mZoomLevel = iZoomLevel;
}

void
ScubaView::Set2DInPlane ( ViewState::Plane iPlane ) {

  mViewState.mInPlane = iPlane;
}

void 
ScubaView::AddLayer ( int iLayerID, int iLevel ) {

  mLevelLayerIDMap[iLevel] = iLayerID;

  // Set the layer's width and height.
  try {
    Layer& layer = Layer::FindByID( iLayerID );
    
    layer.SetWidth( mWidth );
    layer.SetHeight( mHeight );
  }
  catch(...) {
    DebugOutput( << "Couldn't find layer " << iLayerID );
  }
}

void 
ScubaView::RemoveAllLayers () {

  mLevelLayerIDMap.clear();
}

void 
ScubaView::RemoveLayerAtLevel ( int iLevel ) {

  mLevelLayerIDMap.erase( iLevel );
}

void
ScubaView::DoListenToTclCommand( char* isCommand, int iArgc, char** iasArgv ) {

  // SetViewInPlane <viewID> <inPlane>
  if( 0 == strcmp( isCommand, "SetViewInPlane" ) ) {
    if( 3 == iArgc ) {
      int viewID = strtol(iasArgv[1], (char**)NULL, 10);
      if( ERANGE == errno ) {
	sResult = "bad view ID";
	return;
      }

      if( mID == viewID ) {
	
	ViewState::Plane inPlane;
	if( 0 == strcmp( iasArgv[2], "X" )  || 
	    0 == strcmp( iasArgv[2], "x" ) ) {
	  inPlane = ViewState::X;
	} else if( 0 == strcmp( iasArgv[2], "Y" )  || 
	    0 == strcmp( iasArgv[2], "y" ) ) {
	  inPlane = ViewState::Y;
	} else if( 0 == strcmp( iasArgv[2], "Z" )  || 
	    0 == strcmp( iasArgv[2], "z" ) ) {
	  inPlane = ViewState::Z;
	} else {
	  sResult = "bad inPlane \"" + string(iasArgv[2]) + 
	    "\", should be x, y, or z";
	  return;
	}

	Set2DInPlane( inPlane );
      }
    } else {
      sResult = "wrong # args: should be \"SetViewInPlane "
	"viewID inPlane\"";
      DebugOutput( << sResult );
      return;
    }
  }

  // SetViewZoomLevel <viewID> <zoomLevel>
  if( 0 == strcmp( isCommand, "SetViewZoomLevel" ) ) {
    if( 3 == iArgc ) {
      int viewID = strtol(iasArgv[1], (char**)NULL, 10);
      if( ERANGE == errno ) {
	sResult = "bad view ID";
	return;
      }

      if( mID == viewID ) {
	
	float zoomLevel = strtof( iasArgv[2], (char**)NULL );
	if( ERANGE == errno ) {
	  sResult = "bad zoom level";
	  return;
	}
	  
	Set2DZoomLevel( zoomLevel );
      }
    } else {
      sResult = "wrong # args: should be \"SetViewZoomLevel "
	"viewID zoomLevel\"";
      DebugOutput( << sResult );
      return;
    }
  }

  // SetViewRASCenter <viewID> <X> <Y> <Z>
  if( 0 == strcmp( isCommand, "SetViewRASCenter" ) ) {
    if( 5 == iArgc ) {
      int viewID = strtol(iasArgv[1], (char**)NULL, 10);
      if( ERANGE == errno ) {
	sResult = "bad view ID";
	return;
      }

      if( mID == viewID ) {
	
	float x = strtof( iasArgv[2], (char**)NULL );
	if( ERANGE == errno ) {
	  sResult = "bad x coordinate";
	  return;
	}
	float y = strtof( iasArgv[3], (char**)NULL );
	if( ERANGE == errno ) {
	  sResult = "bad y coordinate";
	  return;
	}
	float z = strtof( iasArgv[4], (char**)NULL );
	if( ERANGE == errno ) {
	  sResult = "bad z coordinate";
	  return;
	}
	
	float center[3];
	center[0] = x; center[1] = y; center[2] = z;
	Set2DRASCenter( center );
      }
    } else {
      sResult = "wrong # args: should be \"SetViewRASCenter "
	"viewID X Y Z\"";
      DebugOutput( << sResult );
      return;
    }
  }

  // AddLayerToView <viewID> <layerID> <level>
  if( 0 == strcmp( isCommand, "AddLayerToView" ) ) {
    if( 4 == iArgc ) {
      int viewID = strtol(iasArgv[1], (char**)NULL, 10);
      if( ERANGE == errno ) {
	sResult = "bad view ID";
	return;
      }

      if( mID == viewID ) {
	
	int layerID = strtol( iasArgv[2], (char**)NULL, 10 );
	if( ERANGE == errno ) {
	  sResult = "bad layer ID";
	  return;
	}
	int level = strtol( iasArgv[3], (char**)NULL, 10 );
	if( ERANGE == errno ) {
	  sResult = "bad level";
	  return;
	}

	AddLayer( layerID, level );
      }
    } else {
      sResult = "wrong # args: should be \"AddLayerToView "
	"viewID layerID level\"";
      DebugOutput( << sResult );
      return;
    }
  }

  // RemoveLayerFromViewAtLevel <viewID> <level>
  if( 0 == strcmp( isCommand, "RemoveLayerFromViewAtLevel" ) ) {
    if( 3 == iArgc ) {
      int viewID = strtol(iasArgv[1], (char**)NULL, 10);
      if( ERANGE == errno ) {
	sResult = "bad view ID";
	return;
      }

      if( mID == viewID ) {
	
	int level = strtol( iasArgv[2], (char**)NULL, 10 );
	if( ERANGE == errno ) {
	  sResult = "bad level";
	  return;
	}

	RemoveLayerAtLevel( level );
      }
    } else {
      sResult = "wrong # args: should be \"RemoveLayerFromViewAtLevel "
	"viewID level\"";
      DebugOutput( << sResult );
      return;
    }
  }

  // RemoveAllLayersFromView <viewID>
  if( 0 == strcmp( isCommand, "RemoveAllLayersFromView" ) ) {
    if( 2 == iArgc ) {
      int viewID = strtol(iasArgv[1], (char**)NULL, 10);
      if( ERANGE == errno ) {
	sResult = "bad view ID";
	return;
      }

      if( mID == viewID ) {
	
	RemoveAllLayers();
      }
    } else {
      sResult = "wrong # args: should be \"RemoveAllLayersFromView viewID\"";
      DebugOutput( << sResult );
      return;
    }
  }

}

void
ScubaView::DoDraw() {

  BuildFrameBuffer();
  DrawFrameBuffer();
  DrawOverlay();
}

void
ScubaView::DoReshape( int iWidth, int iHeight ) {

  if( iWidth < 0 || iHeight < 0 ) {
    stringstream sError;
    sError << "Invalid width " << mWidth << " or height " << mHeight;
    DebugOutput( << sError.str() );
    throw new runtime_error( sError.str() );
  }

  // Allocate a new buffer.
  GLubyte* newBuffer = (GLubyte*) malloc( mWidth * mHeight * kBytesPerPixel );
  if( NULL == newBuffer ) {
    stringstream sError;
    sError << "Couldn't allocate buffer for width " << mWidth 
	   << " height " << mHeight;
    DebugOutput( << sError.str() );
    throw new runtime_error( sError.str() );
  }

  // Get rid of the old one.
  if( NULL != mBuffer ) {
    free( mBuffer );
  }
  
  // Save the new one.
  mBuffer = newBuffer;

  // Set the width and height in all the layers.
  map<int,int>::iterator tLevelLayerID;
  for( tLevelLayerID = mLevelLayerIDMap.begin(); 
       tLevelLayerID != mLevelLayerIDMap.end(); ++tLevelLayerID ) {

      int layerID = (*tLevelLayerID).second;

    try {
      Layer& layer = Layer::FindByID( layerID );
      
      layer.SetWidth( mWidth );
      layer.SetHeight( mHeight );
    }
    catch(...) {
      DebugOutput( << "Couldn't find layer " << layerID );
    }
  }

  mbRebuildOverlayDrawList = true;
}

void
ScubaView::DoTimer() {

}
  
void
ScubaView::DoMouseMoved( int inX, int inY, InputState& iState ) {

  float rasX, rasY, rasZ;
  TranslateWindowToRAS( inX, inY, rasX, rasY, rasZ );

  map<string,string> labelValueMap;

  // Go through our draw levels. For each one, get the Layer.
  map<int,int>::iterator tLevelLayerID;
  for( tLevelLayerID = mLevelLayerIDMap.begin(); 
       tLevelLayerID != mLevelLayerIDMap.end(); ++tLevelLayerID ) {

    int layerID = (*tLevelLayerID).second;
    try {
      Layer& layer = Layer::FindByID( layerID );
      
      // Ask the layer for info strings at this point.
      layer.GetInfoAtRAS( rasX, rasY, rasZ, labelValueMap );
    }
    catch(...) {
      DebugOutput( << "Couldn't find layer " << layerID );
    }
  }

  DrawFrameBuffer();
  DrawOverlay();

  // For each one, get the strings.
  int nLine = 0;
  map<string,string>::iterator tLabelValue;
  for( tLabelValue = labelValueMap.begin(); 
       tLabelValue != labelValueMap.end(); ++tLabelValue ) {

    string sLabel = (*tLabelValue).first;
    string sValue = (*tLabelValue).second;

    // Draw them to the screen.
    glRasterPos2i( inX, inY + (nLine * 18));
    glColor3f( 1, 1, 1 );
    for( int nChar = 0; nChar < sLabel.length(); nChar++ ) {
      glutBitmapCharacter( GLUT_BITMAP_8_BY_13, sLabel[nChar] );
    }
    for( int nChar = 0; nChar < sValue.length(); nChar++ ) {
      glutBitmapCharacter( GLUT_BITMAP_8_BY_13, sValue[nChar] );
    }
  }
}

void
ScubaView::DoMouseUp( int inX, int inY, InputState& iState ) {
}

void
ScubaView::DoMouseDown( int inX, int inY, InputState& iState ) {

}

void
ScubaView::DoKeyDown( int inX, int inY, InputState& iState ) {

  string key = iState.Key();
  if( key == msMoveViewLeft ) {
    switch( mViewState.mInPlane ) {
    case ViewState::X: mViewState.mCenterRAS[1] -= 1.0; break;
    case ViewState::Y: mViewState.mCenterRAS[0] -= 1.0; break;
    case ViewState::Z: mViewState.mCenterRAS[0] -= 1.0; break;
    }
    mbRebuildOverlayDrawList = true;
  } else if( key == msMoveViewRight ) {
    switch( mViewState.mInPlane ) {
    case ViewState::X: mViewState.mCenterRAS[1] += 1.0; break;
    case ViewState::Y: mViewState.mCenterRAS[0] += 1.0; break;
    case ViewState::Z: mViewState.mCenterRAS[0] += 1.0; break;
    }
    mbRebuildOverlayDrawList = true;
  } else if( key == msMoveViewDown ) {
    switch( mViewState.mInPlane ) {
    case ViewState::X: mViewState.mCenterRAS[2] -= 1.0; break;
    case ViewState::Y: mViewState.mCenterRAS[2] -= 1.0; break;
    case ViewState::Z: mViewState.mCenterRAS[1] -= 1.0; break;
    }
    mbRebuildOverlayDrawList = true;
  } else if( key == msMoveViewUp ) {
    switch( mViewState.mInPlane ) {
    case ViewState::X: mViewState.mCenterRAS[2] += 1.0; break;
    case ViewState::Y: mViewState.mCenterRAS[2] += 1.0; break;
    case ViewState::Z: mViewState.mCenterRAS[1] += 1.0; break;
    }
    mbRebuildOverlayDrawList = true;
  } else if( key == msMoveViewIn ) {
    switch( mViewState.mInPlane ) {
    case ViewState::X: mViewState.mCenterRAS[0] += 1.0; break;
    case ViewState::Y: mViewState.mCenterRAS[1] += 1.0; break;
    case ViewState::Z: mViewState.mCenterRAS[2] += 1.0; break;
    }
    mbRebuildOverlayDrawList = true;
  } else if( key == msMoveViewOut ) {
    switch( mViewState.mInPlane ) {
    case ViewState::X: mViewState.mCenterRAS[0] -= 1.0; break;
    case ViewState::Y: mViewState.mCenterRAS[1] -= 1.0; break;
    case ViewState::Z: mViewState.mCenterRAS[2] -= 1.0; break;
    }
    mbRebuildOverlayDrawList = true;

  } else if( key == msZoomViewIn ) {
    mViewState.mZoomLevel *= 2.0;
    mbRebuildOverlayDrawList = true;
  } else if( key == msZoomViewOut ) {
    mViewState.mZoomLevel /= 2.0;
    if( mViewState.mZoomLevel < 0.25 ) {
      mViewState.mZoomLevel = 0.25;
    }
    mbRebuildOverlayDrawList = true;

  } else if( key == msInPlaneXKey ) {
    mViewState.mInPlane = ViewState::X;
    mbRebuildOverlayDrawList = true;
  } else if( key == msInPlaneYKey ) {
    mViewState.mInPlane = ViewState::Y;
    mbRebuildOverlayDrawList = true;
  } else if( key == msInPlaneZKey ) {
    mViewState.mInPlane = ViewState::Z;
    mbRebuildOverlayDrawList = true;
  }

  RequestRedisplay();
}

void
ScubaView::DoKeyUp( int inX, int inY, InputState& iState ) {

}

void
ScubaView::TranslateWindowToRAS ( int iXWindow, int iYWindow,
				  float& oXRAS, float& oYRAS, float& oZRAS ) {

  // At zoom level one every pixel is an RAS point, so we're scaled by
  // that and then offset by the RAS center. We find a 3D point by
  // using our mInPlane and the corresponding mCenterRAS in ViewState.

  switch( mViewState.mInPlane ) {
  case ViewState::X:
    oXRAS = mViewState.mCenterRAS[0];
    oYRAS = ConvertWindowToRAS( iXWindow, mViewState.mCenterRAS[1], mWidth );
    oZRAS = ConvertWindowToRAS( iYWindow, mViewState.mCenterRAS[2], mHeight );
    break;
  case ViewState::Y:
    oXRAS = ConvertWindowToRAS( iXWindow, mViewState.mCenterRAS[0], mWidth );
    oYRAS = mViewState.mCenterRAS[1];
    oZRAS = ConvertWindowToRAS( iYWindow, mViewState.mCenterRAS[2], mHeight );
    break;
  case ViewState::Z:
    oXRAS = ConvertWindowToRAS( iXWindow, mViewState.mCenterRAS[0], mWidth );
    oYRAS = ConvertWindowToRAS( iYWindow, mViewState.mCenterRAS[1], mHeight );
    oZRAS = mViewState.mCenterRAS[2];
    break;
  }
}

float
ScubaView::ConvertWindowToRAS ( float iWindow, float iRASCenter, 
				float iWindowDimension ) {

  return (iWindow / mViewState.mZoomLevel) +
    iRASCenter - (iWindowDimension / 2.0 / mViewState.mZoomLevel);
}


void 
ScubaView::BuildFrameBuffer () {

  // Don't draw if our buffer isn't initialized yet.
  if( NULL == mBuffer )
    return;

  // Go through our draw levels. For each one, get the Layer.
  map<int,int>::iterator tLevelLayerID;
  for( tLevelLayerID = mLevelLayerIDMap.begin(); 
       tLevelLayerID != mLevelLayerIDMap.end(); ++tLevelLayerID ) {

    int layerID = (*tLevelLayerID).second;
    try { 
      Layer& layer = Layer::FindByID( layerID );
    
      // tell it to draw into our buffer with our view state information.
      layer.DrawIntoBuffer( mBuffer, mViewState, *this );
    }
    catch(...) {
      DebugOutput( << "Couldn't find layer " << layerID );
    }
  }
}

void
ScubaView::DrawFrameBuffer () {

  glRasterPos2i( 0, 0 );
  glDrawPixels( mWidth, mHeight, GL_RGBA, GL_UNSIGNED_BYTE, mBuffer );
}

void 
ScubaView::BuildOverlay () {

  if( !mbRebuildOverlayDrawList )
    return;

  // Draw the HUD overlay if necessary.
  float left, right, top, bottom, plane;
  switch( mViewState.mInPlane ) {
    case ViewState::X: 
      left =   mViewState.mCenterRAS[1] - (mWidth/2 / mViewState.mZoomLevel);
      right =  mViewState.mCenterRAS[1] + (mWidth/2 / mViewState.mZoomLevel);
      bottom = mViewState.mCenterRAS[2] - (mHeight/2 / mViewState.mZoomLevel);
      top =    mViewState.mCenterRAS[2] + (mHeight/2 / mViewState.mZoomLevel);
      plane =  mViewState.mCenterRAS[0];
      break;
    case ViewState::Y: 
      left =   mViewState.mCenterRAS[0] - (mWidth/2 / mViewState.mZoomLevel);
      right =  mViewState.mCenterRAS[0] + (mWidth/2 / mViewState.mZoomLevel);
      bottom = mViewState.mCenterRAS[2] - (mHeight/2 / mViewState.mZoomLevel);
      top =    mViewState.mCenterRAS[2] + (mHeight/2 / mViewState.mZoomLevel);
      plane =  mViewState.mCenterRAS[1];
      break;
    case ViewState::Z: 
      left =   mViewState.mCenterRAS[0] - (mWidth/2 / mViewState.mZoomLevel);
      right =  mViewState.mCenterRAS[0] + (mWidth/2 / mViewState.mZoomLevel);
      bottom = mViewState.mCenterRAS[1] - (mHeight/2 / mViewState.mZoomLevel);
      top =    mViewState.mCenterRAS[1] + (mHeight/2 / mViewState.mZoomLevel);
      plane =  mViewState.mCenterRAS[2];
      break;
  }


  glNewList( kOverlayDrawListID, GL_COMPILE );

  glColor3f( 1, 1, 1 );

  char sLabel[60];
  sprintf( sLabel, "%.2f", left );
  glRasterPos2i( 0, mHeight / 2 );
  for( int nChar = 0; nChar < strlen( sLabel ); nChar++ ) {
    glutBitmapCharacter( GLUT_BITMAP_8_BY_13, sLabel[nChar] );
  }
  sprintf( sLabel, "%.2f", right );
  glRasterPos2i( mWidth - strlen(sLabel)*8, mHeight / 2 );
  for( int nChar = 0; nChar < strlen( sLabel ); nChar++ ) {
    glutBitmapCharacter( GLUT_BITMAP_8_BY_13, sLabel[nChar] );
  }
  sprintf( sLabel, "%.2f", top );
  glRasterPos2i( mWidth / 2 - (strlen(sLabel)*8 / 2), mHeight-1-13 );
  for( int nChar = 0; nChar < strlen( sLabel ); nChar++ ) {
    glutBitmapCharacter( GLUT_BITMAP_8_BY_13, sLabel[nChar] );
  }
  sprintf( sLabel, "%.2f", bottom );
  glRasterPos2i( mWidth / 2 - (strlen(sLabel)*8 / 2), 4 );
  for( int nChar = 0; nChar < strlen( sLabel ); nChar++ ) {
    glutBitmapCharacter( GLUT_BITMAP_8_BY_13, sLabel[nChar] );
  }
  sprintf( sLabel, "%.2f", plane );
  glRasterPos2i( mWidth - (strlen(sLabel)*8), 1 );
  for( int nChar = 0; nChar < strlen( sLabel ); nChar++ ) {
    glutBitmapCharacter( GLUT_BITMAP_8_BY_13, sLabel[nChar] );
  }
  if( mViewState.mZoomLevel != 1 ) {
    sprintf( sLabel, "%.2fx", mViewState.mZoomLevel );
    glRasterPos2i( 0, mHeight-1-13 );
    for( int nChar = 0; nChar < strlen( sLabel ); nChar++ ) {
      glutBitmapCharacter( GLUT_BITMAP_8_BY_13, sLabel[nChar] );
    }
  }

  glEndList();

  mbRebuildOverlayDrawList = false;
}

void
ScubaView::DrawOverlay () {
  
  if( mbRebuildOverlayDrawList )
    BuildOverlay();

  glCallList( kOverlayDrawListID );
}

