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
  commandMgr.AddCommand( *this, "SetViewInPlane", 2, "viewID inPlane",
			 "Sets the in plane in a view. inPlane should be "
			 "one of the following: x y z" );
  commandMgr.AddCommand( *this, "SetViewZoomLevel", 2, "viewID zoomLevel",
			 "Sets the zoom level in a view. zoomLevel should be "
			 "a float." );
  commandMgr.AddCommand( *this, "SetViewRASCenter", 4, "viewID x y z",
			 "Sets the view center. x, y, and z should be floats "
			 "in world RAS coordinates." );
  commandMgr.AddCommand( *this, "AddLayerToView", 3, "viewID layerID level",
			 "Sets the layer in a view at a given draw level. "
			 "Higher draw levels will draw later." );
  commandMgr.AddCommand( *this, "RemoveAllLayersFromView", 1, "viewID",
			 "Remove all layers from a view." );
  commandMgr.AddCommand( *this, "RemoveLayerFromViewAtLevel", 2, 
			 "viewID layer",
			 "Remove a layer from a view." );
  commandMgr.AddCommand( *this, "GetLabelValuesSet", 2, "viewID setName",
			 "Get a set of label value pairs." );

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

  map<string,string> labelValueMap;
  mLabelValueMaps["cursor"] = labelValueMap;
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

  RequestRedisplay();
}

void 
ScubaView::RemoveAllLayers () {

  mLevelLayerIDMap.clear();
}

void 
ScubaView::RemoveLayerAtLevel ( int iLevel ) {

  mLevelLayerIDMap.erase( iLevel );
}

TclCommandListener::TclCommandResult
ScubaView::DoListenToTclCommand( char* isCommand, int iArgc, char** iasArgv ) {

  // SetViewInPlane <viewID> <inPlane>
  if( 0 == strcmp( isCommand, "SetViewInPlane" ) ) {
    int viewID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad view ID";
      return error;
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
	return error;
      }
      
      Set2DInPlane( inPlane );
    }
  }

  // SetViewZoomLevel <viewID> <zoomLevel>
  if( 0 == strcmp( isCommand, "SetViewZoomLevel" ) ) {
    int viewID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad view ID";
      return error;
    }
    
    if( mID == viewID ) {
      
      float zoomLevel = strtof( iasArgv[2], (char**)NULL );
      if( ERANGE == errno ) {
	sResult = "bad zoom level";
	return error;
      }
      
      Set2DZoomLevel( zoomLevel );
    }
  }

  // SetViewRASCenter <viewID> <X> <Y> <Z>
  if( 0 == strcmp( isCommand, "SetViewRASCenter" ) ) {
    int viewID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad view ID";
      return error;
    }

    if( mID == viewID ) {
      
      float x = strtof( iasArgv[2], (char**)NULL );
      if( ERANGE == errno ) {
	sResult = "bad x coordinate";
	return error;
      }
      float y = strtof( iasArgv[3], (char**)NULL );
      if( ERANGE == errno ) {
	sResult = "bad y coordinate";
	return error;
      }
      float z = strtof( iasArgv[4], (char**)NULL );
      if( ERANGE == errno ) {
	sResult = "bad z coordinate";
	return error;
      }
      
      float center[3];
      center[0] = x; center[1] = y; center[2] = z;
      Set2DRASCenter( center );
    }
  }

  // AddLayerToView <viewID> <layerID> <level>
  if( 0 == strcmp( isCommand, "AddLayerToView" ) ) {
    int viewID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad view ID";
      return error;
    }
    
    if( mID == viewID ) {
      
      int layerID = strtol( iasArgv[2], (char**)NULL, 10 );
      if( ERANGE == errno ) {
	sResult = "bad layer ID";
	return error;
      }
      int level = strtol( iasArgv[3], (char**)NULL, 10 );
      if( ERANGE == errno ) {
	sResult = "bad level";
	return error;
      }
      
      AddLayer( layerID, level );
    }
  }

  // RemoveLayerFromViewAtLevel <viewID> <level>
  if( 0 == strcmp( isCommand, "RemoveLayerFromViewAtLevel" ) ) {
    int viewID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad view ID";
      return error;
    }
    
    if( mID == viewID ) {
      
      int level = strtol( iasArgv[2], (char**)NULL, 10 );
      if( ERANGE == errno ) {
	sResult = "bad level";
	return error;
      }
      
      RemoveLayerAtLevel( level );
    }
  }

  // RemoveAllLayersFromView <viewID>
  if( 0 == strcmp( isCommand, "RemoveAllLayersFromView" ) ) {
    int viewID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad view ID";
      return error;
    }
    
    if( mID == viewID ) {
      
      RemoveAllLayers();
    }
  }

  // GetLabelValuesSet <viewID> <setName>
  if( 0 == strcmp( isCommand, "GetLabelValuesSet" ) ) {
    int viewID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad view ID";
      return error;
    }
    
    if( mID == viewID ) {

      string sSetName = iasArgv[2];

      map<string,map<string,string> >::iterator tMap = 
	mLabelValueMaps.find( sSetName );
      if( tMap == mLabelValueMaps.end() ) {
	stringstream ssResult;
	ssResult << "Set name " << sSetName << " doesn't exist.";
	sResult = ssResult.str();
	return error;
      }

      map<string,string> labelValueMap = mLabelValueMaps[sSetName];
      map<string,string>::iterator tLabelValue;

      stringstream ssFormat;
      stringstream ssResult;
      ssFormat << "L";

      for( tLabelValue = labelValueMap.begin(); 
	   tLabelValue != labelValueMap.end(); ++tLabelValue ) {
	
	ssFormat << "Lssl";
	ssResult << "\"" << (*tLabelValue).first << "\" \"" 
		 << (*tLabelValue).second << "\" ";
      }
      ssFormat << "l";
      
      sReturnFormat = ssFormat.str();
      sReturnValues = ssResult.str();
    }
  }

  return ok;
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

  stringstream sID;
  sID << GetID();
  labelValueMap["View ID"] = sID.str();

  // Get the RAS coords into a string and set that label/value.
  stringstream ssRASCoords;
  ssRASCoords << rasX << " " << rasY << " " << rasZ;
  labelValueMap["RAS"] = ssRASCoords.str();

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

  RequestRedisplay();

  //  DrawFrameBuffer();
  //  DrawOverlay();

  // Set this labelValueMap in the array of label values under the
  // name 'cursor'.
  mLabelValueMaps["cursor"] = labelValueMap;


#if 0
  // For each one, get the strings.
  int nLine = 0;
  map<string,string>::iterator tLabelValue;
  for( tLabelValue = labelValueMap.begin(); 
       tLabelValue != labelValueMap.end(); ++tLabelValue ) {

    string sLabel = (*tLabelValue).first + ": ";
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
    // I don't know why, but the last char in sValue gets overdrawn,
    // so draw this last space to make sure sValue is there.
    glutBitmapCharacter( GLUT_BITMAP_8_BY_13, ' ' );

    nLine++;
  }
#endif
}

void
ScubaView::DoMouseUp( int inX, int inY, InputState& iState ) {

  // No matter what tool we're on, look for ctrl-b{1,2,3} and do some
  // navigation stuff.
  if( iState.IsControlKeyDown() && 
      !iState.IsShiftKeyDown() && !iState.IsAltKeyDown() ) {
    
    // Set the new view center to this point. If they also hit b1 or
    // b3, zoom in or out accordingly.
    float world[3];
    TranslateWindowToRAS( inX, inY, world );
    Set2DRASCenter( world );

    switch( iState.Button() ) {
    case 1:
      mViewState.mZoomLevel *= 2.0;
      break;
    case 3:
      mViewState.mZoomLevel /= 2.0;
      if( mViewState.mZoomLevel < 0.25 ) {
	mViewState.mZoomLevel = 0.25;
      }
      break;
    }

    mbRebuildOverlayDrawList = true;
    RequestRedisplay();
  }

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
      layer.DrawIntoBuffer( mBuffer, mWidth, mHeight, mViewState, *this );
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


  glNewList( kOverlayDrawListID + mID, GL_COMPILE );

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

  glCallList( kOverlayDrawListID + mID );
}

