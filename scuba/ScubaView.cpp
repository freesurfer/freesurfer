#include <stdexcept>
#include "glut.h"
#include "ScubaView.h"
#include "PreferencesManager.h"

using namespace std;

int const ScubaView::kBytesPerPixel = 4;
map<int,bool> ScubaView::mViewIDLinkedList;
int ScubaView::mCurrentBroadcaster = -1;

ScubaView::ScubaView() {
  mBuffer = NULL;
  mbPostRedisplay = false;
  mbRebuildOverlayDrawList = true;
  mViewIDLinkedList[GetID()] = false;
  //  mWorldToView = MatrixAlloc( 4, 4, MATRIX_REAL );
  mWorldRAS = VectorAlloc( 4, MATRIX_REAL );
  mViewRAS = VectorAlloc( 4, MATRIX_REAL );

  // -1 0 0 0
  //  0 0 1 0
  //  0 1 0 0
  //  0 0 0 1

#define FLIP_Z
#define VIEW

  float radians = (float)M_PI;
  MATRIX* rotate = MatrixIdentity( 4, NULL );
#ifdef FLIP_X
  *MATRIX_RELT(rotate,1,1) = 1;
  *MATRIX_RELT(rotate,2,2) = cos(radians);
  *MATRIX_RELT(rotate,2,3) = -sin(radians);
  *MATRIX_RELT(rotate,3,2) = sin(radians);
  *MATRIX_RELT(rotate,3,3) = cos(radians);
  *MATRIX_RELT(rotate,4,4) = 1;
#endif
#ifdef FLIP_Y
  *MATRIX_RELT(rotate,1,1) = cos(radians);
  *MATRIX_RELT(rotate,1,3) = sin(radians);
  *MATRIX_RELT(rotate,2,2) = 1;
  *MATRIX_RELT(rotate,3,1) = -sin(radians);
  *MATRIX_RELT(rotate,3,3) = cos(radians);
  *MATRIX_RELT(rotate,4,4) = 1;
#endif
#ifdef FLIP_Z
  *MATRIX_RELT(rotate,1,1) = cos(radians);
  *MATRIX_RELT(rotate,1,2) = -sin(radians);
  *MATRIX_RELT(rotate,2,1) = sin(radians);
  *MATRIX_RELT(rotate,2,2) = cos(radians);
  *MATRIX_RELT(rotate,3,3) = 1;
  *MATRIX_RELT(rotate,4,4) = 1;
#endif

#ifdef WORLD
  mWorldToView = rotate;
  mViewToWorld = MatrixInverse( mWorldToView, NULL );
#endif
#ifdef VIEW
  mViewToWorld = rotate;
  mWorldToView = MatrixInverse( mViewToWorld, NULL );
#endif

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
  commandMgr.AddCommand( *this, "SetLayerInViewAtLevel", 3, 
			 "viewID layerID level",
			 "Sets the layer in a view at a given draw level. "
			 "Higher draw levels will draw later." );
  commandMgr.AddCommand( *this, "GetLayerInViewAtLevel", 2, "viewID level",
			 "Returns the layer in a view at a given draw level.");
  commandMgr.AddCommand( *this, "RemoveAllLayersFromView", 1, "viewID",
			 "Remove all layers from a view." );
  commandMgr.AddCommand( *this, "RemoveLayerFromViewAtLevel", 2, 
			 "viewID layer",
			 "Remove a layer from a view." );
  commandMgr.AddCommand( *this, "GetLabelValuesSet", 2, "viewID setName",
			 "Get a set of label value pairs." );
  commandMgr.AddCommand( *this, "GetFirstUnusedDrawLevelInView", 1, "viewID",
			 "Returns the first unused draw level." );
  commandMgr.AddCommand( *this, "SetViewLinkedStatus", 2, "viewID linked",
			 "Set the linked status for a view." );
  commandMgr.AddCommand( *this, "GetViewLinkedStatus", 1, "viewID",
			 "Returns the linked status for a view." );
  

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

  // If we're linked, set the center in other linked views.
  if( mViewIDLinkedList[GetID()] && mCurrentBroadcaster == -1 ) {

    mCurrentBroadcaster = GetID();
    map<int,bool>::iterator tIDLinked;
    for( tIDLinked = mViewIDLinkedList.begin();
	 tIDLinked != mViewIDLinkedList.end(); ++tIDLinked ) {
      int viewID = (*tIDLinked).first;
      bool bLinked = (*tIDLinked).second;
      if( bLinked && GetID() != viewID ) {
	try { 
	  View& view = View::FindByID( viewID );
	  // ScubaView& scubaView = dynamic_cast<ScubaView&>(view);
	  ScubaView& scubaView = (ScubaView&)view;
	  scubaView.Set2DRASCenter( iRASCenter );
	  scubaView.RebuildOverlayDrawList();
	  scubaView.RequestRedisplay();
	}
	catch(...) {
	  DebugOutput( << "Couldn't find view " << viewID );
	}
      }
    }
    mCurrentBroadcaster = -1;
  }
  
  // Changed our view center, so we need to rebuild the overlay.
  RebuildOverlayDrawList();
  RequestRedisplay();
}

void
ScubaView::Set2DZoomLevel ( float iZoomLevel ) {

  mViewState.mZoomLevel = iZoomLevel;

  // If we're linked, set the center in other linked views.
  if( mViewIDLinkedList[GetID()] && mCurrentBroadcaster == -1 ) {

    mCurrentBroadcaster = GetID();
    map<int,bool>::iterator tIDLinked;
    for( tIDLinked = mViewIDLinkedList.begin();
	 tIDLinked != mViewIDLinkedList.end(); ++tIDLinked ) {
      int viewID = (*tIDLinked).first;
      bool bLinked = (*tIDLinked).second;
      if( bLinked && GetID() != viewID ) {
	try { 
	  View& view = View::FindByID( viewID );
	  // ScubaView& scubaView = dynamic_cast<ScubaView&>(view);
	  ScubaView& scubaView = (ScubaView&)view;
	  scubaView.Set2DZoomLevel( iZoomLevel );
	  scubaView.RebuildOverlayDrawList();
	  scubaView.RequestRedisplay();
	}
	catch(...) {
	  DebugOutput( << "Couldn't find view " << viewID );
	}
      }
    }
    mCurrentBroadcaster = -1;
  }
  
  // Changed zoom, so we need to rebuild the overlay.
  RebuildOverlayDrawList();
  RequestRedisplay();
}

void
ScubaView::Set2DInPlane ( ViewState::Plane iPlane ) {

  mViewState.mInPlane = iPlane;
}

void 
ScubaView::SetLayerAtLevel ( int iLayerID, int iLevel ) {

  // If we got layer ID -1, remove this layer. Otherwise set it in our
  // map.
  if( iLayerID == -1 ) {

    RemoveLayerAtLevel( iLevel );
    
  } else {
    
    try {
      Layer& layer = Layer::FindByID( iLayerID );
      
      mLevelLayerIDMap[iLevel] = iLayerID;

      // Set the layer's width and height.
      layer.SetWidth( mWidth );
      layer.SetHeight( mHeight );
    }
    catch(...) {
      DebugOutput( << "Couldn't find layer " << iLayerID );
    }
  }

  RequestRedisplay();
}

int
ScubaView::GetLayerAtLevel ( int iLevel ) {

  // Look for this level in the level layer ID map. If we found
  // it, return the layer ID, otherwise return -1.
  map<int,int>::iterator tLevelLayerID = mLevelLayerIDMap.find( iLevel );
  if( tLevelLayerID == mLevelLayerIDMap.end() ) {
    return -1;
  } else {
    return mLevelLayerIDMap[iLevel];
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
ScubaView::CopyLayerSettingsToView ( ScubaView& iView ) {

  map<int,int>::iterator tLevelLayerID;
  for( tLevelLayerID = mLevelLayerIDMap.begin(); 
       tLevelLayerID != mLevelLayerIDMap.end(); ++tLevelLayerID ) {
    
    int level = (*tLevelLayerID).first;
    int layerID = (*tLevelLayerID).second;
    iView.SetLayerAtLevel( level, layerID );

  }
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
      if( 0 == strcmp( iasArgv[2], "X" ) || 
	  0 == strcmp( iasArgv[2], "x" ) ||
	  0 == strcmp( iasArgv[2], "R" ) || 
	  0 == strcmp( iasArgv[2], "r" ) ) {
	inPlane = ViewState::X;
      } else if( 0 == strcmp( iasArgv[2], "Y" ) || 
		 0 == strcmp( iasArgv[2], "y" ) ||
		 0 == strcmp( iasArgv[2], "A" ) || 
		 0 == strcmp( iasArgv[2], "a" ) ) {
	inPlane = ViewState::Y;
      } else if( 0 == strcmp( iasArgv[2], "Z" ) || 
		 0 == strcmp( iasArgv[2], "z" ) ||
		 0 == strcmp( iasArgv[2], "S" ) || 
		 0 == strcmp( iasArgv[2], "s" )) {
	inPlane = ViewState::Z;
      } else {
	sResult = "bad inPlane \"" + string(iasArgv[2]) + 
	  "\", should be x, y, or z, or r, a, or s";
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

  // SetLayerInViewAtLevel <viewID> <layerID> <level>
  if( 0 == strcmp( isCommand, "SetLayerInViewAtLevel" ) ) {
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
      
      SetLayerAtLevel( layerID, level );
    }
  }

  // GetLayerInViewAtLevel <viewID> <level>
  if( 0 == strcmp( isCommand, "GetLayerInViewAtLevel" ) ) {
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
      
      stringstream ssReturn;
      sReturnFormat = "i";
      ssReturn << GetLayerAtLevel( level );
      sReturnValues = ssReturn.str();
      return ok;
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

  // GetFirstUnusedDrawLevelInView <viewID>
  if( 0 == strcmp( isCommand, "GetFirstUnusedDrawLevelInView" ) ) {
    int viewID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad view ID";
      return error;
    }
    
    if( mID == viewID ) {

      sReturnFormat = "i";
      stringstream ssReturnValues;
      ssReturnValues << GetFirstUnusedDrawLevel();
      sReturnValues = ssReturnValues.str();
    }
  }


  // SetViewLinkedStatus <viewID> <linked>
  if( 0 == strcmp( isCommand, "SetViewLinkedStatus" ) ) {
    int viewID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad view ID";
      return error;
    }
    
    if( mID == viewID ) {
      
      if( 0 == strcmp( iasArgv[2], "true" ) || 
	  0 == strcmp( iasArgv[2], "1" )) {
	SetLinkedStatus( true );
      } else if( 0 == strcmp( iasArgv[2], "false" ) ||
		 0 == strcmp( iasArgv[2], "0" ) ) {
	SetLinkedStatus( false );
      } else {
	sResult = "bad linkedStatus \"" + string(iasArgv[2]) +
	  "\", should be true, 1, false, or 0";
	return error;	
      }
    }
  }

  // GetViewLinkedStatus <viewID>
  if( 0 == strcmp( isCommand, "GetViewLinkedStatus" ) ) {
    int viewID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad view ID";
      return error;
    }
    
    if( mID == viewID ) {

      sReturnFormat = "i";
      stringstream ssReturnValues;
      ssReturnValues << (int)GetLinkedStatus();
      sReturnValues = ssReturnValues.str();
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
ScubaView::DoMouseMoved( int iWindow[2], 
			 InputState& iInput, ScubaToolState& iTool ) {

  float ras[3];
  TranslateWindowToRAS( iWindow, ras );

  map<string,string> labelValueMap;

  stringstream sID;
  sID << GetLabel() << " (" << GetID() << ")";
  labelValueMap["View"] = sID.str();

  // Get the RAS coords into a string and set that label/value.
  stringstream ssRASCoords;
  ssRASCoords << ras[0] << " " << ras[1] << " " << ras[2];
  labelValueMap["RAS"] = ssRASCoords.str();

  // Go through our draw levels. For each one, get the Layer.
  map<int,int>::iterator tLevelLayerID;
  for( tLevelLayerID = mLevelLayerIDMap.begin(); 
       tLevelLayerID != mLevelLayerIDMap.end(); ++tLevelLayerID ) {

    int layerID = (*tLevelLayerID).second;
    try {
      Layer& layer = Layer::FindByID( layerID );
      
      // Ask the layer for info strings at this point.
      layer.GetInfoAtRAS( ras, labelValueMap );
    }
    catch(...) {
      DebugOutput( << "Couldn't find layer " << layerID );
    }
  }

  // Set this labelValueMap in the array of label values under the
  // name 'cursor'.
  mLabelValueMaps["cursor"] = labelValueMap;


  // Handle the navigation tool.
  if( iTool.GetMode() == ScubaToolState::navigation ) {

    if( iInput.Button() && !iInput.IsControlKeyDown() ) {
      
      float delta[2];
      delta[0] = (float)(mLastMouseMoved[0] - iWindow[0]) / 
	mViewState.mZoomLevel;
      delta[1] = (float)(mLastMouseMoved[1] - iWindow[1]) /
	mViewState.mZoomLevel;
      
      /* add to the total delta */
      mMouseMoveDelta[0] += delta[0];
      mMouseMoveDelta[1] += delta[1];
      
      /* save this mouse position */
      mLastMouseMoved[0] = iWindow[0];
      mLastMouseMoved[1] = iWindow[1];
      
      float moveLeftRight = 0, moveUpDown = 0, moveInOut = 0, zoomInOut = 0;
      switch( iInput.Button() ) {
      case 1:
	moveLeftRight = mMouseMoveDelta[0];
	moveUpDown    = mMouseMoveDelta[1];
	break;
      case 2:
	moveInOut = mMouseMoveDelta[1];
	break;
      case 3:
	zoomInOut = mMouseMoveDelta[1] / 10.0;
	break;
      }
      
      if( moveLeftRight || moveUpDown || moveInOut ) {
	
	float newRAS[3];
	switch( mViewState.mInPlane ) {
	case ViewState::X:
	  newRAS[0] = mOriginalCenterRAS[0] + moveInOut;
	  newRAS[1] = mOriginalCenterRAS[1] + moveLeftRight;
	  newRAS[2] = mOriginalCenterRAS[2] + moveUpDown;
	  break;
	case ViewState::Y:
	  newRAS[0] = mOriginalCenterRAS[0] + moveLeftRight;
	  newRAS[1] = mOriginalCenterRAS[1] + moveInOut;
	  newRAS[2] = mOriginalCenterRAS[2] + moveUpDown;
	  break;
	case ViewState::Z:
	  newRAS[0] = mOriginalCenterRAS[0] + moveLeftRight;
	  newRAS[1] = mOriginalCenterRAS[1] + moveUpDown;
	  newRAS[2] = mOriginalCenterRAS[2] + moveInOut;
	  break;
	}
	Set2DRASCenter( newRAS );
      }
      
      if( zoomInOut ) {
	
	float newZoom = mOriginalZoom + zoomInOut;
	if( newZoom <= 0.25 ) {
	  newZoom = 0.25;
	}
	Set2DZoomLevel( newZoom );
      }
    }
  }
  
  // Pass this tool to our layers.
  for( tLevelLayerID = mLevelLayerIDMap.begin(); 
       tLevelLayerID != mLevelLayerIDMap.end(); ++tLevelLayerID ) {
    int layerID = (*tLevelLayerID).second;
    try {
      Layer& layer = Layer::FindByID( layerID );
      layer.HandleTool( ras, iTool, iInput );
      if( layer.WantRedisplay() ) {
	RequestRedisplay();
	layer.RedisplayPosted();
      }
    }
    catch(...) {
      DebugOutput( << "Couldn't find layer " << layerID );
    }
  }
}

void
ScubaView::DoMouseUp( int iWindow[2],
		      InputState& iInput, ScubaToolState& iTool ) {

  // No matter what tool we're on, look for ctrl-b{1,2,3} and do some
  // navigation stuff.
  if( iInput.IsControlKeyDown() && 
      !iInput.IsShiftKeyDown() && !iInput.IsAltKeyDown() ) {
    
    // Set the new view center to this point. If they also hit b1 or
    // b3, zoom in or out accordingly.
    float world[3];
    TranslateWindowToRAS( iWindow, world );
    Set2DRASCenter( world );

    switch( iInput.Button() ) {
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


  // Pass this tool to our layers.
  float ras[3];
  TranslateWindowToRAS( iWindow, ras );
  map<int,int>::iterator tLevelLayerID;
  for( tLevelLayerID = mLevelLayerIDMap.begin(); 
       tLevelLayerID != mLevelLayerIDMap.end(); ++tLevelLayerID ) {
    int layerID = (*tLevelLayerID).second;
    try {
      Layer& layer = Layer::FindByID( layerID );
      layer.HandleTool( ras, iTool, iInput );
      if( layer.WantRedisplay() ) {
	RequestRedisplay();
	layer.RedisplayPosted();
      }
    }
    catch(...) {
      DebugOutput( << "Couldn't find layer " << layerID );
    }
  }
}

void
ScubaView::DoMouseDown( int iWindow[2], 
			InputState& iInput, ScubaToolState& iTool ) {

  mLastMouseDown[0] = mLastMouseMoved[0] = iWindow[0];
  mLastMouseDown[1] = mLastMouseMoved[1] = iWindow[1];
  mMouseMoveDelta[0] = 0.0;
  mMouseMoveDelta[1] = 0.0;
  mOriginalCenterRAS[0] = mViewState.mCenterRAS[0];
  mOriginalCenterRAS[1] = mViewState.mCenterRAS[1];
  mOriginalCenterRAS[2] = mViewState.mCenterRAS[2];
  mOriginalZoom = mViewState.mZoomLevel;

  // Pass this tool to our layers.
  float ras[3];
  TranslateWindowToRAS( iWindow, ras );
  map<int,int>::iterator tLevelLayerID;
  for( tLevelLayerID = mLevelLayerIDMap.begin(); 
       tLevelLayerID != mLevelLayerIDMap.end(); ++tLevelLayerID ) {
    int layerID = (*tLevelLayerID).second;
    try {
      Layer& layer = Layer::FindByID( layerID );
      layer.HandleTool( ras, iTool, iInput );
      if( layer.WantRedisplay() ) {
	RequestRedisplay();
	layer.RedisplayPosted();
      }
    }
    catch(...) {
      DebugOutput( << "Couldn't find layer " << layerID );
    }
  }
}

void
ScubaView::DoKeyDown( int iWindow[2], 
		      InputState& iInput, ScubaToolState& iTool ) {

  float moveDistance = 1.0;
  if( iInput.IsControlKeyDown() ) {
    moveDistance = 10.0;
  }

  string key = iInput.Key();
  

  if( key == msMoveViewLeft || key == msMoveViewRight ||
      key == msMoveViewDown || key == msMoveViewUp ||
      key == msMoveViewIn   || key == msMoveViewOut ) {
    
    float newRAS[3];
    newRAS[0] = mViewState.mCenterRAS[0];
    newRAS[1] = mViewState.mCenterRAS[1];
    newRAS[2] = mViewState.mCenterRAS[2];

    if( key == msMoveViewLeft ) {
      switch( mViewState.mInPlane ) {
      case ViewState::X: 
	newRAS[1] = mViewState.mCenterRAS[1] - moveDistance; break;
      case ViewState::Y: 
	newRAS[0] = mViewState.mCenterRAS[0] - moveDistance; break;
      case ViewState::Z: 
	newRAS[0] = mViewState.mCenterRAS[0] - moveDistance; break;
      }
    } else if( key == msMoveViewRight ) {
      switch( mViewState.mInPlane ) {
      case ViewState::X: 
	newRAS[1] = mViewState.mCenterRAS[1] + moveDistance; break;
      case ViewState::Y: 
	newRAS[0] = mViewState.mCenterRAS[0] + moveDistance; break;
      case ViewState::Z: 
	newRAS[0] = mViewState.mCenterRAS[0] + moveDistance; break;
      }
    } else if( key == msMoveViewDown ) {
      switch( mViewState.mInPlane ) {
      case ViewState::X: 
	newRAS[2] = mViewState.mCenterRAS[2] - moveDistance; break;
      case ViewState::Y: 
	newRAS[2] = mViewState.mCenterRAS[2] - moveDistance; break;
      case ViewState::Z: 
	newRAS[1] = mViewState.mCenterRAS[1] - moveDistance; break;
      }
    } else if( key == msMoveViewUp ) {
      switch( mViewState.mInPlane ) {
      case ViewState::X: 
	newRAS[2] = mViewState.mCenterRAS[2] + moveDistance; break;
      case ViewState::Y: 
	newRAS[2] = mViewState.mCenterRAS[2] + moveDistance; break;
      case ViewState::Z: 
	newRAS[1] = mViewState.mCenterRAS[1] + moveDistance; break;
      }
    } else if( key == msMoveViewIn ) {
      switch( mViewState.mInPlane ) {
      case ViewState::X: 
	newRAS[0] = mViewState.mCenterRAS[0] + moveDistance; break;
      case ViewState::Y: 
	newRAS[1] = mViewState.mCenterRAS[1] + moveDistance; break;
      case ViewState::Z: 
	newRAS[2] = mViewState.mCenterRAS[2] + moveDistance; break;
      }
    } else if( key == msMoveViewOut ) {
      switch( mViewState.mInPlane ) {
      case ViewState::X: 
	newRAS[0] = mViewState.mCenterRAS[0] - moveDistance; break;
      case ViewState::Y: 
	newRAS[1] = mViewState.mCenterRAS[1] - moveDistance; break;
      case ViewState::Z: 
	newRAS[2] = mViewState.mCenterRAS[2] - moveDistance; break;
      }
    }
    Set2DRASCenter( newRAS );

  } else if( key == msZoomViewIn || key == msZoomViewOut ) {

    float newZoom;
    if( key == msZoomViewIn ) {
      newZoom = mViewState.mZoomLevel * 2.0;
    } else if( key == msZoomViewOut ) {
      newZoom = mViewState.mZoomLevel / 2.0;
      if( newZoom < 0.25 ) {
	newZoom = 0.25;
      }
    }
    Set2DZoomLevel( newZoom );

  } else if( key == msInPlaneXKey ) {
    mViewState.mInPlane = ViewState::X;
    RebuildOverlayDrawList();
  } else if( key == msInPlaneYKey ) {
    mViewState.mInPlane = ViewState::Y;
    RebuildOverlayDrawList();
  } else if( key == msInPlaneZKey ) {
    mViewState.mInPlane = ViewState::Z;
    RebuildOverlayDrawList();
  }

  RequestRedisplay();
}

void
ScubaView::DoKeyUp( int iWindow[2], 
		    InputState& iInput, ScubaToolState& iTool ) {

}

void
ScubaView::TranslateWindowToRAS ( int iWindow[2], float oRAS[3] ) {

  // At zoom level one every pixel is an RAS point, so we're scaled by
  // that and then offset by the RAS center. We find a 3D point by
  // using our mInPlane and the corresponding mCenterRAS in ViewState.

  float viewRAS[3];
  switch( mViewState.mInPlane ) {
  case ViewState::X:
    viewRAS[0] = mViewState.mCenterRAS[0];
    viewRAS[1] = ConvertWindowToRAS( iWindow[0],mViewState.mCenterRAS[1],mWidth );
    viewRAS[2] = ConvertWindowToRAS(iWindow[1],mViewState.mCenterRAS[2],mHeight );
    break;
  case ViewState::Y:
    viewRAS[0] = ConvertWindowToRAS( iWindow[0],mViewState.mCenterRAS[0],mWidth );
    viewRAS[1] = mViewState.mCenterRAS[1];
    viewRAS[2] = ConvertWindowToRAS(iWindow[1],mViewState.mCenterRAS[2],mHeight );
    break;
  case ViewState::Z:
    viewRAS[0] = ConvertWindowToRAS( iWindow[0],mViewState.mCenterRAS[0],mWidth );
    viewRAS[1] = ConvertWindowToRAS(iWindow[1],mViewState.mCenterRAS[1],mHeight );
    viewRAS[2] = mViewState.mCenterRAS[2];
    break;
  }

  VECTOR_ELT( mViewRAS, 1 ) = viewRAS[0];
  VECTOR_ELT( mViewRAS, 2 ) = viewRAS[1];
  VECTOR_ELT( mViewRAS, 3 ) = viewRAS[2];
  VECTOR_ELT( mViewRAS, 4 ) = 1.0;
  MatrixMultiply( mViewToWorld, mViewRAS, mWorldRAS );
  oRAS[0] = VECTOR_ELT( mWorldRAS, 1 );
  oRAS[1] = VECTOR_ELT( mWorldRAS, 2 );
  oRAS[2] = VECTOR_ELT( mWorldRAS, 3 );
}

float
ScubaView::ConvertWindowToRAS ( float iWindow, float iRASCenter, 
				float iWindowDimension ) {

  return (iWindow / mViewState.mZoomLevel) +
    iRASCenter - (iWindowDimension / 2.0 / mViewState.mZoomLevel);
}

void
ScubaView::TranslateRASToWindow ( float iRAS[3], int oWindow[2] ) {
  
  float viewRAS[3];
  VECTOR_ELT( mWorldRAS, 1 ) = iRAS[0];
  VECTOR_ELT( mWorldRAS, 2 ) = iRAS[1];
  VECTOR_ELT( mWorldRAS, 3 ) = iRAS[2];
  VECTOR_ELT( mWorldRAS, 4 ) = 1.0;
  MatrixMultiply( mWorldToView, mWorldRAS, mViewRAS );
  viewRAS[0] = VECTOR_ELT( mViewRAS, 1 );
  viewRAS[1] = VECTOR_ELT( mViewRAS, 2 );
  viewRAS[2] = VECTOR_ELT( mViewRAS, 3 );

  float xWindow, yWindow;
  switch( mViewState.mInPlane ) {
  case ViewState::X:
    xWindow = ConvertRASToWindow( viewRAS[1],
				  mViewState.mCenterRAS[1], mWidth );
    yWindow = ConvertRASToWindow( viewRAS[2],
				  mViewState.mCenterRAS[2], mHeight );
    break;
  case ViewState::Y:
    xWindow = ConvertRASToWindow( viewRAS[0],
				  mViewState.mCenterRAS[0], mWidth );
    yWindow = ConvertRASToWindow( viewRAS[2],
				  mViewState.mCenterRAS[2], mHeight );
    break;
  case ViewState::Z:
    xWindow = ConvertRASToWindow( viewRAS[0],
				  mViewState.mCenterRAS[0], mWidth );
    yWindow = ConvertRASToWindow( viewRAS[1],
				  mViewState.mCenterRAS[1], mHeight );
    break;
  }

  oWindow[0] = (int) rint( xWindow );
  oWindow[1] = (int) rint( yWindow );
}

float
ScubaView::ConvertRASToWindow ( float iRAS, float iRASCenter, 
				float iWindowDimension ) {

  return ((iRAS - iRASCenter) * mViewState.mZoomLevel) +
    (iWindowDimension / 2.0);
}

int
ScubaView::GetFirstUnusedDrawLevel () {

  int highestLevel = 0;
  map<int,int>::iterator tLevelLayerID;
  for( tLevelLayerID = mLevelLayerIDMap.begin(); 
       tLevelLayerID != mLevelLayerIDMap.end(); ++tLevelLayerID ) {

    int level = (*tLevelLayerID).first;
    int layerID = (*tLevelLayerID).second;
    
    highestLevel = level + 1;
  }

  return highestLevel;
}

void 
ScubaView::BuildFrameBuffer () {

  // Don't draw if our buffer isn't initialized yet.
  if( NULL == mBuffer )
    return;

  // Erase the frame buffer.
  bzero( mBuffer, mHeight * mWidth * kBytesPerPixel );

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
  char sXLabel, sYLabel, sZLabel;
  float left, right, top, bottom, plane;
  switch( mViewState.mInPlane ) {
    case ViewState::X: 
      sXLabel = 'a';
      left =   mViewState.mCenterRAS[1] - (mWidth/2 / mViewState.mZoomLevel);
      right =  mViewState.mCenterRAS[1] + (mWidth/2 / mViewState.mZoomLevel);
      sYLabel = 's';
      bottom = mViewState.mCenterRAS[2] - (mHeight/2 / mViewState.mZoomLevel);
      top =    mViewState.mCenterRAS[2] + (mHeight/2 / mViewState.mZoomLevel);
      sZLabel = 'r';
      plane =  mViewState.mCenterRAS[0];
      break;
    case ViewState::Y: 
      sXLabel = 'r';
      left =   mViewState.mCenterRAS[0] - (mWidth/2 / mViewState.mZoomLevel);
      right =  mViewState.mCenterRAS[0] + (mWidth/2 / mViewState.mZoomLevel);
      sYLabel = 's';
      bottom = mViewState.mCenterRAS[2] - (mHeight/2 / mViewState.mZoomLevel);
      top =    mViewState.mCenterRAS[2] + (mHeight/2 / mViewState.mZoomLevel);
      sZLabel = 'a';
      plane =  mViewState.mCenterRAS[1];
      break;
    case ViewState::Z: 
      sXLabel = 'r';
      left =   mViewState.mCenterRAS[0] - (mWidth/2 / mViewState.mZoomLevel);
      right =  mViewState.mCenterRAS[0] + (mWidth/2 / mViewState.mZoomLevel);
      sYLabel = 'a';
      bottom = mViewState.mCenterRAS[1] - (mHeight/2 / mViewState.mZoomLevel);
      top =    mViewState.mCenterRAS[1] + (mHeight/2 / mViewState.mZoomLevel);
      sZLabel = 's';
      plane =  mViewState.mCenterRAS[2];
      break;
  }


  glNewList( kOverlayDrawListID + mID, GL_COMPILE );

  glColor3f( 1, 1, 1 );

  char sLabel[60];
  sprintf( sLabel, "%c%.2f", sXLabel, left );
  glRasterPos2i( 0, mHeight / 2 );
  for( int nChar = 0; nChar < strlen( sLabel ); nChar++ ) {
    glutBitmapCharacter( GLUT_BITMAP_8_BY_13, sLabel[nChar] );
  }
  sprintf( sLabel, "%c%.2f", sXLabel, right );
  glRasterPos2i( mWidth - strlen(sLabel)*8, mHeight / 2 );
  for( int nChar = 0; nChar < strlen( sLabel ); nChar++ ) {
    glutBitmapCharacter( GLUT_BITMAP_8_BY_13, sLabel[nChar] );
  }
  sprintf( sLabel, "%c%.2f", sYLabel, top );
  glRasterPos2i( mWidth / 2 - (strlen(sLabel)*8 / 2), mHeight-1-13 );
  for( int nChar = 0; nChar < strlen( sLabel ); nChar++ ) {
    glutBitmapCharacter( GLUT_BITMAP_8_BY_13, sLabel[nChar] );
  }
  sprintf( sLabel, "%c%.2f", sYLabel, bottom );
  glRasterPos2i( mWidth / 2 - (strlen(sLabel)*8 / 2), 4 );
  for( int nChar = 0; nChar < strlen( sLabel ); nChar++ ) {
    glutBitmapCharacter( GLUT_BITMAP_8_BY_13, sLabel[nChar] );
  }
  sprintf( sLabel, "%c%.2f", sZLabel, plane );
  glRasterPos2i( mWidth - (strlen(sLabel)*8), 4 );
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
