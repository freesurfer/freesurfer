#include <stdexcept>
#include "glut.h"
#include "ScubaView.h"

using namespace std;

int const ScubaView::kBytesPerPixel = 4;

ScubaView::ScubaView() {
  mBuffer = NULL;

  TclCommandManager& commandMgr = TclCommandManager::GetManager();
  commandMgr.AddCommand( *this, "SetViewInPlane" );
  commandMgr.AddCommand( *this, "SetViewZoomLevel" );
  commandMgr.AddCommand( *this, "SetViewRASCenter" );
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
}

void
ScubaView::DoListenToTclCommand( char* isCommand, int iArgc, char** iasArgv ) {

  // SetViewInPlane <frameID> <inPlane>
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

  // SetViewZoomLevel <frameID> <zoomLevel>
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

  // SetViewRASCenter <frameID> <X> <Y> <Z>
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

}

void
ScubaView::DoDraw() {

  // Don't draw if our buffer isn't initialized yet.
  if( NULL == mBuffer )
    return;

  // Go through our draw levels. For each one, get the Layer.
  map<int,int>::iterator tLevelLayerID;
  for( tLevelLayerID = mLevelLayerIDMap.begin(); 
       tLevelLayerID != mLevelLayerIDMap.end(); ++tLevelLayerID ) {

    int layerID = (*tLevelLayerID).second;
    Layer& layer = Layer::FindByID( layerID );
    
    // tell it to draw into our buffer with our view state information.
    layer.DrawIntoBuffer( mBuffer, mViewState );
  }
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
  GLbyte* newBuffer = (GLbyte*) malloc( mWidth * mHeight * kBytesPerPixel );
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
    Layer& layer = Layer::FindByID( layerID );
    
    // Ask the layer for info strings at this point.
    layer.GetInfoAtRAS( rasX, rasY, rasZ, labelValueMap );
  }

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
      glutBitmapCharacter( GLUT_BITMAP_HELVETICA_18, sLabel[nChar] );
    }
    for( int nChar = 0; nChar < sValue.length(); nChar++ ) {
      glutBitmapCharacter( GLUT_BITMAP_HELVETICA_18, sValue[nChar] );
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
