#include <iostream>
#include <math.h>
#include "ViewState.h"
#include "macros.h"


ViewState::ViewState () {
  mCenterRAS[0] = mCenterRAS[1] = mCenterRAS[2] = 0;
  mZoomLevel = 1;
  mInPlane = X;
  mPlaneNormal[0] = 1; mPlaneNormal[1] = 0; mPlaneNormal[2] = 0;
}

std::ostream& operator << ( std::ostream& os, ViewState& iInput ) { 
  os << "ViewState CenterRAS: " << iInput.mCenterRAS[0] << ", "
     << iInput.mCenterRAS[1] << ", " << iInput.mCenterRAS[2] 
     << " ZoomLevel: " << iInput.mZoomLevel << " InPlane: " 
     << iInput.mInPlane;
  return os;
}

bool
ViewState::IsRASVisibleInPlane ( float iRAS[3], float iRange ) {

  float rasCoord;
  float viewCoord;

  switch( mInPlane ) {
  case ViewState::X:
    rasCoord = iRAS[0];
    viewCoord = mCenterRAS[0];
    break;
  case ViewState::Y:
    rasCoord = iRAS[1];
    viewCoord = mCenterRAS[1];
    break;
  case ViewState::Z:
  default:
    rasCoord = iRAS[2];
    viewCoord = mCenterRAS[2];
    break;
  }

  return( fabs(viewCoord - rasCoord) <= iRange );
}

void
ViewState::ResetUpdateRect () {
  mUpdateRect[0] = -1;
  mUpdateRect[1] = -1;
  mUpdateRect[2] = -1;
  mUpdateRect[3] = -1;
}

void
ViewState::AddUpdateRect ( int iWindowLeft,  int iWindowTop,
			   int iWindowRight, int iWindowBottom ) {

  if( iWindowLeft != -1 ) {
    if( mUpdateRect[0] == -1 ) 
      mUpdateRect[0] = iWindowLeft;
    else
      mUpdateRect[0] = MIN( mUpdateRect[0], iWindowLeft );
  }

  if( iWindowTop != -1 ) {
    if( mUpdateRect[1] == -1 ) 
      mUpdateRect[1] = iWindowTop;
    else
      mUpdateRect[1] = MIN( mUpdateRect[1], iWindowTop );
  }

  if( iWindowRight != -1 ) {
    if( mUpdateRect[2] == -1 ) 
      mUpdateRect[2] = iWindowRight;
    else
      mUpdateRect[2] = MAX( mUpdateRect[2], iWindowRight );
  }

  if( iWindowBottom != -1 ) {
    if( mUpdateRect[3] == -1 ) 
      mUpdateRect[3] = iWindowBottom;
    else
      mUpdateRect[3] = MAX( mUpdateRect[3], iWindowBottom );
  }
}

void
ViewState::UpdateEntireRect () {
  mUpdateRect[0] = 0;
  mUpdateRect[1] = 0;
  mUpdateRect[2] = mBufferWidth;
  mUpdateRect[3] = mBufferHeight;
}
