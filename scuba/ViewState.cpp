#include <iostream>
#include "ViewState.h"


ViewState::ViewState () {
  mCenterRAS[0] = mCenterRAS[1] = mCenterRAS[2] = 0;
  mZoomLevel = 1;
  mInPlane = X;
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
    rasCoord = iRAS[2];
    viewCoord = mCenterRAS[2];
    break;
  }

  return( fabs(viewCoord - rasCoord) <= iRange );
}

