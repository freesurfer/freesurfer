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


