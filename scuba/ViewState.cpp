#include <iostream>
#include "ViewState.h"


ViewState::ViewState () {
  mCenterRAS[0] = mCenterRAS[1] = mCenterRAS[2] = 0;
  mZoomLevel = 1;
  mInPlane = X;
}

std::ostream& operator << ( std::ostream& os, ViewState& iState ) { 
  os << "ViewState CenterRAS: " << iState.mCenterRAS[0] << ", "
     << iState.mCenterRAS[1] << ", " << iState.mCenterRAS[2] 
     << " ZoomLevel: " << iState.mZoomLevel << " InPlane: " 
     << iState.mInPlane;
  return os;
}


