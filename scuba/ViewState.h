#ifndef ViewState_h
#define ViewState_h


// View orientation information for ScubaViews.

class ViewState {
  public:
  enum Plane { X, Y, Z };

  float mCenterRAS[3];
  float mZoomLevel;            
  Plane mInPlane;
};

std::ostream& operator << ( std::ostream& os, ViewState& iState );

#endif
