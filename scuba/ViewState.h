#ifndef ViewState_h
#define ViewState_h


// View orientation information for ScubaViews.

class ViewState {
  public:
  enum Plane { X, Y, Z };

  ViewState();

  // Checks if an RAS coord is visible in the plane within a certain
  // range.
  bool IsRASVisibleInPlane ( float iRAS[3], float iRange );

  float mCenterRAS[3];
  float mZoomLevel;            
  Plane mInPlane;
};

std::ostream& operator << ( std::ostream& os, ViewState& iInput );

#endif
