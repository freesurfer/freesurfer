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
  float mPlaneNormal[3];
  Plane mInPlane;

  int mBufferWidth;
  int mBufferHeight;

  void ResetUpdateRect ();
  void AddUpdateRect ( int iWindowLeft,  int iWindowTop,
		       int iWindowRight, int iWindowBottom );
  void UpdateEntireRect ();
  int mUpdateRect[4];

  void CopyUpdateRect ( int& oWindowLeft,  int& oWindowTop,
			int& oWindowRight, int& oWindowBottom );
  void CopyUpdateRect ( int oWindowUpdate[4] );
};

std::ostream& operator << ( std::ostream& os, ViewState& iInput );

#endif
