#ifndef ViewState_h
#define ViewState_h

// View orientation information for ScubaViews.

class ViewState {
  public:
  enum Plane { X, Y, Z };

  ViewState();

  float*  GetCenterRAS ()    { return mCenterRAS; }
  void    GetCenterRAS ( float oCenterRAS[3] );
  float   GetZoomLevel ()    { return mZoomLevel; }
  float*  GetPlaneNormal ()  { return mPlaneNormal; }
  void    GetPlaneNormal ( float oPlaneNormal[3] );
  Plane   GetInPlane ()      { return mInPlane; }
  int     GetBufferWidth ()  { return mBufferWidth; }
  int     GetBufferHeight () { return mBufferHeight; }

  void SetCenterRAS ( float iCenterRASX, float iCenterRASY,float iCenterRASZ );
  void SetCenterRAS ( float iCenterRAS[3] );
  void SetZoomLevel ( float iZoomLevel );
  void SetPlaneNormal ( float iPlaneNormal[3] );
  void SetPlaneNormal ( float iPlaneNormalX, 
			float iPlaneNormalY, float iPlaneNormalZ );
  void SetInPlane ( Plane iInPlane );
  void SetBufferWidth ( int iWidth );
  void SetBufferHeight ( int iHeight );

  // Checks if an RAS coord is visible in the plane within a certain
  // range.
  bool IsRASVisibleInPlane ( float iRAS[3], float iRange );

  // For saving a view state and comparing it later, e.g. for seeing
  // if a cache associated with a particular view is still valid for
  // the current view.
  void Copy ( ViewState& ioViewState );
  bool IsSameAs ( ViewState& iViewState );

  // Update rectangle management. You can add indiviudal rects to the
  // update rect and it will resize to include it. Rectangles are in
  // buffer space. ResetUpdateRect resets it to 0, and
  // UpdateEntireRect sets it to the entire buffer.
  void ResetUpdateRect ();
  void AddUpdateRect ( int iWindowLeft,  int iWindowTop,
		       int iWindowRight, int iWindowBottom );
  void UpdateEntireRect ();

  void CopyUpdateRect ( int& oWindowLeft,  int& oWindowTop,
			int& oWindowRight, int& oWindowBottom );
  void CopyUpdateRect ( int oWindowUpdate[4] );

 protected:

  int mUpdateRect[4];

  float mCenterRAS[3];
  float mZoomLevel;            
  float mPlaneNormal[3];
  Plane mInPlane;

  int mBufferWidth;
  int mBufferHeight;
};

std::ostream& operator << ( std::ostream& os, ViewState& iInput );

#endif
