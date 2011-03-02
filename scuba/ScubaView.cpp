/**
 * @file  ScubaView.cpp
 * @brief Scuba specific View that manages Layers
 *
 * The ScubaView has draw slots which hold Layers, and composes the
 * graphical contents of those layers in its graphics context with an
 * optional view transform. ScubaViews have their own 2D draw state,
 * and tell Layers what poarts of themselves to draw. ScubaViews also
 * have markers and a shared cursor, and can 'lock on' the cursor to
 * synchronize view states between Views. The View also passes events
 * down to Layers.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:39 $
 *    $Revision: 1.127 $
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
 *
 * Terms and conditions for use, reproduction, distribution and contribution
 * are found in the 'FreeSurfer Software License Agreement' contained
 * in the file 'LICENSE' found in the FreeSurfer distribution, and here:
 *
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
 *
 * Reporting: freesurfer@nmr.mgh.harvard.edu
 *
 */


#include "string_fixed.h"
#include <stdlib.h>
#include <math.h>
#include <stdexcept>
#include <iomanip>
#include "ScubaView.h"
#include "PreferencesManager.h"
#include "ScubaGlobalPreferences.h"
#include "Timer.h"
#include "Point2.h"
#include "Point3.h"
#include "VectorOps.h"
#include "Utilities.h"
#include "PathManager.h"

using namespace std;

int const ScubaView::kBytesPerPixel = 4;
map<int,bool> ScubaView::mViewIDLinkedList;
Point3<float> ScubaView::mCursor( 0, 0, 0 );
int ScubaView::mcMarkers = 25000;
int ScubaView::mNextMarker = 0;
std::map<int,Point3<float> > ScubaView::mMarkerRAS;
std::map<int,bool> ScubaView::mMarkerVisible;
ScubaViewStaticTclListener ScubaView::mStaticListener;
bool const ScubaView::kbDefaultLevelReportInfo = false;

GLenum glError;
#define CheckGLError()  \
  glError = glGetError(); \
  while( glError != GL_NO_ERROR ) { \
    cerr << __LINE__ << " error: " << gluErrorString( glError ) << endl; \
    glError = glGetError(); \
  } \


int const ScubaView::kcInPlaneMarkerColors = 12;
float kInPlaneMarkerColors[ScubaView::kcInPlaneMarkerColors][3] = {
      {
        1, 0, 0
      }
      , { 0, 1, 0 }, { 0, 0, 1 },
      { 1, 1, 0 }, { 1, 0, 1 }, { 0, 1, 1 },
      { 0.25, 0.75, 0 }, { 0.25, 0, 0.75 }, { 0, 0.25, 0.75 },
      { 0.75, 0.25, 0 }, { 0.75, 0, 0.25 }, { 0, 0.75, 0.25 }
    };

ScubaView::ScubaView() {
  mBuffer = NULL;
  mbPostRedisplay = false;
  mbRebuildOverlayDrawList = true;
  mViewIDLinkedList[GetID()] = false;
  mbFlipLeftRightInYZ = true;
  mbLockOnCursor = false;
  int nMarkerColor = GetID() % kcInPlaneMarkerColors;
  mInPlaneMarkerColor[0] = kInPlaneMarkerColors[nMarkerColor][0];
  mInPlaneMarkerColor[1] = kInPlaneMarkerColors[nMarkerColor][1];
  mInPlaneMarkerColor[2] = kInPlaneMarkerColors[nMarkerColor][2];
  mbVisibleInFrame = false;
  mCurrentMovingViewIntersection = -1;
  mThroughPlaneIncrements[0] =
    mThroughPlaneIncrements[1] = mThroughPlaneIncrements[2] = 1.0;
  mLastMouseOver.Set( 0, 0, 0 );

  ScubaGlobalPreferences& globalPrefs =
    ScubaGlobalPreferences::GetPreferences();
  globalPrefs.AddListener( *this );

  // Get the prefs for the lock on cursor status.
  mbLockOnCursor =
    globalPrefs.GetPrefAsBool( ScubaGlobalPreferences::LockOnCursor );

  ScubaViewBroadcaster& broadcaster = ScubaViewBroadcaster::GetBroadcaster();
  broadcaster.AddListener( *this );

  PathManager& pathMgr = PathManager::GetManager();
  pathMgr.AddListener( *this );

  // Try setting our initial transform to the default transform with
  // id 0. If it's not there, create it.
  try {
    mWorldToView = &(ScubaTransform::FindByID( 0 ));
    mWorldToView->AddListener( *this );
  } catch (...) {

    ScubaTransform* transform = new ScubaTransform();
    transform->SetLabel( "Identity" );

    try {
      mWorldToView = &(ScubaTransform::FindByID( 0 ));
      mWorldToView->AddListener( *this );
    } catch (...) {
      DebugOutput( << "Couldn't make default transform!" );
    }
  }

  TclCommandManager& commandMgr = TclCommandManager::GetManager();
  commandMgr.AddCommand( *this, "SetViewInPlane", 2, "viewID inPlane",
                         "Sets the in plane in a view. inPlane should be "
                         "one of the following: x y z" );
  commandMgr.AddCommand( *this, "GetViewInPlane", 1, "viewID",
                         "Returns the in plane in a view." );
  commandMgr.AddCommand( *this, "SetViewZoomLevel", 2, "viewID zoomLevel",
                         "Sets the zoom level in a view. zoomLevel should be "
                         "a float." );
  commandMgr.AddCommand( *this, "GetViewZoomLevel", 1, "viewID",
                         "Returns the zoom level in a view." );
  commandMgr.AddCommand( *this, "SetViewRASCenter", 4, "viewID x y z",
                         "Sets the view center. x, y, and z should be floats "
                         "in world RAS coordinates." );
  commandMgr.AddCommand( *this, "GetViewRASCenter", 1, "viewID",
                         "Returns the view center as a list of x, y, and z "
                         "RAS coordinates." );
  commandMgr.AddCommand( *this, "SetLayerInViewAtLevel", 3,
                         "viewID layerID level",
                         "Sets the layer in a view at a given draw level. "
                         "Higher draw levels will draw later." );
  commandMgr.AddCommand( *this, "GetLayerInViewAtLevel", 2, "viewID level",
                         "Returns the layer in a view at a given draw level.");
  commandMgr.AddCommand( *this, "RemoveAllLayersFromView", 1, "viewID",
                         "Remove all layers from a view." );
  commandMgr.AddCommand( *this, "RemoveLayerFromViewAtLevel", 2,
                         "viewID level",
                         "Remove a layer from a view." );
  commandMgr.AddCommand( *this, "GetListOfLevelsInView", 1, "viewID",
			 "Return a list of level indicies that are empty "
			 "but have been defined, or have layers in them." );
  commandMgr.AddCommand( *this, "SetLevelVisibilityInView", 3,
                         "viewID level visibility",
                         "Sets the visibility for a level in a view." );
  commandMgr.AddCommand( *this, "GetLevelVisibilityInView", 2,
                         "viewID level",
                         "Returns the visibility for a level in a view." );
  commandMgr.AddCommand( *this, "SetLevelReportInfoInView", 3,
                         "viewID level reportInfo",
                         "Sets the flag for reporting info for a level "
                         "in a view." );
  commandMgr.AddCommand( *this, "GetLevelReportInfoInView", 2,
                         "viewID level",
                         "Returns whether a level in a view is reporting "
                         "info." );
  commandMgr.AddCommand( *this, "SetViewStateToLayerBounds", 2,
                         "viewID, layerID", "Sets the view so that the "
                         "layer's data completely fills the view." );
  commandMgr.AddCommand( *this, "GetInfoAtRAS", 2, "viewID setName",
                         "Get an array list of info at an RAS point." );
  commandMgr.AddCommand( *this, "GetFirstUnusedDrawLevelInView", 1, "viewID",
                         "Returns the first unused draw level." );
  commandMgr.AddCommand( *this, "SetViewLinkedStatus", 2, "viewID linked",
                         "Set the linked status for a view." );
  commandMgr.AddCommand( *this, "GetViewLinkedStatus", 1, "viewID",
                         "Returns the linked status for a view." );
  commandMgr.AddCommand( *this, "SetViewLockOnCursor", 2, "viewID lock",
                         "Set a view to keep its view locked on the cursor." );
  commandMgr.AddCommand( *this, "GetViewLockOnCursor", 1, "viewID",
                         "Returns whether or not a view is locked on "
                         "the cursor." );
  commandMgr.AddCommand( *this, "SetViewTransform", 2, "viewID transformID",
                         "Set the view to world transform for a view." );
  commandMgr.AddCommand( *this, "GetViewTransform", 1, "viewID",
                         "Returns the transformID of a view's view to "
                         "world transform." );
  commandMgr.AddCommand( *this, "SetViewFlipLeftRightYZ", 2, "viewID flip",
                         "Set the left-right flip flag for a view." );
  commandMgr.AddCommand( *this, "GetViewFlipLeftRightYZ", 1, "viewID",
                         "Returns the left-right flip flag for a view." );
  commandMgr.AddCommand( *this, "SetViewThroughPlaneIncrement", 3,
                         "viewID throughPlane increment",
                         "Set the amount that using the through plane "
                         "movement keys will increment or decrement the "
                         "through plane RAS value. throughPlane should be "
                         "x, y, or z." );
  commandMgr.AddCommand( *this, "GetViewThroughPlaneIncrement", 2,
                         "viewID throughPlane",
                         "Returns the through plane movement increment "
                         "for throughPlane. throughPlane should be x, "
                         "y, or z." );
  commandMgr.AddCommand( *this, "GetVolumeHistogramInView", 4,
                         "viewID volID roiID numBins",
                         "Returns a histogram of the volume that's visible in "
                         "a view. Returns the format: minBinValue "
                         "binIncrement{binCount0 binCount1 .. binCountN} "
                         "where binCountN is numBins-1." );
  commandMgr.AddCommand( *this, "DoVolumeValueRangeFillInView", 5,
                         "viewID sourceVolID roiID destVolID valueRanges",
                         "Performs a series of value range fills using the "
                         "value ranges in sourceVolID to set values in "
                         "destVolID. valueRanges should be a list of triples: "
			 "the begining of the range, the end of the range, "
			 "and the new value." );
  commandMgr.AddCommand( *this, "ConvertWindowToViewRAS", 3,
                         "viewID windowX windowY",
                         "Returns the RAS coordinates of the input window "
                         "coordinate." );

  // Get some prefs values
  msMoveViewLeft  = ScubaKeyCombo::MakeKeyCombo();
  msMoveViewRight = ScubaKeyCombo::MakeKeyCombo();
  msMoveViewUp    = ScubaKeyCombo::MakeKeyCombo();
  msMoveViewDown  = ScubaKeyCombo::MakeKeyCombo();
  msMoveViewIn    = ScubaKeyCombo::MakeKeyCombo();
  msMoveViewOut   = ScubaKeyCombo::MakeKeyCombo();
  msZoomViewIn    = ScubaKeyCombo::MakeKeyCombo();
  msZoomViewOut   = ScubaKeyCombo::MakeKeyCombo();

  ScubaGlobalPreferences& prefs = ScubaGlobalPreferences::GetPreferences();
  msMoveViewLeft->SetFromString( prefs.GetPrefAsString( ScubaGlobalPreferences::KeyMoveViewLeft ) );
  msMoveViewRight->SetFromString( prefs.GetPrefAsString( ScubaGlobalPreferences::KeyMoveViewRight ) );
  msMoveViewUp->SetFromString( prefs.GetPrefAsString( ScubaGlobalPreferences::KeyMoveViewUp ) );
  msMoveViewDown->SetFromString( prefs.GetPrefAsString( ScubaGlobalPreferences::KeyMoveViewDown ) );
  msMoveViewIn->SetFromString( prefs.GetPrefAsString( ScubaGlobalPreferences::KeyMoveViewIn ) );
  msMoveViewOut->SetFromString( prefs.GetPrefAsString( ScubaGlobalPreferences::KeyMoveViewOut ) );
  msZoomViewIn->SetFromString( prefs.GetPrefAsString( ScubaGlobalPreferences::KeyZoomViewIn ) );
  msZoomViewOut->SetFromString( prefs.GetPrefAsString( ScubaGlobalPreferences::KeyZoomViewOut ) );
  mbFlipLeftRightInYZ =
    prefs.GetPrefAsBool( ScubaGlobalPreferences::ViewFlipLeftRight );

  list<Layer::InfoAtRAS> lInfo;
  mInfoAtRASMap["mouse"] = lInfo;
  mInfoAtRASMap["cursor"] = lInfo;

  // Inits label value lists to something decent, to create empty
  // filler space for display.
  float origin[3] = {0, 0, 0};
  RebuildLabelValueInfo( origin, "cursor" );
  RebuildLabelValueInfo( origin, "mouse" );

  // For some reason we can't use glGenLists yet because it will
  // return 0, an invalid index. So just set our draw list ID to 0 now
  // and we'll make one later.
  mDrawListID = 0;

  mViewState.ResetUpdateRect();
}

ScubaView::~ScubaView() {
  if ( NULL != mBuffer ) {
    free( mBuffer );
  }

  delete msMoveViewLeft;
  delete msMoveViewRight;
  delete msMoveViewUp;
  delete msMoveViewDown;
  delete msMoveViewIn;
  delete msMoveViewOut;
  delete msZoomViewIn;
  delete msZoomViewOut;

  // Stop listening.
  ScubaGlobalPreferences& globalPrefs =
    ScubaGlobalPreferences::GetPreferences();
  globalPrefs.RemoveListener( *this );

  ScubaViewBroadcaster& broadcaster = ScubaViewBroadcaster::GetBroadcaster();
  broadcaster.RemoveListener( *this );

  PathManager& pathMgr = PathManager::GetManager();
  pathMgr.RemoveListener( *this );

  if ( mWorldToView )
    mWorldToView->RemoveListener( *this );

  map<int,int>::iterator tLevelLayerID;
  for ( tLevelLayerID = mLevelLayerIDMap.begin();
        tLevelLayerID != mLevelLayerIDMap.end(); ++tLevelLayerID ) {
    int layerID = (*tLevelLayerID).second;
    try {
      Layer& layer = Layer::FindByID( layerID );
      layer.RemoveListener( *this );
    } catch (...) {}
  }
}

void
ScubaView::Set2DRASCenter ( float iRASCenter[3] ) {

  mViewState.SetCenterRAS( iRASCenter );

  // Broadcast this change.
  ScubaViewBroadcaster& broadcaster = ScubaViewBroadcaster::GetBroadcaster();
  broadcaster.SendBroadcast( "2DRASCenterChanged", (void*)&mID );

  CalcViewToWindowTransform();

  // Changed our view center, so we need to rebuild the overlay.
  RebuildOverlayDrawList();
  RequestRedisplay();
}

void
ScubaView::Set2DZoomLevel ( float iZoomLevel ) {

  mViewState.SetZoomLevel( iZoomLevel );

  // Broadcast this change.
  ScubaViewBroadcaster& broadcaster = ScubaViewBroadcaster::GetBroadcaster();
  broadcaster.SendBroadcast( "2DZoomLevelChanged", (void*)&mID );

  // Changed zoom, so we need to rebuild the overlay.
  RebuildOverlayDrawList();
  RequestRedisplay();
}

void
ScubaView::Set2DInPlane ( ViewState::Plane iPlane ) {

  // If we are going to a new plane, reset our plane normal.
  if ( mViewState.GetInPlane() != iPlane ) {
    switch ( iPlane ) {
    case ViewState::X:
      mViewState.SetPlaneNormal( 1, 0, 0 );
      break;
    case ViewState::Y:
      mViewState.SetPlaneNormal( 0, 1, 0 );
      break;
    case ViewState::Z:
      mViewState.SetPlaneNormal( 0, 0, 1 );
      break;
    }
  }
  mViewState.SetInPlane( iPlane );

  // Broadcast this change.
  ScubaViewBroadcaster& broadcaster = ScubaViewBroadcaster::GetBroadcaster();
  broadcaster.SendBroadcast( "2DInPlaneChanged", (void*)&mID );

  // If we're centering around the cursor and we can no longer see it,
  // change our plane to the cursor.
  if ( mbLockOnCursor &&
       !mViewState.IsRASVisibleInPlane( mCursor.xyz(),
                                        GetThroughPlaneIncrement(mViewState.GetInPlane()))) {

    float newCenter[3];
    mViewState.GetCenterRAS( newCenter );

    switch ( mViewState.GetInPlane() ) {
    case ViewState::X:
      newCenter[0] = mCursor[0];
      break;
    case ViewState::Y:
      newCenter[1] = mCursor[1];
      break;
    case ViewState::Z:
      newCenter[2] = mCursor[2];
      break;
    }
    Set2DRASCenter( newCenter );
  }

  CalcViewToWindowTransform();

  // We changed our orientation so we must recalc all the other views'
  // intersection points.
  CalcAllViewIntersectionPoints();

  // Changed in plane, so we need to rebuild the overlay.
  RebuildOverlayDrawList();
  RequestRedisplay();
}

void
ScubaView::Set2DPlaneNormalOrthogonalToInPlane () {

  // Reset our plane normal.
  float plane[3];
  switch ( mViewState.GetInPlane() ) {
  case ViewState::X:
    plane[0] = 1;
    plane[1] = 0;
    plane[2] = 0;
    break;
  case ViewState::Y:
    plane[0] = 0;
    plane[1] = 1;
    plane[2] = 0;
    break;
  case ViewState::Z:
    plane[0] = 0;
    plane[1] = 0;
    plane[2] = 1;
    break;
  }

  Set2DPlaneNormal( plane );
}

void
ScubaView::Set2DPlaneNormal ( float iNormal[3] ) {

  switch ( mViewState.GetInPlane() ) {
  case ViewState::X:
    if ( fabs(iNormal[0]) >= 0.5 ) {
      mViewState.SetPlaneNormal( iNormal );
    }
    break;
  case ViewState::Y:
    if ( fabs(iNormal[1]) >= 0.5 ) {
      mViewState.SetPlaneNormal( iNormal );
    }
    break;
  case ViewState::Z:
    if ( fabs(iNormal[2]) >= 0.5 ) {
      mViewState.SetPlaneNormal( iNormal );
    }
    break;
  }

  // Broadcast this change.
  ScubaViewBroadcaster& broadcaster = ScubaViewBroadcaster::GetBroadcaster();
  broadcaster.SendBroadcast( "2DPlaneNormalChanged", (void*)&mID );

  CalcViewToWindowTransform();

  // We changed our orientation so we must recalc all the other views'
  // intersection points.
  CalcAllViewIntersectionPoints();

  // Changed in plane, so we need to rebuild the overlay.
  RebuildOverlayDrawList();
  RequestRedisplay();
}

void
ScubaView::Get2DRASCenter ( float oRASCenter[3] ) {

  mViewState.GetCenterRAS( oRASCenter );
}

float
ScubaView::Get2DZoomLevel () {

  return mViewState.GetZoomLevel();
}

ViewState::Plane
ScubaView::Get2DInPlane () {

  return mViewState.GetInPlane();
}

string
ScubaView::Get2DInPlaneAsString () {
  switch ( mViewState.GetInPlane() ) {
  case ViewState::X:
    return "x";
    break;
  case ViewState::Y:
    return "y";
    break;
  case ViewState::Z:
    return "z";
    break;
  }
  return "Unknown";
}

void
ScubaView::Get2DPlaneNormal ( float oNormal[3] ) {

  mViewState.GetPlaneNormal( oNormal );
}


void
ScubaView::SetLayerAtLevel ( int iLayerID, int iLevel ) {

  // If we got layer ID -1, remove this layer. Otherwise set it in our
  // map.
  if ( iLayerID == -1 ) {

    RemoveLayerAtLevel( iLevel );

  } else {

    try {
      Layer& layer = Layer::FindByID( iLayerID );

      mLevelLayerIDMap[iLevel] = iLayerID;

      // Set the layer's width and height.
      layer.SetWidth( mWidth);
      layer.SetHeight( mHeight );

      // Set pixel size.
      layer.SetBytesPerPixel( kBytesPerPixel );

      // Start out visible.
      if ( mLevelVisibilityMap.find( iLevel ) == mLevelVisibilityMap.end() ) {
        mLevelVisibilityMap[iLevel] = true;
      }

      // Default reporting info value.
      if ( mLevelReportInfoMap.find( iLevel ) == mLevelReportInfoMap.end() ) {
        mLevelReportInfoMap[iLevel] = kbDefaultLevelReportInfo;
      }

      // Make a new display list for this level if necessary.
      if ( mLevelGLListIDMap.find( iLevel ) == mLevelGLListIDMap.end() ) {
        mLevelGLListIDMap[iLevel] = glGenLists( 1 );
      }

      // Listen to it.
      layer.AddListener( *this );

      // Rebuild our label value lists.
      float origin[3] = {0, 0, 0};
      RebuildLabelValueInfo( mCursor.xyz(), "cursor" );
      RebuildLabelValueInfo( origin, "mouse" );
    } catch (...) {
      DebugOutput( << "Couldn't find layer " << iLayerID );
    }
  }

  // Get the increments from the top level.
  int highestLevel = 0;
  map<int,int>::iterator tLevelLayerID;
  for ( tLevelLayerID = mLevelLayerIDMap.begin();
        tLevelLayerID != mLevelLayerIDMap.end(); ++tLevelLayerID ) {
    int level = (*tLevelLayerID).first;
    if ( level > highestLevel )
      highestLevel = level;
  }
  Layer& layer = Layer::FindByID( mLevelLayerIDMap[highestLevel] );
  layer.GetPreferredThroughPlaneIncrements( mThroughPlaneIncrements );

  RequestRedisplay();
}

int
ScubaView::GetLayerAtLevel ( int iLevel ) {

  // Look for this level in the level layer ID map. If we found
  // it, return the layer ID, otherwise return -1.
  map<int,int>::iterator tLevelLayerID = mLevelLayerIDMap.find( iLevel );
  if ( tLevelLayerID == mLevelLayerIDMap.end() ) {
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

int
ScubaView::GetFirstUnusedDrawLevel () {

  // This is actually the lowest empty level,
  int lowestEmptyLevel = 0;
  map<int,int>::iterator tLevelLayerID;
  for ( tLevelLayerID = mLevelLayerIDMap.begin();
        tLevelLayerID != mLevelLayerIDMap.end(); ++tLevelLayerID ) {

    int level = (*tLevelLayerID).first;

    // If the level above this one is empty,
    if ( mLevelLayerIDMap.find(level+1) == mLevelLayerIDMap.end () ||
         mLevelLayerIDMap[level+1] == -1 ) {

      // And if the level above this one is the lowest empty level,
      if ( level + 1 > lowestEmptyLevel ) {

        // Save this level.
        lowestEmptyLevel = level + 1;
      }
    }
  }

  return lowestEmptyLevel;
}

void
ScubaView::SetDrawLevelVisibility ( int inLevel, bool ibVisible ) {

  if ( mLevelVisibilityMap[inLevel] != ibVisible ) {

    mLevelVisibilityMap[inLevel] = ibVisible;

    // This affects what label information we display, so rebuild the
    // label value info.
    RebuildLabelValueInfo( mCursor.xyz(), "cursor" );
    RebuildLabelValueInfo( mCursor.xyz(), "mouse" );

    RequestRedisplay();
  }
}

bool
ScubaView::GetDrawLevelVisibility ( int inLevel ) {

  if ( mLevelVisibilityMap.find( inLevel ) == mLevelVisibilityMap.end() )
    mLevelVisibilityMap[inLevel] = true;

  return mLevelVisibilityMap[inLevel];
}

void
ScubaView::SetDrawLevelReportInfo ( int inLevel, bool ibReportInfo ) {

  if ( mLevelReportInfoMap[inLevel] != ibReportInfo ) {

    mLevelReportInfoMap[inLevel] = ibReportInfo;

    // This affects what label information we display, so rebuild the
    // label value info.
    RebuildLabelValueInfo( mCursor.xyz(), "cursor" );
    RebuildLabelValueInfo( mCursor.xyz(), "mouse" );

    RequestRedisplay();
  }
}

bool
ScubaView::GetDrawLevelReportInfo ( int inLevel ) {

  // If we don't have a setting for this level yet, set it to our
  // default.
  if ( mLevelReportInfoMap.find( inLevel ) == mLevelReportInfoMap.end() )
    mLevelReportInfoMap[inLevel] = kbDefaultLevelReportInfo;

  return mLevelReportInfoMap[inLevel];
}

void
ScubaView::SetViewStateToLayerBounds ( int iLayerID ) {

  try {
    // Get the layer.
    Layer& layer = Layer::FindByID( iLayerID );

    // Get the main data collection and get its RAS bounds.
    DataCollection* col = layer.GetMainDataCollection();
    if ( NULL != col ) {

      // Get the bounds.
      float RAS[6];
      col->GetDataRASBounds( RAS );

      // Calc the center.
      float centerRAS[3];
      centerRAS[0] = ((RAS[1] - RAS[0]) / 2.0) + RAS[0];
      centerRAS[1] = ((RAS[3] - RAS[2]) / 2.0) + RAS[2];
      centerRAS[2] = ((RAS[5] - RAS[4]) / 2.0) + RAS[4];

      // If this is a volume layer, this unfortunately could put us
      // right at the edge of a slice. This will mess up some of our
      // intersection stuff because now two slices of voxels intersect
      // with our plane. So adjust the coords by half a voxel.
      if ( col->GetTypeDescription() == "Volume" ) {
        VolumeCollection* vol = (VolumeCollection*) col;
        float voxelZero[3], voxelOne[3];
        voxelZero[0] = voxelZero[1] = voxelZero[2] = 0;
        voxelOne[0] = voxelOne[1] = voxelOne[2] = 1;
        float RASVoxelZero[3], RASVoxelOne[3];
        vol->MRIIndexToRAS( voxelZero, RASVoxelZero );
        vol->MRIIndexToRAS( voxelOne, RASVoxelOne );
        centerRAS[0] += (RASVoxelOne[0] - RASVoxelZero[0]) / 2.0;
        centerRAS[1] += (RASVoxelOne[1] - RASVoxelZero[1]) / 2.0;
        centerRAS[2] += (RASVoxelOne[2] - RASVoxelZero[2]) / 2.0;
      }

      // Set the center.
      Set2DRASCenter( centerRAS );

      // Calculate the zoom necessary.
      float dataWidth = 0, dataHeight = 0;
      switch ( mViewState.GetInPlane() ) {
      case ViewState::X:
        dataWidth  = RAS[3] - RAS[2];
        dataHeight = RAS[5] - RAS[4];
        break;
      case ViewState::Y:
        dataWidth  = RAS[1] - RAS[0];
        dataHeight = RAS[5] - RAS[4];
        break;
      case ViewState::Z:
        dataWidth  = RAS[1] - RAS[0];
        dataHeight = RAS[3] - RAS[2];
        break;
      }

      float xZoomLevel, yZoomLevel;
      xZoomLevel = (float)mViewState.GetBufferWidth() / (float)dataWidth;
      yZoomLevel = (float)mViewState.GetBufferHeight() / (float)dataHeight;
      // Scale it by 0.9 so there's a little bit of border around the
      // view.
      Set2DZoomLevel( 0.9 * (xZoomLevel<yZoomLevel?xZoomLevel:yZoomLevel) );
    }

  } catch (...) {
    DebugOutput( << "Couldn't find layer " << iLayerID );
  }
}

void
ScubaView::CopyLayerSettingsToView ( ScubaView& iView ) {

  map<int,int>::iterator tLevelLayerID;
  for ( tLevelLayerID = mLevelLayerIDMap.begin();
        tLevelLayerID != mLevelLayerIDMap.end(); ++tLevelLayerID ) {

    int level = (*tLevelLayerID).first;
    int layerID = (*tLevelLayerID).second;
    iView.SetLayerAtLevel( level, layerID );
  }
}

void
ScubaView::SetWorldToViewTransform ( int iTransformID ) {

  try {
    mWorldToView->RemoveListener( *this );
    mWorldToView = &(ScubaTransform::FindByID( iTransformID ));
    mWorldToView->AddListener( *this );
    CalcWorldToWindowTransform();
    RequestRedisplay();
  } catch (...) {
    DebugOutput( << "Couldn't find transform " << iTransformID );
  }
}

int
ScubaView::GetWorldToViewTransform () {

  return mWorldToView->GetID();
}

void
ScubaView::SetThroughPlaneIncrement ( ViewState::Plane iInPlane,
                                      float iIncrement ) {
  mThroughPlaneIncrements[iInPlane] = iIncrement;
}

float
ScubaView::GetThroughPlaneIncrement ( ViewState::Plane iInPlane ) {

  return mThroughPlaneIncrements[iInPlane];
}

list<Layer::InfoAtRAS>&
ScubaView::GetInfoAtRASList ( string isSet ) {

  map<string,list<Layer::InfoAtRAS> >::iterator tMap =
    mInfoAtRASMap.find( isSet );

  if ( tMap != mInfoAtRASMap.end() ) {
    return mInfoAtRASMap[isSet];
  } else {
    throw runtime_error( "Couldn't find that set." );
  }
}

TclCommandListener::TclCommandResult
ScubaView::DoListenToTclCommand( char* isCommand,
                                 int, char** iasArgv ) {

  // SetViewInPlane <viewID> <inPlane>
  if ( 0 == strcmp( isCommand, "SetViewInPlane" ) ) {
    int viewID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad view ID";
      return error;
    }

    if ( mID == viewID ) {

      ViewState::Plane inPlane;
      if ( 0 == strcmp( iasArgv[2], "X" ) ||
           0 == strcmp( iasArgv[2], "x" ) ||
           0 == strcmp( iasArgv[2], "R" ) ||
           0 == strcmp( iasArgv[2], "r" ) ) {
        inPlane = ViewState::X;
      } else if ( 0 == strcmp( iasArgv[2], "Y" ) ||
                  0 == strcmp( iasArgv[2], "y" ) ||
                  0 == strcmp( iasArgv[2], "A" ) ||
                  0 == strcmp( iasArgv[2], "a" ) ) {
        inPlane = ViewState::Y;
      } else if ( 0 == strcmp( iasArgv[2], "Z" ) ||
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

  // GetViewInPlane <viewID>
  if ( 0 == strcmp( isCommand, "GetViewInPlane" ) ) {
    int viewID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad view ID";
      return error;
    }

    if ( mID == viewID ) {

      sReturnFormat = "s";
      switch ( mViewState.GetInPlane() ) {
      case ViewState::X:
        sReturnValues = "x";
        break;
      case ViewState::Y:
        sReturnValues = "y";
        break;
      case ViewState::Z:
        sReturnValues = "z";
        break;
      }
    }
  }

  // SetViewZoomLevel <viewID> <zoomLevel>
  if ( 0 == strcmp( isCommand, "SetViewZoomLevel" ) ) {
    int viewID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad view ID";
      return error;
    }

    if ( mID == viewID ) {

      float zoomLevel = (float) strtod( iasArgv[2], (char**)NULL );
      if ( ERANGE == errno ) {
        sResult = "bad zoom level";
        return error;
      }

      Set2DZoomLevel( zoomLevel );
    }
  }

  // GetViewZoomLevel <viewID>
  if ( 0 == strcmp( isCommand, "GetViewZoomLevel" ) ) {
    int viewID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad view ID";
      return error;
    }

    if ( mID == viewID ) {

      stringstream ssReturn;
      sReturnFormat = "f";
      ssReturn << Get2DZoomLevel();
      sReturnValues = ssReturn.str();
      return ok;
    }
  }

  // SetViewRASCenter <viewID> <X> <Y> <Z>
  if ( 0 == strcmp( isCommand, "SetViewRASCenter" ) ) {
    int viewID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad view ID";
      return error;
    }

    if ( mID == viewID ) {

      float x = (float) strtod( iasArgv[2], (char**)NULL );
      if ( ERANGE == errno ) {
        sResult = "bad x coordinate";
        return error;
      }
      float y = (float) strtod( iasArgv[3], (char**)NULL );
      if ( ERANGE == errno ) {
        sResult = "bad y coordinate";
        return error;
      }
      float z = (float) strtod( iasArgv[4], (char**)NULL );
      if ( ERANGE == errno ) {
        sResult = "bad z coordinate";
        return error;
      }

      float center[3];
      center[0] = x;
      center[1] = y;
      center[2] = z;
      Set2DRASCenter( center );
    }
  }

  // GetViewRASCenter <viewID> <X> <Y> <Z>
  if ( 0 == strcmp( isCommand, "GetViewRASCenter" ) ) {
    int viewID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad view ID";
      return error;
    }

    if ( mID == viewID ) {

      float center[3];
      Get2DRASCenter( center );

      stringstream ssReturn;
      sReturnFormat = "Lfffl";
      ssReturn << center[0] << " " << center[1] << " " << center[2];
      sReturnValues = ssReturn.str();
      return ok;
    }
  }

  // SetLayerInViewAtLevel <viewID> <layerID> <level>
  if ( 0 == strcmp( isCommand, "SetLayerInViewAtLevel" ) ) {
    int viewID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad view ID";
      return error;
    }

    if ( mID == viewID ) {

      int layerID = strtol( iasArgv[2], (char**)NULL, 10 );
      if ( ERANGE == errno ) {
        sResult = "bad layer ID";
        return error;
      }
      int level = strtol( iasArgv[3], (char**)NULL, 10 );
      if ( ERANGE == errno ) {
        sResult = "bad level";
        return error;
      }

      SetLayerAtLevel( layerID, level );
    }
  }

  // GetLayerInViewAtLevel <viewID> <level>
  if ( 0 == strcmp( isCommand, "GetLayerInViewAtLevel" ) ) {
    int viewID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad view ID";
      return error;
    }

    if ( mID == viewID ) {

      int level = strtol( iasArgv[2], (char**)NULL, 10 );
      if ( ERANGE == errno ) {
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
  if ( 0 == strcmp( isCommand, "RemoveLayerFromViewAtLevel" ) ) {
    int viewID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad view ID";
      return error;
    }

    if ( mID == viewID ) {

      int level = strtol( iasArgv[2], (char**)NULL, 10 );
      if ( ERANGE == errno ) {
        sResult = "bad level";
        return error;
      }

      RemoveLayerAtLevel( level );
    }
  }

  // GetListOfLevelsInView <viewID>
  if ( 0 == strcmp( isCommand, "GetListOfLevelsInView" ) ) {
    int viewID;
    try {
      viewID = TclCommandManager::ConvertArgumentToInt( iasArgv[1] );
    } catch ( runtime_error& e ) {
      sResult = string("bad viewID: ") + e.what();
      return error;
    }

    if ( mID == viewID ) {
      stringstream ssFormat;
      stringstream ssResult;
      ssFormat << "L";
     
      for( map<int,int>::iterator tLayerID = mLevelLayerIDMap.begin();
	   tLayerID != mLevelLayerIDMap.end(); ++tLayerID ) {
	ssFormat << "i";
	ssResult << tLayerID->first << " ";
      }

      ssFormat << "l"; 

      sReturnFormat = ssFormat.str();
      sReturnValues = ssResult.str();
    }
  }

  // SetLevelVisibilityInView <viewID> <level> <visibility>
  if ( 0 == strcmp( isCommand, "SetLevelVisibilityInView" ) ) {

    int viewID;
    try {
      viewID = TclCommandManager::ConvertArgumentToInt( iasArgv[1] );
    } catch ( runtime_error& e ) {
      sResult = string("bad viewID: ") + e.what();
      return error;
    }

    if ( mID == viewID ) {

      int nLevel;
      try {
        nLevel = TclCommandManager::ConvertArgumentToInt( iasArgv[2] );
      } catch ( runtime_error& e ) {
        sResult = string("bad level: ") + e.what();
        return error;
      }

      bool bVisible;
      try {
        bVisible = TclCommandManager::ConvertArgumentToBoolean( iasArgv[3] );
      } catch ( runtime_error& e ) {
        sResult = string("bad visibility: ") + e.what();
        return error;
      }

      SetDrawLevelVisibility( nLevel, bVisible );
    }
  }

  // GetLevelVisibilityInView <viewID> <level>
  if ( 0 == strcmp( isCommand, "GetLevelVisibilityInView" ) ) {
    int viewID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad view ID";
      return error;
    }

    if ( mID == viewID ) {

      int nLevel;
      try {
        nLevel = TclCommandManager::ConvertArgumentToInt( iasArgv[2] );
      } catch ( runtime_error& e ) {
        sResult = string("bad level: ") + e.what();
        return error;
      }

      sReturnValues = TclCommandManager::ConvertBooleanToReturnValue( GetDrawLevelVisibility( nLevel ) );
      sReturnFormat = "i";
    }
  }

  // SetLevelReportInfoInView <viewID> <level> <reportInfo>
  if ( 0 == strcmp( isCommand, "SetLevelReportInfoInView" ) ) {

    int viewID;
    try {
      viewID = TclCommandManager::ConvertArgumentToInt( iasArgv[1] );
    } catch ( runtime_error& e ) {
      sResult = string("bad viewID: ") + e.what();
      return error;
    }

    if ( mID == viewID ) {

      int nLevel;
      try {
        nLevel = TclCommandManager::ConvertArgumentToInt( iasArgv[2] );
      } catch ( runtime_error& e ) {
        sResult = string("bad level: ") + e.what();
        return error;
      }

      bool bReportInfo;
      try {
        bReportInfo =
          TclCommandManager::ConvertArgumentToBoolean( iasArgv[3] );
      } catch ( runtime_error& e ) {
        sResult = string("bad reportInfo: ") + e.what();
        return error;
      }

      SetDrawLevelReportInfo( nLevel, bReportInfo );
    }
  }

  // GetLevelReportInfoInView <viewID> <level>
  if ( 0 == strcmp( isCommand, "GetLevelReportInfoInView" ) ) {
    int viewID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad view ID";
      return error;
    }

    if ( mID == viewID ) {

      int nLevel;
      try {
        nLevel = TclCommandManager::ConvertArgumentToInt( iasArgv[2] );
      } catch ( runtime_error& e ) {
        sResult = string("bad level: ") + e.what();
        return error;
      }

      sReturnValues = TclCommandManager::ConvertBooleanToReturnValue( GetDrawLevelReportInfo( nLevel ) );
      sReturnFormat = "i";
    }
  }

  // SetViewStateToLayerBounds <viewID> <layerID>
  if ( 0 == strcmp( isCommand, "SetViewStateToLayerBounds" ) ) {
    int viewID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad view ID";
      return error;
    }

    if ( mID == viewID ) {

      int layerID;
      try {
        layerID = TclCommandManager::ConvertArgumentToInt( iasArgv[2] );
      } catch ( runtime_error& e ) {
        sResult = string("bad layer ID: ") + e.what();
        return error;
      }

      SetViewStateToLayerBounds( layerID );
    }
  }

  // GetInfoAtRAS <viewID> <setName>
  if ( 0 == strcmp( isCommand, "GetInfoAtRAS" ) ) {
    int viewID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad view ID";
      return error;
    }

    if ( mID == viewID ) {

      string sSetName = iasArgv[2];

      try {
        list<Layer::InfoAtRAS> lInfo = GetInfoAtRASList( sSetName );
        list<Layer::InfoAtRAS>::iterator tInfo;

        stringstream ssFormat;
        stringstream ssResult;
        ssFormat << "L";

        for ( tInfo = lInfo.begin(); tInfo != lInfo.end(); ++tInfo ) {

          ssFormat << "Lssssssssssl";
          ssResult << "\"label\" \"" << (*tInfo).GetLabel() << "\" "
          << "\"value\" \"" << (*tInfo).GetValue() << "\" "
          << "\"callback\" \"" << (*tInfo).GetTclCallback() << "\" "
          << "\"filter\" \"" << (*tInfo).GetInputFilter() << "\" "
          << "\"shortenHint\" \"" << (*tInfo).GetShortenHint()<<"\" ";
        }
        ssFormat << "l";

        sReturnFormat = ssFormat.str();
        sReturnValues = ssResult.str();
      } catch ( exception& e ) {
        sResult = e.what();
        return error;
      }

    }
  }

  // GetFirstUnusedDrawLevelInView <viewID>
  if ( 0 == strcmp( isCommand, "GetFirstUnusedDrawLevelInView" ) ) {
    int viewID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad view ID";
      return error;
    }

    if ( mID == viewID ) {

      sReturnFormat = "i";
      stringstream ssReturnValues;
      ssReturnValues << GetFirstUnusedDrawLevel();
      sReturnValues = ssReturnValues.str();
    }
  }


  // SetViewLinkedStatus <viewID> <linked>
  if ( 0 == strcmp( isCommand, "SetViewLinkedStatus" ) ) {
    int viewID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad view ID";
      return error;
    }

    if ( mID == viewID ) {

      if ( 0 == strcmp( iasArgv[2], "true" ) ||
           0 == strcmp( iasArgv[2], "1" )) {
        SetLinkedStatus( true );
      } else if ( 0 == strcmp( iasArgv[2], "false" ) ||
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
  if ( 0 == strcmp( isCommand, "GetViewLinkedStatus" ) ) {
    int viewID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad view ID";
      return error;
    }

    if ( mID == viewID ) {

      sReturnFormat = "i";
      stringstream ssReturnValues;
      ssReturnValues << (int)GetLinkedStatus();
      sReturnValues = ssReturnValues.str();
    }
  }

  // SetViewLockOnCursor <viewID> <lock>
  if ( 0 == strcmp( isCommand, "SetViewLockOnCursor" ) ) {
    int viewID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad view ID";
      return error;
    }

    if ( mID == viewID ) {

      if ( 0 == strcmp( iasArgv[2], "true" ) ||
           0 == strcmp( iasArgv[2], "1" )) {
        SetLockOnCursor( true );
      } else if ( 0 == strcmp( iasArgv[2], "false" ) ||
                  0 == strcmp( iasArgv[2], "0" ) ) {
        SetLockOnCursor( false );
      } else {
        sResult = "bad lock \"" + string(iasArgv[2]) +
                  "\", should be true, 1, false, or 0";
        return error;
      }
    }
  }

  // GetViewLockOnCursor <viewID>
  if ( 0 == strcmp( isCommand, "GetViewLockOnCursor" ) ) {
    int viewID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad view ID";
      return error;
    }

    if ( mID == viewID ) {

      sReturnFormat = "i";
      stringstream ssReturnValues;
      ssReturnValues << (int)GetLockOnCursor();
      sReturnValues = ssReturnValues.str();
    }
  }

  // SetViewTransform <viewID> <transformID>
  if ( 0 == strcmp( isCommand, "SetViewTransform" ) ) {
    int viewID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad view ID";
      return error;
    }

    if ( mID == viewID ) {

      int transformID = strtol( iasArgv[2], (char**)NULL, 10 );
      if ( ERANGE == errno ) {
        sResult = "bad transformID";
        return error;
      }

      SetWorldToViewTransform( transformID );
    }
  }

  // GetViewTransform <viewID>
  if ( 0 == strcmp( isCommand, "GetViewTransform" ) ) {
    int viewID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad view ID";
      return error;
    }

    if ( mID == viewID ) {

      sReturnFormat = "i";
      stringstream ssReturnValues;
      ssReturnValues << (int)GetWorldToViewTransform();
      sReturnValues = ssReturnValues.str();
    }
  }

  // SetViewFlipLeftRightYZ <viewID> <flip>
  if ( 0 == strcmp( isCommand, "SetViewFlipLeftRightYZ" ) ) {
    int viewID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad view ID";
      return error;
    }

    if ( mID == viewID ) {

      if ( 0 == strcmp( iasArgv[2], "true" ) ||
           0 == strcmp( iasArgv[2], "1" )) {
        SetFlipLeftRightYZ( true );
      } else if ( 0 == strcmp( iasArgv[2], "false" ) ||
                  0 == strcmp( iasArgv[2], "0" ) ) {
        SetFlipLeftRightYZ( false );
      } else {
        sResult = "bad flip value \"" + string(iasArgv[2]) +
                  "\", should be true, 1, false, or 0";
        return error;
      }
    }
  }

  // GetViewFlipLeftRightYZ <viewID>
  if ( 0 == strcmp( isCommand, "GetViewFlipLeftRightYZ" ) ) {
    int viewID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad view ID";
      return error;
    }

    if ( mID == viewID ) {

      sReturnFormat = "i";
      stringstream ssReturnValues;
      ssReturnValues << (int)GetFlipLeftRightYZ();
      sReturnValues = ssReturnValues.str();
    }
  }

  // SetViewThroughPlaneIncrement <viewID> <inPlane> <increment>
  if ( 0 == strcmp( isCommand, "SetViewThroughPlaneIncrement" ) ) {
    int viewID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad view ID";
      return error;
    }

    if ( mID == viewID ) {

      float increment = (float) strtod( iasArgv[3], (char**)NULL );
      if ( ERANGE == errno ) {
        sResult = "bad increment";
        return error;
      }

      if ( 0 == strcmp( iasArgv[2], "x" ) ||
           0 == strcmp( iasArgv[2], "X" ) ) {

        SetThroughPlaneIncrement( ViewState::X, increment );

      } else if ( 0 == strcmp( iasArgv[2], "y" ) ||
                  0 == strcmp( iasArgv[2], "Y" ) ) {

        SetThroughPlaneIncrement( ViewState::Y, increment );

      } else if ( 0 == strcmp( iasArgv[2], "z" ) ||
                  0 == strcmp( iasArgv[2], "Z" ) ) {

        SetThroughPlaneIncrement( ViewState::Z, increment );

      } else {
        stringstream ssResult;
        ssResult << "bad throughPlane \"" << iasArgv[1] << "\", should be "
        << "x, y, or z.";
        sResult = ssResult.str();
        return error;
      }
    }
  }

  // GetViewThroughPlaneIncrement <viewID> <inPlane>
  if ( 0 == strcmp( isCommand, "GetViewThroughPlaneIncrement" ) ) {
    int viewID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad view ID";
      return error;
    }

    if ( mID == viewID ) {

      float increment = 0;
      if ( 0 == strcmp( iasArgv[2], "x" ) ||
           0 == strcmp( iasArgv[2], "X" ) ) {

        increment = GetThroughPlaneIncrement( ViewState::X );

      } else if ( 0 == strcmp( iasArgv[2], "y" ) ||
                  0 == strcmp( iasArgv[2], "Y" ) ) {

        increment = GetThroughPlaneIncrement( ViewState::Y );

      } else if ( 0 == strcmp( iasArgv[2], "z" ) ||
                  0 == strcmp( iasArgv[2], "Z" ) ) {

        increment = GetThroughPlaneIncrement( ViewState::Z );

      } else {
        stringstream ssResult;
        ssResult << "bad throughPlane \"" << iasArgv[2] << "\", should be "
        << "x, y, or z.";
        sResult = ssResult.str();
        return error;
      }

      sReturnFormat = "f";
      stringstream ssReturnValues;
      ssReturnValues << increment;
      sReturnValues = ssReturnValues.str();
    }
  }

  // GetVolumeHistogramInView <viewID> <volID> <roiID> <numBins>
  if ( 0 == strcmp( isCommand, "GetVolumeHistogramInView" ) ) {

    int viewID;
    try {
      viewID = TclCommandManager::ConvertArgumentToInt( iasArgv[1] );
    } catch ( runtime_error& e ) {
      sResult = string("bad viewID: ") + e.what();
      return error;
    }

    if ( mID == viewID ) {

      int collectionID;
      try {
        collectionID = TclCommandManager::ConvertArgumentToInt( iasArgv[2] );
      } catch ( runtime_error& e ) {
        sResult = string("bad collectionID: ") + e.what();
        return error;
      }

      // Find a volume collection.
      VolumeCollection* vol = NULL;
      try {
        DataCollection* col = &DataCollection::FindByID( collectionID );
        //    VolumeCollection* vol = dynamic_cast<VolumeCollection*>(col);
        vol = (VolumeCollection*)col;
      } catch (...) {
        throw runtime_error( "Couldn't find volume or data collection "
                             "wasn't a volume.." );
      }

      // Find an ROI if they want.
      int roiID;
      try {
        roiID = TclCommandManager::ConvertArgumentToInt( iasArgv[3] );
      } catch ( runtime_error& e ) {
        sResult = string("bad roiID: ") + e.what();
        return error;
      }

      ScubaROIVolume* volROI = NULL;
      if ( roiID != -1 ) {
        try {
          ScubaROI* roi = &ScubaROI::FindByID( roiID );
          volROI = (ScubaROIVolume*)roi;
        } catch (...) {
          throw runtime_error( "Couldn't find ROI or ROI "
                               "wasn't a volume ROI.." );
        }
      }

      int cBins;
      try {
        cBins = TclCommandManager::ConvertArgumentToInt( iasArgv[4] );
      } catch ( runtime_error& e ) {
        sResult = string("bad numBins: ") + e.what();
        return error;
      }

      float minBinValue, binIncrement;
      map<int,int> binCounts;
      GetVolumeHistogramInView( *vol, volROI, cBins,
                                minBinValue, binIncrement, binCounts );

      // Build the output.
      stringstream ssValues;
      stringstream ssFormat;

      ssFormat << "Lff";
      ssValues << minBinValue << " " << binIncrement << " ";

      ssFormat << "L";
      for ( int nBin = 0; nBin < (int)binCounts.size(); nBin++ ) {
        ssFormat << "i";
        ssValues << binCounts[nBin] << " ";
      }
      ssFormat << "ll";

      sReturnFormat = ssFormat.str();
      sReturnValues = ssValues.str();
    }

  }

  // DoVolumeValueRangeFillInView <viewID> <sourceVolID> <roiID> <destVolID> <valueRanges>
  if ( 0 == strcmp( isCommand, "DoVolumeValueRangeFillInView" ) ) {

    int viewID;
    try {
      viewID = TclCommandManager::ConvertArgumentToInt( iasArgv[1] );
    } catch ( runtime_error& e ) {
      sResult = string("bad viewID: ") + e.what();
      return error;
    }

    if ( mID == viewID ) {

      int sourceVolID;
      try {
        sourceVolID = TclCommandManager::ConvertArgumentToInt( iasArgv[2] );
      } catch ( runtime_error& e ) {
        sResult = string("bad sourceVolID: ") + e.what();
        return error;
      }

      // Find a volume collection.
      VolumeCollection* sourceVol = NULL;
      try {
        DataCollection* col = &DataCollection::FindByID( sourceVolID );
        sourceVol = dynamic_cast<VolumeCollection*>(col);
      } catch (...) {
        throw runtime_error( "Couldn't find source volume or data collection "
                             "wasn't a volume." );
      }

      // Find an ROI if they want.
      int roiID;
      try {
        roiID = TclCommandManager::ConvertArgumentToInt( iasArgv[3] );
      } catch ( runtime_error& e ) {
        sResult = string("bad roiID: ") + e.what();
        return error;
      }

      ScubaROIVolume* volROI = NULL;
      if ( roiID != -1 ) {
        try {
          ScubaROI* roi = &ScubaROI::FindByID( roiID );
          volROI = (ScubaROIVolume*)roi;
        } catch (...) {
          throw runtime_error( "Couldn't find ROI or ROI "
			       "wasn't a volume ROI." );
        }
      }

      int destVolID;
      try {
        destVolID = TclCommandManager::ConvertArgumentToInt( iasArgv[4] );
      } catch ( runtime_error& e ) {
        sResult = string("bad destVolID: ") + e.what();
        return error;
      }

      // Find a volume collection.
      VolumeCollection* destVol = NULL;
      try {
        DataCollection* col = &DataCollection::FindByID( destVolID );
        //    VolumeCollection* vol = dynamic_cast<VolumeCollection*>(col);
        destVol = (VolumeCollection*)col;
      } catch (...) {
        throw runtime_error( "Couldn't find dest volume or data collection "
                             "wasn't a volume." );
      }

      // Parse out the range elements.
      vector<ValueRangeFillElement> lElements;
      try {
	stringstream ssElements( iasArgv[5] );
	while( !ssElements.eof() ) {
	  float beginValue, endValue, newValue;
	  ssElements >> beginValue;
	  ssElements >> endValue;
	  ssElements >> newValue;
	  ValueRangeFillElement element( beginValue, endValue, newValue );
	  lElements.push_back( element );
	}
      } catch (...) {
        throw runtime_error( "Invalid value range list." );
      }

      DoVolumeValueRangeFill( *sourceVol, volROI, *destVol, lElements );
    }
  }

  // ConvertWindowToViewRAS <viewID> <X> <Y>
  if ( 0 == strcmp( isCommand, "ConvertWindowToViewRAS" ) ) {
    int viewID;
    try {
      viewID = TclCommandManager::ConvertArgumentToInt( iasArgv[1] );
    } catch ( runtime_error& e ) {
      sResult = string("bad viewID: ") + e.what();
      return error;
    }

    if ( mID == viewID ) {

      int window[2];
      float ras[3];
      try {
        window[0] = TclCommandManager::ConvertArgumentToInt( iasArgv[2] );
        window[1] = TclCommandManager::ConvertArgumentToInt( iasArgv[3] );
      } catch ( runtime_error& e ) {
        sResult = string("bad window coordinate: ") + e.what();
        return error;
      }

      TranslateWindowToRAS( window, ras );

      stringstream ssReturn;
      sReturnFormat = "Lfffl";
      ssReturn << ras[0] << " " << ras[1] << " " << ras[2];
      sReturnValues = ssReturn.str();
      return ok;
    }
  }


  return ok;
}

void
ScubaView::DoListenToMessage ( string isMessage, void* iData ) {

  if ( isMessage == "transformChanged" ) {

    // World to view transform changed. Our overlay coords are
    // different now. Also recalc the overall transform.
    RebuildOverlayDrawList();
    CalcWorldToWindowTransform();
    RequestRedisplay();
  }

  if ( isMessage == "layerChanged" ) {

    // Might have changed the values or labels.
    RebuildLabelValueInfo( mCursor.xyz(), "cursor" );
    RebuildLabelValueInfo( mLastMouseOver.xyz(), "mouse" );

    RequestRedisplay();

    SendBroadcast( "viewChanged", NULL );
  }

  if ( isMessage == "layerDeleted" ) {

    // A layer is being deleted. Grab the ID, figure out if it's in
    // any of oour levels, and if so, clear them.
    int deleteLayerID = *(int*)iData;
    map<int,int>::iterator tLevelLayerID;
    for ( tLevelLayerID = mLevelLayerIDMap.begin();
          tLevelLayerID != mLevelLayerIDMap.end(); ++tLevelLayerID ) {
      int layerID = (*tLevelLayerID).second;
      int nLevel = (*tLevelLayerID).first;
      if ( layerID == deleteLayerID ) {
        RemoveLayerAtLevel( nLevel );
      }
    }
  }

  if ( isMessage == "layerInfoSettingsChanged" ) {
    RebuildLabelValueInfo( mCursor.xyz(), "cursor" );
    RebuildLabelValueInfo( mLastMouseOver.xyz(), "mouse" );
  }

  if ( isMessage == "cursorChanged" ) {

    // Rebuild the cursor info.
    RebuildLabelValueInfo( mCursor.xyz(), "cursor" );

    // If we're locked on the cursor, and the cursor is no longer in
    // our view plane, set our view now. But only focus on the
    // cursor's new plane to minimize the view jumping around. Unless
    // the cursor is no longer visible in the view, then recenter.
    if ( mbLockOnCursor ) {

      float newCenter[3];
      mViewState.GetCenterRAS( newCenter );

      // Make sure we stay in plane.
      if ( !mViewState.IsRASVisibleInPlane( mCursor.xyz(),
                                            GetThroughPlaneIncrement( mViewState.GetInPlane() ))) {

        switch ( mViewState.GetInPlane() ) {
        case ViewState::X:
          newCenter[0] = mCursor[0];
          break;
        case ViewState::Y:
          newCenter[1] = mCursor[1];
          break;
        case ViewState::Z:
          newCenter[2] = mCursor[2];
          break;
        }
      }

      // Make sure this is in the window view. Convert to window
      // coords and if any are outside of the view...
      int window[2];
      TranslateRASToWindow( mCursor.xyz(), window );
      if ( window[0] < 0 || window[0] > mWidth ||
           window[1] < 0 || window[1] > mHeight ) {

        // Just center around cursor.
        newCenter[0] = mCursor[0];
        newCenter[1] = mCursor[1];
        newCenter[2] = mCursor[2];
      }

      Set2DRASCenter( newCenter );
    }

    RebuildOverlayDrawList();
    RequestRedisplay();
  }

  if ( isMessage == "markerChanged" ) {
    RebuildOverlayDrawList();
    RequestRedisplay();
  }

  if ( isMessage == "DrawCoordinateOverlay" ||
       isMessage == "DrawPlaneIntersections" ||
       isMessage == "DrawMarkers" ||
       isMessage == "DrawPaths") {
    RebuildOverlayDrawList(); // our overlay will be different
    RequestRedisplay();
  }

  if ( isMessage == "2DRASCenterChanged" ||
       isMessage == "2DZoomLevelChanged" ||
       isMessage == "2DInPlaneChanged" ||
       isMessage == "2DPlaneNormalChanged" ) {

    int viewID = *(int*)iData;

    // If somebody else's RAS center or inplane changed, we need to
    // recalc our inplane intersection points.
    CalcViewIntersectionPoints( viewID );

    // If we're linked, we need to change our view.
    if ( mViewIDLinkedList[GetID()] && mViewIDLinkedList[viewID] ) {
      View& view = View::FindByID( viewID );
      // ScubaView& scubaView = dynamic_cast<ScubaView&>(view);
      ScubaView& scubaView = (ScubaView&)view;

      // Change center or zoom level if linked. Don't link the in plane.
      if ( isMessage == "2DRASCenterChanged" ) {
        float RASCenter[3];
        scubaView.Get2DRASCenter( RASCenter );
        Set2DRASCenter( RASCenter );
      } else if ( isMessage == "2DZoomLevelChanged" ) {
        float zoomLevel;
        zoomLevel = scubaView.Get2DZoomLevel();
        Set2DZoomLevel( zoomLevel );
      }
    }

    // If we're visible request redisplays.
    if ( IsVisibleInFrame() ) {
      RebuildOverlayDrawList();
      RequestRedisplay();
    }
  }

  // We cache these values but down have to act on them right away.
  if ( isMessage == "KeyMoveViewLeft" )
    msMoveViewLeft->SetFromString( *(string*)iData );
  if ( isMessage == "KeyMoveViewRight" )
    msMoveViewRight->SetFromString( *(string*)iData );
  if ( isMessage == "KeyMoveViewUp" )
    msMoveViewUp->SetFromString( *(string*)iData );
  if ( isMessage == "KeyMoveViewDown" )
    msMoveViewDown->SetFromString( *(string*)iData );
  if ( isMessage == "KeyMoveViewIn" )
    msMoveViewIn->SetFromString( *(string*)iData );
  if ( isMessage == "KeyMoveViewOut" )
    msMoveViewOut->SetFromString( *(string*)iData );
  if ( isMessage == "KeyZoomViewIn" )
    msZoomViewIn->SetFromString( *(string*)iData );
  if ( isMessage == "KeyZoomViewOut" )
    msZoomViewOut->SetFromString( *(string*)iData );

  // New view, get some info about it.
  if ( isMessage == "NewView" ) {
    int viewID = *(int*)iData;
    if ( viewID != GetID() ) {
      CalcViewIntersectionPoints( viewID );
    }
  }

  if ( isMessage == "pathChanged" ||
       isMessage == "pathVertexAdded" ) {
    RebuildOverlayDrawList();
    RequestRedisplay();
  }

  View::DoListenToMessage( isMessage, iData );
}

void
ScubaView::DoDraw() {

  ScubaGlobalPreferences& prefs = ScubaGlobalPreferences::GetPreferences();
  ::Timer timer;
  timer.Start();

  BuildFrameBuffer();
  DrawFrameBuffer();
  DrawOverlay();

  // Go through our draw levels. For each one, call the draw list.
  map<int,int>::iterator tLevelLayerID;
  for ( tLevelLayerID = mLevelLayerIDMap.begin();
        tLevelLayerID != mLevelLayerIDMap.end(); ++tLevelLayerID ) {
    int nLevel = (*tLevelLayerID).first;
    glCallList( mLevelGLListIDMap[nLevel] );
  }

  int msec = timer.TimeNow();

  if ( prefs.GetPrefAsBool( ScubaGlobalPreferences::ShowFPS )) {
    float fps = 1.0 / ((float)msec/1000.0);

    stringstream ssCommand;
    ssCommand << "SetStatusBarText \"" << fps << " fps\"";

    TclCommandManager& mgr = TclCommandManager::GetManager();
    mgr.SendCommand( ssCommand.str() );
  }

}

void
ScubaView::DoReshape( int iWidth, int iHeight ) {

  if ( iWidth < 0 || iHeight < 0 ) {
    stringstream sError;
    sError << "Invalid width " << mWidth << " or height " << mHeight;
    DebugOutput( << sError.str() );
    throw runtime_error( sError.str() );
  }

  // Set the view state buffer height and width.
  mViewState.SetBufferWidth( iWidth );
  mViewState.SetBufferHeight( iHeight );

  // Allocate a new buffer.
  GLubyte* newBuffer = (GLubyte*) malloc( mWidth * mHeight * kBytesPerPixel );
  if ( NULL == newBuffer ) {
    stringstream sError;
    sError << "Couldn't allocate buffer for width " << mWidth
    << " height " << mHeight;
    DebugOutput( << sError.str() );
    throw runtime_error( sError.str() );
  }

  // Get rid of the old one.
  if ( NULL != mBuffer ) {
    free( mBuffer );
  }

  // Save the new one.
  mBuffer = newBuffer;

  // Set the width and height in all the layers.
  map<int,int>::iterator tLevelLayerID;
  for ( tLevelLayerID = mLevelLayerIDMap.begin();
        tLevelLayerID != mLevelLayerIDMap.end(); ++tLevelLayerID ) {

    int layerID = (*tLevelLayerID).second;

    try {
      Layer& layer = Layer::FindByID( layerID );

      layer.SetWidth( mWidth );
      layer.SetHeight( mHeight );
    } catch (...) {
      DebugOutput( << "Couldn't find layer " << layerID );
    }
  }

  RebuildOverlayDrawList();
}

void
ScubaView::DoTimer() {

  map<int,int>::iterator tLevelLayerID;
  for ( tLevelLayerID = mLevelLayerIDMap.begin();
        tLevelLayerID != mLevelLayerIDMap.end(); ++tLevelLayerID ) {
    int layerID = (*tLevelLayerID).second;
    try {
      Layer& layer = Layer::FindByID( layerID );
      layer.Timer();
      if ( layer.WantRedisplay() ) {
        RequestRedisplay();
        layer.RedisplayPosted();
      }
    } catch (...) {
      DebugOutput( << "Couldn't find layer " << layerID );
    }
  }
}

void
ScubaView::DoMouseMoved( int iWindow[2],
                         InputState& iInput, ScubaToolState& iTool ) {

  float ras[3];
  TranslateWindowToRAS( iWindow, ras );

  // Rebuild our label value info because the mouse has moved.
  RebuildLabelValueInfo( ras, "mouse" );
  mLastMouseOver.Set( ras );

  // Handle the navigation tool and plane tool.
  if ( iTool.GetMode() == ScubaToolState::navigation ||
       iTool.GetMode() == ScubaToolState::plane ) {

    if ( iInput.Button() &&
         !iInput.IsControlKeyDown() && !iInput.IsShiftKeyDown() ) {

      float delta[2];
      delta[0] = (float)(mLastMouseMoved[0] - iWindow[0]) /
                 mViewState.GetZoomLevel();
      delta[1] = (float)(mLastMouseMoved[1] - iWindow[1]) /
                 mViewState.GetZoomLevel();

      /* add to the total delta */
      mMouseMoveDelta[0] += delta[0];
      mMouseMoveDelta[1] += delta[1];

      /* save this mouse position */
      mLastMouseMoved[0] = iWindow[0];
      mLastMouseMoved[1] = iWindow[1];

      float moveLeftRight = 0, moveUpDown = 0, moveInOut = 0, zoomInOut = 0;
      switch ( iInput.Button() ) {
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

      if ( mbFlipLeftRightInYZ &&
           (mViewState.GetInPlane() == ViewState::Y ||
            mViewState.GetInPlane() == ViewState::Z) ) {
        moveLeftRight = -moveLeftRight;
      }

      if ( moveLeftRight || moveUpDown || moveInOut ) {

        if ( iTool.GetMode() == ScubaToolState::navigation ) {

          // Put the window relative translations into window relative.
          Point3<float> move=NULL;
          switch ( mViewState.GetInPlane() ) {
          case ViewState::X:
            move.Set( moveInOut, moveLeftRight, moveUpDown );
            break;
          case ViewState::Y:
            move.Set( moveLeftRight, moveInOut, moveUpDown );
            break;
          case ViewState::Z:
            move.Set( moveLeftRight, moveUpDown, moveInOut );
            break;
          }

          // Do the move.
          Point3<float> newCenterRAS;
          TranslateRASInWindowSpace( mOriginalCenterRAS, move.xyz(),
                                     newCenterRAS.xyz() );

          // Set the new center.
          Set2DRASCenter( newCenterRAS.xyz() );

        } else if ( iTool.GetMode() == ScubaToolState::plane ) {

          // If we have a view to move...
          if ( mCurrentMovingViewIntersection != -1 ) {

            // Get the view..
            View& view = View::FindByID( mCurrentMovingViewIntersection );
            // ScubaView& scubaView = dynamic_cast<ScubaView&>(view);
            ScubaView& scubaView = (ScubaView&)view;

            switch ( iInput.Button() ) {

              // Button 1, move that view's center RAS.
            case 1: {

              // Put the window relative translations into window
              // relative. Note we don't move in or out in this tool.
              Point3<float> move=NULL;
              switch ( mViewState.GetInPlane() ) {
              case ViewState::X:
                move.Set( 0, -moveLeftRight, -moveUpDown );
                break;
              case ViewState::Y:
                move.Set( -moveLeftRight, 0, -moveUpDown );
                break;
              case ViewState::Z:
                move.Set( -moveLeftRight, -moveUpDown, 0 );
                break;
              }

              // Do the move.
              Point3<float> newCenterRAS;
              TranslateRASInWindowSpace( mOriginalCenterRAS, move.xyz(),
                                         newCenterRAS.xyz() );

              // Set the new center in the view.
              scubaView.Set2DRASCenter( newCenterRAS.xyz() );

            }
            break;

            // Button 2, rotate.
            case 2: {

              // Find the biggest delta and divide it by the screen
              // dimension.
              float delta = (fabs(mMouseMoveDelta[0]) >
                             fabs(mMouseMoveDelta[1]) ?
                             mMouseMoveDelta[0] / (float)mWidth :
                             mMouseMoveDelta[1] / (float)mHeight);

              // Multiple by ~2pi to get radians.
              float rads = delta * 6.3;

              // We're going to rotate around our plane normal. With
              // an origin of 0.
              Point3<float> axis( mViewState.GetPlaneNormal() );
              Point3<float> origin( 0,0,0 );
              Matrix44 rotate;
              rotate.MakeRotation( origin.xyz(), axis.xyz(), rads );

              // Now rotate the original normal to get a new normal.
              Point3<float> newPlaneNormal;
              rotate.MultiplyVector3( mOriginalPlaneNormal.xyz(),
                                      newPlaneNormal.xyz() );

              // Normalize it and set it in the view.
              VectorOps::Normalize( newPlaneNormal );
              scubaView.Set2DPlaneNormal( newPlaneNormal.xyz() );
            }
            break;

            default:
              break;
            }
          }
        }
      }

      if ( zoomInOut ) {

        float newZoom = mOriginalZoom + zoomInOut;
        if ( newZoom <= 0.25 ) {
          newZoom = 0.25;
        }
        Set2DZoomLevel( newZoom );
      }
    }
  }

  // If not a straight control-click, pass this tool to our layers.
  if ( !(iInput.IsControlKeyDown() &&
         !iInput.IsShiftKeyDown() && !iInput.IsAltKeyDown()) ) {
    map<int,int>::iterator tLevelLayerID;
    for ( tLevelLayerID = mLevelLayerIDMap.begin();
          tLevelLayerID != mLevelLayerIDMap.end(); ++tLevelLayerID ) {
      int layerID = (*tLevelLayerID).second;
      try {
        Layer& layer = Layer::FindByID( layerID );
        layer.HandleTool( ras, mViewState, *this, iTool, iInput );
        if ( layer.WantRedisplay() ) {
          RequestRedisplay();
          layer.RedisplayPosted();
        }
      } catch ( runtime_error& e) {
        throw( e );
      } catch (...) {
        DebugOutput( << "Couldn't find layer " << layerID );
      }
    }
  }
}

void
ScubaView::DoMouseUp( int iWindow[2],
                      InputState& iInput, ScubaToolState& iTool ) {

  // No matter what tool we're on, look for ctrl-b{1,2,3} and do some
  // navigation stuff.
  if ( iInput.IsControlKeyDown() &&
       !iInput.IsShiftKeyDown() && !iInput.IsAltKeyDown() ) {

    // Set the new view center to this point. If they also hit b1 or
    // b3, zoom in or out accordingly.
    float world[3];
    TranslateWindowToRAS( iWindow, world );
    Set2DRASCenter( world );

    switch ( iInput.Button() ) {
    case 1:
      mViewState.SetZoomLevel( mViewState.GetZoomLevel() * 2.0 );
      break;
    case 3:
      mViewState.SetZoomLevel( mViewState.GetZoomLevel() / 2.0 );
      if ( mViewState.GetZoomLevel() < 0.25 ) {
        mViewState.SetZoomLevel( 0.25 );
      }
      break;
    }

    RebuildOverlayDrawList();
    RequestRedisplay();
  }

  // Always set cursor on mouse up button except for nav tool and plane tool.
  if ( iTool.GetMode() != ScubaToolState::navigation &&
       iTool.GetMode() != ScubaToolState::plane &&
       !iInput.IsShiftKeyDown() &&
       !iInput.IsControlKeyDown() &&
       iInput.Button() == 1 ) {

    float world[3];
    TranslateWindowToRAS( iWindow, world );
    SetCursor( world );
    RebuildLabelValueInfo( mCursor.xyz(), "cursor" );
  }

  // Handle marker tool.
  if ( iTool.GetMode() == ScubaToolState::marker &&
       !iInput.IsControlKeyDown() ) {

    float world[3];
    TranslateWindowToRAS( iWindow, world );

    switch ( iInput.Button() ) {
    case 2:
      SetNextMarker( world );
      break;
    case 3:
      HideNearestMarker( world );
      break;
    }
  }

  // Unselect view intersection line if plane tool.
  if ( iTool.GetMode() == ScubaToolState::plane ) {

    // If we clicked with button 3, set the plane to its orthogonal
    // position.
    if ( mCurrentMovingViewIntersection != -1 &&
         iInput.Button() == 3 ) {

      // Get the view..
      View& view = View::FindByID( mCurrentMovingViewIntersection );
      // ScubaView& scubaView = dynamic_cast<ScubaView&>(view);
      ScubaView& scubaView = (ScubaView&)view;

      scubaView.Set2DPlaneNormalOrthogonalToInPlane();
    }

    mCurrentMovingViewIntersection = -1;
    RebuildOverlayDrawList();
    RequestRedisplay();
  }

  // If not a straight control-click, pass this tool to our layers.
  if ( !(iInput.IsControlKeyDown() &&
         !iInput.IsShiftKeyDown() && !iInput.IsAltKeyDown()) ) {
    float ras[3];
    TranslateWindowToRAS( iWindow, ras );
    map<int,int>::iterator tLevelLayerID;
    for ( tLevelLayerID = mLevelLayerIDMap.begin();
          tLevelLayerID != mLevelLayerIDMap.end(); ++tLevelLayerID ) {
      int layerID = (*tLevelLayerID).second;
      try {
        Layer& layer = Layer::FindByID( layerID );
        layer.HandleTool( ras, mViewState, *this, iTool, iInput );
        if ( layer.WantRedisplay() ) {
          RequestRedisplay();
          layer.RedisplayPosted();
        }
      } catch ( runtime_error& e) {
        throw( e );
      } catch (...) {
        DebugOutput( << "Couldn't find layer " << layerID );
      }
    }
  }
}

void
ScubaView::DoMouseDown( int iWindow[2],
                        InputState& iInput, ScubaToolState& iTool ) {

  Point3<float> ras;
  TranslateWindowToRAS( iWindow, ras.xyz() );

  mLastMouseDown[0] = mLastMouseMoved[0] = iWindow[0];
  mLastMouseDown[1] = mLastMouseMoved[1] = iWindow[1];
  mMouseMoveDelta[0] = 0.0;
  mMouseMoveDelta[1] = 0.0;
  mViewState.GetCenterRAS( mOriginalCenterRAS );
  mOriginalZoom = mViewState.GetZoomLevel();

  // If this is the plane tool, find the nearest plane line.
  if ( iTool.GetMode() == ScubaToolState::plane ) {

    // Find the closest inplane line from another view.
    float minDistance = 999999;
    int closestViewID = -1;
    list<int> viewIDs;
    GetIDList( viewIDs );
    list<int>::iterator tViewID;
    for ( tViewID = viewIDs.begin(); tViewID != viewIDs.end(); ++tViewID ) {

      int viewID = *tViewID;
      if ( viewID != GetID () ) {

        View& view = View::FindByID( viewID );
        // ScubaView& scubaView = dynamic_cast<ScubaView&>(view);
        ScubaView& scubaView = (ScubaView&)view;

        if ( scubaView.IsVisibleInFrame() ) {

          Point3<float> P1, P2;
          P1.Set( mViewIDViewIntersectionPointMap[viewID][0] );
          P2.Set( mViewIDViewIntersectionPointMap[viewID][1] );

          float distance =
            Utilities::DistanceFromLineToPoint3f( P1, P2, ras );
          if ( distance <= minDistance ) {
            minDistance = distance;
            closestViewID = viewID;
          }
        }
      }
    }

    mCurrentMovingViewIntersection = closestViewID;
    if ( mCurrentMovingViewIntersection != -1 ) {
      // Get the views current RAS center.
      View& view = View::FindByID( mCurrentMovingViewIntersection );
      // ScubaView& scubaView = dynamic_cast<ScubaView&>(view);
      ScubaView& scubaView = (ScubaView&)view;
      scubaView.Get2DRASCenter( mOriginalCenterRAS );
      scubaView.Get2DPlaneNormal( mOriginalPlaneNormal.xyz() );

      RebuildOverlayDrawList();
      RequestRedisplay();
    }
  }

  // If not a straight control-click, pass this tool to our layers.
  if ( !(iInput.IsControlKeyDown() &&
         !iInput.IsShiftKeyDown() && !iInput.IsAltKeyDown()) ) {
    map<int,int>::iterator tLevelLayerID;
    for ( tLevelLayerID = mLevelLayerIDMap.begin();
          tLevelLayerID != mLevelLayerIDMap.end(); ++tLevelLayerID ) {
      int layerID = (*tLevelLayerID).second;
      try {
        Layer& layer = Layer::FindByID( layerID );
        layer.HandleTool( ras.xyz(), mViewState, *this, iTool, iInput );
        if ( layer.WantRedisplay() ) {
          RequestRedisplay();
          layer.RedisplayPosted();
        }
      } catch ( runtime_error& e) {
        throw( e );
      } catch (...) {
        DebugOutput( << "Couldn't find layer " << layerID );
      }
    }
  }
}

void
ScubaView::DoKeyDown( int iWindow[2],
                      InputState& iInput, ScubaToolState& iTool ) {

  ScubaKeyCombo const& key = iInput.Key();

  float ras[3];
  TranslateWindowToRAS( iWindow, ras );

  // Start with a move distance of 1. If we're moving in plane, set
  // that to the specific in plane increment. If control is down,
  // multiplay that value by 10.
  float moveDistance = 1.0;
  if ( key.IsSameAs(msMoveViewIn) || key.IsSameAs(msMoveViewOut) ) {
    moveDistance = GetThroughPlaneIncrement( mViewState.GetInPlane() );
  }
  if ( iInput.IsControlKeyDown() ) {
    moveDistance = 10.0;
  }


  if ( key.IsSameAs(msMoveViewLeft) || key.IsSameAs(msMoveViewRight) ||
       key.IsSameAs(msMoveViewDown) || key.IsSameAs(msMoveViewUp) ||
       key.IsSameAs(msMoveViewIn)   || key.IsSameAs(msMoveViewOut) ) {

    float move[3] = {0, 0, 0};

    if ( key.IsSameAs(msMoveViewLeft) ) {
      switch ( mViewState.GetInPlane() ) {
      case ViewState::X:
        move[1] = -moveDistance;
        break;
      case ViewState::Y:
        move[0] = moveDistance;
        break;
      case ViewState::Z:
        move[0] = moveDistance;
        break;
      }
    } else if ( key.IsSameAs(msMoveViewRight) ) {
      switch ( mViewState.GetInPlane() ) {
      case ViewState::X:
        move[1] = moveDistance;
        break;
      case ViewState::Y:
        move[0] = -moveDistance;
        break;
      case ViewState::Z:
        move[0] = -moveDistance;
        break;
      }
    } else if ( key.IsSameAs(msMoveViewDown) ) {
      switch ( mViewState.GetInPlane() ) {
      case ViewState::X:
        move[2] = -moveDistance;
        break;
      case ViewState::Y:
        move[2] = -moveDistance;
        break;
      case ViewState::Z:
        move[1] = -moveDistance;
        break;
      }
    } else if ( key.IsSameAs(msMoveViewUp) ) {
      switch ( mViewState.GetInPlane() ) {
      case ViewState::X:
        move[2] = moveDistance;
        break;
      case ViewState::Y:
        move[2] = moveDistance;
        break;
      case ViewState::Z:
        move[1] = moveDistance;
        break;
      }
    } else if ( key.IsSameAs(msMoveViewIn) ) {
      switch ( mViewState.GetInPlane() ) {
      case ViewState::X:
        move[0] = moveDistance;
        break;
      case ViewState::Y:
        move[1] = moveDistance;
        break;
      case ViewState::Z:
        move[2] = moveDistance;
        break;
      }
    } else if ( key.IsSameAs(msMoveViewOut) ) {
      switch ( mViewState.GetInPlane() ) {
      case ViewState::X:
        move[0] = -moveDistance;
        break;
      case ViewState::Y:
        move[1] = -moveDistance;
        break;
      case ViewState::Z:
        move[2] = -moveDistance;
        break;
      }
    }

    // Do the move.
    float centerRAS[3];
    float newCenterRAS[3];
    Get2DRASCenter( centerRAS );
    TranslateRASInWindowSpace( centerRAS, move, newCenterRAS );

    Set2DRASCenter( newCenterRAS );

    // Rebuild our label value info because the view has moved.
    RebuildLabelValueInfo( newCenterRAS, "mouse" );

  } else if ( key.IsSameAs(msZoomViewIn) || key.IsSameAs(msZoomViewOut) ) {

    float newZoom = mViewState.GetZoomLevel();
    if ( key.IsSameAs(msZoomViewIn) ) {
      newZoom = mViewState.GetZoomLevel() * 2.0;
    } else if ( key.IsSameAs(msZoomViewOut) ) {
      newZoom = mViewState.GetZoomLevel() / 2.0;
      if ( newZoom < 0.25 ) {
        newZoom = 0.25;
      }
    }
    Set2DZoomLevel( newZoom );

    // Rebuild our label value info because the view has moved.
    float ras[3];
    TranslateWindowToRAS( iWindow, ras );
    RebuildLabelValueInfo( ras, "mouse" );
  }

  // Pass this tool to our layers.
  map<int,int>::iterator tLevelLayerID;
  for ( tLevelLayerID = mLevelLayerIDMap.begin();
        tLevelLayerID != mLevelLayerIDMap.end(); ++tLevelLayerID ) {
    int layerID = (*tLevelLayerID).second;
    try {
      Layer& layer = Layer::FindByID( layerID );
      layer.HandleTool( ras, mViewState, *this, iTool, iInput );
      if ( layer.WantRedisplay() ) {
        RequestRedisplay();
        layer.RedisplayPosted();
      }
    } catch ( runtime_error& e) {
      throw( e );
    } catch (...) {
      DebugOutput( << "Couldn't find layer " << layerID );
    }
  }

  RequestRedisplay();
}

void
ScubaView::DoKeyUp( int[2],
                    InputState&, ScubaToolState& ) {}

void
ScubaView::CalcViewToWindowTransform () {

  Point3<float> N( mViewState.GetPlaneNormal() );
  N = VectorOps::Normalize( N );

  Point3<float> D;
  switch ( mViewState.GetInPlane() ) {
  case ViewState::X:
    D.Set( 1, 0, 0 );
    break;
  case ViewState::Y:
    D.Set( 0, 1, 0 );
    break;
  case ViewState::Z:
    D.Set( 0, 0, 1 );
    break;
  }
  D = VectorOps::Normalize( D );

  double rads = VectorOps::RadsBetweenVectors( N, D );
  if ( mViewState.GetInPlane() == ViewState::X ) {
    rads = -rads;
  }

  Point3<float> axis = VectorOps::Cross( N, D );
  mViewToWindow.MakeRotation( mViewState.GetCenterRAS(),
                              axis.xyz(), rads );

  CalcWorldToWindowTransform();
}

void
ScubaView::CalcWorldToWindowTransform () {

  Transform44 viewToWorld = mWorldToView->Inverse();
  Transform44 tmp = mViewToWindow * viewToWorld;
  mWorldToWindow = tmp;
}

void
ScubaView::TranslateWindowToRAS ( int const iWindow[2], float oRAS[3] ) {

  // At zoom level one every pixel is an RAS point, so we're scaled by
  // that and then offset by the RAS center. We find a 3D point by
  // using our mInPlane and the corresponding mCenterRAS in ViewState.

  int xWindow = iWindow[0];
  float windowRAS[3];
  switch ( mViewState.GetInPlane() ) {
  case ViewState::X:
    windowRAS[0] = mViewState.GetCenterRAS()[0];
    windowRAS[1] = ConvertWindowToRAS( xWindow,
                                       mViewState.GetCenterRAS()[1], mWidth );
    windowRAS[2] = ConvertWindowToRAS(iWindow[1],
                                      mViewState.GetCenterRAS()[2], mHeight );
    break;
  case ViewState::Y:
    if ( mbFlipLeftRightInYZ ) {
      xWindow = mWidth - xWindow;
    }
    windowRAS[0] = ConvertWindowToRAS( xWindow,
                                       mViewState.GetCenterRAS()[0], mWidth );
    windowRAS[1] = mViewState.GetCenterRAS()[1];
    windowRAS[2] = ConvertWindowToRAS( iWindow[1],
                                       mViewState.GetCenterRAS()[2], mHeight );
    break;
  case ViewState::Z:
    if ( mbFlipLeftRightInYZ ) {
      xWindow = mWidth - xWindow;
    }
    windowRAS[0] = ConvertWindowToRAS( xWindow,
                                       mViewState.GetCenterRAS()[0], mWidth );
    windowRAS[1] = ConvertWindowToRAS( iWindow[1],
                                       mViewState.GetCenterRAS()[1], mHeight );
    windowRAS[2] = mViewState.GetCenterRAS()[2];
    break;
  }

  mWorldToWindow.InvMultiplyVector3( windowRAS, oRAS );
}

float
ScubaView::ConvertWindowToRAS ( float iWindow, float iRASCenter,
                                float iWindowDimension ) {

  return ((iWindow - iWindowDimension / 2.0) / mViewState.GetZoomLevel()) +
         iRASCenter;
}

void
ScubaView::TranslateRASToWindow ( float const iRAS[3], int oWindow[2] ) {

  float windowRAS[3];
  mWorldToWindow.MultiplyVector3( iRAS, windowRAS );

  float xWindow = 0, yWindow = 0;
  switch ( mViewState.GetInPlane() ) {
  case ViewState::X:
    xWindow = ConvertRASToWindow( windowRAS[1],
                                  mViewState.GetCenterRAS()[1], mWidth );
    yWindow = ConvertRASToWindow( windowRAS[2],
                                  mViewState.GetCenterRAS()[2], mHeight );
    break;
  case ViewState::Y:
    xWindow = ConvertRASToWindow( windowRAS[0],
                                  mViewState.GetCenterRAS()[0], mWidth );
    yWindow = ConvertRASToWindow( windowRAS[2],
                                  mViewState.GetCenterRAS()[2], mHeight );
    if ( mbFlipLeftRightInYZ ) {
      xWindow = mWidth - xWindow;
    }
    break;
  case ViewState::Z:
    xWindow = ConvertRASToWindow( windowRAS[0],
                                  mViewState.GetCenterRAS()[0], mWidth );
    yWindow = ConvertRASToWindow( windowRAS[1],
                                  mViewState.GetCenterRAS()[1], mHeight );
    if ( mbFlipLeftRightInYZ ) {
      xWindow = mWidth - xWindow;
    }
    break;
  }

  oWindow[0] = (int)floor( xWindow + 0.5 );
  oWindow[1] = (int)floor( yWindow + 0.5 );
}

float
ScubaView::ConvertRASToWindow ( float iRAS, float iRASCenter,
                                float iWindowDimension ) {

  return ((iRAS - iRASCenter) * mViewState.GetZoomLevel()) +
         (iWindowDimension / 2.0);
}

void
ScubaView::TranslateRASInWindowSpace ( float iRAS[3], float iMove[3],
                                       float oRAS[3] ) {

  // Translate the view point into a window point.
  Point3<float> window;
  mViewToWindow.MultiplyVector3( iRAS, window.xyz() );

  // Apply the move.
  window[0] += iMove[0];
  window[1] += iMove[1];
  window[2] += iMove[2];

  // Convert back into view.
  mViewToWindow.InvMultiplyVector3( window.xyz(), oRAS );
}

void
ScubaView::CalcAllViewIntersectionPoints () {

  list<int> viewIDs;
  GetIDList( viewIDs );
  list<int>::iterator tViewID;
  for ( tViewID = viewIDs.begin(); tViewID != viewIDs.end(); ++tViewID ) {
    int viewID = *tViewID;
    if ( viewID != GetID () ) {
      CalcViewIntersectionPoints( viewID );
    }
  }
}

void
ScubaView::CalcViewIntersectionPoints ( int iViewID ) {

  View& view = View::FindByID( iViewID );
  // ScubaView& scubaView = dynamic_cast<ScubaView&>(view);
  ScubaView& scubaView = (ScubaView&)view;

  // Get the dot of our plane normal and its plane normal. If not
  // zero...
  Point3<float> n1( mViewState.GetPlaneNormal() );
  Point3<float> n2;
  scubaView.Get2DPlaneNormal( n2.xyz() );
  if ( !VectorOps::AreVectorsParallel( n1, n2 ) ) {

    // Get p1 and p2, the center RAS points for our plane
    // and their plane.
    Point3<float> p1( mViewState.GetCenterRAS() );
    Point3<float> p2;
    scubaView.Get2DRASCenter( p2.xyz() );

    // p3 is our RAS point for the edge of the window, and
    // n3 is the normal for the 'side' of the viewing 'box',
    // pointing into the middle of the window. This is
    // definitely perpendicular to n1, but we need to check
    // against n2, because if it is, then we need to change
    // the normal to try a plane in the other
    // orientation. i.e. if we're checking the left side
    // plane, we'll change to the top side plane.
    Point2<int> windowTopLeft( 0, 0 );
    Point3<float> p3;
    TranslateWindowToRAS( windowTopLeft.xy(), p3.xyz() );
    Point3<float> n3;
    switch ( Get2DInPlane() ) {
    case ViewState::X:
      n3.Set( 0, 1, 0 );
      break;
    case ViewState::Y:
      n3.Set( 1, 0, 0 );
      break;
    case ViewState::Z:
      n3.Set( 1, 0, 0 );
      break;
    }

    if ( VectorOps::AreVectorsParallel( n2, n3 ) ) {
      switch ( Get2DInPlane() ) {
      case ViewState::X:
        n3.Set( 0, 0, 1 );
        break;
      case ViewState::Y:
        n3.Set( 0, 0, 1 );
        break;
      case ViewState::Z:
        n3.Set( 0, 1, 0 );
        break;
      }
    }

    // Intersect the three planes. This gives us an RAS
    // interesction.
    float p1dn1 = VectorOps::Dot( p1, n1 );
    float p2dn2 = VectorOps::Dot( p2, n2 );
    float p3dn3 = VectorOps::Dot( p3, n3 );
    Point3<float> n2xn3 = VectorOps::Cross( n2, n3 );
    Point3<float> n3xn1 = VectorOps::Cross( n3, n1 );
    Point3<float> n1xn2 = VectorOps::Cross( n1, n2 );
    Point3<float> P_1( p1dn1 * n2xn3 +
                       p2dn2 * n3xn1 +
                       p3dn3 * n1xn2);
    Point3<float> P1 = P_1 / VectorOps::TripleScalar( n1, n2, n3 );

    // Now do the right or bottom plane.
    Point2<int> windowBottomRight( mWidth-1, mHeight-1 );
    TranslateWindowToRAS( windowBottomRight.xy(), p3.xyz() );
    p3dn3 = VectorOps::Dot( p3, n3 );
    Point3<float> P_2( p1dn1 * n2xn3 +
                       p2dn2 * n3xn1 +
                       p3dn3 * n1xn2);
    Point3<float> P2 = P_2 / VectorOps::TripleScalar( n1, n2, n3 );

    // Save the results.
    mViewIDViewIntersectionPointMap[iViewID][0].Set( P1 );
    mViewIDViewIntersectionPointMap[iViewID][1].Set( P2 );
  }
}


void
ScubaView::SetCursor ( float iRAS[3] ) {

  // Set the cursor;
  mCursor.Set( iRAS );

  ScubaViewBroadcaster& broadcaster = ScubaViewBroadcaster::GetBroadcaster();
  broadcaster.SendBroadcast( "cursorChanged", NULL );
}

void
ScubaView::GetCursor ( float oRAS[3] ) {

  // Return the cursor;
  oRAS[0] = mCursor[0];
  oRAS[1] = mCursor[1];
  oRAS[2] = mCursor[2];
}

void
ScubaView::SetNextMarker ( float iRAS[3] ) {

  if ( mcMarkers > 0 ) {

    if ( mNextMarker >= mcMarkers )
      mNextMarker = 0;

    mMarkerRAS[mNextMarker] = iRAS;
    mMarkerVisible[mNextMarker] = true;

    mNextMarker++;

    ScubaViewBroadcaster& broadcaster = ScubaViewBroadcaster::GetBroadcaster();
    broadcaster.SendBroadcast( "markerChanged", NULL );
  }
}

void
ScubaView::HideNearestMarker ( float iRAS[3] ) {

  float closestDistance = 10000;
  int nClosestMarker = -1;
  for ( int nMarker = 0; nMarker < mcMarkers; nMarker++ ) {
    float distance =
      sqrt( (mMarkerRAS[nMarker].x() - iRAS[0]) *
            (mMarkerRAS[nMarker].x() - iRAS[0]) +
            (mMarkerRAS[nMarker].y() - iRAS[1]) *
            (mMarkerRAS[nMarker].y() - iRAS[1]) +
            (mMarkerRAS[nMarker].z() - iRAS[2]) *
            (mMarkerRAS[nMarker].z() - iRAS[2]) );
    if ( distance < closestDistance ) {
      closestDistance = distance;
      nClosestMarker = nMarker;
    }
  }

  if ( -1 != nClosestMarker ) {
    mMarkerVisible[nClosestMarker] = false;

    ScubaViewBroadcaster& broadcaster = ScubaViewBroadcaster::GetBroadcaster();
    broadcaster.SendBroadcast( "markerChanged", NULL );
  }
}

void
ScubaView::SetNumberOfMarkers ( int icMarkers ) {

  if ( icMarkers < mcMarkers ) {
    for ( int nMarker = icMarkers; nMarker < mcMarkers; nMarker++ ) {
      mMarkerVisible[nMarker] = false;
    }
  }

  mcMarkers = icMarkers;

  if ( mNextMarker >= mcMarkers || mNextMarker < 0 )
    mNextMarker = 0;
}

bool
ScubaView::IsNthMarkerVisible ( int inMarker ) {

  if ( inMarker >= 0 && inMarker < mcMarkers ) {
    return mMarkerVisible[inMarker];
  } else {
    throw runtime_error( "Marker index is out of bounds" );
  }
}

void
ScubaView::GetNthMarker ( int inMarker, float oMarkerRAS[3] ) {

  if ( inMarker >= 0 && inMarker < mcMarkers ) {
    oMarkerRAS[0] = mMarkerRAS[inMarker].x();
    oMarkerRAS[1] = mMarkerRAS[inMarker].y();
    oMarkerRAS[2] = mMarkerRAS[inMarker].z();
  } else {
    throw runtime_error( "Marker index is out of bounds" );
  }
}

void
ScubaView::ExportMarkersToControlPointsForVolume ( string ifnControlPoints,
    VolumeCollection& iVolume ) {

  // Make a list out of our visible control points.
  list<Point3<float> > lMarkers;
  int cMarkers = ScubaView::GetNumberOfMarkers();
  for ( int nMarker = 0; nMarker < cMarkers; nMarker++ ) {
    if ( ScubaView::IsNthMarkerVisible( nMarker ) ) {
      float markerRAS[3];
      ScubaView::GetNthMarker( nMarker, markerRAS );
      lMarkers.push_back( Point3<float>( markerRAS ) );
    }
  }

  iVolume.ExportControlPoints( ifnControlPoints, lMarkers );
}

void
ScubaView::ImportMarkersFromControlPointsForVolume ( string ifnControlPoints,
    VolumeCollection& iVolume ) {

  list<Point3<float> > lControlPoints;
  iVolume.ImportControlPoints( ifnControlPoints, lControlPoints );

  // Set the number of markers if we don't have enough.
  int cControlPoints = lControlPoints.size();
  int cMarkers = ScubaView::GetNumberOfMarkers();
  if ( cMarkers < cControlPoints ) {
    ScubaView::SetNumberOfMarkers( cControlPoints );
  }

  list<Point3<float> >::iterator tControlPoint;
  for ( tControlPoint = lControlPoints.begin();
        tControlPoint != lControlPoints.end();
        ++tControlPoint ) {

    ScubaView::SetNextMarker( (*tControlPoint).xyz() );
  }
}


void
ScubaView::SetFlipLeftRightYZ ( bool iFlip ) {

  mbFlipLeftRightInYZ = iFlip;
  RequestRedisplay();
}


void
ScubaView::GetVolumeHistogramInView ( VolumeCollection& iSourceVol,
                                      ScubaROIVolume* iROI,
                                      int icBins,
                                      float& oMinBinValue,
                                      float& oBinIncrement,
                                      map<int,int>& oBinCounts ) {

  // First we need to generate a list of RAS points in the view. We
  // step through pixel by pixel but as we only want to do voxel
  // resolution for the volume, we neeed to only add RAS points that
  // have a unique voxel index. So make a volume the same size as the
  // volume that we use to track which voxels have been added. Then go
  // through all the RAS points on screen, convert to voxel index, and
  // check if they've been added. If not, add them to the RAS list.
  int range[3];
  iSourceVol.GetMRIIndexRange( range );
  Volume3<bool> bAdded( range[0], range[1], range[2], false );

  int window[2];
  Point3<float> RAS;
  list<Point3<float> > RASPoints;
  Point3<int> index;
  for ( window[1] = 0; window[1] < mHeight; window[1]++ ) {
    for ( window[0] = 0; window[0] < mWidth; window[0]++ ) {
      TranslateWindowToRAS( window, RAS.xyz() );
      VolumeLocation loc( iSourceVol.MakeVolumeLocationFromRAS( RAS.xyz() ) );

      if ( iSourceVol.IsInBounds( loc ) ) {
        iSourceVol.RASToMRIIndex( RAS.xyz(), index.xyz() );

        // If not already added, add it.
        if ( !bAdded.Get( index[0], index[1], index[2] ) ) {
          bAdded.Set( index[0], index[1], index[2], true );
          RASPoints.push_back( RAS );
        }
      }
    }
  }

  // Make sure we had some points in the view.
  if ( RASPoints.size() == 0 ) {
    oMinBinValue = 0;
    oBinIncrement = 0;
    oBinCounts.clear();
    return;
  }

  // Now just get a histogram for those points.
  iSourceVol.MakeHistogram ( RASPoints, iROI, icBins,
                             oMinBinValue, oBinIncrement, oBinCounts );
}


void
ScubaView::DoVolumeValueRangeFill ( VolumeCollection& iSourceVol,
				    ScubaROIVolume* iROI,
				    VolumeCollection& iDestVol,
				    vector<ValueRangeFillElement>& lElements ){

  // Do the same walk through of voxels that we do in the get histogram.
  int range[3];
  iSourceVol.GetMRIIndexRange( range );

  int window[2];
  Point3<float> RAS;
  list<Point3<float> > RASPoints;
  Point3<int> index;
  for ( window[1] = 0; window[1] < mHeight; window[1]++ ) {
    for ( window[0] = 0; window[0] < mWidth; window[0]++ ) {
      TranslateWindowToRAS( window, RAS.xyz() );
      VolumeLocation loc( iSourceVol.MakeVolumeLocationFromRAS( RAS.xyz() ) );
      if ( iSourceVol.IsInBounds( loc ) ) {
        iSourceVol.RASToMRIIndex( RAS.xyz(), index.xyz() );

        // If they gave us an ROI to use, make sure this is in the
        // ROI.
        if ( NULL != iROI &&
	     !iROI->IsVoxelSelected( index.xyz() ) )
	  continue;
	
        // We need to see if we should fill this voxel. Go through our
        // list and see if it falls into a range. If so, edit the
        // value in the dest.
        float value = iSourceVol.GetMRINearestValue( loc );
        vector<ValueRangeFillElement>::iterator tElement;
        for ( tElement = lElements.begin();
	      tElement != lElements.end(); ++tElement ) {
          if ( value > tElement->mBegin && value < tElement->mEnd ) {
            iDestVol.SetMRIValue( loc, tElement->mValue );
            break;
          }
        }
      }
    }
  }
}



void
ScubaView::BuildFrameBuffer () {

  // Don't draw if our buffer isn't initialized yet.
  if ( NULL == mBuffer )
    return;

  // Erase the frame buffer.
  memset( mBuffer, 0, mHeight * mWidth * kBytesPerPixel );

  // Go through our draw levels. For each one, get the Layer.
  map<int,int>::iterator tLevelLayerID;
  for ( tLevelLayerID = mLevelLayerIDMap.begin();
        tLevelLayerID != mLevelLayerIDMap.end(); ++tLevelLayerID ) {

    int nLevel = (*tLevelLayerID).first;
    if ( !mLevelVisibilityMap[nLevel] ) {
      // This clears the draw list that we might already have saved
      // for this layer.
      glNewList( mLevelGLListIDMap[nLevel], GL_COMPILE );
      glBegin( GL_POINTS );
      glVertex2d( 0, 0 );
      glEnd();
      glEndList();
      continue;
    }

    int layerID = (*tLevelLayerID).second;
    try {
      Layer& layer = Layer::FindByID( layerID );

      // tell it to draw into our buffer with our view state information.
      layer.DrawIntoBuffer( mBuffer, mWidth, mHeight, mViewState, *this );

      // Start a draw list, tell the layer to draw it's 3d stuff, then
      // close the list.
      glNewList( mLevelGLListIDMap[nLevel], GL_COMPILE );
      layer.DrawIntoGL( mViewState, *this );
      glEndList();

    } catch (...) {
      cerr << "Couldn't find layer " << layerID << endl;
    }
  }
}

void
ScubaView::DrawFrameBuffer () {

#if 0
  glEnable( GL_TEXTURE_2D );
  CheckGLError();

  glTexImage2D( GL_TEXTURE_2D, 0, GL_RGBA,
                mWidth, mHeight, 0,
                GL_RGBA, GL_UNSIGNED_BYTE,
                mBuffer );
  CheckGLError();

  glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP );
  CheckGLError();
  glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP );
  CheckGLError();
  glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST );
  CheckGLError();
  glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST );
  CheckGLError();
  glTexEnvf( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL );
  CheckGLError();

  glBegin( GL_QUADS );
  CheckGLError();

  glTexCoord2f( 0.0f, 0.0f );
  glVertex3f  ( 0.0f, 0.0f, 0.0f );

  glTexCoord2f( 0.0f, 1.0f );
  glVertex3f  ( 0.0f, (float)mHeight, 0.0f );

  glTexCoord2f( 1.0f, 1.0f );
  glVertex3f  ( (float)mWidth, (float)mHeight, 0.0f );

  glTexCoord2f( 1.0f, 0.0f );
  glVertex3f  ( (float)mWidth, 0.0f, 0.0f );

  glEnd();
  glGetError(); // clear error

#endif

  // Get the update bounds.
  int windowUpdateBounds[4];
  mViewState.CopyUpdateRect( windowUpdateBounds );

  // Configure OpenGL so it only draws the subimage defined by the
  // bounds.
  glPixelStorei( GL_UNPACK_ROW_LENGTH, mWidth );
  glPixelStorei( GL_UNPACK_SKIP_PIXELS, windowUpdateBounds[0] );
  glPixelStorei( GL_UNPACK_SKIP_ROWS, windowUpdateBounds[1] );

  // Draw.
  glRasterPos2i( windowUpdateBounds[0], windowUpdateBounds[1] );
  glDrawPixels( windowUpdateBounds[2] - windowUpdateBounds[0],
                windowUpdateBounds[3] - windowUpdateBounds[1],
                GL_RGBA, GL_UNSIGNED_BYTE, mBuffer );

#if 0
  cerr << "Rect : ("
  << windowUpdateBounds[0] << ", " << windowUpdateBounds[1] << ") ("
  << windowUpdateBounds[2] << ", " << windowUpdateBounds[3] << ") "
  << "width " << windowUpdateBounds[2] - windowUpdateBounds[0]
  << ", height " << windowUpdateBounds[3] - windowUpdateBounds[1]
  << endl;

  glLineWidth( 1 );
  glColor3f( 0, 1, 0 );
  glBegin( GL_LINE_STRIP );
  glVertex2d( mViewState.mUpdateRect[0], mViewState.mUpdateRect[1] );
  glVertex2d( mViewState.mUpdateRect[2], mViewState.mUpdateRect[1] );
  glVertex2d( mViewState.mUpdateRect[2], mViewState.mUpdateRect[3] );
  glVertex2d( mViewState.mUpdateRect[0], mViewState.mUpdateRect[3] );
  glVertex2d( mViewState.mUpdateRect[0], mViewState.mUpdateRect[1] );
  glEnd();
#endif

  // Clear the update rect.
  mViewState.ResetUpdateRect();
}

void
ScubaView::BuildOverlay () {

  if ( !mbRebuildOverlayDrawList )
    return;

  // Create a draw list ID if we don't have one yet. Open the overlay
  // display list.
  if ( 0 == mDrawListID ) {
    mDrawListID = glGenLists( 1 );
  }
  glNewList( mDrawListID, GL_COMPILE );


  // Draw the HUD overlay if necessary. We need to take our edge
  // window coords and translate them to RAS coords and draw them on
  // the sides of the screens. Note w'ere only really using one of the
  // coords in each left/right/bottom/top/plane calculation because we
  // don't care about the other dimensions.
  ScubaGlobalPreferences& prefs = ScubaGlobalPreferences::GetPreferences();
  if ( prefs.GetPrefAsBool( ScubaGlobalPreferences::DrawCoordinateOverlay )) {

    int window[2] = {0, 0};
    float ras[3];
    char sXLabel, sYLabel, sZLabel;
    float left, right, top, bottom, plane;
    switch ( mViewState.GetInPlane() ) {
    case ViewState::X:
      sXLabel = 'a';
      sYLabel = 's';
      sZLabel = 'r';
      window[0] = 0;
      TranslateWindowToRAS( window, ras );
      left  = ras[1];
      window[0] = mWidth;
      TranslateWindowToRAS( window, ras );
      right = ras[1];
      window[1] = 0;
      TranslateWindowToRAS( window, ras );
      bottom= ras[2];
      window[1] = mHeight;
      TranslateWindowToRAS( window, ras );
      top   = ras[2];
      plane = ras[0];
      break;
    case ViewState::Y:
      sXLabel = 'r';
      sYLabel = 's';
      sZLabel = 'a';
      window[0] = 0;
      TranslateWindowToRAS( window, ras );
      left  = ras[0];
      window[0] = mWidth;
      TranslateWindowToRAS( window, ras );
      right = ras[0];
      window[1] = 0;
      TranslateWindowToRAS( window, ras );
      bottom= ras[2];
      window[1] = mHeight;
      TranslateWindowToRAS( window, ras );
      top   = ras[2];
      plane = ras[1];
      break;
    case ViewState::Z:
    default:
      sXLabel = 'r';
      sYLabel = 'a';
      sZLabel = 's';
      window[0] = 0;
      TranslateWindowToRAS( window, ras );
      left  = ras[0];
      window[0] = mWidth;
      TranslateWindowToRAS( window, ras );
      right = ras[0];
      window[1] = 0;
      TranslateWindowToRAS( window, ras );
      bottom= ras[1];
      window[1] = mHeight;
      TranslateWindowToRAS( window, ras );
      top   = ras[1];
      plane = ras[2];
      break;
    }

    glColor3f( 1, 1, 1 );

    // Now draw the labels we calc'd before.
    char sLabel[60];
    sprintf( sLabel, "%c%.2f", sXLabel, left );
    glRasterPos2i( 0, mHeight / 2 + 7 );
    for ( int nChar = 0; nChar < (int)strlen( sLabel ); nChar++ ) {
      glutBitmapCharacter( GLUT_BITMAP_8_BY_13, sLabel[nChar] );
    }
    sprintf( sLabel, "%c%.2f", sXLabel, right );
    glRasterPos2i( mWidth - strlen(sLabel)*8, mHeight / 2 + 7 );
    for ( int nChar = 0; nChar < (int)strlen( sLabel ); nChar++ ) {
      glutBitmapCharacter( GLUT_BITMAP_8_BY_13, sLabel[nChar] );
    }
    sprintf( sLabel, "%c%.2f", sYLabel, top );
    glRasterPos2i( mWidth / 2 - (strlen(sLabel)*8 / 2), mHeight-1-13 );
    for ( int nChar = 0; nChar < (int)strlen( sLabel ); nChar++ ) {
      glutBitmapCharacter( GLUT_BITMAP_8_BY_13, sLabel[nChar] );
    }
    sprintf( sLabel, "%c%.2f", sYLabel, bottom );
    glRasterPos2i( mWidth / 2 - (strlen(sLabel)*8 / 2), 4 );
    for ( int nChar = 0; nChar < (int)strlen( sLabel ); nChar++ ) {
      glutBitmapCharacter( GLUT_BITMAP_8_BY_13, sLabel[nChar] );
    }
    sprintf( sLabel, "%c%.2f", sZLabel, plane );
    glRasterPos2i( mWidth - (strlen(sLabel)*8), 4 );
    for ( int nChar = 0; nChar < (int)strlen( sLabel ); nChar++ ) {
      glutBitmapCharacter( GLUT_BITMAP_8_BY_13, sLabel[nChar] );
    }
    if ( mViewState.GetZoomLevel() != 1 ) {
      sprintf( sLabel, "%.2fx", mViewState.GetZoomLevel() );
      glRasterPos2i( 0, mHeight-1-13 );
      for ( int nChar = 0; nChar < (int)strlen( sLabel ); nChar++ ) {
        glutBitmapCharacter( GLUT_BITMAP_8_BY_13, sLabel[nChar] );
      }
    }
  }

  if ( prefs.GetPrefAsBool( ScubaGlobalPreferences::DrawPlaneIntersections )) {

    // Draw our marker color around us.
    glColor3f( mInPlaneMarkerColor[0],
               mInPlaneMarkerColor[1], mInPlaneMarkerColor[2] );
    glLineWidth( 1 );
    glBegin( GL_LINE_STRIP );
    glVertex2d( 2, 2 );
    glVertex2d( mWidth-3, 2 );
    glVertex2d( mWidth-3, mHeight-3 );
    glVertex2d( 2, mHeight-3 );
    glVertex2d( 2, 2 );
    glEnd();


    // For each other visible view...
    list<int> viewIDs;
    GetIDList( viewIDs );
    list<int>::iterator tViewID;
    for ( tViewID = viewIDs.begin(); tViewID != viewIDs.end(); ++tViewID ) {

      int viewID = *tViewID;

      if ( viewID != GetID () ) {

        View& view = View::FindByID( viewID );
        // ScubaView& scubaView = dynamic_cast<ScubaView&>(view);
        ScubaView& scubaView = (ScubaView&)view;

        try {
          if ( scubaView.IsVisibleInFrame() ) {

            Point3<float> P1, P2;
            P1.Set( mViewIDViewIntersectionPointMap[viewID][0] );
            P2.Set( mViewIDViewIntersectionPointMap[viewID][1] );

            // Transform to window points. These are the two points
            // to connect to draw a line.
            Point2<int> drawPoint1;
            Point2<int> drawPoint2;
            TranslateRASToWindow( P1.xyz(), drawPoint1.xy() );
            TranslateRASToWindow( P2.xyz(), drawPoint2.xy() );

            // Get its marker color.
            float color[3];
            scubaView.GetInPlaneMarkerColor( color );

            // If this is the one we're currently moving, draw it thicker.
            if ( viewID == mCurrentMovingViewIntersection ) {
              glLineWidth( 3 );
            } else {
              glLineWidth( 1 );
            }

            // Draw the line.
            glColor3f( color[0], color[1], color[2] );
            glBegin( GL_LINES );
            glVertex2d( drawPoint1.x(), drawPoint1.y() );
            glVertex2d( drawPoint2.x(), drawPoint2.y() );
            glEnd();

            // Now just get the RAS center and draw a little circle there.
            Point3<float> centerRAS;
            Point2<int> centerWindow;
            scubaView.Get2DRASCenter( centerRAS.xyz() );
            TranslateRASToWindow( centerRAS.xyz(), centerWindow.xy() );

            glColor3f( color[0], color[1], color[2] );
            glBegin( GL_LINE_STRIP );
            glVertex2d( centerWindow.x(), centerWindow.y()-4 );
            glVertex2d( centerWindow.x()-4, centerWindow.y() );
            glVertex2d( centerWindow.x(), centerWindow.y()+4 );
            glVertex2d( centerWindow.x()+4, centerWindow.y() );
            glVertex2d( centerWindow.x(), centerWindow.y()-4 );
            glEnd();


          }
        } catch (...) {}
      }
    }
  }

  if ( prefs.GetPrefAsBool( ScubaGlobalPreferences::DrawMarkers )) {

    // Draw our markers.
    float range = GetThroughPlaneIncrement( mViewState.GetInPlane() ) / 2.0;

    if ( mViewState.IsRASVisibleInPlane( mCursor.xyz(), range ) ) {

      int cursorWindow[2];
      TranslateRASToWindow( mCursor.xyz(), cursorWindow );
      glLineWidth( 1 );
      glColor3f( 1,0,0 );
      glBegin( GL_LINES );
      glVertex2d( cursorWindow[0] - 5, cursorWindow[1] );
      glVertex2d( cursorWindow[0] + 6, cursorWindow[1] );
      glVertex2d( cursorWindow[0], cursorWindow[1] - 5 );
      glVertex2d( cursorWindow[0], cursorWindow[1] + 6 );
      glEnd();
    }

    for ( int nMarker = 0; nMarker < mcMarkers; nMarker++ ) {
      if ( mMarkerVisible[nMarker] &&
           mViewState.IsRASVisibleInPlane( mMarkerRAS[nMarker].xyz(), range ) ) {
        int markerWindow[2];
        TranslateRASToWindow( mMarkerRAS[nMarker].xyz(), markerWindow );
        glLineWidth( 1 );
        glColor3f( 0,1,0 );
        glBegin( GL_LINES );
        glVertex2d( markerWindow[0] - 5, markerWindow[1] );
        glVertex2d( markerWindow[0] + 6, markerWindow[1] );
        glVertex2d( markerWindow[0], markerWindow[1] - 5 );
        glVertex2d( markerWindow[0], markerWindow[1] + 6 );
        glEnd();
      }
    }
  }



  // Line range.
  float range = GetThroughPlaneIncrement( mViewState.GetInPlane() ) / 2.0;

  // Drawing paths.
  if ( prefs.GetPrefAsBool( ScubaGlobalPreferences::DrawPaths )) {

    vector<Path<float>*>::const_iterator tPath;
    PathManager& pathMgr = PathManager::GetManager();
    vector<Path<float>*> const& pathList = pathMgr.GetPathList();
    for ( tPath = pathList.begin(); tPath != pathList.end(); ++tPath ) {
      Path<float> const* path = *tPath;
      if ( path->GetNumVertices() > 0 ) {
        Point3<float> beginRAS = path->GetVertexAtIndex( 0 );
        if ( mViewState.IsRASVisibleInPlane( beginRAS.xyz(), range ) ) {

          if ( path->IsSelected() ) {
            glColor3f( 0, 1, 0 );
          } else                     {
            glColor3f( 1, 0, 0 );
          }

          int cVertices = path->GetNumVertices();
          for ( int nCurVertex = 1; nCurVertex < cVertices; nCurVertex++ ) {

            int nBackVertex = nCurVertex - 1;

            Point3<float> curVertex  = path->GetVertexAtIndex( nCurVertex );
            Point3<float> backVertex = path->GetVertexAtIndex( nBackVertex );

            int curWindow[2], backWindow[2];
            TranslateRASToWindow( curVertex.xyz(), curWindow );
            TranslateRASToWindow( backVertex.xyz(), backWindow );

            glLineWidth( 1 );
            glBegin( GL_LINES );
            glVertex2d( backWindow[0], backWindow[1] );
            glVertex2d( curWindow[0], curWindow[1] );
            glEnd();
          }
        }
      }
    }
  }

  glEndList();

  mbRebuildOverlayDrawList = false;
}

void
ScubaView::RebuildLabelValueInfo ( float  iRAS[3],
                                   string isLabel) {

  // Clear existing list.
  mInfoAtRASMap[isLabel].clear();

  // Get the RAS coords into a string and set that label/value.
  stringstream ssRASCoords;
  ssRASCoords << setprecision(2) << iRAS[0] 
	      << " " << iRAS[1] << " " << iRAS[2];

  Layer::InfoAtRAS info;
  info.SetLabel( "RAS" );
  info.SetValue( ssRASCoords.str() );
  info.SetTclCallback( "SetViewRASCursor" );
  info.SetInputFilter( "3sf" );
  mInfoAtRASMap[isLabel].push_back( info );
  info.Clear();

#if 0
  int window[2];
  TranslateRASToWindow( iRAS, window );
  stringstream ssWindowCoords;
  ssWindowCoords << window[0] << " " << window[1];
  info.SetLabel( "Window" );
  info.SetValue( ssWindowCoords.str() );
  mInfoAtRASMap[isLabel].push_back( info );
  info.Clear();
#endif

  // Go through our draw levels. For each one, get the Layer.
  map<int,int>::iterator tLevelLayerID;
  for ( tLevelLayerID = mLevelLayerIDMap.begin();
        tLevelLayerID != mLevelLayerIDMap.end(); ++tLevelLayerID ) {

    int nLevel = (*tLevelLayerID).first;
    int layerID = (*tLevelLayerID).second;

    // If this level is not reporting info, skip it.
    if ( !GetDrawLevelReportInfo(nLevel) )
      continue;

    try {
      Layer& layer = Layer::FindByID( layerID );

      // Ask the layer for info strings at this point.
      layer.GetInfoAtRAS( iRAS, mInfoAtRASMap[isLabel] );
    } catch (...) {
      DebugOutput( << "Couldn't find layer " << layerID );
    }
  }
}

void
ScubaView::GetInPlaneMarkerColor ( float oColor[3] ) {

  oColor[0] = mInPlaneMarkerColor[0];
  oColor[1] = mInPlaneMarkerColor[1];
  oColor[2] = mInPlaneMarkerColor[2];
}

void
ScubaView::DrawOverlay () {

  if ( mbRebuildOverlayDrawList ) {
    BuildOverlay();
  }

  glCallList( mDrawListID );
}

ScubaViewBroadcaster::ScubaViewBroadcaster () :
    Broadcaster( "ScubaViewBroadcaster" ),
    mCurrentBroadcaster( -1 ) {}

View*
ScubaViewFactory::NewView() {
  ScubaView* view = new ScubaView();

  // Notify other views of this new one.
  int id = view->GetID();
  ScubaViewBroadcaster& broadcaster = ScubaViewBroadcaster::GetBroadcaster();
  broadcaster.SendBroadcast( "NewView", (void*)&id );

  return view;
}


ScubaViewBroadcaster&
ScubaViewBroadcaster::GetBroadcaster () {
  static ScubaViewBroadcaster* sBroadcaster = NULL;
  if ( NULL == sBroadcaster ) {
    sBroadcaster = new ScubaViewBroadcaster();
  }
  return *sBroadcaster;
}

void
ScubaViewBroadcaster::SendBroadcast ( std::string isMessage, void* iData ) {

  // If this is a message for linked views, we want to basically put a
  // mutex on broadcasting.
  if ( isMessage == "2DRASCenterChanged" ||
       isMessage == "2DZoomLevelChanged" ||
       isMessage == "2DInPlaneChanged" ) {

    // If -1, no current broadcaster. Save this one's ID and pass the
    // message along. If it is -1, don't pass the broadcast.
    if ( mCurrentBroadcaster == -1 ) {

      int viewID = *(int*)iData;
      mCurrentBroadcaster = viewID;

      Broadcaster::SendBroadcast( isMessage, iData );

      mCurrentBroadcaster = -1;
    }


    // Pass on other kinds of messages.
  } else {

    Broadcaster::SendBroadcast( isMessage, iData );
  }
}

ScubaViewStaticTclListener::ScubaViewStaticTclListener () {

  TclCommandManager& commandMgr = TclCommandManager::GetManager();
  commandMgr.AddCommand( *this, "SetNumberOfViewMarkers", 1, "numMarkers",
                         "Sets the number of view markers." );
  commandMgr.AddCommand( *this, "GetNumberOfViewMarkers", 0, "",
                         "Returns the number of view markers." );
  commandMgr.AddCommand( *this, "ExportMarkersToControlPoints", 2,
                         "collectionID fileName",
                         "Writes the markers to a control.dat file using "
                         "a volume collection to transform them." );
  commandMgr.AddCommand( *this, "ImportMarkersFromControlPoints", 2,
                         "collectionID fileName",
                         "Imports markers from a control.dat file using "
                         "a volume collection to transform them.." );
  commandMgr.AddCommand( *this, "SetViewRASCursor", 3, "x y z",
                         "Sets the cursor in RAS coords." );
  commandMgr.AddCommand( *this, "GetViewRASCursor", 0, "",
                         "Returns the cursor in RAS coords in a list of "
                         "x y z coords." );
}

ScubaViewStaticTclListener::~ScubaViewStaticTclListener () {}

TclCommandListener::TclCommandResult
ScubaViewStaticTclListener::DoListenToTclCommand ( char* isCommand,
    int, char** iasArgv ) {

  // SetNumberOfViewMarkers <number>
  if ( 0 == strcmp( isCommand, "SetNumberOfViewMarkers" ) ) {
    int cMarkers = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad number of markers";
      return error;
    }

    ScubaView::SetNumberOfMarkers( cMarkers );
  }

  // GetNumberOfViewMarkers
  if ( 0 == strcmp( isCommand, "GetNumberOfViewMarkers" ) ) {

    sReturnFormat = "i";
    stringstream ssReturnValues;
    ssReturnValues << ScubaView::GetNumberOfMarkers();
    sReturnValues = ssReturnValues.str();
  }

  // ImportMarkersFromControlPoints <collectionID> <fileName>
  if ( 0 == strcmp( isCommand, "ImportMarkersFromControlPoints" ) ) {
    int collectionID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad collection ID";
      return error;
    }

    // Find a volume collection.
    VolumeCollection* vol = NULL;
    try {
      DataCollection* col = &DataCollection::FindByID( collectionID );
      //    VolumeCollection* vol = dynamic_cast<VolumeCollection*>(col);
      vol = (VolumeCollection*)col;
    } catch (...) {
      throw runtime_error( "Couldn't find volume or data collection "
                           "wasn't a volume.." );
    }

    string fnControlPoints = iasArgv[2];
    ScubaView::ImportMarkersFromControlPointsForVolume( fnControlPoints,
        *vol );

  }

  // ExportMarkersToControlPoints <collectionID> <fileName>
  if ( 0 == strcmp( isCommand, "ExportMarkersToControlPoints" ) ) {
    int collectionID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad collection ID";
      return error;
    }

    // Find a volume collection.
    VolumeCollection* vol = NULL;
    try {
      DataCollection* col = &DataCollection::FindByID( collectionID );
      //    VolumeCollection* vol = dynamic_cast<VolumeCollection*>(col);
      vol = (VolumeCollection*)col;
    } catch (...) {
      throw runtime_error( "Couldn't find volume or data collection "
                           "wasn't a volume.." );
    }

    string fnControlPoints = iasArgv[2];
    ScubaView::ExportMarkersToControlPointsForVolume( fnControlPoints, *vol );
  }

  // SetViewRASCursor <X> <Y> <Z>
  if ( 0 == strcmp( isCommand, "SetViewRASCursor" ) ) {

    float x = (float) strtod( iasArgv[1], (char**)NULL );
    if ( ERANGE == errno ) {
      sResult = "bad x coordinate";
      return error;
    }
    float y = (float) strtod( iasArgv[2], (char**)NULL );
    if ( ERANGE == errno ) {
      sResult = "bad y coordinate";
      return error;
    }
    float z = (float) strtod( iasArgv[3], (char**)NULL );
    if ( ERANGE == errno ) {
      sResult = "bad z coordinate";
      return error;
    }

    float ras[3];
    ras[0] = x;
    ras[1] = y;
    ras[2] = z;
    ScubaView::SetCursor( ras );

  }

  // GetViewRASCursor <viewID> <X> <Y> <Z>
  if ( 0 == strcmp( isCommand, "GetViewRASCursor" ) ) {

    float cursor[3];
    ScubaView::GetCursor( cursor );

    stringstream ssReturn;
    sReturnFormat = "Lfffl";
    ssReturn << cursor[0] << " " << cursor[1] << " " << cursor[2];
    sReturnValues = ssReturn.str();
    return ok;
  }

  return ok;
}



