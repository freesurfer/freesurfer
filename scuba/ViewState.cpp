/**
 * @file  ViewState.cpp
 * @brief View orientation information for ScubaViews
 *
 * Settable and accesible information for a
 * View. IsRASVisibleInPlane() is used to tell when an RAS point is in
 * the current view given a range (usually half the voxel size in the
 * inplane direction). Also manages update rects.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:39 $
 *    $Revision: 1.16 $
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


#include <iostream>
#include <math.h>
#include "ViewState.h"
#include "macros.h"

using namespace std;

ViewState::ViewState () {
  mCenterRAS[0] = mCenterRAS[1] = mCenterRAS[2] = 0;
  mZoomLevel = 1;
  mInPlane = X;
  mPlaneNormal[0] = 1;
  mPlaneNormal[1] = 0;
  mPlaneNormal[2] = 0;
}

std::ostream& operator << ( std::ostream& os, ViewState& iInput ) {
  os << "ViewState CenterRAS: " << iInput.GetCenterRAS()[0] << ", "
  << iInput.GetCenterRAS()[1] << ", " << iInput.GetCenterRAS()[2]
  << " ZoomLevel: " << iInput.GetZoomLevel() << " InPlane: "
  << iInput.GetInPlane();
  return os;
}

void
ViewState::GetCenterRAS ( float oCenterRAS[3] ) {

  oCenterRAS[0] = mCenterRAS[0];
  oCenterRAS[1] = mCenterRAS[1];
  oCenterRAS[2] = mCenterRAS[2];
}

void
ViewState::GetPlaneNormal ( float oPlaneNormal[3] ) {

  oPlaneNormal[0] = mPlaneNormal[0];
  oPlaneNormal[1] = mPlaneNormal[1];
  oPlaneNormal[2] = mPlaneNormal[2];
}

void
ViewState::SetCenterRAS ( float iCenterRASX,
                          float iCenterRASY, float iCenterRASZ ) {

  mCenterRAS[0] = iCenterRASX;
  mCenterRAS[1] = iCenterRASY;
  mCenterRAS[2] = iCenterRASZ;
}

void
ViewState::SetCenterRAS ( float iCenterRAS[3] ) {

  mCenterRAS[0] = iCenterRAS[0];
  mCenterRAS[1] = iCenterRAS[1];
  mCenterRAS[2] = iCenterRAS[2];
}

void
ViewState::SetZoomLevel ( float iZoomLevel ) {

  mZoomLevel = iZoomLevel;
}

void
ViewState::SetPlaneNormal ( float iPlaneX, float iPlaneY, float iPlaneZ ) {

  mPlaneNormal[0] = iPlaneX;
  mPlaneNormal[1] = iPlaneY;
  mPlaneNormal[2] = iPlaneZ;
}

void
ViewState::SetPlaneNormal ( float iPlaneNormal[3] ) {

  mPlaneNormal[0] = iPlaneNormal[0];
  mPlaneNormal[1] = iPlaneNormal[1];
  mPlaneNormal[2] = iPlaneNormal[2];
}

void
ViewState::SetInPlane ( Plane iInPlane ) {

  mInPlane = iInPlane;
}

void
ViewState::SetBufferWidth ( int iWidth ) {

  mBufferWidth = iWidth;
}

void
ViewState::SetBufferHeight ( int iHeight ) {

  mBufferHeight = iHeight;
}


bool
ViewState::IsRASVisibleInPlane ( float iRAS[3], float iRange ) {

  float rasCoord;
  float viewCoord;

  switch ( mInPlane ) {
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
ViewState::SetFrom ( ViewState& ioViewState ) {

  SetCenterRAS( ioViewState.GetCenterRAS() );
  SetZoomLevel( ioViewState.GetZoomLevel() );
  SetPlaneNormal( ioViewState.GetPlaneNormal() );
  SetInPlane( ioViewState.GetInPlane() );
  SetBufferWidth( ioViewState.GetBufferWidth() );
  SetBufferHeight( ioViewState.GetBufferHeight() );
  int windowUpdateBounds[4];
  ioViewState.CopyUpdateRect( windowUpdateBounds );
  mUpdateRect[0] = windowUpdateBounds[0];
  mUpdateRect[1] = windowUpdateBounds[1];
  mUpdateRect[2] = windowUpdateBounds[2];
  mUpdateRect[3] = windowUpdateBounds[3];
}

bool
ViewState::IsSameAs ( ViewState& iViewState ) {

  if ( iViewState.GetCenterRAS()[0] == GetCenterRAS()[0] &&
       iViewState.GetCenterRAS()[1] == GetCenterRAS()[1] &&
       iViewState.GetCenterRAS()[2] == GetCenterRAS()[2] &&

       iViewState.GetZoomLevel() == GetZoomLevel() &&

       iViewState.GetPlaneNormal()[0] == GetPlaneNormal()[0] &&
       iViewState.GetPlaneNormal()[1] == GetPlaneNormal()[1] &&
       iViewState.GetPlaneNormal()[2] == GetPlaneNormal()[2] &&

       iViewState.GetInPlane() == GetInPlane() &&

       iViewState.GetBufferWidth() == GetBufferWidth() &&

       iViewState.GetBufferHeight() == GetBufferHeight() ) {

    return true;

  } else {

    return false;
  }
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

  if ( iWindowLeft != -1 ) {
    if ( mUpdateRect[0] == -1 )
      mUpdateRect[0] = MAX( 0, iWindowLeft );
    else
      mUpdateRect[0] = MAX( 0, MIN( mUpdateRect[0], iWindowLeft ) );
  }

  if ( iWindowTop != -1 ) {
    if ( mUpdateRect[1] == -1 )
      mUpdateRect[1] = MAX( 0, iWindowTop );
    else
      mUpdateRect[1] = MAX( 0, MIN( mUpdateRect[1], iWindowTop ) );
  }

  if ( iWindowRight != -1 ) {
    if ( mUpdateRect[2] == -1 )
      mUpdateRect[2] = MIN( mBufferWidth-1, iWindowRight );
    else
      mUpdateRect[2] = MIN( mBufferWidth-1, MAX( mUpdateRect[2], iWindowRight ));
  }

  if ( iWindowBottom != -1 ) {
    if ( mUpdateRect[3] == -1 )
      mUpdateRect[3] = MIN( mBufferHeight-1, iWindowBottom );
    else
      mUpdateRect[3] = MIN( mBufferHeight-1, MAX( mUpdateRect[3],iWindowBottom));
  }
}

void
ViewState::UpdateEntireRect () {
  mUpdateRect[0] = 0;
  mUpdateRect[1] = 0;
  mUpdateRect[2] = mBufferWidth-1;
  mUpdateRect[3] = mBufferHeight-1;
}

void
ViewState::CopyUpdateRect ( int& oWindowLeft,  int& oWindowTop,
                            int& oWindowRight, int& oWindowBottom ) {

  if ( mUpdateRect[0] == -1 )
    oWindowLeft = 0;
  else
    oWindowLeft = mUpdateRect[0];

  if ( mUpdateRect[1] == -1 )
    oWindowTop = 0;
  else
    oWindowTop = mUpdateRect[1];

  if ( mUpdateRect[2] == -1 )
    oWindowRight = mBufferWidth-1;
  else
    oWindowRight = mUpdateRect[2];

  if ( mUpdateRect[3] == -1 )
    oWindowBottom = mBufferHeight-1;
  else
    oWindowBottom = mUpdateRect[3];
}

void
ViewState::CopyUpdateRect ( int oWindowUpdate[4] ) {

  CopyUpdateRect( oWindowUpdate[0], oWindowUpdate[1],
                  oWindowUpdate[2], oWindowUpdate[3] );
}
