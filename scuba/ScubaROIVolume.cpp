/**
 * @file  ScubaROIVolume.cpp
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 02:09:14 $
 *    $Revision: 1.9 $
 *
 * Copyright (C) 2002-2007,
 * The General Hospital Corporation (Boston, MA). 
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
 *
 */


#include <errno.h>
#include "string_fixed.h"
#include <stdexcept>
#include "ScubaROIVolume.h"

using namespace std;

ScubaROIVolume::ScubaROIVolume () {

  mBounds[0] = mBounds[1] = mBounds[2] = 0;
  mVoxels = NULL;
  mcSelectedVoxels = 0;
  mbDirtyList = true;

  // Init selected bounds cache.
  mSelectedBounds[0] = 99999;
  mSelectedBounds[1] = -99999;
  mSelectedBounds[2] = 99999;
  mSelectedBounds[3] = -99999;
  mSelectedBounds[4] = 99999;
  mSelectedBounds[5] = -99999;
}

ScubaROIVolume::~ScubaROIVolume () {

  if ( NULL != mVoxels ) {
    for ( int nZ = 0; nZ < mBounds[2]; nZ++ ) {
      for ( int nY = 0; nY < mBounds[1]; nY++ ) {
        free( mVoxels[nZ][nY] );
      }
      free( mVoxels[nZ] );
    }
    free( mVoxels );
  }
}

void
ScubaROIVolume::SetROIBounds ( int const iBounds[3] ) {

  if ( iBounds[0] <= 0 || iBounds[1] <= 0 || iBounds[2] <= 0 )
    throw runtime_error( "out of bounds" );

  // Delete existing bounds.
  if ( NULL != mVoxels ) {
    for ( int nZ = 0; nZ < mBounds[2]; nZ++ ) {
      for ( int nY = 0; nY < mBounds[1]; nY++ ) {
        free( mVoxels[nZ][nY] );
      }
      free( mVoxels[nZ] );
    }
    free( mVoxels );
  }

  // Save the new bounds.
  mBounds[0] = iBounds[0];
  mBounds[1] = iBounds[1];
  mBounds[2] = iBounds[2];

  // Init selected bounds cache.
  mSelectedBounds[0] = 99999;
  mSelectedBounds[1] = -99999;
  mSelectedBounds[2] = 99999;
  mSelectedBounds[3] = -99999;
  mSelectedBounds[4] = 99999;
  mSelectedBounds[5] = -99999;

  // Allocate new volume.
  mVoxels = (bool***) calloc( mBounds[2], sizeof(bool**) );
  for ( int nZ = 0; nZ < mBounds[2]; nZ++ ) {
    mVoxels[nZ] = (bool**) calloc( mBounds[1], sizeof(bool*) );
    for ( int nY = 0; nY < mBounds[1]; nY++ ) {
      mVoxels[nZ][nY] = (bool*) calloc( mBounds[0], sizeof(bool) );
    }
  }

  // No selected voxels.
  mcSelectedVoxels = 0;

  // ROI is changed.
  ROIChanged();
}

void
ScubaROIVolume::GetROIBounds ( int oBounds[3] ) const {

  oBounds[0] = mBounds[0];
  oBounds[1] = mBounds[1];
  oBounds[2] = mBounds[2];
}

void
ScubaROIVolume::SelectVoxel ( int const iVoxel[3] ) {

  // Check the bounds.
  if ( iVoxel[0] < 0 || iVoxel[0] >= mBounds[0] ||
       iVoxel[1] < 0 || iVoxel[1] >= mBounds[1] ||
       iVoxel[2] < 0 || iVoxel[2] >= mBounds[2] ) {
    throw runtime_error( "out of bounds" );
  }

  // If not selected...
  if ( !mVoxels[iVoxel[2]][iVoxel[1]][iVoxel[0]] ) {

    // Inc the count and set the flag. Call our changed function.
    mcSelectedVoxels++;
    mVoxels[iVoxel[2]][iVoxel[1]][iVoxel[0]] = true;
    ROIChanged();

    // Change the size of our selected bounds cache if we need to.
    if ( iVoxel[0] < mSelectedBounds[0] ) mSelectedBounds[0] = iVoxel[0];
    if ( iVoxel[0] > mSelectedBounds[1] ) mSelectedBounds[1] = iVoxel[0];
    if ( iVoxel[1] < mSelectedBounds[2] ) mSelectedBounds[2] = iVoxel[1];
    if ( iVoxel[1] > mSelectedBounds[3] ) mSelectedBounds[3] = iVoxel[1];
    if ( iVoxel[2] < mSelectedBounds[4] ) mSelectedBounds[4] = iVoxel[2];
    if ( iVoxel[2] > mSelectedBounds[5] ) mSelectedBounds[5] = iVoxel[2];
  }
}

void
ScubaROIVolume::UnselectVoxel ( int const iVoxel[3] ) {

  // Check the bounds.
  if ( iVoxel[0] < 0 || iVoxel[0] >= mBounds[0] ||
       iVoxel[1] < 0 || iVoxel[1] >= mBounds[1] ||
       iVoxel[2] < 0 || iVoxel[2] >= mBounds[2] ) {
    throw runtime_error( "out of bounds" );
  }

  // If selected...
  if ( mVoxels[iVoxel[2]][iVoxel[1]][iVoxel[0]] ) {

    // Dec the count and set the flag. Call our changed function.
    mcSelectedVoxels--;
    mVoxels[iVoxel[2]][iVoxel[1]][iVoxel[0]] = false;

    ROIChanged();
  }
}

bool
ScubaROIVolume::IsVoxelSelected ( int const iVoxel[3] ) const {

  // Check the bounds.
  if ( iVoxel[0] < 0 || iVoxel[0] >= mBounds[0] ||
       iVoxel[1] < 0 || iVoxel[1] >= mBounds[1] ||
       iVoxel[2] < 0 || iVoxel[2] >= mBounds[2] ) {
    stringstream ssErr;
    ssErr << "Out of bounds voxel " << iVoxel[0] << ", "
    << iVoxel[1] << ", " << iVoxel[2] << endl;
    throw runtime_error( ssErr.str() );
  }

  // Return the voxel of the flag.
  return mVoxels[iVoxel[2]][iVoxel[1]][iVoxel[0]];
}

list<Point3<int> >
ScubaROIVolume::GetSelectedVoxelList () {

  // If the flag is dirty, iterate over the selected bounds and
  // rebuild the list.
  if ( mbDirtyList ) {

    mlSelectedVoxels.clear();

    for ( int nZ = mSelectedBounds[4]; nZ <= mSelectedBounds[5]; nZ++ ) {
      for ( int nY = mSelectedBounds[2]; nY <= mSelectedBounds[3]; nY++ ) {
        for ( int nX = mSelectedBounds[0]; nX <= mSelectedBounds[1]; nX++ ) {
          if ( mVoxels[nZ][nY][nX] ) {
            Point3<int> voxel( nX, nY, nZ );
            mlSelectedVoxels.push_back( voxel );
          }
        }
      }
    }

    mbDirtyList = false;
  }

  return mlSelectedVoxels;
}

void
ScubaROIVolume::ROIChanged () {

  // Set our flag so we rebuild the list when GetSelectedVoxelList()
  // is called.
  mbDirtyList = true;

  // Call superclass function.
  ScubaROI::ROIChanged();
}
