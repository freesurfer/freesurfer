/**
 * @file  ScubaROIVolume.cpp
 * @brief Implementation of ScubaROI representing a volume
 *
 * This is an ROI whose elements can be selected and unselected via
 * voxel coordinates, used by VolumeCollection.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:38 $
 *    $Revision: 1.12 $
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


#include <limits>
#include <errno.h>
#include "string_fixed.h"
#include <stdexcept>
#include "ScubaROIVolume.h"

using namespace std;

ScubaROIVolume::ScubaROIVolume () :
  mVoxels( NULL ),
  mcSelectedVoxels( 0 ),
  mbDirtyList( true ) {

  // Init selected bounds cache.
  mSelectedBounds[0] = numeric_limits<int>::max();
  mSelectedBounds[1] = numeric_limits<int>::min();
  mSelectedBounds[2] = numeric_limits<int>::max();
  mSelectedBounds[3] = numeric_limits<int>::min();
  mSelectedBounds[4] = numeric_limits<int>::max();
  mSelectedBounds[5] = numeric_limits<int>::min();
}

ScubaROIVolume::~ScubaROIVolume () {

  delete mVoxels;
}

void
ScubaROIVolume::SetROIBounds ( int const iBounds[3] ) {

  if ( iBounds[0] <= 0 || iBounds[1] <= 0 || iBounds[2] <= 0 )
    throw runtime_error( "out of bounds" );

  // Delete existing volume and allocate a new one..
  delete mVoxels;
  mVoxels = new Volume3<bool>( iBounds[0], iBounds[1], iBounds[2], false );

  // Init selected bounds cache.
  mSelectedBounds[0] = numeric_limits<int>::max();
  mSelectedBounds[1] = numeric_limits<int>::min();
  mSelectedBounds[2] = numeric_limits<int>::max();
  mSelectedBounds[3] = numeric_limits<int>::min();
  mSelectedBounds[4] = numeric_limits<int>::max();
  mSelectedBounds[5] = numeric_limits<int>::min();

  // No selected voxels.
  mcSelectedVoxels = 0;

  // ROI is changed.
  ROIChanged();
}

void
ScubaROIVolume::GetROIBounds ( int oBounds[3] ) const {

  mVoxels->GetBounds( oBounds[0], oBounds[1], oBounds[2] );
}

void
ScubaROIVolume::SelectVoxel ( int const iVoxel[3] ) {

  // If not selected...
  if ( !mVoxels->Get( iVoxel[0], iVoxel[1], iVoxel[2] ) ) {

    // Inc the count and set the flag. Call our changed function.
    mcSelectedVoxels++;
    mVoxels->Set( iVoxel[0], iVoxel[1], iVoxel[2], true );
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

  // If selected...
  if ( mVoxels->Get( iVoxel[0], iVoxel[1], iVoxel[2] ) ) {

    // Dec the count and set the flag. Call our changed function.
    mcSelectedVoxels--;
    mVoxels->Set( iVoxel[0], iVoxel[1], iVoxel[2], false );

    ROIChanged();
  }
}

bool
ScubaROIVolume::IsVoxelSelected ( int const iVoxel[3] ) const {

  // Return the voxel of the flag.
  return mVoxels->Get( iVoxel[0], iVoxel[1], iVoxel[2] );
}

list<Point3<int> >
ScubaROIVolume::GetSelectedVoxelList () const {

  // If the flag is dirty, iterate over the selected bounds and
  // rebuild the list.
  if ( mbDirtyList ) {

    mlSelectedVoxels.clear();

    for ( int nZ = mSelectedBounds[4]; nZ <= mSelectedBounds[5]; nZ++ )
      for ( int nY = mSelectedBounds[2]; nY <= mSelectedBounds[3]; nY++ )
        for ( int nX = mSelectedBounds[0]; nX <= mSelectedBounds[1]; nX++ )
          if ( mVoxels->Get( nX, nY, nZ ) ) {
            Point3<int> voxel( nX, nY, nZ );
            mlSelectedVoxels.push_back( voxel );
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
