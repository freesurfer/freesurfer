#include <errno.h>
#include "string_fixed.h"
#include <stdexcept>
#include "ScubaROIVolume.h"

using namespace std;

ScubaROIVolume::ScubaROIVolume () {

  mBounds[0] = mBounds[1] = mBounds[2] = 0;
  mVoxels = NULL;
  mcSelectedVoxels = 0;
}

ScubaROIVolume::~ScubaROIVolume () {

  if( NULL != mVoxels ) {
    for( int nZ = 0; nZ < mBounds[2]; nZ++ ) {
      for( int nY = 0; nY < mBounds[1]; nY++ ) {
	free( mVoxels[nZ][nY] );
      }
      free( mVoxels[nZ] );
    }
    free( mVoxels );
  }
}

void 
ScubaROIVolume::SetROIBounds ( int const iBounds[3] ) {

  if( iBounds[0] <= 0 || iBounds[1] <= 0 || iBounds[2] <= 0 ) 
    throw runtime_error( "out of bounds" );

  if( NULL != mVoxels ) {
    for( int nZ = 0; nZ < mBounds[2]; nZ++ ) {
      for( int nY = 0; nY < mBounds[1]; nY++ ) {
	free( mVoxels[nZ][nY] );
      }
      free( mVoxels[nZ] );
    }
    free( mVoxels );
  }

  mBounds[0] = iBounds[0];
  mBounds[1] = iBounds[1];
  mBounds[2] = iBounds[2];

  mVoxels = (bool***) calloc( mBounds[2], sizeof(bool**) );
  for( int nZ = 0; nZ < mBounds[2]; nZ++ ) {
    mVoxels[nZ] = (bool**) calloc( mBounds[1], sizeof(bool*) );
    for( int nY = 0; nY < mBounds[1]; nY++ ) {
      mVoxels[nZ][nY] = (bool*) calloc( mBounds[0], sizeof(bool) );
    }
  }

  mcSelectedVoxels = 0;
}

void 
ScubaROIVolume::GetROIBounds ( int oBounds[3] ) const {

  oBounds[0] = mBounds[0];
  oBounds[1] = mBounds[1];
  oBounds[2] = mBounds[2];
}

void
ScubaROIVolume::SelectVoxel ( int const iVoxel[3] ) {

  if( iVoxel[0] < 0 || iVoxel[0] >= mBounds[0] ||
      iVoxel[1] < 0 || iVoxel[1] >= mBounds[1] ||
      iVoxel[2] < 0 || iVoxel[2] >= mBounds[2] ) {
    throw runtime_error( "out of bounds" );
  }    
    
  if( !mVoxels[iVoxel[2]][iVoxel[1]][iVoxel[0]] ) {

    mcSelectedVoxels++;
    mVoxels[iVoxel[2]][iVoxel[1]][iVoxel[0]] = true;

    ROIChanged();
  }
}

void 
ScubaROIVolume::UnselectVoxel ( int const iVoxel[3] ) {

  if( iVoxel[0] < 0 || iVoxel[0] >= mBounds[0] ||
      iVoxel[1] < 0 || iVoxel[1] >= mBounds[1] ||
      iVoxel[2] < 0 || iVoxel[2] >= mBounds[2] ) {
    throw runtime_error( "out of bounds" );
  }    

  if( mVoxels[iVoxel[2]][iVoxel[1]][iVoxel[0]] ) {

    mcSelectedVoxels--;
    mVoxels[iVoxel[2]][iVoxel[1]][iVoxel[0]] = false;

    ROIChanged();
  }
}

bool
ScubaROIVolume::IsVoxelSelected ( int const iVoxel[3] ) const {

  if( iVoxel[0] < 0 || iVoxel[0] >= mBounds[0] ||
      iVoxel[1] < 0 || iVoxel[1] >= mBounds[1] ||
      iVoxel[2] < 0 || iVoxel[2] >= mBounds[2] ) {
    stringstream ssErr;
    ssErr << "Out of bounds voxel " << iVoxel[0] << ", "
	  << iVoxel[1] << ", " << iVoxel[2] << endl;
    throw runtime_error( ssErr.str() );
  }    

  return mVoxels[iVoxel[2]][iVoxel[1]][iVoxel[0]];
}

