#include <stdexcept>
#include "ScubaROIVolume.h"

using namespace std;

ScubaROIVolume::ScubaROIVolume () {

  mBounds[0] = mBounds[1] = mBounds[2] = 0;
  mVoxels = NULL;
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
ScubaROIVolume::SetROIBounds ( int iBounds[3] ) {

  if( iBounds[0] <= 0 || iBounds[1] <= 0 || iBounds[2] <= 0 ) 
    throw runtime_error( "out of bounds" );

  mBounds[0] = iBounds[0];
  mBounds[1] = iBounds[1];
  mBounds[2] = iBounds[2];

  if( NULL != mVoxels ) {
    for( int nZ = 0; nZ < mBounds[2]; nZ++ ) {
      for( int nY = 0; nY < mBounds[1]; nY++ ) {
	free( mVoxels[nZ][nY] );
      }
      free( mVoxels[nZ] );
    }
    free( mVoxels );
  }

  mVoxels = (bool***) calloc( mBounds[2], sizeof(bool**) );
  for( int nZ = 0; nZ < mBounds[2]; nZ++ ) {
    mVoxels[nZ] = (bool**) calloc( mBounds[1], sizeof(bool*) );
    for( int nY = 0; nY < mBounds[1]; nY++ ) {
      mVoxels[nZ][nY] = (bool*) calloc( mBounds[0], sizeof(bool) );
    }
  }

}

void
ScubaROIVolume::SelectVoxel ( int iVoxel[3] ) {

  if( iVoxel[0] < 0 || iVoxel[0] >= mBounds[0] ||
      iVoxel[1] < 0 || iVoxel[1] >= mBounds[1] ||
      iVoxel[2] < 0 || iVoxel[2] >= mBounds[2] ) {
    throw runtime_error( "out of bounds" );
  }    
    
  mVoxels[iVoxel[2]][iVoxel[1]][iVoxel[0]] = true;
}

void 
ScubaROIVolume::UnselectVoxel ( int iVoxel[3] ) {

  if( iVoxel[0] < 0 || iVoxel[0] >= mBounds[0] ||
      iVoxel[1] < 0 || iVoxel[1] >= mBounds[1] ||
      iVoxel[2] < 0 || iVoxel[2] >= mBounds[2] ) {
    throw runtime_error( "out of bounds" );
  }    

  mVoxels[iVoxel[2]][iVoxel[1]][iVoxel[0]] = false;
}

bool
ScubaROIVolume::IsVoxelSelected ( int iVoxel[3] ) {

  if( iVoxel[0] < 0 || iVoxel[0] >= mBounds[0] ||
      iVoxel[1] < 0 || iVoxel[1] >= mBounds[1] ||
      iVoxel[2] < 0 || iVoxel[2] >= mBounds[2] ) {
    throw runtime_error( "out of bounds" );
  }    

  return mVoxels[iVoxel[2]][iVoxel[1]][iVoxel[0]];
}

