#ifndef ScubaROIVolume_h
#define ScubaROIVolume_h

#include "ScubaROI.h"

class ScubaROIVolume : public ScubaROI {
  
  friend class ScubaROIVolumeTester;

 public:

  ScubaROIVolume ();
  ~ScubaROIVolume ();

  void SetROIBounds ( int iBounds[3] );
  void SelectVoxel ( int iVoxel[3] );
  void UnselectVoxel ( int iVoxel[3] );

  bool IsVoxelSelected ( int iVoxel[3] );

 protected:

  int mBounds[3];
  bool*** mVoxels;
};



#endif
