#ifndef ScubaROIVolume_h
#define ScubaROIVolume_h

#include "ScubaROI.h"

class ScubaROIVolume : public ScubaROI {
  
  friend class ScubaROIVolumeTester;

 public:

  ScubaROIVolume ();
  ~ScubaROIVolume ();

  void SetROIBounds ( int const iBounds[3] );
  void GetROIBounds ( int oBounds[3] ) const;

  void SelectVoxel ( int const iVoxel[3] );
  void UnselectVoxel ( int const iVoxel[3] );

  bool IsVoxelSelected ( int const iVoxel[3] ) const;

  int NumSelectedVoxels () const { return mcSelectedVoxels; }

 protected:

  int mBounds[3];
  bool*** mVoxels;
  int mcSelectedVoxels;
};



#endif
