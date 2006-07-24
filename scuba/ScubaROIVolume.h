#ifndef ScubaROIVolume_h
#define ScubaROIVolume_h

#include "Point3.h"
#include "ScubaROI.h"

class ScubaROIVolume : public ScubaROI {
  
  friend class ScubaROIVolumeTester;

 public:

  ScubaROIVolume ();
  ~ScubaROIVolume ();

  // Set the bounds of the volume. Deletes the existing voxel volume
  // if it exists.
  void SetROIBounds ( int const iBounds[3] );

  // Get the bounds of the volume.
  void GetROIBounds ( int oBounds[3] ) const;

  // Select or unselect a voxel.
  void SelectVoxel ( int const iVoxel[3] );
  void UnselectVoxel ( int const iVoxel[3] );
  
  // Whether a specific voxel is selected.
  bool IsVoxelSelected ( int const iVoxel[3] ) const;

  // The number of selected voxels.
  int NumSelectedVoxels () const { return mcSelectedVoxels; }

  // Get a list of selected voxels.
  std::list<Point3<int> > GetSelectedVoxelList ();

 protected:

  // Sets the mbDirtyList flag then calls the superclass function.
  virtual void ROIChanged ();

  // A volume of booleans for which voxels are selected, and the
  // bounds.
  int mBounds[3];
  bool*** mVoxels;

  // The number of selected voxels.
  int mcSelectedVoxels;

  // A rough cache of the bounds of the selected voxels.
  int mSelectedBounds[6]; // 0 = xmin, 1 = xmax, 2 = ymin, etc

  // A cache of the selected voxels. Invalidated when a voxel is
  // selected or unselected.
  bool mbDirtyList;
  std::list<Point3<int> > mlSelectedVoxels;
};

#endif
