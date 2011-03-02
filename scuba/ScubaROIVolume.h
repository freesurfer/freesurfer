/**
 * @file  ScubaROIVolume.h
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
 *    $Revision: 1.6 $
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


#ifndef ScubaROIVolume_h
#define ScubaROIVolume_h

#include "Point3.h"
#include "ScubaROI.h"
#include "Volume3.h"

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
  int NumSelectedVoxels () const {
    return mcSelectedVoxels;
  }

  // Get a list of selected voxels.
  std::list<Point3<int> > GetSelectedVoxelList () const;

protected:

  // Sets the mbDirtyList flag then calls the superclass function.
  virtual void ROIChanged ();

  // A volume of booleans for which voxels are selected.
  Volume3<bool>* mVoxels;

  // The number of selected voxels.
  int mcSelectedVoxels;

  // A rough cache of the bounds of the selected voxels.
  int mSelectedBounds[6]; // 0 = xmin, 1 = xmax, 2 = ymin, etc

  // A cache of the selected voxels. Invalidated when a voxel is
  // selected or unselected.
  mutable bool mbDirtyList;
  mutable std::list<Point3<int> > mlSelectedVoxels;
};

#endif
