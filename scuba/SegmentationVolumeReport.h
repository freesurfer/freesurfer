/**
 * @file  SegmentationVolumeReport.h
 * @brief Calculates information based on segmentation groups
 *
 * Given a segmentation volume and one or more intensity volumes, will
 * make reports about the average intensity per segmentation value in
 * a variety of formats.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:39 $
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


#ifndef SegmentationVolumeReport_H
#define SegmentationVolumeReport_H

#include <list>
#include <map>
#include "TclCommandManager.h"
#include "VolumeCollection.h"
#include "ScubaROIVolume.h"
#include "ScubaColorLUT.h"

class SegmentationVolumeReport : public TclCommandListener {

  friend class SegmentationVolumeReportTester;

public:

  SegmentationVolumeReport();

  static SegmentationVolumeReport& GetReport ();

  // Clear all info about the report.
  void Clear ();

  // Declare the volume to use as the main segmentation. This will be
  // the volume whose segmentation structures are counted.
  void SetSegmentation ( VolumeCollection& iSeg );

  // Add an additional volume. Average intensities will be calculated
  // for all structures in the report. If the additional intensity
  // report is requested, all intensities in the structures will be
  // individually listed for each volume.
  void AddIntensityVolume ( VolumeCollection& iVol );

  // Limit considered voxels to an ROI. Only voxels that are in an
  // ROI will be counted for each strucutre.
  void DontUseROI ();
  void UseVolumeForROI ( VolumeCollection& iVol );
  void UseROI ( ScubaROIVolume& iROI );

  // Specify the LUT to use, and which structures. The LUT basically
  // just defines the labels in the report.
  void SetColorLUT ( ScubaColorLUT& iLUT );
  void AddSegmentationStructure ( int inStructure );

  // Generate the tab-delimited reports.
  void MakeVolumeReport ( std::string ifnReport );

  void MakeIntensityReport ( std::string ifnReport );

  virtual TclCommandResult DoListenToTclCommand ( char* isCommand, int,
      char** iasArgv );

protected:

  VolumeCollection* mSegVol;
  std::list<VolumeCollection*> mlIntVols;

  bool mbUseROI;
  VolumeCollection* mROIVol;
  ScubaROIVolume* mROI;

  ScubaColorLUT* mLUT;
  std::list<int> mlStructures;

  bool mbReportDirty;
  void MakeVolumeReport();

  // Generated after report is ready.
  std::map<int,float> mStructureToVolumeMap;

  typedef std::map<int,float> tStructureToIntensityAverageMap;
  std::map<VolumeCollection*,tStructureToIntensityAverageMap>
  mVolumeToIntensityAverageMap;

  // List of voxels in each structure.
  typedef std::map<VolumeCollection*,std::list<VolumeLocation> >
  tVolumeToVoxelListMap;
  std::map<int,tVolumeToVoxelListMap> mStructureToVolumeVoxelListMap;
};

#endif
