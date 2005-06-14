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

  void MakeVolumeReport();

  // Generated after report is ready.
  std::map<int,float> mStructureToVolumeMap;
  
  typedef std::map<int,float> tStructureToIntensityAverageMap;
  std::map<VolumeCollection*,tStructureToIntensityAverageMap> 
    mVolumeToIntensityAverageMap;
};

#endif
