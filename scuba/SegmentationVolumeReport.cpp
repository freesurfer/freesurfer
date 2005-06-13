#include "SegmentationVolumeReport.h"

SegmentationVolumeReport& 
SegmentationVolumeReport::GetReport () {

  static SegmentationVolumeReport* sReport = NULL;
  if( NULL == sReport ) {

    sReport = new SegmentationVolumeReport();

    TclCommandManager& commandMgr = TclCommandManager::GetManager();
    commandMgr.AddCommand( *sReport, "DoStuff", 1, "arg",
			   "Calc stuff." );
  }

  return *sReport;
}

SegmentationVolumeReport::SegmentationVolumeReport() {

  mSegVol = NULL;
  mbUseROI = false;
  mROIVol = NULL;
  mROI = NULL;
  mLUT = NULL;
}

void  
SegmentationVolumeReport::SetSegmentation ( VolumeCollection& iSeg ) {

  mSegVol = &iSeg;
}

void  
SegmentationVolumeReport::AddIntensityVolume ( VolumeCollection& iVol ) {

  mlIntVols.push_back( &iVol );
}

void  
SegmentationVolumeReport::DontUseROI () {

  mbUseROI = false;
}

void  
SegmentationVolumeReport::UseVolumeForROI ( VolumeCollection& iVol ) {

  mbUseROI = true;
  mROIVol = &iVol;
}

void  
SegmentationVolumeReport::UseROI ( ScubaROIVolume& iROI ) {

  mbUseROI = true;
  mROI = &iROI;
}

void  
SegmentationVolumeReport::SetColorLUT ( ScubaColorLUT& iLUT ) {

  mLUT = &iLUT;
}

void  
SegmentationVolumeReport::AddSegmentationStructure ( int inStructure ) {

  mlStructures.push_back( inStructure );
}

void  
SegmentationVolumeReport::MakeVolumeReport ( std::string ifnReport ) {

  MakeVolumeReport();
}

void  
SegmentationVolumeReport::MakeIntensityReport ( std::string ifnReport ) {

}

TclCommandManager::TclCommandResult
SegmentationVolumeReport::DoListenToTclCommand ( char* isCommand, int, 
						 char** iasArgv ) {
}

void
SegmentationVolumeReport::MakeVolumeReport () {

  // For each structure...
  list<int>::iterator tStructure;
  for( tStructure = mlStructures.first(); tStructure != mlStructures.end();
       ++tStructure ) {

    int nStructure = *tStructure;

    // Get a list of voxels with this structure value from our
    // segmentation volume.
    list<VolumeLocation> lVoxels;
    list<VolumeLocation>::iterator tVoxel;
    lVoxels = mSegVol.GetVoxelsInStructure( nStructure );

    // If we're using an ROI, go through all the voxels and if they
    // are not in the ROI, remove them from our list.
    if( mbUseROI && NULL != mROIVol && NULL != mROI ) {
      for( tVoxel = lVoxels.first(); tVoxel != lVoxels.end(); ++tVoxel ) {

	// We have a location from the seg vol, which may not be the
	// same volume as the ROI volume. So make a location for the
	// ROI volume from the RAS of the seg volume location.
	VolumeLocation& loc = *tLoc;
	VolumeLocation roiLoc = mROIVol->MakeLocationFromRAS( loc.RAS() );
	if( !mROIVol->IsSelected( roiLoc ) ) {
	  lVoxels.erase( tVoxel );
	}
      }
    }
    
    // Get the RAS volume of this many voxels from our segmentation
    // volume.
    mStructureToVolumeMap[nStructure] = 
      mSegVol.GetRASVolumeOfNVoxels( lVoxels.size() );

    // For each intensity volume...
    list<VolumeCollection*>::iterator tVol;
    for( tVol = mlIntVols.begin(); tVol != mlIntVols.end(); ++tVol ) {

      VolumeCollection& vol = *(*tVol);

      // Get the average intensity for this list of voxels.
      mVolumeToIntensityAverageMap[vol][nStructure] = 
	vol.GetAverageIntensity( lVoxels );
    }
  }
}
