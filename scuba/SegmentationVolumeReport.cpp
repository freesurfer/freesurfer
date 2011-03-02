/**
 * @file  SegmentationVolumeReport.cpp
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
 *    $Revision: 1.9 $
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


#include <fstream>
#include "SegmentationVolumeReport.h"

using namespace std;

SegmentationVolumeReport&
SegmentationVolumeReport::GetReport () {

  static SegmentationVolumeReport* sReport = NULL;
  if ( NULL == sReport ) {

    sReport = new SegmentationVolumeReport();

    TclCommandManager& commandMgr = TclCommandManager::GetManager();
    commandMgr.AddCommand( *sReport, "ClearSegVolReport", 0, "",
                           "Clears all info in the segmentation "
                           "volume report." );
    commandMgr.AddCommand( *sReport, "SetSegVolReportSegmentation", 1, "segID",
                           "Set the segmentation volume in the segmentation "
                           "volume report." );
    commandMgr.AddCommand( *sReport, "AddSegVolReportIntensityVolume", 1,
                           "volID", "Add an intensity volume to the "
                           "segmentation volume report." );
    commandMgr.AddCommand( *sReport, "SetROIForSegVolReport", 2, "volID roiID",
                           "Set the volume and ROI to use in the "
                           "segmentation volume report." );
    commandMgr.AddCommand( *sReport, "DontUseROIInSegVolReport", 0, "",
                           "Don't use an ROI in the segmentation volume "
                           "report." );
    commandMgr.AddCommand( *sReport, "SetSegVolReportLUT", 1, "lutID",
                           "Set the LUT in the segmentation volume report." );
    commandMgr.AddCommand( *sReport, "AddSegmentationToSegVolReport", 1,
                           "structure", "Add a segmentation to the "
                           "segmentation volume report." );
    commandMgr.AddCommand( *sReport, "MakeSegVolReport", 1, "fnReport",
                           "Make the segmentation volume report." );
    commandMgr.AddCommand( *sReport, "MakeSegVolIntensityReport", 1,
                           "fnReport", "Make the intensity report from the "
                           "segmentation volume report." );
  }

  return *sReport;
}

SegmentationVolumeReport::SegmentationVolumeReport() {

  mSegVol = NULL;
  mbUseROI = false;
  mROIVol = NULL;
  mROI = NULL;
  mLUT = NULL;
  mbReportDirty = true;
}

void
SegmentationVolumeReport::Clear () {

  mSegVol = NULL;
  mlIntVols.clear();
  mbUseROI = false;
  mROIVol = NULL;
  mROI = NULL;
  mLUT = NULL;
  mlStructures.clear();
  mbReportDirty = true;
  mStructureToVolumeMap.clear();
  mVolumeToIntensityAverageMap.clear();
  mStructureToVolumeVoxelListMap.clear();
}

void
SegmentationVolumeReport::SetSegmentation ( VolumeCollection& iSeg ) {

  mSegVol = &iSeg;
  mbReportDirty = true;
}

void
SegmentationVolumeReport::AddIntensityVolume ( VolumeCollection& iVol ) {

  mlIntVols.push_back( &iVol );
  mbReportDirty = true;
}

void
SegmentationVolumeReport::DontUseROI () {

  mbUseROI = false;
  mbReportDirty = true;
}

void
SegmentationVolumeReport::UseVolumeForROI ( VolumeCollection& iVol ) {

  mbUseROI = true;
  mROIVol = &iVol;
  mbReportDirty = true;
}

void
SegmentationVolumeReport::UseROI ( ScubaROIVolume& iROI ) {

  mbUseROI = true;
  mROI = &iROI;
  mbReportDirty = true;
}

void
SegmentationVolumeReport::SetColorLUT ( ScubaColorLUT& iLUT ) {

  mLUT = &iLUT;
  mbReportDirty = true;
}

void
SegmentationVolumeReport::AddSegmentationStructure ( int inStructure ) {

  mlStructures.push_back( inStructure );
  mbReportDirty = true;
}

void
SegmentationVolumeReport::MakeVolumeReport ( string ifnReport ) {

  try {
    // Check file name first.
    ofstream fReport( ifnReport.c_str(), ios::out );

    // Generate the report data.
    MakeVolumeReport();

    // Write out the file.

    // SUBJ Struct Struct-intVol-intensity...
    fReport << "SUBJ";
    list<int>::iterator tStructure;
    for ( tStructure = mlStructures.begin(); tStructure != mlStructures.end();
          ++tStructure ) {
      int nStructure = *tStructure;
      string sStructure = mLUT->GetLabelAtIndex( nStructure );
      fReport << "\t" << sStructure;

      list<VolumeCollection*>::iterator tIntVol;
      for ( tIntVol = mlIntVols.begin(); tIntVol != mlIntVols.end();
            ++tIntVol ) {
        VolumeCollection* intVol = *tIntVol;
        string sIntVol = intVol->GetLabel();
        fReport << "\t" << sStructure << "-" << sIntVol
        << "-intensity";
      }
    }
    fReport << endl;


    // SubjName volume intensityAverage...
    fReport << mSegVol->GetLabel();
    for ( tStructure = mlStructures.begin(); tStructure != mlStructures.end();
          ++tStructure ) {
      int nStructure = *tStructure;
      fReport << "\t" << mStructureToVolumeMap[nStructure];

      list<VolumeCollection*>::iterator tIntVol;
      for ( tIntVol = mlIntVols.begin(); tIntVol != mlIntVols.end();
            ++tIntVol ) {
        VolumeCollection* intVol = *tIntVol;
        fReport << "\t" << mVolumeToIntensityAverageMap[intVol][nStructure];
      }
    }
    fReport << endl;

  } catch ( exception& e ) {
    throw runtime_error( "Error writing " + ifnReport + ": " + e.what() );
  } catch ( ... ) {
    throw runtime_error( "Error writing " + ifnReport );
  }
}


void
SegmentationVolumeReport::MakeIntensityReport ( std::string ifnReport ) {

  try {
    // Check file name first.
    ofstream fReport( ifnReport.c_str(), ios::out );

    // Generate the report data.
    MakeVolumeReport();

    // Struct-intVol...
    bool bFirst = true;
    list<int>::iterator tStructure;
    list<VolumeCollection*>::iterator tIntVol;
    for ( tStructure = mlStructures.begin(); tStructure != mlStructures.end();
          ++tStructure ) {
      int nStructure = *tStructure;
      string sStructure = mLUT->GetLabelAtIndex( nStructure );
      for ( tIntVol = mlIntVols.begin(); tIntVol != mlIntVols.end();
            ++tIntVol ) {
        VolumeCollection* intVol = *tIntVol;
        string sIntVol = intVol->GetLabel();

        if ( !bFirst ) {
          fReport << "\t";
        }
        fReport << sIntVol << "-" << sStructure;
        bFirst = false;
      }
    }
    fReport << endl;

    // Because our table output has lists of voxels in columns, but we
    // need to write it in rows, we need to put our voxel data in a
    // random access structure. So we'll use a volume with structures
    // indices (0-based) on the x axis, vol indices on the y, and
    // voxel lists in the z direction.

    // Find the max length of the list of voxels in each structure.
    int cVoxels = 0;
    for ( tIntVol = mlIntVols.begin(); tIntVol != mlIntVols.end();
          ++tIntVol ) {
      VolumeCollection* intVol = *tIntVol;
      for ( tStructure = mlStructures.begin(); tStructure != mlStructures.end();
            ++tStructure ) {
        int nStructure = *tStructure;

        if ((int)mStructureToVolumeVoxelListMap[nStructure][intVol].size() > cVoxels)
          cVoxels = mStructureToVolumeVoxelListMap[nStructure][intVol].size();
      }
    }

    int cStructures = mlStructures.size();
    int cVolumes = mlIntVols.size();

    // Allocate a volume of cStructures x cIntVols x cVoxels to hold our
    // VolumeLocation pointers.
    Volume3<VolumeLocation*> volLocs( cStructures, cVolumes, cVoxels+1, NULL );

    // Go through everything and fill it out.
    int nVolume = 0;
    int nZeroBasedStructure = 0;
    for ( tIntVol = mlIntVols.begin(); tIntVol != mlIntVols.end();
          ++tIntVol ) {

      nZeroBasedStructure = 0;
      VolumeCollection* intVol = *tIntVol;
      for ( tStructure = mlStructures.begin(); tStructure != mlStructures.end();
            ++tStructure ) {
        int nStructure = *tStructure;

        list<VolumeLocation>::iterator tLocation;
        int nLocation = 0;
        for ( tLocation =
                mStructureToVolumeVoxelListMap[nStructure][intVol].begin();
              tLocation !=
              mStructureToVolumeVoxelListMap[nStructure][intVol].end();
              ++tLocation ) {
          volLocs.Set( nZeroBasedStructure, nVolume, nLocation, &(*tLocation));
          nLocation++;
        }
        nZeroBasedStructure++;
      }
      nVolume++;
    }

    // Now fill them out by row.
    for ( int nVoxel = 0; nVoxel < cVoxels; nVoxel++ ) {
      for ( nZeroBasedStructure = 0; nZeroBasedStructure < cStructures; nZeroBasedStructure++ ) {
        for ( int nVolume = 0; nVolume < cVolumes; nVolume++ ) {

          VolumeLocation* loc =
            volLocs.Get( nZeroBasedStructure, nVolume, nVoxel );
          if ( NULL != loc ) {
            fReport << Point3<int>(loc->Index());
          }
          fReport << "\t";
        }
      }
      fReport << endl;
    }
  } catch ( exception& e ) {
    throw runtime_error( "Error writing " + ifnReport + ": " + e.what() );
  } catch ( ... ) {
    throw runtime_error( "Error writing " + ifnReport );
  }
}

TclCommandManager::TclCommandResult
SegmentationVolumeReport::DoListenToTclCommand ( char* isCommand, int,
    char** iasArgv ) {

  // ClearSegVolReport
  if ( 0 == strcmp( isCommand, "ClearSegVolReport" ) ) {

    Clear();
    return ok;
  }

  // SetSegVolReportSegmentation <segID>
  if ( 0 == strcmp( isCommand, "SetSegVolReportSegmentation" ) ) {

    int segID;
    try {
      segID = TclCommandManager::ConvertArgumentToInt( iasArgv[1] );
      DataCollection::FindByID( segID );
    } catch ( runtime_error& e ) {
      sResult = string("bad segID: ") + e.what();
      return error;
    }

    VolumeCollection& seg =
      (VolumeCollection&) DataCollection::FindByID( segID );
    SetSegmentation( seg );

    return ok;
  }

  // AddSegVolReportIntensityVolume <volID>
  if ( 0 == strcmp( isCommand, "AddSegVolReportIntensityVolume" ) ) {

    int volID;
    try {
      volID = TclCommandManager::ConvertArgumentToInt( iasArgv[1] );
      DataCollection::FindByID( volID );
    } catch ( runtime_error& e ) {
      sResult = string("bad volID: ") + e.what();
      return error;
    }

    VolumeCollection& intVol =
      (VolumeCollection&) DataCollection::FindByID( volID );
    AddIntensityVolume( intVol );

    return ok;
  }

  // SetROIForSegVolReport <volID> <roiID>
  if ( 0 == strcmp( isCommand, "SetROIForSegVolReport" ) ) {

    int volID;
    try {
      volID = TclCommandManager::ConvertArgumentToInt( iasArgv[1] );
      DataCollection::FindByID( volID );
    } catch ( runtime_error& e ) {
      sResult = string("bad volID: ") + e.what();
      return error;
    }

    int roiID;
    try {
      roiID = TclCommandManager::ConvertArgumentToInt( iasArgv[2] );
      ScubaROI::FindByID( volID );
    } catch ( runtime_error& e ) {
      sResult = string("bad roiID: ") + e.what();
      return error;
    }

    VolumeCollection& vol =
      (VolumeCollection&) DataCollection::FindByID( volID );
    ScubaROIVolume& roi = (ScubaROIVolume&) ScubaROI::FindByID( roiID );

    UseVolumeForROI( vol );
    UseROI( roi );

    return ok;
  }

  // DontUseROIInSegVolReportk
  if ( 0 == strcmp( isCommand, "DontUseROIInSegVolReportk" ) ) {

    DontUseROI();

    return ok;
  }

  // SetSegVolReportLUT <lutID>
  if ( 0 == strcmp( isCommand, "SetSegVolReportLUT" ) ) {

    int lutID;
    try {
      lutID = TclCommandManager::ConvertArgumentToInt( iasArgv[1] );
      ScubaColorLUT::FindByID( lutID );
    } catch ( runtime_error& e ) {
      sResult = string("bad lutID: ") + e.what();
      return error;
    }

    ScubaColorLUT& lut = (ScubaColorLUT&) ScubaColorLUT::FindByID( lutID );
    SetColorLUT( lut );

    return ok;
  }

  // AddSegmentationToSegVolReport <structure>
  if ( 0 == strcmp( isCommand, "AddSegmentationToSegVolReport" ) ) {

    int nStructure;
    try {
      nStructure = TclCommandManager::ConvertArgumentToInt( iasArgv[1] );
    } catch ( runtime_error& e ) {
      sResult = string("bad structure: ") + e.what();
      return error;
    }

    AddSegmentationStructure( nStructure );

    return ok;
  }

  // MakeSegVolReport <fnReport>
  if ( 0 == strcmp( isCommand, "MakeSegVolReport" ) ) {

    string fnReport( iasArgv[1] );
    MakeVolumeReport( fnReport );

    return ok;
  }

  // MakeSegVolIntensityReport <fnReport>
  if ( 0 == strcmp( isCommand, "MakeSegVolIntensityReport" ) ) {

    string fnReport( iasArgv[1] );
    MakeIntensityReport( fnReport );

    return ok;
  }

  return ok;
}

void
SegmentationVolumeReport::MakeVolumeReport () {

  if ( !mbReportDirty )
    return;

  // For each structure...
  list<int>::iterator tStructure;
  for ( tStructure = mlStructures.begin(); tStructure != mlStructures.end();
        ++tStructure ) {

    int nStructure = *tStructure;

    // Get a list of voxels with this structure value from our
    // segmentation volume.
    list<VolumeLocation> lLocations;
    list<VolumeLocation>::iterator tLoc;
    mSegVol->GetVoxelsWithValue( nStructure, lLocations );

    // If we're using an ROI, go through all the voxels and if they
    // are not in the ROI, remove them from our list.
    if ( mbUseROI && NULL != mROIVol && NULL != mROI ) {
      for ( tLoc = lLocations.begin(); tLoc != lLocations.end(); ++tLoc ) {

        // We have a location from the seg vol, which may not be the
        // same volume as the ROI volume. So make a location for the
        // ROI volume from the RAS of the seg volume location.
        VolumeLocation& loc = *tLoc;
        VolumeLocation roiLoc( mROIVol->MakeVolumeLocationFromRAS( loc.RAS() ) );
        if ( !mROIVol->IsSelected( roiLoc ) ) {
          lLocations.erase( tLoc );
        }
      }
    }

    // Get the RAS volume of this many voxels from our segmentation
    // volume.
    mStructureToVolumeMap[nStructure] =
      mSegVol->GetRASVolumeOfNVoxels( lLocations.size() );

    // For each intensity volume...
    list<VolumeCollection*>::iterator tVol;
    for ( tVol = mlIntVols.begin(); tVol != mlIntVols.end(); ++tVol ) {

      VolumeCollection* vol = (*tVol);

      // Have to make a list of locations for this volume from the
      // other one.
      mStructureToVolumeVoxelListMap[nStructure][vol].clear();
      for ( tLoc = lLocations.begin(); tLoc != lLocations.end(); ++tLoc ) {
        VolumeLocation& loc = *tLoc;
        VolumeLocation intLoc( vol->MakeVolumeLocationFromRAS( loc.RAS() ) );
        mStructureToVolumeVoxelListMap[nStructure][vol].push_back( intLoc );
      }

      // Get the average intensity for this list of voxels.
      if ( mStructureToVolumeVoxelListMap[nStructure][vol].size() > 0 ) {
        mVolumeToIntensityAverageMap[vol][nStructure] =
          vol->GetAverageValue( mStructureToVolumeVoxelListMap[nStructure][vol] );
      } else {
        mVolumeToIntensityAverageMap[vol][nStructure] = 0;
      }
    }
  }

  mbReportDirty = false;
}
