/**
 * @file  ScubaVolumeROIIntensityChart.cpp
 * @brief Given a volume and an ROI, uses a ChartWindow to plot intenity data
 *
 * Uses a ChartWindow to plot the voxel index on the x axes and the
 * intensity on the y axis. Voxel index can be sorted by x, y, or z.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:39 $
 *    $Revision: 1.7 $
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


#include <queue>
#include "ScubaVolumeROIIntensityChart.h"

using namespace std;

bool ScubaVolumeROIIntensityChartFactory::mbAddedTclCommands = false;


ScubaVolumeROIIntensityChart::
ScubaVolumeROIIntensityChart( VolumeCollection& iVolume,
                              ScubaROIVolume& iROI, SortOrder iOrder ) :
    Listener("ScubaVolumeROIIntensityChart"),
    mVolume(iVolume),
    mROI(iROI),
    mChart(NULL),
    mSortOrder(iOrder) {

  mChart = ChartWindow::NewChartWindow();
  if ( NULL == mChart ) {
    throw runtime_error( "Couldn't create chart." );
  }

  // Listen to some stuff.
  mChart->AddListener( *this );
  mVolume.AddListener( *this );
  mROI.AddListener( *this );

  Draw();
}

ScubaVolumeROIIntensityChart::~ScubaVolumeROIIntensityChart () {

  // Stop listening.
  mVolume.RemoveListener( *this );
  mROI.RemoveListener( *this );
}

void
ScubaVolumeROIIntensityChart::DoListenToMessage ( string isMessage,
    void*  iData ) {

  if ( isMessage == "dataChanged" ) {
    Draw();
  }
  if ( isMessage == "roiChanged" ) {
    Draw();
  }
  if ( isMessage == "chartDeleted" ) {
    delete this;
  }
}

void
ScubaVolumeROIIntensityChart::Draw() {

  // Switch on our sort order and call the second Draw function
  // passing in the class of the proper comparator.
  switch ( mSortOrder ) {
  case x:
    DrawWithSort<VolumeLocationRASXComparatorGT>();
    break;
  case y:
    DrawWithSort<VolumeLocationRASYComparatorGT>();
    break;
  case z:
    DrawWithSort<VolumeLocationRASZComparatorGT>();
    break;
  }
}

template <typename T> void
ScubaVolumeROIIntensityChart::DrawWithSort() {

  // Clear it.
  mChart->ClearData();

  list<ChartWindow::PointData> lData;
  ChartWindow::PointData data;

  // Get the bounds of the ROI to iterate over.
  int bounds[3];
  mROI.GetROIBounds( bounds );

  // Check against the bounds of the volume.
  int mriIdxRange[3];
  mVolume.GetMRIIndexRange( mriIdxRange );
  if ( bounds[0] != mriIdxRange[0] ||
       bounds[1] != mriIdxRange[1] ||
       bounds[2] != mriIdxRange[2] ) {
    throw runtime_error( "Bounds of volume don't match bounds of ROI." );
  }

  // Create the queue with the comparator class that we got.
  priority_queue<VolumeLocation,vector<VolumeLocation>,T> lLocs;

  // Iterate over the ROI...
  int idx[3];
  for ( idx[2] = 0; idx[2] < bounds[2]; idx[2]++ ) {
    for ( idx[1] = 0; idx[1] < bounds[1]; idx[1]++ ) {
      for ( idx[0] = 0; idx[0] < bounds[0]; idx[0]++ ) {

        // If this voxel is selected...
        if ( mROI.IsVoxelSelected( idx ) ) {

          // Make a location from this index.
          VolumeLocation loc( mVolume.MakeVolumeLocationFromIndex( idx ) );

          // Add it our sorted queue of selected points.
          lLocs.push( loc );
        }
      }
    }
  }

  // At this point our queue is sorted, so we just go through the
  // queue picking stuff off the top.
  int cSelected = 0;
  int nPoint = 0;
  while ( !lLocs.empty() ) {

    // Get the top element and pop it.
    VolumeLocation loc = lLocs.top();
    lLocs.pop();

    // Copy the data into the point data structure. The x value is an
    // incrementing index and the y is the intensity.
    data.mX = nPoint;
    data.mY = mVolume.GetMRINearestValue( loc );

    // Our label is the RAS coordinate.
    stringstream ssLabel;
    ssLabel << "(" << loc.RAS(0) << ", " << loc.RAS(1)
    << ", " << loc.RAS(2) << ")";
    data.msLabel = ssLabel.str();

    // Add the data to the list of points we'll send to the chart.
    lData.push_back( data );

    nPoint++;
    cSelected++;
  }

  // Build the title and info string for the chart.
  stringstream sTitle, sInfo;
  sTitle << "Volume ROI Intensity Plot: Vol " << mVolume.GetID ()
  << " ROI " << mROI.GetID();

  sInfo << "Volume: \\\"" << mVolume.GetLabel()
  << "\\\" ROI: \\\"" << mROI.GetLabel()
  << "\\\" Num voxels: " << cSelected;

  // Set all the data chart and tell it to draw.
  mChart->SetTitle( sTitle.str() );
  mChart->SetInfo( sInfo.str() );
  mChart->SetXAxisLabel( "Voxels" );
  mChart->SetYAxisLabel( "Intensity" );
  mChart->SetPointData( lData );
  mChart->Draw();
}

bool
ScubaVolumeROIIntensityChart::VolumeLocationRASXComparatorGT::operator()
  ( VolumeLocation const& v1, VolumeLocation const& v2 ) const {

  // Compare the first VolumeLocation to the second and return true if
  // the x coordinate is greater in the first one. If they are equal,
  // move to the y coordinate, and then to z.
  if ( fabs(v1.RAS(0) - v2.RAS(0)) > 0.00001 ) { //if v1.RAS(0) != v2.RAS(0)
    return (v1.RAS(0) > v2.RAS(0));
  } else if ( fabs(v1.RAS(1) - v2.RAS(1)) > 0.00001 ) {
    return (v1.RAS(1) > v2.RAS(1));
  } else {
    return (v1.RAS(2) > v2.RAS(2));
  }
}

bool
ScubaVolumeROIIntensityChart::VolumeLocationRASYComparatorGT::operator()
  ( VolumeLocation const& v1, VolumeLocation const& v2 ) const {
  
  // Same, but y first, then z, then x.
  if ( fabs(v1.RAS(1) - v2.RAS(1)) > 0.00001 ) {
    return (v1.RAS(1) > v2.RAS(1));
  } else if ( fabs(v1.RAS(2) - v2.RAS(2)) > 0.00001 ) {
    return (v1.RAS(2) > v2.RAS(2));
  } else {
    return (v1.RAS(0) > v2.RAS(0));
  }
}

bool
ScubaVolumeROIIntensityChart::VolumeLocationRASZComparatorGT::operator()
  ( VolumeLocation const& v1, VolumeLocation const& v2 ) const {

  // Same, but z first, then x, then y.
  if ( fabs(v1.RAS(2) - v2.RAS(2)) > 0.00001 ) {
    return (v1.RAS(2) > v2.RAS(2));
  } else if ( fabs(v1.RAS(0) - v2.RAS(0)) > 0.00001 ) {
    return (v1.RAS(0) > v2.RAS(0));
  } else {
    return (v1.RAS(1) > v2.RAS(1));
  }
}


ScubaVolumeROIIntensityChartFactory&
ScubaVolumeROIIntensityChartFactory::GetFactory() {

  static ScubaVolumeROIIntensityChartFactory sFactory;

  if ( !mbAddedTclCommands ) {

    mbAddedTclCommands = true;
    TclCommandManager& commandMgr = TclCommandManager::GetManager();
    commandMgr.AddCommand( sFactory, "MakeNewVolumeROIIntensityChart", 3,
                           "volumeID roiID sortOrder", "Makes a new chart "
                           "window plotting the intensities of voxels in a "
                           "volume ROI. sortOrder should be x, y, or z." );
  }

  return sFactory;
}

TclCommandListener::TclCommandResult
ScubaVolumeROIIntensityChartFactory::DoListenToTclCommand( char* isCommand,
    int, char** iasArgv ) {

  // MakeVolumeROIIntensityChart <volumeID> <roiID>
  if ( 0 == strcmp( isCommand, "MakeNewVolumeROIIntensityChart" ) ) {

    // Get the collection ID.
    int collectionID;
    try {
      collectionID = TclCommandManager::ConvertArgumentToInt( iasArgv[1] );
    } catch ( runtime_error& e ) {
      sResult = string("bad collectionID: ") + e.what();
      return error;
    }

    // Try to find the object.
    DataCollection* col;
    try {
      col = &DataCollection::FindByID( collectionID );
    } catch ( runtime_error& e ) {
      sResult = string("Couldn't find collection: ") + e.what();
      return error;
    }

    // Make sure it's a volume.
    if ( col->GetTypeDescription() != "Volume" ) {
      throw runtime_error( "Collection wasn't a volume." );
    }

    // We can do this because we checked the type.
    VolumeCollection* vol = (VolumeCollection*) col;


    // Now try to get the ROI.
    int roiID;
    try {
      roiID = TclCommandManager::ConvertArgumentToInt( iasArgv[2] );
    } catch ( runtime_error& e ) {
      sResult = string("bad roiID: ") + e.what();
      return error;
    }

    // If it's not in this volume, return an error.
    if ( !vol->IsROIInThisCollection( roiID ) ) {
      sResult = "That ROI is not in that volume.";
      return error;
    }

    // Get the object.
    ScubaROI* roi;
    try {
      roi = &ScubaROI::FindByID( roiID );
    } catch ( runtime_error& e ) {
      sResult = string("Couldn't find roi: ") + e.what();
      return error;
    }

    // This conversion is safe because if it's in the volume, it's a
    // ScubaROIVolume.
    ScubaROIVolume* roiVolume = (ScubaROIVolume*)roi;

    // Try to parse the sortOrder.
    ScubaVolumeROIIntensityChart::SortOrder sortOrder;
    if ( iasArgv[3][0] == 'x' || iasArgv[3][0] == 'X' ) {
      sortOrder = ScubaVolumeROIIntensityChart::x;
    } else if ( iasArgv[3][0] == 'y' || iasArgv[3][0] == 'Y' ) {
      sortOrder = ScubaVolumeROIIntensityChart::y;
    } else if ( iasArgv[3][0] == 'z' || iasArgv[3][0] == 'Z' ) {
      sortOrder = ScubaVolumeROIIntensityChart::z;
    } else {
      sResult = string("bad sortOrder: ") + iasArgv[3][0];
      return error;
    }

    // Create the chart.
    new ScubaVolumeROIIntensityChart( *vol, *roiVolume, sortOrder );
  }

  return ok;
}
