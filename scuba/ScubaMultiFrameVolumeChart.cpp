/**
 * @file  ScubaMultiFrameVolumeChart.cpp
 * @brief Creates charts for values in multiframe volumes
 *
 * Uses the ChartWindow class to create charts based on volumes that
 * will get the current cursor and draw the volume's values at that
 * cursor over multiple frames. If there is time metadata available,
 * it will use that to format the graph into conditions and times.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:38 $
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
#include "ScubaMultiFrameVolumeChart.h"
#include "ScubaView.h"
#include "ScubaLayer2DMRI.h"

using namespace std;

bool ScubaMultiFrameVolumeChartFactory::mbAddedTclCommands = false;

ScubaMultiFrameVolumeChart::
ScubaMultiFrameVolumeChart ( VolumeCollection& iVolume ) :
    Listener( "ScubaMultiFrameVolumeChart" ),
    mVolume(iVolume),
    mROI(NULL),
    mChart(NULL) {

  // Create the chart and listen to it.
  mChart = ChartWindow::NewChartWindow();
  if ( NULL == mChart ) {
    throw runtime_error( "Couldn't create chart." );
  }

  // Listen to some stuff.
  mChart->AddListener( *this );
  mVolume.AddListener( *this );

  // Listen to the view boradcaster for cursor changes.
  ScubaViewBroadcaster::GetBroadcaster().AddListener( *this );

  // We want to listen to any layers that use this volume, so we can
  // update the current timepoint and condition when the layer changes
  // it. Find all layers, and listen to them.
  list<int> lIDs;
  Layer::GetIDList( lIDs );
  list<int>::iterator tID;
  for ( tID = lIDs.begin(); tID != lIDs.end(); ++tID ) {
    int id = *tID;
    Layer& layer = Layer::FindByID( id );
    DataCollection* col = layer.GetMainDataCollection();
    if ( NULL != col &&
         col->GetID() == mVolume.GetID() ) {
      layer.AddListener( *this );
    }
  }

  Draw();
}

ScubaMultiFrameVolumeChart::
ScubaMultiFrameVolumeChart ( VolumeCollection& iVolume,
                             ScubaROIVolume&   iROI ) :
    Listener( "ScubaMultiFrameVolumeChart" ),
    mVolume(iVolume),
    mROI(&iROI),
    mChart(NULL) {

  // Create the chart and listen to it.
  mChart = ChartWindow::NewChartWindow();
  if ( NULL == mChart ) {
    throw runtime_error( "Couldn't create chart." );
  }

  // Listen to some stuff.
  mChart->AddListener( *this );
  mVolume.AddListener( *this );
  mROI->AddListener( *this );

  // Listen to the view boradcaster for cursor changes.
  ScubaViewBroadcaster::GetBroadcaster().AddListener( *this );

  // We want to listen to any layers that use this volume, so we can
  // update the current timepoint and condition when the layer changes
  // it. Find all layers, and listen to them.
  list<int> lIDs;
  Layer::GetIDList( lIDs );
  list<int>::iterator tID;
  for ( tID = lIDs.begin(); tID != lIDs.end(); ++tID ) {
    int id = *tID;
    Layer& layer = Layer::FindByID( id );
    DataCollection* col = layer.GetMainDataCollection();
    if ( NULL != col &&
         col->GetID() == mVolume.GetID() ) {
      layer.AddListener( *this );
    }
  }

  Draw();
}

ScubaMultiFrameVolumeChart::~ScubaMultiFrameVolumeChart () {

  // Stop listening.
  //  mChart->RemoveListener( *this );
  mVolume.RemoveListener( *this );

  // Stop listening to all those listeners.
  list<int> lIDs;
  Layer::GetIDList( lIDs );
  list<int>::iterator tID;
  for ( tID = lIDs.begin(); tID != lIDs.end(); ++tID ) {
    int id = *tID;
    Layer& layer = Layer::FindByID( id );
    DataCollection* col = layer.GetMainDataCollection();
    if ( NULL != col &&
         col->GetID() == mVolume.GetID() ) {
      layer.RemoveListener( *this );
    }
  }
}

void
ScubaMultiFrameVolumeChart::DoListenToMessage ( string isMessage,
    void*  iData ) {

  if ( isMessage == "dataChanged" ) {
    Draw();
  }
  if ( isMessage == "roiChanged" ) {
    Draw();
  }
  if ( isMessage == "layerChanged" ) {
    Draw();
  }
  if ( isMessage == "cursorChanged" ) {
    Draw();
  }
  if ( isMessage == "chartDeleted" ) {
    delete this;
  }
}

void
ScubaMultiFrameVolumeChart::Draw() {

  // Clear it.
  mChart->ClearData();

  // We're going to draw stuff differently according to whether this
  // volume is just a series of frame, or whether we have meta info
  // about how to divide those frames up into time points and
  // conditions. A tp/cond volume may have more than one group, and
  // will need a legend. We can also call different frames time
  // instead of just frame numbers.

  // Set our x axis label according to whether this is a volume with
  // timepoint data.
  if ( mVolume.InterpretFramesAsTimePoints() ) {
    mChart->SetXAxisLabel( "Time" );
  } else {
    mChart->SetXAxisLabel( "Frame" );
  }

  // Y axis label is just "value".
  mChart->SetYAxisLabel( "Value" );

  // If we have time info, we maybe have multiple conditions, which
  // we'll draw as difference groups. Otherwise, we just have one
  // group.
  int cGroup = 1;
  if ( mVolume.InterpretFramesAsTimePoints() ) {
    cGroup = mVolume.GetNumberOfConditions();
  } else {
    cGroup = 1;
  }

  // If we have more than one group, show the legend.
  if ( cGroup > 1 ) {
    mChart->SetShowLegend( true );
  }

  // For each group...
  for ( int nGroup = 0; nGroup < cGroup; nGroup++ ) {

    // If we have time info, we'll add a label for this condition,
    // as we'll need it in the legend.
    stringstream ssGroupLabel;
    if ( mVolume.InterpretFramesAsTimePoints() ) {
      ssGroupLabel << "Condition " << nGroup;
    }

    // Set the group info for this group.
    mChart->SetGroupLabel( nGroup, ssGroupLabel.str() );
    mChart->SetGroupConnected( nGroup, 1 );

    // Find out how many elements on the x axis we have. For a
    // volume with time info, this is the number of time points,
    // otherwise it's just the number of frames.
    int cX;
    if ( mVolume.InterpretFramesAsTimePoints() ) {
      cX = mVolume.GetNumberOfTimePoints();
    } else {
      cX = mVolume.GetNumberOfFrames();
    }

    // For each frame or time point...
    for ( int nX = 0; nX < cX; nX++ ) {

      // We need to find the proper frame. For a volume with time
      // info, we use the condition and time point to calculate a
      // frame, otherwise it's just the x index here.
      int nFrame;
      if ( mVolume.InterpretFramesAsTimePoints() ) {
        nFrame = mVolume.ConvertConditionAndTimePointToFrame( nGroup, nX );
      } else {
        nFrame = nX;
      }

      // We're going to find the value and the label for this
      // point. If this chart doesn't have an ROI, it's just the value
      // at the cursor. If it does have an ROI, we need to find the
      // value at each ROI location and average it.
      float value = 0;
      stringstream ssPointLabel;

      // No ROI...
      if ( NULL == mROI ) {

        // Get the cursor.
        float cursorRAS[3];
        ScubaView::GetCursor( cursorRAS );

        // Make a location in the volume.
        VolumeLocation 
	  loc( mVolume.MakeVolumeLocationFromRAS( cursorRAS, nFrame ) );

        // Check against the bounds of the volume.
        if ( mVolume.IsInBounds( loc ) ) {

          // Just value at the cursor.
          value = mVolume.GetMRINearestValue( loc );
          ssPointLabel << "Frame " << nX;
        }

      } else {

        // As long as this is a non-empty ROI...
        if ( mROI->NumSelectedVoxels() > 0 ) {

          // Iterate over the ROI and sum the values.
          list<Point3<int> > lSelected = mROI->GetSelectedVoxelList();
          list<Point3<int> >::iterator tSelected;
          for ( tSelected = lSelected.begin();
                tSelected != lSelected.end(); ++tSelected ) {
            Point3<int> voxel = *tSelected;
            VolumeLocation 
	      loc( mVolume.MakeVolumeLocationFromIndex(voxel.xyz(),
						       nFrame) );
            value += mVolume.GetMRINearestValue( loc );
          }

          // Average the value.
          value /= mROI->NumSelectedVoxels();
        }
      }

      // Set data in this point and add it to the list.
      ChartWindow::PointData data;
      data.mX = nX;
      data.mY = value;
      data.msLabel = ssPointLabel.str();
      mChart->AddPointData( nGroup, data );
    }
  }

  // Go through all the layers. If they are showing this volume, get
  // the current time point and condition from the layer. Draw a time
  // point marker.
  list<int> lIDs;
  Layer::GetIDList( lIDs );
  list<int>::iterator tID;
  for ( tID = lIDs.begin(); tID != lIDs.end(); ++tID ) {
    int id = *tID;
    Layer& layer = Layer::FindByID( id );
    DataCollection* col = layer.GetMainDataCollection();
    if ( NULL != col &&
         col->GetID() == mVolume.GetID() ) {

      ScubaLayer2DMRI& mriLayer = (ScubaLayer2DMRI&) layer;

      ChartWindow::MarkerData marker;

      if ( mVolume.InterpretFramesAsTimePoints() ) {

        marker.mValue = mriLayer.GetCurrentTimePoint();

        // Get the color of the group matching the current condition
        // and use that for the marker color.
        int condition = mriLayer.GetCurrentCondition();
        int color[3];
        mChart->GetGroupColor( condition, color );
        marker.mColorRGBi[0] = color[0];
        marker.mColorRGBi[1] = color[1];
        marker.mColorRGBi[2] = color[2];

      } else {

        marker.mValue = mriLayer.GetCurrentFrame();
      }

      marker.msLabel = col->GetLabel();

      mChart->AddXAxisMarker( marker );
    }
  }

  // Build the title and info string for the chart.
  stringstream sTitle, sInfo;

  if ( NULL == mROI ) {

    // Get the cursor and put that in the strings.
    float cursorRAS[3];
    ScubaView::GetCursor( cursorRAS );

    sTitle << "Multiframe Volume Chart: Vol " << mVolume.GetID ()
    << " RAS " << Point3<float>(cursorRAS);
    sInfo << "Volume: \\\"" << mVolume.GetLabel()
    << "\\\" RAS: \\\"" << Point3<float>(cursorRAS);

  } else {

    // Get the ROI ID and number of voxels in the ROI.
    sTitle << "Multiframe Volume Chart: Vol " << mVolume.GetID ()
    << " ROI " << mROI->GetID();
    sInfo << "Volume: \\\"" << mVolume.GetLabel()
    << "\\\" ROI: \\\"" << mROI->GetLabel()
    << "\\\" (" << mROI->NumSelectedVoxels() << " voxels)";
  }

  mChart->SetTitle( sTitle.str() );
  mChart->SetInfo( sInfo.str() );

  // Tell the chart to draw.
  mChart->Draw();
}

ScubaMultiFrameVolumeChartFactory&
ScubaMultiFrameVolumeChartFactory::GetFactory() {

  static ScubaMultiFrameVolumeChartFactory sFactory;

  if ( !mbAddedTclCommands ) {

    mbAddedTclCommands = true;
    TclCommandManager& commandMgr = TclCommandManager::GetManager();
    commandMgr.AddCommand( sFactory, "MakeNewMultiFrameVolumeChart", 1,
                           "volumeID", "Makes a new chart window "
                           "plotting the voxel values over multiple frames." );
    commandMgr.AddCommand( sFactory, "MakeNewMultiFrameVolumeChartWithROI", 2,
                           "volumeID roiID", "Makes a new chart window "
                           "plotting the average of an ROI over multiple "
                           "frames." );
  }

  return sFactory;
}

TclCommandListener::TclCommandResult
ScubaMultiFrameVolumeChartFactory::DoListenToTclCommand ( char*  isCommand,
    int    iArgc,
    char** iasArgv ) {

  // MakeNewMultiFrameVolumeChart <volumeID>
  if ( 0 == strcmp( isCommand, "MakeNewMultiFrameVolumeChart" ) ) {

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

    // Create the chart.
    new ScubaMultiFrameVolumeChart( *vol );
  }

  // MakeNewMultiFrameVolumeChartWithROI <volumeID> <roiID>
  if ( 0 == strcmp( isCommand, "MakeNewMultiFrameVolumeChartWithROI" ) ) {

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


    // Create the chart.
    new ScubaMultiFrameVolumeChart( *vol, *roiVolume );
  }

  return ok;
}
