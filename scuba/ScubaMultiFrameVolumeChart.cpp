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
  mChart(NULL) {

  // Create the chart and listen to it.
  mChart = ChartWindow::NewChartWindow();
  if( NULL == mChart ) {
    throw runtime_error( "Couldn't create chart." );
  }

  // Listen to some stuff.
  mChart->AddListener( this );
  mVolume.AddListener( this );

  // Listen to the view boradcaster for cursor changes.
  ScubaViewBroadcaster::GetBroadcaster().AddListener( this );

  // We want to listen to any layers that use this volume, so we can
  // update the current timepoint and condition when the layer changes
  // it. Find all layers, and listen to them.
  list<int> lIDs;
  Layer::GetIDList( lIDs );
  list<int>::iterator tID;
  for( tID = lIDs.begin(); tID != lIDs.end(); ++tID ) {
    int id = *tID;
    Layer& layer = Layer::FindByID( id );
    DataCollection* col = layer.GetMainDataCollection();
    if( NULL != col &&
	col->GetID() == mVolume.GetID() ) {
      layer.AddListener( this );
    }
  }

  Draw();
}

ScubaMultiFrameVolumeChart::~ScubaMultiFrameVolumeChart () {

  // Stop listening.
  //  mChart->RemoveListener( this );
  mVolume.RemoveListener( this );

  // Stop listening to all those listeners.
  list<int> lIDs;
  Layer::GetIDList( lIDs );
  list<int>::iterator tID;
  for( tID = lIDs.begin(); tID != lIDs.end(); ++tID ) {
    int id = *tID;
    Layer& layer = Layer::FindByID( id );
    DataCollection* col = layer.GetMainDataCollection();
    if( NULL != col &&
	col->GetID() == mVolume.GetID() ) {
      layer.RemoveListener( this );
    }
  }
}

void
ScubaMultiFrameVolumeChart::DoListenToMessage ( string isMessage,
						void*  iData ) {

  if( isMessage == "dataChanged" ) {
    Draw();
  }
  if( isMessage == "layerChanged" ) {
    Draw();
  }
  if( isMessage == "cursorChanged" ) {
    Draw();
  }
  if( isMessage == "chartDeleted" ) {
    delete this;
  }
}

void
ScubaMultiFrameVolumeChart::Draw() {

  // Clear it.
  mChart->ClearData();

  list<ChartWindow::PointData> lData;
  ChartWindow::PointData data;

  // Get the cursor.
  float cursorRAS[3];
  ScubaView::GetCursor( cursorRAS );

  // Make a location in the volume.
  VolumeLocation& loc = 
    (VolumeLocation&) mVolume.MakeLocationFromRAS( cursorRAS );

  // Check against the bounds of the volume.
  if( mVolume.IsInBounds( loc ) ) {

    if( mVolume.InterpretFramesAsTimePoints() ) {

      mChart->SetXAxisLabel( "Frame" );
      mChart->SetShowLegend( true );

      // We're going to make a group for each condition, and go
      // through the time points adding points in that group.
      for( int nCond = 0; nCond < mVolume.GetNumberOfConditions(); nCond++ ) {

	stringstream ssLabel;
	ssLabel << "Condition " << nCond;
	mChart->SetGroupLabel( nCond, ssLabel.str() );
	mChart->SetGroupConnected( nCond, 1 );

	for( int nTP = 0; nTP < mVolume.GetNumberOfTimePoints(); nTP++ ) {

	  int nFrame = 
	    mVolume.ConvertConditionAndTimePointToFrame( nCond, nTP );
	  loc.SetFrame( nFrame );
	  float value = mVolume.GetMRINearestValue( loc );

	  data.mX = nTP;
	  data.mY = value;
	  
	  mChart->AddPointData( nCond, data );
	}
      }

    } else {
      
      // No conditions, so just loop through all the frames and get
      // the value at each one. Plot those points.
      mChart->SetXAxisLabel( "Frame" );
      mChart->SetShowLegend( true );

      mChart->SetGroupConnected( 0, 1 );

      for( int nFrame = 0; nFrame < mVolume.GetNumberOfFrames(); nFrame++ ) {

	loc.SetFrame( nFrame );
	float value = mVolume.GetMRINearestValue( loc );

	data.mX = nFrame;
	data.mY = value;
	stringstream ssLabel;
	ssLabel << "Frame " << nFrame;
	data.msLabel = ssLabel.str();

	mChart->AddPointData( 0, data );
      }
    }
  }

  delete &loc;


  // Go through all the layers. If they are showing this volume, get
  // the current time point and condition from the layer. Draw a time
  // point marker.
  list<int> lIDs;
  Layer::GetIDList( lIDs );
  list<int>::iterator tID;
  for( tID = lIDs.begin(); tID != lIDs.end(); ++tID ) {
    int id = *tID;
    Layer& layer = Layer::FindByID( id );
    DataCollection* col = layer.GetMainDataCollection();
    if( NULL != col &&
	col->GetID() == mVolume.GetID() ) {

      ScubaLayer2DMRI& mriLayer = (ScubaLayer2DMRI&) layer;

      ChartWindow::MarkerData marker;

      if( mVolume.InterpretFramesAsTimePoints() ) {
	
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
  sTitle << "Multiframe Volume Chart: Vol " << mVolume.GetID ()
	 << " RAS " << Point3<float>(cursorRAS);

  sInfo << "Volume: \\\"" << mVolume.GetLabel() 
	<< "\\\" RAS: \\\"" << Point3<float>(cursorRAS);

  // Set all the data chart and tell it to draw.
  mChart->SetTitle( sTitle.str() );
  mChart->SetInfo( sInfo.str() );
  mChart->SetYAxisLabel( "Value" );
  mChart->Draw();
}


ScubaMultiFrameVolumeChartFactory& 
ScubaMultiFrameVolumeChartFactory::GetFactory() {

  static ScubaMultiFrameVolumeChartFactory sFactory;

  if( !mbAddedTclCommands ) {

    mbAddedTclCommands = true;
    TclCommandManager& commandMgr = TclCommandManager::GetManager();
    commandMgr.AddCommand( sFactory, "MakeNewMultiFrameVolumeChart", 1,
			   "volumeID", "Makes a new chart window "
			   "plotting the voxel values over multiple frames." );
  }

  return sFactory;
}

TclCommandListener::TclCommandResult
ScubaMultiFrameVolumeChartFactory::DoListenToTclCommand ( char*  isCommand,
							  int    iArgc,
							  char** iasArgv ) {
  
  // MakeNewMultiFrameVolumeChart <volumeID>
  if( 0 == strcmp( isCommand, "MakeNewMultiFrameVolumeChart" ) ) {

    // Get the collection ID.
    int collectionID;
    try {
      collectionID = TclCommandManager::ConvertArgumentToInt( iasArgv[1] );
    }
    catch( runtime_error& e ) {
      sResult = string("bad collectionID: ") + e.what();
      return error;
    }

    // Try to find the object.
    DataCollection* col;
    try { 
      col = &DataCollection::FindByID( collectionID );
    }
    catch( runtime_error& e ) {
      sResult = string("Couldn't find collection: ") + e.what();
      return error;
    }
    
    // Make sure it's a volume.
    if( col->GetTypeDescription() != "Volume" ) {
      throw runtime_error( "Collection wasn't a volume." );
    }

    // We can do this because we checked the type.
    VolumeCollection* vol = (VolumeCollection*) col;

    // Create the chart.
    new ScubaMultiFrameVolumeChart( *vol );
  }

  return ok;
}
