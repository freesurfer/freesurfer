/**
 * @file  vtkKWScubaLayerCollection.cxx
 * @brief A collection of layers for different vtkKWScubaView display modes
 *
 * A class that collects multiples layers designed to render one kind
 * of data in multiple display modes. The vtkKWScubaView will ask the
 * collection for a layer for a certain display mode, and the
 * collection is responsible for returning the right kind of layer,
 * and instantiating it and setting it up with initial data if
 * necessary. Subclasses should include function necessary to set
 * filenames for data to be loaded, so they can later create instances
 * of layers and load them. The collection can also keep common data 
 * objects that all layers will use.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:40 $
 *    $Revision: 1.4 $
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

#include "vtkKWScubaLayerCollection.h"
#include "vtkObjectFactory.h"
#include "vtkKWWidget.h"
#include "vtkKWEntry.h"
#include "vtkKWScaleWithEntry.h"
#include "vtkKWFrame.h"
#include "vtkKWEntryWithLabel.h"

using namespace std;

vtkStandardNewMacro( vtkKWScubaLayerCollection );
vtkCxxRevisionMacro( vtkKWScubaLayerCollection, "$Revision: 1.4 $" );

vtkKWScubaLayerCollection::vtkKWScubaLayerCollection () :
  Broadcaster( "vtkKWScubaLayerCollection" ),
  Listener( "vtkKWScubaLayerCollection" ),
  msLabel(""),
  mOpacity( 1.0 ),
  mbFastMode( false ),
  m2DRASZ( 0 ),
  mInPlane( 0 ),
  m3DRASX( 0 ),
  m3DRASY( 0 ),
  m3DRASZ( 0 ),
  mEntryLabel( NULL ),
  mScaleOpacity( NULL ),
  mLayerSettingsPanel(NULL), 
  mDisplayMode(vtkKWScubaView::TwoDee),
  mCurrentView( NULL ) {

}

vtkKWScubaLayerCollection::~vtkKWScubaLayerCollection () {

  ModeViewLayerMapType::iterator tMap;
  for( tMap = maLayers.begin(); tMap != maLayers.end(); ++tMap ) {
    ViewLayerMapType& map = tMap->second;
    ViewLayerMapType::iterator tLayer;
    for( tLayer = map.begin(); tLayer != map.end(); ++tLayer ) {
      vtkKWScubaLayer* layer = tLayer->second;
      if( layer )
	layer->Delete();
    }
  }
}

char const*
vtkKWScubaLayerCollection::GetLabel () const {

  return msLabel.c_str();
}

vtkKWScubaLayer* 
vtkKWScubaLayerCollection::GetCurrentLayer () const {

  ModeViewLayerMapType::const_iterator tMap;
  tMap = maLayers.find( mDisplayMode );
  if( tMap != maLayers.end() ) {
    const ViewLayerMapType& map = tMap->second;
    ViewLayerMapType::const_iterator tLayer;
    tLayer = map.find( mCurrentView );
    if( tLayer != map.end() )
      return tLayer->second;
  }

  return NULL;
}

void
vtkKWScubaLayerCollection::SetDisplayModeAndView ( vtkKWScubaView::DisplayMode iMode, vtkKWScubaView* iView ) {

  // This function needs to do two things. First, make sure we have a
  // layer for this display mode and view, and if not, we create one,
  // initiate it, and stash it in our map.
  ModeViewLayerMapType::const_iterator tMap;
  ViewLayerMapType::const_iterator tLayer;
  bool bCreate = false;

  // Try to find an entry for this display mode. If we do, also try to
  // find an entry for this view.
  tMap = maLayers.find( iMode );
  if( tMap != maLayers.end() ) {
    const ViewLayerMapType& map = tMap->second;
    tLayer = map.find( iView );
    if( tLayer == map.end () ||
	NULL == tLayer->second )
      bCreate = true;
  } else {
    bCreate = true;
  }

  // If we didn't find either...
  if( bCreate ) {

    // Use MakeLayerForDisplayMode to make a layer.
    try { 
      vtkKWScubaLayer* layer = this->MakeLayerForDisplayMode( iMode );
      if( NULL == layer )
	throw runtime_error( "MakeLayerForDisplayMode returned NULL" );
      
      // Set up the layer.
      layer->SetCollectionProperties( this );
      layer->SetViewProperties( iView );
      layer->SetApplication( this->GetApplication() );
      layer->Create();
      layer->LoadDataFromProperties();
 
      // Link us as listeners.
      AddListener( layer );
      layer->AddListener( this );

      // Now save a reference to the new layer.
      maLayers[iMode][iView] = layer;
    }
    catch( exception& e ) {
      string sError = string( "Error in SetDisplayMode: " ) + e.what();
      maLayers[iMode][iView] = NULL;
      throw runtime_error( sError );
    }    
  }
  
  // Next, we need to make sure that this layer is showing its
  // settings. If this is a different mode or view, we need to change
  // settings. We keep layer settings in a subpanel, so if it exists
  // (means we're also the current layer collection and we're showing
  // our settings), unpack the current layer and pack the new one.
  if( NULL != mLayerSettingsPanel &&
      (mDisplayMode != iMode || mCurrentView != iView) ) {

    // Try to unpack the current layer.
    vtkKWScubaLayer* oldLayer = this->GetCurrentLayer();
    if( NULL != oldLayer ) {
      mLayerSettingsPanel->UnpackChildren();
      oldLayer->DepopulateControlPage();
    }

    // Populate with the new layer.
    vtkKWScubaLayer* newLayer = maLayers[iMode][iView];
    newLayer->PopulateControlPage( mLayerSettingsPanel );
  }
  
  // Save the new mode and view.
  mDisplayMode = iMode;
  mCurrentView = iView;
}

void
vtkKWScubaLayerCollection::PopulateControlPage ( vtkKWWidget* iPanel ) {

  vtkKWEntryWithLabel* labeledEntry = vtkKWEntryWithLabel::New();
  labeledEntry->SetParent( iPanel );
  labeledEntry->Create();
  labeledEntry->SetLabelText( "Name: " );

  mEntryLabel = labeledEntry->GetWidget();
  mEntryLabel->SetValue( this->GetLabel() );
  mEntryLabel->SetCommand( this, "SetLabel" );

  mScaleOpacity = vtkKWScaleWithEntry::New();
  mScaleOpacity->SetParent( iPanel );
  mScaleOpacity->SetOrientationToHorizontal();
  mScaleOpacity->Create();
  mScaleOpacity->SetLabelText( "Opacity: " );
  mScaleOpacity->SetRange( 0, 1 );
  mScaleOpacity->SetResolution( 0.1 );
  mScaleOpacity->SetEntryWidth( 3 );
  mScaleOpacity->SetCommand( this, "SetOpacity" );
  mScaleOpacity->SetValue( mOpacity );
  
  // First pack our common collection controls.
  this->Script( "pack %s %s -side top -fill x -anchor nw",
                labeledEntry->GetWidgetName(),
                mScaleOpacity->GetWidgetName() );

  labeledEntry->Delete();

  // Add the collection subclass's controls.
  this->AddControls( iPanel );

  // Make a settings panel for our layer.
  mLayerSettingsPanel = vtkKWFrame::New();
  mLayerSettingsPanel->SetParent( iPanel );
  mLayerSettingsPanel->Create();

  // Pack the panel for the layer.
  this->Script( "pack %s -side top -fill x -anchor nw",
		mLayerSettingsPanel->GetWidgetName() );

  // If we have a current layer, tell the layer to populate this panel.
  if( NULL != this->GetCurrentLayer() )
    this->GetCurrentLayer()->PopulateControlPage( mLayerSettingsPanel );
}

void
vtkKWScubaLayerCollection::DepopulateControlPage () {


  // Tell the layer to remove its controls.
  if( this->GetCurrentLayer() )
    this->GetCurrentLayer()->DepopulateControlPage();
  
  if( mEntryLabel ) {
    mEntryLabel->Delete();
    mEntryLabel = NULL;
  }

  if( mScaleOpacity ) {
    mScaleOpacity->Delete();
    mScaleOpacity = NULL;
  }

  if( mLayerSettingsPanel ) {
    mLayerSettingsPanel->Delete();
    mLayerSettingsPanel = NULL;
  }

  // Tell the subclass to remove its controls.
  this->RemoveControls();
}

void
vtkKWScubaLayerCollection::DoListenToMessage ( std::string const isMessage,
					       void* const iData ) {
  // Rebroadcast messages to our layers.
  this->SendBroadcast( isMessage, iData );
}

void
vtkKWScubaLayerCollection::SetOpacity ( float iOpacity ) {
  mOpacity = iOpacity;

  // Set the scale control's value.
  if ( mScaleOpacity )
    mScaleOpacity->SetValue( iOpacity );

  this->OpacityChanged();
}

float
vtkKWScubaLayerCollection::GetOpacity () const {
  return mOpacity;
}

void
vtkKWScubaLayerCollection::OpacityChanged () {
  this->SendBroadcast( "OpacityChanged", NULL );
}

void
vtkKWScubaLayerCollection::SetFastMode ( bool ibFastMode ) {
  mbFastMode = ibFastMode;
  this->FastModeChanged();
}

bool
vtkKWScubaLayerCollection::GetFastMode () const {
  return mbFastMode;
}

void
vtkKWScubaLayerCollection::FastModeChanged () {
  this->SendBroadcast( "FastModeChanged" );
}

void
vtkKWScubaLayerCollection::Set2DInPlane ( int iPlane ) {
  if ( iPlane != mInPlane ) {
    mInPlane = iPlane;
    this->Layer2DInfoChanged();
  }
}

int
vtkKWScubaLayerCollection::Get2DInPlane () const {
  return mInPlane;
}

void
vtkKWScubaLayerCollection::Set2DRASZ ( float iRASZ ) {
  if ( iRASZ != m2DRASZ ) {
    m2DRASZ = iRASZ;
    this->Layer2DInfoChanged();
  }
}

float
vtkKWScubaLayerCollection::Get2DRASZ () const {
  return m2DRASZ;
}

void
vtkKWScubaLayerCollection::Layer2DInfoChanged () {
  this->SendBroadcast( "Layer2DInfoChanged" );
}

void
vtkKWScubaLayerCollection::Set3DRASX ( float i3DRASX ) {
  if ( i3DRASX != m3DRASX ) {
    m3DRASX = i3DRASX;
    this->Layer3DInfoChanged();
  }
}

float
vtkKWScubaLayerCollection::Get3DRASX () const {
  return m3DRASX;
}

void
vtkKWScubaLayerCollection::Set3DRASY ( float i3DRASY ) {
  if ( i3DRASY != m3DRASY ) {
    m3DRASY = i3DRASY;
    this->Layer3DInfoChanged();
  }
}

float
vtkKWScubaLayerCollection::Get3DRASY () const {
  return m3DRASY;
}

void
vtkKWScubaLayerCollection::Set3DRASZ ( float i3DRASZ ) {
  if ( i3DRASZ != m3DRASZ ) {
    m3DRASZ = i3DRASZ;
    this->Layer3DInfoChanged();
  }
}

float
vtkKWScubaLayerCollection::Get3DRASZ () const {
  return m3DRASZ;
}

void
vtkKWScubaLayerCollection::Layer3DInfoChanged () {
  this->SendBroadcast( "Layer3DInfoChanged" );
}

void
vtkKWScubaLayerCollection::AddControls ( vtkKWWidget* iPanel ) {

}

void
vtkKWScubaLayerCollection::RemoveControls () {
  
  if( mEntryLabel ) {
    mEntryLabel->Delete();
    mEntryLabel = NULL;
  }

  if( mScaleOpacity ) {
    mScaleOpacity->Delete();
    mScaleOpacity = NULL;
  }
}

void
vtkKWScubaLayerCollection::SetLabel ( const char* isLabel ) {

  // If we have another collection with the same label, append a
  // number to our name.
  char sLabel[1024];
  strncpy( sLabel, isLabel, sizeof(sLabel) );
  int nName = 2;
  bool bMatch = false;
  list<vtkKWScubaLayerCollection*> lCol;
  list<vtkKWScubaLayerCollection*>::const_iterator tCol;
  vtkKWScubaLayerCollection::GetPointerList( lCol );
  do { 
    bMatch = false;
    for( tCol = lCol.begin(); tCol != lCol.end(); ++tCol ) {
      if( *tCol ) {
	vtkKWScubaLayerCollection const* testCol = *tCol;
	if( testCol != this ) {
	  char const* sTestLabel = testCol->GetLabel();
	  if( string(sLabel) == sTestLabel ) {
	    snprintf( sLabel, sizeof(sLabel), "%s (%d)", isLabel, nName++ );
	    bMatch = true;
	  }
	}
      }
    }
  } while( bMatch );
  
  msLabel = sLabel;
}

vtkKWScubaLayer*
vtkKWScubaLayerCollection::MakeLayerForDisplayMode ( vtkKWScubaView::DisplayMode iMode ) {
  
  throw runtime_error( "MakeLayerForDisplayMode not impelemented" );
}
