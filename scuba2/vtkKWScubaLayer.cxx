/**
 * @file  vtkKWScubaLayer.cxx
 * @brief Base class for layers
 *
 * Provides basic event stubs for layers. Subclasses
 * should load data, create vtk objects, use AddProp() to add them to
 * the view, and handle events.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:40 $
 *    $Revision: 1.5 $
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


#include <stdexcept>
#include "vtkKWScubaLayer.h"
#include "vtkObjectFactory.h"
#include "vtkKWScubaLayerCollection.h"
#include "vtkKWWidget.h"
#include "vtkProp.h"
#include "vtkPropCollection.h"

using namespace std;

vtkStandardNewMacro( vtkKWScubaLayer );
vtkCxxRevisionMacro( vtkKWScubaLayer, "$Revision: 1.5 $" );

map<vtkProp*,vtkKWScubaLayer*> vtkKWScubaLayer::mPropToLayerMap;

vtkKWScubaLayer::vtkKWScubaLayer () :
    Broadcaster( "vtkKWScubaLayer" ),
    Listener( "vtkKWScubaLayer" ),
    mProperties( NULL ),
    mViewProperties( NULL ),
    mProps( NULL ),
    mbInfoChanged( false ) {

  // Create our prop collection.
  mProps = vtkPropCollection::New();
}

vtkKWScubaLayer::~vtkKWScubaLayer () {

  // Remove our props from the map.
  mProps->InitTraversal();
  for ( vtkProp* prop = mProps->GetNextProp();
        NULL != prop;
        prop = mProps->GetNextProp() ) {
    mPropToLayerMap[prop] = NULL;
  }
  
  // Delete the props.
  mProps->Delete();
}

void
vtkKWScubaLayer::SetCollectionProperties ( ScubaCollectionProperties const* iProperties ) {
  mProperties = iProperties;
}

void
vtkKWScubaLayer::SetViewProperties ( ScubaViewProperties const* iProperties ) {
  mViewProperties = iProperties;
}

void
vtkKWScubaLayer::DoListenToMessage ( string const isMessage,
				     void* const ) {

  if( isMessage == "DataAvailabilityChanged" ) {
    this->LoadDataFromProperties();
  }
}

const char*
vtkKWScubaLayer::GetLabel () const {

  if( mProperties )
    return mProperties->GetLabel();
  else
    return NULL;
}

vtkPropCollection*
vtkKWScubaLayer::GetPropCollection () const {
  return mProps;
}

void
vtkKWScubaLayer::PopulateControlPage ( vtkKWWidget* iPanel ) {
  // Add the subclass's controls.
  this->AddControls( iPanel );
}

void
vtkKWScubaLayer::DepopulateControlPage () {
  // Tell subclassses to remove controls.
  this->RemoveControls();
}

void
vtkKWScubaLayer::GetRASBounds ( float ioBounds[6] ) const {
  ioBounds[0] = ioBounds[1] = ioBounds[2] =
    ioBounds[3] = ioBounds[4] = ioBounds[5] = 0;
}


void
vtkKWScubaLayer::Get2DRASZIncrementHint ( float ioHint[3] ) const {
  ioHint[0] = ioHint[1] = ioHint[2] = 1;
}

bool
vtkKWScubaLayer::IsInfoChanged () const {
  return mbInfoChanged;
}

void
vtkKWScubaLayer::InfoUpdated () {
  mbInfoChanged = false;
}

void
vtkKWScubaLayer::AddProp ( vtkProp* iProp ) {

  // Add it to our list to render. Link it to us in the map.
  mProps->AddItem( iProp );
  mPropToLayerMap[iProp] = this;

  // Notify the view so we can draw it.
  int ID = this->GetID();
  this->SendBroadcast( "PropListChanged", (void*)&ID );
}

void
vtkKWScubaLayer::RemoveProp ( vtkProp* iProp ) {

  // Remove it from our list to render. Unlink it to us in the map.
  if( mProps->IsItemPresent( iProp ) ) {
    mProps->RemoveItem( iProp );
    mPropToLayerMap.erase( iProp );

    // Notify the view so it can remove it.
    int ID = this->GetID();
    this->SendBroadcast( "PropListChanged", (void*)&ID );
  }
}

void
vtkKWScubaLayer::RemoveAllProps () {

  // Remove all items in the collection. First go through them and
  // unlink them to us.
  mProps->InitTraversal();
  for ( vtkProp* prop = mProps->GetNextProp();
	NULL != prop;
	prop = mProps->GetNextProp() ) {
    mPropToLayerMap.erase( prop );
  }

  mProps->RemoveAllItems();
}

void
vtkKWScubaLayer::AddControls ( vtkKWWidget* iPanel ) {}

void
vtkKWScubaLayer::RemoveControls () {}

void
vtkKWScubaLayer::InfoChanged () {
  mbInfoChanged = true;
}

void
vtkKWScubaLayer::PipelineChanged () {
  int ID = this->GetID();
  this->SendBroadcast( "PipelineChanged", (void*)&ID );
}

vtkKWScubaLayer*
vtkKWScubaLayer::GetLayerFromProp( vtkProp* iProp ) {

  vtkKWScubaLayer* layer = mPropToLayerMap[iProp];
  if ( NULL == layer ) {
    throw runtime_error( "Couldn't find layer." );
  }
  return layer;
}

