/**
 * @file  vtkKWScubaView.cxx
 * @brief A VTK view containing ScubaLayers
 *
 * Implementation of a vtkKWRenderWidget that can have multiple
 * ScubaLayers. Also handles passing events down to them and managing
 * their viewing state.
 *
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:40 $
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

#include <limits>
#include <stdexcept>
#include <assert.h>
#include "IconLoader.h"
#include "ScubaInfoItem.h"
#include "vtkAxes.h"
#include "vtkCamera.h"
#include "vtkCellPicker.h"
#include "vtkInteractorObserver.h"
#include "vtkKWComboBox.h"
#include "vtkKWComboBoxWithLabel.h"
#include "vtkKWFrameWithLabel.h"
#include "vtkKWIcon.h"
#include "vtkKWScubaLayer.h"
#include "vtkKWScubaLayerCollection.h"
#include "vtkKWScubaTool.h"
#include "vtkKWScubaView.h"
#include "vtkKWScubaWindow.h"
#include "vtkKWPushButton.h"
#include "vtkKWRadioButton.h"
#include "vtkKWRadioButtonSet.h"
#include "vtkKWScaleWithEntry.h"
#include "vtkKWToolbar.h"
#include "vtkLight.h"
#include "vtkObjectFactory.h"
#include "vtkPolyDataMapper.h"
#include "vtkProp3DCollection.h"
#include "vtkPropCollection.h"
#include "vtkProperty.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkRenderer.h"
#include "vtkTransform.h"
#include "vtkTransformPolyDataFilter.h"
#include "vtkTubeFilter.h"

using namespace std;

vtkStandardNewMacro( vtkKWScubaView );
vtkCxxRevisionMacro( vtkKWScubaView, "$Revision: 1.9 $" );

map<vtkRenderWindow*,vtkSmartPointer<vtkKWScubaView> > vtkKWScubaView::mRenderWindowToViewMap;

vtkKWScubaView::vtkKWScubaView () :
    vtkKWRenderWidget(),
    Listener( "vtkKWScubaView" ),
    Broadcaster( "vtkKWScubaView" ),
    msLabel( "" ),
    mbDirty( false ),
    mbInfoChanged( false ),
    mDisplayMode( TwoDee ),
    mbFastMode( false ),
    m2DRASZ( 0 ),
    m2DInPlane( -1 ),
    m3DRASX( 0 ),
    m3DRASY( 0 ),
    m3DRASZ( 0 ) {

  mCursorRASCoords[0] = 0;
  mCursorRASCoords[1] = 0;
  mCursorRASCoords[2] = 0;

  mMouseOverRASCoords[0] = 0;
  mMouseOverRASCoords[1] = 0;
  mMouseOverRASCoords[2] = 0;

}

vtkKWScubaView::~vtkKWScubaView () {

  // Clear our entry in the static map.
  if ( NULL != this->GetRenderWindow() )
    mRenderWindowToViewMap.erase( this->GetRenderWindow() );
}

void
vtkKWScubaView::CreateWidget () {

  // Create the superclass widget. This also sets up our
  // vtkRenderWindow stuff.
  this->Superclass::CreateWidget();

  // Make our cameras.
  m2DCamera = vtkSmartPointer<vtkCamera>::New();
  m3DCamera = vtkSmartPointer<vtkCamera>::New();

  this->GetRenderer()->SetActiveCamera( m2DCamera );
  
  // set the background to something neutral
  this->GetRenderer()->SetBackground( 0.7, 0.7, 0.9 );
  
  // Sanity check.
  assert( this->GetRenderWindow() );
  assert( this->GetRenderer() );
  
  // Make some head lights. This still needs tweaking.
  double color[] = { 0.1, 0.1, 0.1 };
  vtkSmartPointer<vtkLight> light = 
    vtkSmartPointer<vtkLight>::New();
  light->SetPosition( 1, 0, 0 );
  light->SetLightTypeToHeadlight();  
  light->SetAmbientColor( color );
  light->SetDiffuseColor( color );
  light->SetSpecularColor( color );
  this->GetRenderer()->AddLight( light );

  // Set the picker.
  vtkSmartPointer<vtkCellPicker> picker = 
    vtkSmartPointer<vtkCellPicker>::New();
  this->GetRenderWindowInteractor()->SetPicker( picker );

  // Associate us with the window.
  mRenderWindowToViewMap[this->GetRenderWindow()] = this;

#if 0
  // Axes pipeline.
  vtkSmartPointer<vtkAxes> axes = 
    vtkSmartPointer<vtkAxes>::New();
  axes->SetOrigin( 0, 0, 0 );
  axes->SetScaleFactor( 20 );

  vtkSmartPointer<vtkTubeFilter> axesTubes = 
    vtkSmartPointer<vtkTubeFilter>::New();
  axesTubes->SetInputConnection( axes->GetOutputPort() );
  axesTubes->SetRadius( 1 );
  axesTubes->SetNumberOfSides( 6 );

  vtkSmartPointer<vtkPolyDataMapper> axesMapper = 
    vtkSmartPointer<vtkPolyDataMapper>::New();
  axesMapper->SetInputConnection( axesTubes->GetOutputPort() );

  vtkSmartPointer<vtkActor> axesActor = 
    vtkSmartPointer<vtkActor>::New();
  axesActor->SetMapper( axesMapper );
  axesActor->SetPickable( false );
  axesActor->GetProperty()->SetOpacity( 0.5 );

  this->GetRenderer()->AddActor( axesActor );
#endif

  // Reset the camera to show all the actors.
  this->GetRenderer()->ResetCamera();

  // Set up our initial view.
  this->Set2DInPlane( 0 );
}

void
vtkKWScubaView::SetLabel ( const char* isLabel ) {
  msLabel = isLabel;
}

const char*
vtkKWScubaView::GetLabel () const {
  return msLabel.c_str();
}

void
vtkKWScubaView::PopulateToolbar ( vtkKWToolbar* iToolbar ) {

  // Zoom Out
  mBtnZoomOut = vtkSmartPointer<vtkKWPushButton>::New();
  mBtnZoomOut->SetParent( iToolbar->GetFrame() );
  mBtnZoomOut->Create();
  mBtnZoomOut->SetText( "Zoom Out" );
  mBtnZoomOut->SetBalloonHelpString( "Zoom Out" );
  mBtnZoomOut->SetImageToPredefinedIcon( vtkKWIcon::IconMagGlass );
  mBtnZoomOut->SetCommand( this, "ZoomOut" );
  try { IconLoader::SetPushButtonIcon( "ZoomOut", mBtnZoomOut ); }
  catch( exception& e ) { cerr << "Error loading icon: " << e.what() << endl; }
  iToolbar->AddWidget( mBtnZoomOut );

  // Zoom In
  mBtnZoomIn = vtkSmartPointer<vtkKWPushButton>::New();
  mBtnZoomIn->SetParent( iToolbar->GetFrame() );
  mBtnZoomIn->Create();
  mBtnZoomIn->SetText( "Zoom In" );
  mBtnZoomIn->SetBalloonHelpString( "Zoom In" );
  mBtnZoomIn->SetImageToPredefinedIcon( vtkKWIcon::IconMagGlass );
  mBtnZoomIn->SetCommand( this, "ZoomIn" );
  try { IconLoader::SetPushButtonIcon( "ZoomIn", mBtnZoomIn ); }
  catch( exception& e ) { cerr << "Error loading icon: " << e.what() << endl; }
  iToolbar->AddWidget( mBtnZoomIn );

  // Rotate X Pos
  mBtnRotateXPos = vtkSmartPointer<vtkKWPushButton>::New();
  mBtnRotateXPos->SetParent( iToolbar->GetFrame() );
  mBtnRotateXPos->Create();
  iToolbar->AddWidget( mBtnRotateXPos );
  mBtnRotateXPos->SetText( "Rotate" );
  mBtnRotateXPos->SetBalloonHelpString( "Rotate" );
  mBtnRotateXPos->SetCommand( this, "AnimateCameraElevatePositive" );
  try { IconLoader::SetPushButtonIcon( "RotateXPos", mBtnRotateXPos); }
  catch (...) {}

  // Rotate X Neg
  mBtnRotateXNeg = vtkSmartPointer<vtkKWPushButton>::New();
  mBtnRotateXNeg->SetParent( iToolbar->GetFrame() );
  mBtnRotateXNeg->Create();
  iToolbar->AddWidget( mBtnRotateXNeg );
  mBtnRotateXNeg->SetText( "Rotate" );
  mBtnRotateXNeg->SetBalloonHelpString( "Rotate" );
  mBtnRotateXNeg->SetCommand( this, "AnimateCameraElevateNegative" );
  try { IconLoader::SetPushButtonIcon( "RotateXNeg", mBtnRotateXNeg); }
  catch (...) {}

  // Rotate Y Pos
  mBtnRotateYPos = vtkSmartPointer<vtkKWPushButton>::New();
  mBtnRotateYPos->SetParent( iToolbar->GetFrame() );
  mBtnRotateYPos->Create();
  iToolbar->AddWidget( mBtnRotateYPos );
  mBtnRotateYPos->SetText( "Rotate" );
  mBtnRotateYPos->SetBalloonHelpString( "Rotate" );
  mBtnRotateYPos->SetCommand( this, "AnimateCameraAzimuthNegative" );
  try { IconLoader::SetPushButtonIcon( "RotateYPos", mBtnRotateYPos); }
  catch (...) {}

  // Rotate Y Neg
  mBtnRotateYNeg = vtkSmartPointer<vtkKWPushButton>::New();
  mBtnRotateYNeg->SetParent( iToolbar->GetFrame() );
  mBtnRotateYNeg->Create();
  iToolbar->AddWidget( mBtnRotateYNeg );
  mBtnRotateYNeg->SetText( "Rotate" );
  mBtnRotateYNeg->SetBalloonHelpString( "Rotate" );
  mBtnRotateYNeg->SetCommand( this, "AnimateCameraAzimuthPositive" );
  try { IconLoader::SetPushButtonIcon( "RotateYNeg", mBtnRotateYNeg); }
  catch (...) {}

  // Rotate Z Pos
  mBtnRotateZPos = vtkSmartPointer<vtkKWPushButton>::New();
  mBtnRotateZPos->SetParent( iToolbar->GetFrame() );
  mBtnRotateZPos->Create();
  iToolbar->AddWidget( mBtnRotateZPos );
  mBtnRotateZPos->SetText( "Rotate" );
  mBtnRotateZPos->SetBalloonHelpString( "Rotate" );
  mBtnRotateZPos->SetCommand( this, "AnimateCameraRollNegative" );
  try { IconLoader::SetPushButtonIcon( "RotateZPos", mBtnRotateZPos); }
  catch (...) {}

  // Rotate Z Neg
  mBtnRotateZNeg = vtkSmartPointer<vtkKWPushButton>::New();
  mBtnRotateZNeg->SetParent( iToolbar->GetFrame() );
  mBtnRotateZNeg->Create();
  iToolbar->AddWidget( mBtnRotateZNeg );
  mBtnRotateZNeg->SetText( "Rotate" );
  mBtnRotateZNeg->SetBalloonHelpString( "Rotate" );
  mBtnRotateZNeg->SetCommand( this, "AnimateCameraRollPositive" );
  try { IconLoader::SetPushButtonIcon( "RotateZNeg", mBtnRotateZNeg); }
  catch (...) {}

}

void
vtkKWScubaView::DepopulateToolbar () {

  mBtnZoomOut = NULL;
  mBtnZoomIn = NULL;
  mBtnRotateXPos = NULL;
  mBtnRotateXNeg = NULL;
  mBtnRotateYPos = NULL;
  mBtnRotateYNeg = NULL;
  mBtnRotateZPos = NULL;
  mBtnRotateZNeg = NULL;
}


void
vtkKWScubaView::PopulateControlPage ( vtkKWWidget* iPanel ) {

  // Display mode radio buttons ------------------------------------------
  vtkSmartPointer<vtkKWRadioButtonSet> radBtnSetDisplayMode =
    vtkSmartPointer<vtkKWRadioButtonSet>::New();
  radBtnSetDisplayMode->SetParent( iPanel );
  radBtnSetDisplayMode->Create();
  radBtnSetDisplayMode->PackHorizontallyOn();

  mRadBtnDisplayMode2D.
    TakeReference( radBtnSetDisplayMode->AddWidget( (int)TwoDee ) );
  mRadBtnDisplayMode2D->SetText( "2D" );
  mRadBtnDisplayMode2D->SetCommand( this, "SetDisplayModeTo2D" );
  if( TwoDee == mDisplayMode )
    mRadBtnDisplayMode2D->SelectedStateOn();

  mRadBtnDisplayMode3D.
    TakeReference( radBtnSetDisplayMode->AddWidget( (int)ThreeDee ) );
  mRadBtnDisplayMode3D->SetText( "3D" );
  mRadBtnDisplayMode3D->SetCommand( this, "SetDisplayModeTo3D" );
  if( ThreeDee == mDisplayMode )
    mRadBtnDisplayMode3D->SelectedStateOn();
  
  this->Script( "pack %s -side top -fill x -anchor nw",
		radBtnSetDisplayMode->GetWidgetName() );

  // Frame for the menus ------------------------------------------------
  // We'll pack the individual controls in the PackDisplayModeControls()
  // function.
  vtkSmartPointer<vtkKWFrameWithLabel> labeledFrame = 
    vtkSmartPointer<vtkKWFrameWithLabel>::New();
  labeledFrame->SetParent( iPanel );
  labeledFrame->Create();
  labeledFrame->SetLabelText( "Display Mode controls" );

  mFrameDisplayModeControls = labeledFrame->GetFrame();
  
  this->Script( "pack %s -side top -fill x -anchor nw",
		labeledFrame->GetWidgetName() );

  // 2D in plane radio buttons -------------------------------------------
  mRadBtnSet2DInPlane = vtkSmartPointer<vtkKWRadioButtonSet>::New();
  mRadBtnSet2DInPlane->SetParent( mFrameDisplayModeControls );
  mRadBtnSet2DInPlane->Create();
  mRadBtnSet2DInPlane->PackHorizontallyOn();

  mRadBtn2DInPlaneX.
    TakeReference( mRadBtnSet2DInPlane->AddWidget( 0 ) );
  mRadBtn2DInPlaneX->SetText( "X" );
  mRadBtn2DInPlaneX->SetCommand( this, "Set2DInPlane 0" );
  mRadBtn2DInPlaneX->SetBalloonHelpString( "X In Plane" );
  mRadBtn2DInPlaneX->SelectedStateOn();

  mRadBtn2DInPlaneY.
    TakeReference( mRadBtnSet2DInPlane->AddWidget( 1 ) );
  mRadBtn2DInPlaneY->SetText( "Y" );
  mRadBtn2DInPlaneY->SetCommand( this, "Set2DInPlane 1" );
  mRadBtn2DInPlaneY->SetBalloonHelpString( "Y In Plane" );

  mRadBtn2DInPlaneZ.
    TakeReference( mRadBtnSet2DInPlane->AddWidget( 2 ) );
  mRadBtn2DInPlaneZ->SetText( "Z" );
  mRadBtn2DInPlaneZ->SetCommand( this, "Set2DInPlane 2" );
  mRadBtn2DInPlaneZ->SetBalloonHelpString( "Z In Plane" );

  vtkSmartPointer<vtkKWRadioButton> btn =
    mRadBtnSet2DInPlane->GetWidget( this->Get2DInPlane() );
  if ( btn.GetPointer() ) {
    btn->SelectedStateOn();
  }
  
  // 2dRASZ scale -------------------------------------------------------
  mScale2DRASZ = vtkSmartPointer<vtkKWScaleWithEntry>::New();
  mScale2DRASZ->SetParent( mFrameDisplayModeControls );
  mScale2DRASZ->SetLabelText( "X RAS" );
  mScale2DRASZ->SetOrientationToHorizontal();
  mScale2DRASZ->SetCommand( this, "Set2DRASZ " );
  mScale2DRASZ->Create();

  float min = 0, max = 0;
  this->Get2DRASZRange( min, max );
  mScale2DRASZ->SetRange( min, max );
  mScale2DRASZ->SetResolution( this->Get2DRASZIncrementHint() );
  mScale2DRASZ->SetValue( this->Get2DRASZ() );

  // 3DRAS scales -------------------------------------------------------
  mScale3DRASX = vtkSmartPointer<vtkKWScaleWithEntry>::New();
  mScale3DRASX->SetParent( mFrameDisplayModeControls );
  mScale3DRASX->SetLabelText( "X RAS" );
  mScale3DRASX->SetOrientationToHorizontal();
  mScale3DRASX->SetCommand( this, "Set3DRASX " );
  mScale3DRASX->Create();
  mScale3DRASX->SetRange( min, max );
  mScale3DRASX->SetResolution( this->Get2DRASZIncrementHint() );
  mScale3DRASX->SetValue( this->Get3DRASX() );

  mScale3DRASY = vtkSmartPointer<vtkKWScaleWithEntry>::New();
  mScale3DRASY->SetParent( mFrameDisplayModeControls );
  mScale3DRASY->SetLabelText( "Y RAS" );
  mScale3DRASY->SetOrientationToHorizontal();
  mScale3DRASY->SetCommand( this, "Set3DRASY " );
  mScale3DRASY->Create();
  mScale3DRASY->SetRange( min, max );
  mScale3DRASY->SetResolution( this->Get2DRASZIncrementHint() );
  mScale3DRASY->SetValue( this->Get3DRASY() );

  mScale3DRASZ = vtkSmartPointer<vtkKWScaleWithEntry>::New();
  mScale3DRASZ->SetParent( mFrameDisplayModeControls );
  mScale3DRASZ->SetLabelText( "Z RAS" );
  mScale3DRASZ->SetOrientationToHorizontal();
  mScale3DRASZ->SetCommand( this, "Set3DRASZ " );
  mScale3DRASZ->Create();
  mScale3DRASZ->SetRange( min, max );
  mScale3DRASZ->SetResolution( this->Get2DRASZIncrementHint() );
  mScale3DRASZ->SetValue( this->Get3DRASZ() );

  // Call this to pack the right display mode controls.
  this->PackDisplayModeControls();

  // Layer combo boxes --------------------------------------------------

  // Make the frame for the menus.
  mFrameSlotMenus = vtkSmartPointer<vtkKWFrame>::New();
  mFrameSlotMenus->SetParent( iPanel );
  mFrameSlotMenus->Create();

  this->Script( "pack %s -side top -fill x -anchor nw",
		mFrameSlotMenus->GetWidgetName() );

  // We make a menu for every existing slot that we have, even if that
  // slot is empty. This may skip slots that aren't being used.
  SlotCollectionMapType::iterator tCol;
  for( tCol = maCol.begin(); tCol != maCol.end(); ++tCol )
    this->MakeMenuForSlot( tCol->first );

  // Populate the menus we just created and set their current values.
  this->UpdateLayerComboBoxes();

}

void
vtkKWScubaView::DepopulateControlPage () {

  mRadBtnDisplayMode2D = NULL;
  mRadBtnDisplayMode3D = NULL;

  // Delete all menus and clear the map.
  maSlotMenu.clear();

  mFrameSlotMenus = NULL;
  mRadBtnSet2DInPlane = NULL;
  mRadBtn2DInPlaneX = NULL;
  mRadBtn2DInPlaneY = NULL;
  mRadBtn2DInPlaneZ = NULL;
  mScale2DRASZ = NULL;
  mScale3DRASX = NULL;
  mScale3DRASY = NULL;
  mScale3DRASZ = NULL;
  
  mFrameDisplayModeControls = NULL;
}

vtkKWScubaLayerCollection*
vtkKWScubaView::GetCollectionAtSlot ( int inSlot ) const {

  // Try and find this slot number and return the contents if found
  // (could be NULL, meaning the slot is empty). If there is no slot
  // defined, erturn NULL.
  SlotCollectionMapType::const_iterator tCol = maCol.find( inSlot );
  if( tCol != maCol.end() )
    return tCol->second;
  else 
    return NULL;
}

int
vtkKWScubaView::GetFirstUnusedLayerSlot () const {

  // Go through the map looking for two things: 1) A skipped slot. If
  // a slot is skipped, return the index of the first skipped slot. 2)
  // An empty slot.
  int nLastSlot = -1;
  SlotCollectionMapType::const_iterator tCol;
  for( tCol = maCol.begin(); tCol != maCol.end(); ++tCol ) {

    int nSlot = tCol->first;
    vtkKWScubaLayerCollection* col = tCol->second;

    if( nSlot != nLastSlot + 1 )
      return nLastSlot + 1;	// Skipped a slot.

    if( NULL == col )
      return nSlot;		// Empty slot.

    nLastSlot = nSlot;
  }
  
  // Just return the next slot number.
  return nLastSlot + 1;
}

int
vtkKWScubaView::GetHighestFilledLayerSlot () const {

  // Go through all the slots and find the highest filled one.
  int nHighestSlot = -1;
  SlotCollectionMapType::const_iterator tCol;
  for( tCol = maCol.begin(); tCol != maCol.end(); ++tCol )
    if( NULL != tCol->second.GetPointer() && // Slot is filled
	tCol->first > nHighestSlot ) // Slot number is higher
      nHighestSlot = tCol->first;
  
  return nHighestSlot;
}

void
vtkKWScubaView::SetLayerCollectionAtSlot ( int inSlot, 
					   vtkKWScubaLayerCollection* iCol ) {

  if( inSlot < 0 )
    throw runtime_error( "Invalid slot number" );

  // Make sure we don't already have this layer in any slot, or if
  // we're setting a slot to NULL.
  if ( NULL == iCol ||
       (iCol && -1 == this->GetSlotOfLayerCollection( iCol )) ) {

    // Set the slot.
    maCol[inSlot] = iCol;

    // If we're not just setting it to NULL...
    if ( NULL != maCol[inSlot].GetPointer() ) {

      // Listen to this layer.
      maCol[inSlot]->AddListener( this );

      // Set the current display mode.
      maCol[inSlot]->SetDisplayModeAndView( mDisplayMode, this );

      // Set the 2D info in the layer collection from our info.
      maCol[inSlot]->Set2DInPlane( m2DInPlane );
      maCol[inSlot]->Set2DRASZ( m2DRASZ );
    }

    // We need to add the props for our new layer setup.
    this->AddLayerPropsToView();
  }

  // If there's no slot defined here, there's no menu for it yet, so
  // we have to make one.
  if( NULL != mFrameSlotMenus.GetPointer() ) {

    if( maSlotMenu.end() == maSlotMenu.find( inSlot ) )
      this->MakeMenuForSlot( inSlot );
    
    // Set the proper name in our menu.
    if( NULL != maSlotMenu[inSlot].GetPointer() ) {
      if ( NULL != maCol[inSlot].GetPointer() )
	maSlotMenu[inSlot]->SetValue( maCol[inSlot]->GetLabel() );
      else
	maSlotMenu[inSlot]->SetValue( "None" );
    }
  }

  // Redraw.
  this->GetRenderWindow()->Render();

  // If we have a z slider, adjust the range now, since our top layer
  // could have changed.
  this->AdjustViewPlaneSliders();
}

void
vtkKWScubaView::SetLayerCollectionAtSlotByComboBoxValue ( int inSlot,
						  const char* isValue ) {

  if ( maSlotMenu[inSlot].GetPointer() ) {

    // If we can find this string in the combo box...
    if ( maSlotMenu[inSlot]->HasValue( isValue ) ) {

      // Try to get an entry. If it's in our map...
      int nEntry = maSlotMenu[inSlot]->GetValueIndex( isValue );
      if ( maMenuIndexCollection.end() !=
           maMenuIndexCollection.find( nEntry ) ) {

        // Set the layer at this slot from the pointer our map.
        this->SetLayerCollectionAtSlot( inSlot,maMenuIndexCollection[nEntry] );
      } else {
        throw runtime_error( "Tool table doens't have that entry" );
      }
    } else {
      throw runtime_error( "Tool menu doens't have that entry" );
    }
  }
}

int
vtkKWScubaView::GetSlotOfLayerCollection ( const vtkKWScubaLayerCollection* iCol ) const {

  // Go through our collection map and look for this collection.
  SlotCollectionMapType::const_iterator tCol;
  for( tCol = maCol.begin(); tCol != maCol.end(); ++tCol )
    if ( tCol->second.GetPointer() == iCol )
      return tCol->first;
      
  return -1;
}

int
vtkKWScubaView::GetSlotOfLayer ( const vtkKWScubaLayer* iLayer ) const {

  // Go through our collection map, get the layer collection from each
  // one, and compare the current layer.
  SlotCollectionMapType::const_iterator tCol;
  for( tCol = maCol.begin(); tCol != maCol.end(); ++tCol )
    if ( NULL != tCol->second.GetPointer() )
      if ( tCol->second->GetCurrentLayer() == iLayer )
	return tCol->first;

  return -1;
}

void
vtkKWScubaView::MakeMenuForSlot ( int inSlot ) {

  if( NULL == mFrameSlotMenus )
    throw runtime_error( "Trying to make menus with no frame" );
  
  char sLabel[1024];
  char sCmd[1024];

  // Make a menu.
  vtkSmartPointer<vtkKWComboBoxWithLabel> labeledComboBox = 
    vtkSmartPointer<vtkKWComboBoxWithLabel>::New();
  labeledComboBox->SetParent( mFrameSlotMenus );
  labeledComboBox->Create();
  sprintf( sLabel, "Slot %d:", inSlot );
  labeledComboBox->SetLabelText( sLabel );
  
  // Assign the menu to this slot.
  maSlotMenu[inSlot] = labeledComboBox->GetWidget();
  sprintf( sCmd, "SetLayerCollectionAtSlotByComboBoxValue %d", inSlot );
  maSlotMenu[inSlot]->SetCommand( this, sCmd );
  
  this->Script( "pack %s -side top -fill x -anchor nw",
		labeledComboBox->GetWidgetName() );
}

void
vtkKWScubaView::LayerListChanged () {
  this->UpdateLayerComboBoxes();
}

void
vtkKWScubaView::ResetAllCameras () {

  vtkSmartPointer<vtkCamera> savedCamera = 
    this->GetRenderer()->GetActiveCamera();

  // Set each of our cameras in the render and let the renderer reset
  // the camera, sizing it to fit the new data.
  this->GetRenderer()->SetActiveCamera( m2DCamera );
  this->GetRenderer()->ResetCamera();

  this->GetRenderer()->SetActiveCamera( m3DCamera );
  this->GetRenderer()->ResetCamera();
  m3DCamera->Azimuth( 45 );
  m3DCamera->Elevation( 45 );

  // Go back to the original camera (now reset) and render.
  this->GetRenderer()->SetActiveCamera( savedCamera );
  this->GetRenderer()->Render();


}

void
vtkKWScubaView::AdjustViewPlaneSliders () {

  // Adjust our 2D slider if we have it.
  if ( mScale2DRASZ.GetPointer() ) {

    float min = 0, max = 0;
    this->Get2DRASZRange( min, max );

    mScale2DRASZ->SetRange( min, max );
    mScale2DRASZ->SetResolution( this->Get2DRASZIncrementHint() );
  }

  // Adjust our 3D sliders if we have them.
  if ( mScale3DRASX.GetPointer() ||
       mScale3DRASY.GetPointer() ||
       mScale3DRASZ.GetPointer() ) {

    float min[3], max[3];
    float hint[3];
    this->Get3DRASRange( min, max );
    this->Get3DRASIncrementHint( hint );

    if ( mScale3DRASX.GetPointer() ) {
      mScale3DRASX->SetRange( min[0], max[0] );
      mScale3DRASX->SetResolution( hint[0] );
    }
    if ( mScale3DRASY.GetPointer() ) {
      mScale3DRASY->SetRange( min[1], max[1] );
      mScale3DRASY->SetResolution( hint[1] );
    }
    if ( mScale3DRASZ.GetPointer() ) {
      mScale3DRASZ->SetRange( min[2], max[2] );
      mScale3DRASZ->SetResolution( hint[2] );
    }
  }
}

void
vtkKWScubaView::SetDisplayMode ( vtkKWScubaView::DisplayMode iMode ) {
  mDisplayMode = iMode;

  this->DisplayModeChanged();
}

void
vtkKWScubaView::SetDisplayModeTo2D () {
  this->SetDisplayMode( vtkKWScubaView::TwoDee );
}

void
vtkKWScubaView::SetDisplayModeTo3D () {
  this->SetDisplayMode( vtkKWScubaView::ThreeDee );
}

vtkKWScubaView::DisplayMode
vtkKWScubaView::GetDisplayMode () const {
  return mDisplayMode;
}

void
vtkKWScubaView::DisplayModeChanged () {

  // Pack the right controls.
  this->PackDisplayModeControls();

  // Switch out the cameras.
  if( TwoDee == mDisplayMode )
    this->GetRenderer()->SetActiveCamera( m2DCamera );
  else if( ThreeDee == mDisplayMode )
    this->GetRenderer()->SetActiveCamera( m3DCamera );

  // That changed our layers, so add the props.
  this->AddLayerPropsToView();

  // Update our buttons.
  if( TwoDee == mDisplayMode ) {

    if( mRadBtnDisplayMode2D.GetPointer() && 
	!mRadBtnDisplayMode2D->GetSelectedState() ) 
      mRadBtnDisplayMode2D->SelectedStateOn();

    if( mBtnRotateXPos.GetPointer() ) mBtnRotateXPos->SetStateToDisabled();
    if( mBtnRotateXNeg.GetPointer() ) mBtnRotateXNeg->SetStateToDisabled();
    if( mBtnRotateYPos.GetPointer() ) mBtnRotateYPos->SetStateToDisabled();
    if( mBtnRotateYNeg.GetPointer() ) mBtnRotateYNeg->SetStateToDisabled();
    if( mBtnRotateZPos.GetPointer() ) mBtnRotateZPos->SetStateToDisabled();
    if( mBtnRotateZNeg.GetPointer() ) mBtnRotateZNeg->SetStateToDisabled();

  } else if( ThreeDee == mDisplayMode ) {

    if( mRadBtnDisplayMode3D.GetPointer() && 
	!mRadBtnDisplayMode3D->GetSelectedState() ) 
      mRadBtnDisplayMode3D->SelectedStateOn();

    if( mBtnRotateXPos.GetPointer() ) mBtnRotateXPos->SetStateToNormal();
    if( mBtnRotateXNeg.GetPointer() ) mBtnRotateXNeg->SetStateToNormal();
    if( mBtnRotateYPos.GetPointer() ) mBtnRotateYPos->SetStateToNormal();
    if( mBtnRotateYNeg.GetPointer() ) mBtnRotateYNeg->SetStateToNormal();
    if( mBtnRotateZPos.GetPointer() ) mBtnRotateZPos->SetStateToNormal();
    if( mBtnRotateZNeg.GetPointer() ) mBtnRotateZNeg->SetStateToNormal();

  } 

  // Render the new view.
  this->GetRenderWindow()->Render();
}

void
vtkKWScubaView::StartFastMode () {

  mbFastMode = true;

  // Turn on fast mode on all layers.
  SlotCollectionMapType::const_iterator tCol;
  for ( tCol = maCol.begin(); tCol != maCol.end(); ++tCol )
    if ( NULL != tCol->second.GetPointer() ) {
      try { 
	tCol->second->SetFastMode( true );
      }
      catch( exception& e ) {
	cerr << e.what() << endl;
      }
    }
  
  // This will work on LOD actors.
  this->GetRenderWindow()->SetDesiredUpdateRate( 30.0 );
}

void
vtkKWScubaView::StopFastMode () {

  mbFastMode = false;

  // Turn off fast mode on all layers.
  SlotCollectionMapType::const_iterator tCol;
  for ( tCol = maCol.begin(); tCol != maCol.end(); ++tCol )
    if ( NULL != tCol->second.GetPointer() ) {
      try { 
	tCol->second->SetFastMode( false );
      }
      catch( exception& e ) {
	cerr << e.what() << endl;
      }
    }
  
  // This will work on LOD actors.
  this->GetRenderWindow()->SetDesiredUpdateRate( 0 );

  // Render the new view.
  this->GetRenderWindow()->Render();
}

bool 
vtkKWScubaView::GetFastMode () const {
  return mbFastMode;
}

bool
vtkKWScubaView::IsDataDirty () const {
  return mbDirty;
}

void
vtkKWScubaView::Set2DInPlane ( int iPlane ) {

  int oldInPlane = m2DInPlane;
  m2DInPlane = iPlane;

  // We want to switch orientations, but keep the same relative
  // position, distance, and focal point. So get those now relative to
  // the current orientation.
  double oldFocalPoint[3], oldDistance, oldPosition[3];
  vtkSmartPointer<vtkCamera> camera = 
    this->GetRenderer()->GetActiveCamera();
  camera->GetFocalPoint( oldFocalPoint );
  oldDistance = camera->GetDistance();
  camera->GetPosition( oldPosition );

  // We have 3D position and focal point, but we need to get the 2D
  // version, based on what we display on our screen.
  double oldRelativePosX = 0, oldRelativePosY = 0;
  double oldRelativeFocX = 0, oldRelativeFocY = 0;
  switch ( oldInPlane ) {
  case 0:
    oldRelativePosX = oldPosition[1];
    oldRelativePosY = oldPosition[2];
    oldRelativeFocX = oldFocalPoint[1];
    oldRelativeFocY = oldFocalPoint[2];
    break;
  case 1:
    oldRelativePosX = oldPosition[0];
    oldRelativePosY = oldPosition[2];
    oldRelativeFocX = oldFocalPoint[0];
    oldRelativeFocY = oldFocalPoint[2];
    break;
  case 2:
    oldRelativePosX = oldPosition[0];
    oldRelativePosY = oldPosition[1];
    oldRelativeFocX = oldFocalPoint[0];
    oldRelativeFocY = oldFocalPoint[1];
    break;
  }

  // Set the new camera position using the relative position and focal
  // point from before.
  switch ( m2DInPlane ) {
  case 0:
    camera->SetFocalPoint( 0, oldRelativeFocX, oldRelativeFocY );
    camera->SetViewUp( 0, 0, 1 );
    camera->SetPosition( oldDistance, oldRelativePosX, oldRelativePosY );
    break;
  case 1:
    camera->SetFocalPoint( oldRelativeFocX, 0, oldRelativeFocY );
    camera->SetViewUp( 0, 0, 1 );
    camera->SetPosition( oldRelativePosX, oldDistance, oldRelativePosY );
    break;
  case 2:
    camera->SetFocalPoint( oldRelativeFocX, oldRelativeFocY, 0 );
    camera->SetViewUp( 0, 1, 0 );
    camera->SetPosition( oldRelativePosX, oldRelativePosY, -oldDistance );
    break;
  }

  // Tell all the layers to set their inplane.
  SlotCollectionMapType::iterator tCol;
  for ( tCol = maCol.begin(); tCol != maCol.end(); ++tCol )
    if ( NULL != tCol->second ) {
      try {
	tCol->second->Set2DInPlane( m2DInPlane );
      }
      catch( exception& e  ) {
	cerr << e.what() << endl; 
      }
    }

  // Render the new view.
  if ( this->GetRenderWindow() )
    this->GetRenderWindow()->Render();

  // Adjust the sliders
  this->AdjustViewPlaneSliders();

  // The x/y/z inplane radio buttons.
  if( mRadBtnSet2DInPlane.GetPointer() ) {
    vtkSmartPointer<vtkKWRadioButton> btn =
      mRadBtnSet2DInPlane->GetWidget( this->Get2DInPlane() );
    if ( btn.GetPointer() ) {
      btn->SelectedStateOn();
    }
  }

  // Update our scale text.
  if( mScale2DRASZ.GetPointer() ) {
    switch ( iPlane ) {
    case 0:
      mScale2DRASZ->SetLabelText( "X RAS" );
      break;
    case 1:
      mScale2DRASZ->SetLabelText( "Y RAS" );
      break;
    case 2:
      mScale2DRASZ->SetLabelText( "Z RAS" );
      break;
    }
  }

  // Update the plane radio buttons.
  if ( mRadBtnSet2DInPlane.GetPointer() ) {
    vtkSmartPointer<vtkKWRadioButton> btn = 
      mRadBtnSet2DInPlane->GetWidget( iPlane );
    if ( btn.GetPointer() ) {
      btn->SelectedStateOn();
    }
  }
}

int
vtkKWScubaView::Get2DInPlane () const {
  return m2DInPlane;
}

void
vtkKWScubaView::Set2DRASZ ( float iRASZ ) {
  m2DRASZ = iRASZ;

  // Tell all the layers to set RASZ.
  SlotCollectionMapType::iterator tCol;
  for ( tCol = maCol.begin(); tCol != maCol.end(); ++tCol )
    if ( NULL != tCol->second ) {
      try {
	tCol->second->Set2DRASZ( m2DRASZ );
      }
      catch( exception& e  ) {
	cerr << e.what() << endl; 
      }
    }

  // Our info is now changed.
  this->InfoChanged();

  // Render the new view.
  if ( this->GetRenderWindow() )
    this->GetRenderWindow()->Render();

  // Update our scale.
  if ( mScale2DRASZ )
    mScale2DRASZ->SetValue( this->Get2DRASZ() );
}

float
vtkKWScubaView::Get2DRASZ () const {
  return m2DRASZ;
}

void
vtkKWScubaView::Get2DRASZRange ( float& oMin, float& oMax ) const {

  // Get the range from all layers. Find the least widest range of
  // all.
  float min = 99999;
  float max = -99999;
  SlotCollectionMapType::const_iterator tCol;
  for( tCol = maCol.begin(); tCol != maCol.end(); ++tCol )
    if( NULL != tCol->second ) {
      try {
	float bounds[6] = { 0,0,0,0,0,0 };
	vtkKWScubaLayer* layer = tCol->second->GetCurrentLayer();
	if( NULL != layer ) {
	  layer->GetRASBounds( bounds );
	  switch ( m2DInPlane ) {
	  case 0:
	    if ( bounds[0] < min ) min = bounds[0];
	    if ( bounds[1] > max ) max = bounds[1];
	    break;
	  case 1:
	    if ( bounds[2] < min ) min = bounds[2];
	    if ( bounds[3] > max ) max = bounds[3];
	    break;
	  case 2:
	    if ( bounds[4] < min ) min = bounds[4];
	    if ( bounds[5] > max ) max = bounds[5];
	    break;
	  }
	}
      }
      catch( exception& e ) {
	cerr << e.what() << endl;
      }
    }

  // Set the min and max if we got them, or else return 0-1.
  if ( min != 99999 )
    oMin = min;
  else
    oMin = 0;
  if ( max != -99999 )
    oMax = max;
  else
    max = 1;
}

float
vtkKWScubaView::Get2DRASZIncrementHint () const {

  // Try to get the increment for the highest layer, the one that's
  // being drawn on top. Return the dimension of our current inplane.
  float hint = 1.0;
  try {
    int nSlot = GetHighestFilledLayerSlot();
    SlotCollectionMapType::const_iterator tCol;
    if ( (tCol = maCol.find(nSlot)) != maCol.end() &&
	 NULL != tCol->second ) {

      vtkKWScubaLayer* layer = tCol->second->GetCurrentLayer();
      if( NULL != layer ) {

	float hints[3] = {0,0,0};
	layer->Get2DRASZIncrementHint( hints );
	
	if ( hints[m2DInPlane] > 0 && hints[m2DInPlane] < hint )
	  hint = hints[m2DInPlane];
      }
    }
  } 
  catch( exception& e ) {
    // No filled slots, so just use 1.0.
  }

  return hint;
}

void
vtkKWScubaView::Get3DRASRange ( float oMin[3], float oMax[3] ) const {

  // Get the range from all layers. Find the least widest range of
  // all.
  float min[3] = { numeric_limits<float>::max(), 
		   numeric_limits<float>::max(),
		   numeric_limits<float>::max() };
  float max[3] = { numeric_limits<float>::min(),
		   numeric_limits<float>::min(),
		   numeric_limits<float>::min() };
  SlotCollectionMapType::const_iterator tCol;
  for( tCol = maCol.begin(); tCol != maCol.end(); ++tCol )
    if( NULL != tCol->second ) {
      try {
	float bounds[6] = { 0,0,0,0,0,0 };
	vtkKWScubaLayer* layer = tCol->second->GetCurrentLayer();
	if( NULL != layer ) {
	  layer->GetRASBounds( bounds );
	  if ( bounds[0] < min[0] ) min[0] = bounds[0];
	  if ( bounds[1] > max[0] ) max[0] = bounds[1];
	  if ( bounds[2] < min[1] ) min[1] = bounds[2];
	  if ( bounds[3] > max[1] ) max[1] = bounds[3];
	  if ( bounds[4] < min[2] ) min[2] = bounds[4];
	  if ( bounds[5] > max[2] ) max[2] = bounds[5];
	}
      }
      catch( exception& e ) {
	cerr << e.what() << endl;
      }
    }

  // Set the min and max if we got them, or else return 0-1.
  for( int nOrientation = 0; nOrientation < 3; nOrientation++ ) {
    if ( min[nOrientation] != numeric_limits<float>::max() )
      oMin[nOrientation] = min[nOrientation];
    else
      oMin[nOrientation] = 0;
    if ( max[nOrientation] != numeric_limits<float>::min() )
      oMax[nOrientation] = max[nOrientation];
    else
      max[nOrientation] = 1;
  }
}

void
vtkKWScubaView::Get3DRASIncrementHint ( float oInc[3] ) const {

  // Try to get the increment for the highest layer, the one that's
  // being drawn on top.
  float bestHints[3] = { 1.0, 1.0, 1.0 };
  try {
    int nSlot = GetHighestFilledLayerSlot();
    SlotCollectionMapType::const_iterator tCol;
    if ( (tCol = maCol.find(nSlot)) != maCol.end() &&
	 NULL != tCol->second ) {

      vtkKWScubaLayer* layer = tCol->second->GetCurrentLayer();
      if( NULL != layer ) {

	float hints[3] = {0,0,0};
	layer->Get2DRASZIncrementHint( hints );
	
	for( int nOrientation = 0; nOrientation < 3; nOrientation++ )
	  if ( hints[nOrientation] > 0 && 
	       hints[nOrientation] < bestHints[nOrientation] )
	    bestHints[nOrientation] = hints[nOrientation];
      }
    }
  } 
  catch( exception& e ) {
    // No filled slots, so just use 1.0.
  }

  for( int nOrientation = 0; nOrientation < 3; nOrientation++ )
    oInc[nOrientation] = bestHints[nOrientation];
      
}

void
vtkKWScubaView::Set2DZoomLevel ( float iZoomlevel ) {

  if ( this->GetRenderWindow() ) {

    // This is implemented using the view angle. 30 is the default
    // starting view angle, so calc a new angle based on our zoom
    // level as a factor of that.
    double newAngle = 30.0 / iZoomlevel;

    vtkCamera *camera = this->GetRenderer()->GetActiveCamera();
    camera->SetViewAngle( newAngle );
    this->GetRenderWindowInteractor()->Render();

  }
}

float
vtkKWScubaView::Get2DZoomLevel () {

  if ( this->GetRenderWindow() ) {

    // This is implemented using the view angle. 30 is the default
    // starting view angle, so return our zoom level as a factor of
    // that.
    vtkCamera *camera = this->GetRenderer()->GetActiveCamera();
    return 30.0 / camera->GetViewAngle();

  } else {

    return 0;

  }
}

void
vtkKWScubaView::ZoomOut () {

  this->Set2DZoomLevel( this->Get2DZoomLevel() / 2.0 );
}

void
vtkKWScubaView::ZoomIn () {

  this->Set2DZoomLevel( this->Get2DZoomLevel() * 2.0 );
}

void 
vtkKWScubaView::AnimateCameraElevatePositive () {
  
  this->StartFastMode();
  for( int nStep = 0; nStep < 90; nStep++ ) {
    this->GetRenderer()->GetActiveCamera()->Elevation( 1.0 );
    this->GetRenderer()->GetActiveCamera()->OrthogonalizeViewUp();
    this->Render();
  }
  this->StopFastMode();
}

void 
vtkKWScubaView::AnimateCameraElevateNegative () {

  this->StartFastMode();
  for( int nStep = 0; nStep < 90; nStep++ ) {
    this->GetRenderer()->GetActiveCamera()->Elevation( -1.0 );
    this->GetRenderer()->GetActiveCamera()->OrthogonalizeViewUp();
    this->Render();
  }
  this->StopFastMode();
}

void 
vtkKWScubaView::AnimateCameraAzimuthNegative () {

  this->StartFastMode();
  for( int nStep = 0; nStep < 90; nStep++ ) {
    this->GetRenderer()->GetActiveCamera()->Azimuth( -1.0 );
    this->GetRenderer()->GetActiveCamera()->OrthogonalizeViewUp();
    this->Render();
  }
  this->StopFastMode();
}

void 
vtkKWScubaView::AnimateCameraAzimuthPositive () {

  this->StartFastMode();
  for( int nStep = 0; nStep < 90; nStep++ ) {
    this->GetRenderer()->GetActiveCamera()->Azimuth( 1.0 );
    this->GetRenderer()->GetActiveCamera()->OrthogonalizeViewUp();
    this->Render();
  }
  this->StopFastMode();
}

void 
vtkKWScubaView::AnimateCameraRollNegative () {

  this->StartFastMode();
  for( int nStep = 0; nStep < 90; nStep++ ) {
    this->GetRenderer()->GetActiveCamera()->Roll( -1.0 );
    this->Render();
  }
  this->StopFastMode();
}

void 
vtkKWScubaView::AnimateCameraRollPositive () {

  this->StartFastMode();
  for( int nStep = 0; nStep < 90; nStep++ ) {
    this->GetRenderer()->GetActiveCamera()->Roll( 1.0 );
    this->Render();
  }
  this->StopFastMode();
}

void
vtkKWScubaView::Set3DRASX ( float i3DRASX ) {
  m3DRASX = i3DRASX;

  // Tell all the layers to set 3DRASX.
  SlotCollectionMapType::iterator tCol;
  for ( tCol = maCol.begin(); tCol != maCol.end(); ++tCol )
    if ( NULL != tCol->second.GetPointer() ) {
      try {
	tCol->second->Set3DRASX( m3DRASX );
      }
      catch( exception& e  ) {
	cerr << e.what() << endl; 
      }
    }

  // Render the new view.
  if ( this->GetRenderWindow() )
    this->GetRenderWindow()->Render();

  // Update our scale.
  if ( mScale3DRASX.GetPointer() )
    mScale3DRASX->SetValue( this->Get3DRASX() );
}

float
vtkKWScubaView::Get3DRASX () const {
  return m3DRASX;
}

void
vtkKWScubaView::Set3DRASY ( float i3DRASY ) {
  m3DRASY = i3DRASY;

  // Tell all the layers to set 3DRASY.
  SlotCollectionMapType::iterator tCol;
  for ( tCol = maCol.begin(); tCol != maCol.end(); ++tCol )
    if ( NULL != tCol->second.GetPointer() ) {
      try {
	tCol->second->Set3DRASY( m3DRASY );
      }
      catch( exception& e  ) {
	cerr << e.what() << endl; 
      }
    }

  // Render the new view.
  if ( this->GetRenderWindow() )
    this->GetRenderWindow()->Render();

  // Update our scale.
  if ( mScale3DRASY.GetPointer() )
    mScale3DRASY->SetValue( this->Get3DRASY() );
}

float
vtkKWScubaView::Get3DRASY () const {
  return m3DRASY;
}

void
vtkKWScubaView::Set3DRASZ ( float i3DRASZ ) {
  m3DRASZ = i3DRASZ;

  // Tell all the layers to set 3DRASZ.
  SlotCollectionMapType::iterator tCol;
  for ( tCol = maCol.begin(); tCol != maCol.end(); ++tCol )
    if ( NULL != tCol->second.GetPointer() ) {
      try {
	tCol->second->Set3DRASZ( m3DRASZ );
      }
      catch( exception& e  ) {
	cerr << e.what() << endl; 
      }
    }

  // Render the new view.
  if ( this->GetRenderWindow() )
    this->GetRenderWindow()->Render();

  // Update our scale.
  if ( mScale3DRASZ.GetPointer() )
    mScale3DRASZ->SetValue( this->Get3DRASZ() );
}

float
vtkKWScubaView::Get3DRASZ () const {
  return m3DRASZ;
}

void
vtkKWScubaView::Get3DRAS ( float& o3DRASX, float& o3DRASY, float& o3DRASZ ) const {
  o3DRASX = m3DRASX;
  o3DRASY = m3DRASY;
  o3DRASZ = m3DRASZ;
}


void
vtkKWScubaView::MouseMoveEvent ( vtkKWScubaWindow* iWindow, int iCoords[2] ) {

  // Call the MouseMoveEvent function on the tool and get the world
  // coords back.
  double worldCoords[3] = { 0,0,0 };
  this->PickPointAndCallToolEvent( iWindow, iCoords,
                                   &vtkKWScubaTool::MouseMoveEvent,
                                   worldCoords );

  // Save the mouseover coords and note our info has changed. We
  // reconvert becase we could be on a different RASZ by now.
  float ras[3] = { 0,0,0 };
  this->ConvertWorldCoordsToRAS( worldCoords, ras );
  mMouseOverRASCoords[0] = ras[0];
  mMouseOverRASCoords[1] = ras[1];
  mMouseOverRASCoords[2] = ras[2];

  this->InfoChanged();
}

void
vtkKWScubaView::LeftButtonDownEvent ( vtkKWScubaWindow* iWindow,
                                      int iCoords[2] ) {

  // Call the LeftMouseDownEvent function on the tool.
  this->PickPointAndCallToolEvent( iWindow, iCoords,
                                   &vtkKWScubaTool::LeftMouseDownEvent, NULL );
}

void
vtkKWScubaView::LeftButtonUpEvent ( vtkKWScubaWindow* iWindow,
                                    int iCoords[2] ) {

  // Call the LeftMouseUpEvent function on the tool.
  this->PickPointAndCallToolEvent( iWindow, iCoords,
                                   &vtkKWScubaTool::LeftMouseUpEvent, NULL );
}

void
vtkKWScubaView::MiddleButtonDownEvent ( vtkKWScubaWindow* iWindow,
                                        int iCoords[2] ) {

  // Call the MiddleMouseDownEvent function on the tool.
  this->PickPointAndCallToolEvent( iWindow, iCoords,
                                   &vtkKWScubaTool::MiddleMouseDownEvent,NULL);
}

void
vtkKWScubaView::MiddleButtonUpEvent ( vtkKWScubaWindow* iWindow,
                                      int iCoords[2] ) {

  // Call the MiddleMouseUpEvent function on the tool.
  this->PickPointAndCallToolEvent( iWindow, iCoords,
                                   &vtkKWScubaTool::MiddleMouseUpEvent, NULL );
}

void
vtkKWScubaView::RightButtonDownEvent ( vtkKWScubaWindow* iWindow,
                                       int iCoords[2] ) {

  // Call the RightMouseDownEvent function on the tool.
  this->PickPointAndCallToolEvent( iWindow, iCoords,
                                   &vtkKWScubaTool::RightMouseDownEvent, NULL);
}

void
vtkKWScubaView::RightButtonUpEvent ( vtkKWScubaWindow* iWindow,
                                     int iCoords[2] ) {

  // Call the RightMouseUpEvent function on the tool.
  this->PickPointAndCallToolEvent( iWindow, iCoords,
                                   &vtkKWScubaTool::RightMouseUpEvent, NULL );
}

void
vtkKWScubaView::KeyDownEvent ( vtkKWScubaWindow* iWindow,
                               int iCoords[2] ) {

  // Call the KeyDownEvent function on the tool.
  this->PickPointAndCallToolEvent( iWindow, iCoords,
                                   &vtkKWScubaTool::KeyDownEvent, NULL );
}

void
vtkKWScubaView::KeyUpEvent ( vtkKWScubaWindow* iWindow,
                             int iCoords[2] ) {}

void
vtkKWScubaView::EnterEvent ( vtkKWScubaWindow* iWindow,
                             int iCoords[2] ) {

  // Call the EnterEvent function on the tool.
  this->PickPointAndCallToolEvent( iWindow, iCoords,
                                   &vtkKWScubaTool::EnterEvent, NULL );
}

void
vtkKWScubaView::LeaveEvent ( vtkKWScubaWindow* iWindow,
                             int iCoords[2] ) {

  // Call the LeaveEvent function on the tool.
  this->PickPointAndCallToolEvent( iWindow, iCoords,
                                   &vtkKWScubaTool::LeaveEvent, NULL );
}

void
vtkKWScubaView::PickPointAndCallToolEvent ( vtkKWScubaWindow* iWindow,
    int iCoords[2],
    void(vtkKWScubaTool::*iToolFunc)(vtkKWScubaWindow*,vtkKWScubaView*,vtkKWScubaLayer*,float[3]),
    double* ioWorldCoords ) {

  // Get the current tool.
  vtkKWScubaTool* tool = iWindow->GetCurrentTool();

  double worldCoords[3] = { 0,0,0 };
  float ras[3] = { 0,0,0 };
  vtkKWScubaLayer* highestLayer = NULL;

  // If the tool wants us to, use the picker to pick and find the
  // world coords of the point picked.
  if( !tool->SuspendPickEvents() ) {
    
    // Get this RWI's picker and downcast it to the cell picker we
    // created before.
    vtkAbstractPicker* abstractPicker =
      this->GetRenderWindowInteractor()->GetPicker();
    vtkCellPicker* picker = vtkCellPicker::SafeDownCast( abstractPicker );
    if ( NULL == picker ) {
    throw runtime_error( "Incorrect picker!" );
  }
    
  // Pick the coords.
  picker->Pick( iCoords[0], iCoords[1], 0, this->GetRenderer() );
  picker->GetPickPosition( worldCoords );
  
  // Convert those to RAS coords.
  this->ConvertWorldCoordsToRAS( worldCoords, ras );
  
  // Get the list of all actors hit by the picker. For each one, find
  // a corresponding layer, and find the slot number of that
  // layer. Search for the highest slot.
  int highestSlot = -1;
  vtkProp3DCollection* props = picker->GetProp3Ds();
  if ( props->GetNumberOfItems() > 0 ) {
    props->InitTraversal();
    vtkProp* prop;
    while ( (prop = props->GetNextProp()) ) {
      vtkKWScubaLayer* layer = vtkKWScubaLayer::GetLayerFromProp( prop );
      if ( layer ) {
        int slot = this->GetSlotOfLayer( layer );
        if ( slot > highestSlot ) {
          highestLayer = layer;
          highestSlot = slot;
        }
      }
    }
  }
  
  // If we had a layer, find out if its info changed.
  if ( NULL != highestLayer )
    if ( highestLayer->IsInfoChanged() )
    this->InfoChanged();
  
    // Return the world coords picked if they want.
    if ( NULL != ioWorldCoords ) {
      ioWorldCoords[0] = worldCoords[0];
      ioWorldCoords[1] = worldCoords[1];
      ioWorldCoords[2] = worldCoords[2];
    }
  }
  
  // Call the event function on the tool.
  (tool->*iToolFunc)( iWindow, this, highestLayer, ras );
}

void
vtkKWScubaView::ConvertWorldCoordsToRAS ( double const iWorldCoords[3],
					  float  oRAS[3] ) const {

  // Our world coords are similar to RAS coords, but in 2d mode we
  // just need to make sure the Z component is set to our 2D RASZ
  // value.
  oRAS[0] = iWorldCoords[0];
  oRAS[1] = iWorldCoords[1];
  oRAS[2] = iWorldCoords[2];

  if( mDisplayMode == TwoDee ) 
    oRAS[m2DInPlane] = m2DRASZ;
}

void
vtkKWScubaView::SetCursorFromWorld ( double const iWorldCoords[3] ) {

  // Convert the world coord to an RAS and set the cursor.
  float ras[3] = {0,0,0};
  this->ConvertWorldCoordsToRAS( iWorldCoords, ras );
  this->SetCursorFromRAS( ras );
}

void
vtkKWScubaView::SetCursorFromRAS ( float const iRASCoords[3] ) {

  // Set the cursor.
  mCursorRASCoords[0] = iRASCoords[0];
  mCursorRASCoords[1] = iRASCoords[1];
  mCursorRASCoords[2] = iRASCoords[2];

  // Need to redisplay our info.
  this->InfoChanged();
}

void
vtkKWScubaView::GetCursorInfoItems ( list<ScubaInfoItem>& ilInfo ) {

  // Get info items using our saved cursor location.
  this->GetInfoItemsAtRAS( mCursorRASCoords, ilInfo );
}

void
vtkKWScubaView::GetMouseOverInfoItems ( list<ScubaInfoItem>& ilInfo ) {

  // Get info items using the last mouseover location.
  this->GetInfoItemsAtRAS( mMouseOverRASCoords, ilInfo );
}

void
vtkKWScubaView::GetInfoItemsAtRAS ( float iRAS[3],
                                    list<ScubaInfoItem>& ilInfo ) {

  // Build our RAS info item.
  ScubaInfoItem info;
  info.SetLabel( "RAS" );
  char sRAS[1024];
  snprintf( sRAS, sizeof(sRAS), "%.2f %.2f %.2f", iRAS[0], iRAS[1], iRAS[2] );
  info.SetValue( sRAS );
  info.SetShortenHint( false );

  ilInfo.push_back( info );

  // Ask all of the layers for their info items at this point.
  SlotCollectionMapType::iterator tCol;
  for ( tCol = maCol.begin(); tCol != maCol.end(); ++tCol ) {
    if ( NULL != tCol->second.GetPointer() ) {
      try {
	vtkKWScubaLayer* layer = tCol->second->GetCurrentLayer();	
	if( NULL != layer )
	  layer->GetInfoItems( iRAS, ilInfo );
      }
      catch( exception& e ) {
	cerr << e.what() << endl;
      }
    }
  }
}

bool
vtkKWScubaView::IsInfoChanged () {
  return mbInfoChanged;
}

void
vtkKWScubaView::InfoUpdated () {

  // Clear our flag.
  mbInfoChanged = false;

  // Clear all the layers' flags.
  SlotCollectionMapType::iterator tCol;
  for ( tCol = maCol.begin(); tCol != maCol.end(); ++tCol ) {
    if ( NULL != tCol->second.GetPointer() ) {
      try {
	vtkKWScubaLayer* layer = tCol->second->GetCurrentLayer();
	if( NULL != layer )
	  layer->InfoUpdated();
      }
      catch( exception& e ) {
	cerr << e.what() << endl;
      }
    }
  }
}

void
vtkKWScubaView::InfoChanged () {
  mbInfoChanged = true;
}

void
vtkKWScubaView::AddLayerPropsToView () {

  // Because we can't sort props, and we want to display them with
  // the highest slots on top, we need to remove all our props and
  // add them back in the correct order.
  this->GetRenderer()->RemoveAllViewProps();
  
  // Get the props from the proper view of each layer collection.
  SlotCollectionMapType::iterator tCol;
  for ( tCol = maCol.begin(); tCol != maCol.end(); ++tCol ) {
    if ( NULL != tCol->second.GetPointer() ) {

      try {

	vtkKWScubaLayerCollection* col = tCol->second;

	// Set the current display mode in the collection.
	col->SetDisplayModeAndView( mDisplayMode, this );

	// Get the props from the current layer.
	vtkKWScubaLayer* layer = col->GetCurrentLayer();
	if( NULL != layer ) {
	  vtkPropCollection* props = layer->GetPropCollection();
	  
	  // Add all props to our renderer.
	  props->InitTraversal();
	  for ( vtkProp* prop = props->GetNextProp();
		NULL != prop;
		prop = props->GetNextProp() ) {
	    this->GetRenderer()->AddViewProp(prop);
	  }
	}
      }
      catch( exception& e ) {
	cerr << e.what() << endl;
      }
    }
  }
  
}

void
vtkKWScubaView::DoListenToMessage ( string const isMessage, 
				    void* const iData ) {

  // If the pipeline changed, just redraw.
  if ( isMessage == "PipelineChanged" ) {
    this->GetRenderWindow()->Render();
  }

  // If the layer label changed, pass the message along to the window.
  if ( isMessage == "LayerLabelChanged" ) {
    this->SendBroadcast( isMessage, iData );
  }

  // If a prop is added or removed, we need to rebuild the list.
  if( isMessage == "PropListChanged" ) {
    this->AddLayerPropsToView();
    this->GetRenderWindow()->Render();
  }
}

vtkKWScubaView*
vtkKWScubaView::GetViewFromRenderWindow( vtkRenderWindow* iWindow ) {

  // Static function to find a view from a render window, glueing the
  // vtk framework to our framework.
  vtkKWScubaView* view = mRenderWindowToViewMap[iWindow];
  if ( NULL == view ) {
    throw runtime_error( "Couldn't find view." );
  }
  return view;
}

void
vtkKWScubaView::PanBetweenWindowCoords ( int iDelta[2] ) {

  double viewFocus[4] = { 0,0,0,0 }, 
      viewPoint[3] = { 0,0,0 };
  double start[4] = { 0,0,0,0 }, 
    end[4] = { 0,0,0,0 },
      motionVector[3] = { 0,0,0 };

  // Start at a 0,0,0 point in world coords.
  vtkInteractorObserver::
    ComputeDisplayToWorld( this->GetRenderer(), 0, 0, 0, start );

  // Get a deltax,deltay,0 point in world coords.
  vtkInteractorObserver::
    ComputeDisplayToWorld( this->GetRenderer(),
			   (double)iDelta[0], (double)iDelta[1], 0, end );
  
  // Get the vector between the two points.
  motionVector[0] = start[0] - end[0];
  motionVector[1] = start[1] - end[1];
  motionVector[2] = start[2] - end[2];

  // Add the motion vector to the view focus.
  vtkCamera *camera = this->GetRenderer()->GetActiveCamera();
  camera->GetFocalPoint( viewFocus );
  camera->SetFocalPoint( motionVector[0] + viewFocus[0],
                         motionVector[1] + viewFocus[1],
                         motionVector[2] + viewFocus[2] );

  // Add the motion vector to the view point.
  camera->GetPosition( viewPoint );
  camera->SetPosition( motionVector[0] + viewPoint[0],
                       motionVector[1] + viewPoint[1],
                       motionVector[2] + viewPoint[2] );

  this->GetRenderer()->ResetCameraClippingRange();
  this->GetRenderer()->UpdateLightsGeometryToFollowCamera();

  this->Render();
}

void
vtkKWScubaView::ScrollBetweenWindowCoords ( int iDelta[2] ) {

  // Calcuate a new Z by multiplying the current delta with the
  // suggested Z increment.
  float newZ = this->Get2DRASZ() +
               (this->Get2DRASZIncrementHint() * (float)(iDelta[1]));

  // Get the current RAS range.
  float min = 0, max = 0;
  this->Get2DRASZRange( min, max );

  // If it's in range, set the new RASZ.
  if ( newZ >= min && newZ <= max ) {
    this->Set2DRASZ( newZ );
    this->Render();
  }
}

void
vtkKWScubaView::ZoomBetweenWindowCoords ( int iDelta[2] ) {

  // Get a y factor based on the size of the screen and the y delta.
  double *center = this->GetRenderer()->GetCenter();
  double dyf = 10.0 * (double)(iDelta[1]) / (double)(center[1]);
  double factor = pow((double)1.1, dyf);

  // Zoom by that amount.
  this->GetRenderer()->GetActiveCamera()->Zoom( factor );

  this->GetRenderer()->ResetCameraClippingRange();
  this->GetRenderer()->UpdateLightsGeometryToFollowCamera();

  this->Render();
}

void
vtkKWScubaView::RotateBetweenWindowCoords ( int iDelta[2] ) {

  int *size = this->GetRenderer()->GetRenderWindow()->GetSize();

  double delta_elevation = -20.0 / (double)size[1];
  double delta_azimuth = -20.0 / (double)size[0];
  
  double rxf = (double)iDelta[0] * delta_azimuth * 10.0;
  double ryf = (double)iDelta[1] * delta_elevation * 10.0;
  
  this->GetRenderer()->GetActiveCamera()->Azimuth(rxf);
  this->GetRenderer()->GetActiveCamera()->Elevation(ryf);
  this->GetRenderer()->GetActiveCamera()->OrthogonalizeViewUp();

  this->GetRenderer()->ResetCameraClippingRange();
  this->GetRenderer()->UpdateLightsGeometryToFollowCamera();

  this->Render();
}

void
vtkKWScubaView::PackDisplayModeControls () {
  
  if( NULL == mFrameDisplayModeControls )
    return;

  // Unpack our controls frame and pack the controls for our new
  // mode. These are already created in PopulateControlPage(), we're
  // just going to pack them here.
  mFrameDisplayModeControls->UnpackChildren();
  mFrameDisplayModeControls->SetHeight( 1 );

  if( TwoDee == mDisplayMode ) {
    
    if( mScale2DRASZ.GetPointer() && 
	mRadBtnSet2DInPlane.GetPointer() ) {
      this->Script( "pack %s %s -side top -fill x",
		    mScale2DRASZ->GetWidgetName(),
		    mRadBtnSet2DInPlane->GetWidgetName() );
    }

  } else if( ThreeDee == mDisplayMode ) {

    if( mScale3DRASX.GetPointer() && 
	mScale3DRASY.GetPointer() && 
	mScale3DRASZ.GetPointer() ) {
      this->Script( "pack %s %s %s -side top -fill x",
		    mScale3DRASX->GetWidgetName(),
		    mScale3DRASY->GetWidgetName(),
		    mScale3DRASZ->GetWidgetName() );
    }

  }
}

void
vtkKWScubaView::UpdateLayerComboBoxes () {

  // If we don't have any menus, bail.
  if( 0 == maSlotMenu.size() )
    return;

  // Clear the map of menu items to collections
  maMenuIndexCollection.clear();

  // Add "None" entry with a null layer pointer to all menus.
  SlotMenuMapType::iterator tMenu;
  for( tMenu = maSlotMenu.begin(); tMenu != maSlotMenu.end(); ++tMenu )
    tMenu->second->AddValue( "None" );

  // Ad a menu item entry for None, meaning that selecting the 0th
  // item in any menu will set that slot's layer collection to NULL.
  maMenuIndexCollection[0] = NULL;

  // Get a list of layer collections.
  list<vtkKWScubaLayerCollection*> lCols;
  vtkKWScubaLayerCollection::GetPointerList( lCols );
  list<vtkKWScubaLayerCollection*>::iterator tCol;

  // For each one...
  int nEntry = 1;
  for ( tCol = lCols.begin(); tCol != lCols.end(); ++tCol ) {
    vtkKWScubaLayerCollection* col = *tCol;

    // Add this entry to all menus boxes, and associate this entry with
    // the pointer to the collection.
    for( tMenu = maSlotMenu.begin(); tMenu != maSlotMenu.end(); ++tMenu )
      tMenu->second->AddValue( col->GetLabel() );
    maMenuIndexCollection[nEntry++] = col;
  }

  // Set the current layer names in the menus.
  for( tMenu = maSlotMenu.begin(); tMenu != maSlotMenu.end(); ++tMenu ) {

    // Get the collection for this slot.
    vtkKWScubaLayerCollection* col = maCol[tMenu->first];

    // Set it to the collection label if we have one, or just none.
    if ( NULL != col )
      tMenu->second->SetValue( col->GetLabel() );
    else
      tMenu->second->SetValue( "None" );
  }
}

