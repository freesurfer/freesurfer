/**
 * @file  vtkKWOrientMRIWindow.cxx
 * @brief Loads and works on data, handles UI commands
 *
 * Populates menus and toolbars with commands. Handles dialog boxes
 * for loading and saving data. Owns the data objects. Calculates the
 * new transform based on the camera orientation and writes it to the
 * volume.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: kteich $
 *    $Date: 2007/03/27 22:40:28 $
 *    $Revision: 1.8 $
 *
 * Copyright (C) 2002-2007,
 * The General Hospital Corporation (Boston, MA). 
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
 *
 */


#include <stdexcept>
#include "vtkKWOrientMRIWindow.h"

#include "IconLoader.h"
#include "vtkCamera.h"
#include "vtkFSVolumeSource.h"
#include "vtkKWApplication.h"
#include "vtkKWIcon.h"
#include "vtkKWLoadSaveDialog.h"
#include "vtkKWMenu.h"
#include "vtkKWMessageDialog.h"
#include "vtkKWPushButton.h"
#include "vtkKWToolbar.h"
#include "vtkKWToolbarSet.h"
#include "vtkLookupTable.h"
#include "vtkMatrix4x4.h"
#include "vtkObjectFactory.h"
#include "vtkRenderer.h"

using namespace std;

vtkStandardNewMacro( vtkKWOrientMRIWindow );
vtkCxxRevisionMacro( vtkKWOrientMRIWindow, "$Revision: 1.8 $" );

vtkKWOrientMRIWindow::vtkKWOrientMRIWindow () :
    vtkKWWindow(),
    mView( NULL ),
    mbDirty( false ),
    mVolume( NULL ),
    mLUT( NULL ),
    mOriginalVoxelToRASMatrix( NULL ),
    mOriginalView( NULL ),
    mOriginalViewI( NULL ) {
  for( int nCmd = 0; nCmd < kcCommands; nCmd++ ) {
    maCommandEnabled[nCmd] = false;
    maMenuItems[nCmd].menu = NULL;
    maMenuItems[nCmd].nItem = -1;
    maPushButtons[nCmd] = NULL;
  }
}

vtkKWOrientMRIWindow::~vtkKWOrientMRIWindow () {

  if( mView )
    mView->Delete();
  if( mVolume )
    mVolume->Delete();
  if( mLUT )
    mLUT->Delete();
  if( mOriginalVoxelToRASMatrix )
    mOriginalVoxelToRASMatrix->Delete();
  if( mOriginalView )
    mOriginalView->Delete();
  if( mOriginalViewI )
    mOriginalViewI->Delete();
}

void
vtkKWOrientMRIWindow::Create () {

  this->SupportHelpOn();

  this->Superclass::Create();

  this->SetPanelLayoutToSecondaryBelowView();
  this->SecondaryPanelVisibilityOff();
  this->MainPanelVisibilityOff();

  // Create our interior view.
  mView = vtkKWOrientMRIView::New();
  mView->SetParent( GetViewFrame() );
  mView->Create();
  this->Script( "pack %s -expand yes -fill both -anchor c",
                mView->GetWidgetName() );


  // Make a toolbar.
  vtkKWToolbar* toolbar = vtkKWToolbar::New();
  toolbar->SetName( "Main" );
  toolbar->SetParent( GetMainToolbarSet()->GetToolbarsFrame() );
  toolbar->Create();
  this->GetMainToolbarSet()->AddToolbar( toolbar );

  // Make our commands. For each one, we make a toolbar button,
  // associate the command with the button, and then make a menu
  // entry.
  int nFilePos = GetFileMenuInsertPosition();
  int nEditPos = 0;
  int nViewPos = GetViewMenuInsertPosition();

  memset( maPushButtons, 0, sizeof(vtkKWPushButton*) * kcCommands );
  memset( maMenuItems, 0, sizeof(MenuItem) * kcCommands );
  memset( maCommandEnabled, 0, sizeof(bool) * kcCommands );

  // Init our toolbar buttons with common settings.
  for ( int nCmd = 0; nCmd < kcCommands; nCmd++ ) {
    switch( nCmd ) {
    case CmdLoadVolume:
    case CmdSaveVolume:
    case CmdTransformVolume:
    case CmdRevertVolume:
    case CmdRestoreView:
    case CmdZoomOut:
    case CmdZoomIn:
    case CmdRotateXPos:
    case CmdRotateXNeg:
    case CmdRotateYPos:
    case CmdRotateYNeg:
    case CmdRotateZPos:
    case CmdRotateZNeg:
      maPushButtons[(Command)nCmd] = vtkKWPushButton::New();
      maPushButtons[(Command)nCmd]->SetParent( toolbar->GetFrame() );
      maPushButtons[(Command)nCmd]->Create();
      maPushButtons[(Command)nCmd]->SetWidth( 30 );
      maPushButtons[(Command)nCmd]->SetHeight( 30 );
      maPushButtons[(Command)nCmd]->
	SetImageToPredefinedIcon( vtkKWIcon::IconFileOpen ); // Placeholder.
      toolbar->AddWidget( maPushButtons[(Command)nCmd] );
      maPushButtons[(Command)nCmd]->Delete();
    default:
      break;
    }
  }

  // Load Volume
  maPushButtons[CmdLoadVolume]->SetText( "Load Volume" );
  maPushButtons[CmdLoadVolume]->SetBalloonHelpString( "Load Volume" );
  maPushButtons[CmdLoadVolume]->
    SetImageToPredefinedIcon( vtkKWIcon::IconFileOpen );
  maPushButtons[CmdLoadVolume]->
    SetImageToPredefinedIcon( vtkKWIcon::IconFileOpen );
  maPushButtons[CmdLoadVolume]->SetCommand( this, "LoadVolumeFromDlog" );
  try {
    IconLoader::SetPushButtonIcon( "LoadVolume", 
				   maPushButtons[CmdLoadVolume] );
  } catch (...) {}

  // Save Volume
  maPushButtons[CmdSaveVolume]->SetText( "Save Volume" );
  maPushButtons[CmdSaveVolume]->SetBalloonHelpString( "Save Volume" );
  maPushButtons[CmdSaveVolume]->
    SetImageToPredefinedIcon( vtkKWIcon::IconFloppy );
  maPushButtons[CmdSaveVolume]->SetCommand( this, "SaveVolumeWithConfirm" );
  try {
    IconLoader::SetPushButtonIcon( "SaveVolume", 
				   maPushButtons[CmdSaveVolume] );
  } catch (...) {}

  // Transform Volume
  maPushButtons[CmdTransformVolume]->SetText( "Transform Volume" );
  maPushButtons[CmdTransformVolume]->SetBalloonHelpString( "Transform Volume");
  maPushButtons[CmdTransformVolume]->
    SetImageToPredefinedIcon( vtkKWIcon::IconGridLinear );
  maPushButtons[CmdTransformVolume]->SetCommand( this, "TransformVolume" );
  try {
    IconLoader::SetPushButtonIcon( "TransformVolume",
				   maPushButtons[CmdTransformVolume] );
  } catch (...) {}

  // Revert.
  maPushButtons[CmdRevertVolume]->SetText( "Revert Volume" );
  maPushButtons[CmdRevertVolume]->SetBalloonHelpString( "Revert Volume" );
  maPushButtons[CmdRevertVolume]->
    SetImageToPredefinedIcon( vtkKWIcon::IconBrowserBack );
  maPushButtons[CmdRevertVolume]->SetCommand( this, "RevertToSavedTransform" );
  try {
    IconLoader::SetPushButtonIcon( "RevertVolume",
				   maPushButtons[CmdRevertVolume] );
  } catch (...) {}

  // Restore View
  maPushButtons[CmdRestoreView]->SetText( "Restore View" );
  maPushButtons[CmdRestoreView]->SetBalloonHelpString( "Restore View" );
  maPushButtons[CmdRestoreView]->
    SetImageToPredefinedIcon( vtkKWIcon::IconReload );
  maPushButtons[CmdRestoreView]->SetCommand( this, "RestoreView" );
  try {
    IconLoader::SetPushButtonIcon( "RestoreView", 
				   maPushButtons[CmdRestoreView] );
  } catch (...) {}

  // Zoom Out
  maPushButtons[CmdZoomOut]->SetText( "Zoom Out" );
  maPushButtons[CmdZoomOut]->SetBalloonHelpString( "Zoom Out" );
  maPushButtons[CmdZoomOut]->
    SetImageToPredefinedIcon( vtkKWIcon::IconMagGlass );
  maPushButtons[CmdZoomOut]->SetCommand( this, "ZoomOut" );
  try {
    IconLoader::SetPushButtonIcon( "ZoomOut", maPushButtons[CmdZoomOut] );
  } catch (...) {}

  // Zoom In
  maPushButtons[CmdZoomIn]->SetText( "Zoom In" );
  maPushButtons[CmdZoomIn]->SetBalloonHelpString( "Zoom In" );
  maPushButtons[CmdZoomIn]->SetCommand( this, "ZoomIn" );
  try {
    IconLoader::SetPushButtonIcon( "ZoomIn", maPushButtons[CmdZoomIn] );
  } catch (...) {}

  // Rotate X Pos
  maPushButtons[CmdRotateXPos]->SetText( "Rotate" );
  maPushButtons[CmdRotateXPos]->SetBalloonHelpString( "Rotate" );
  maPushButtons[CmdRotateXPos]->
    SetCommand( mView, "AnimateCameraElevatePositive" );
  try {
    IconLoader::SetPushButtonIcon( "RotateXPos", maPushButtons[CmdRotateXPos]);
  } catch (...) {}

  // Rotate X Neg
  maPushButtons[CmdRotateXNeg]->SetText( "Rotate" );
  maPushButtons[CmdRotateXNeg]->SetBalloonHelpString( "Rotate" );
  maPushButtons[CmdRotateXNeg]->
    SetCommand( mView, "AnimateCameraElevateNegative" );
  try {
    IconLoader::SetPushButtonIcon( "RotateXNeg", maPushButtons[CmdRotateXNeg]);
  } catch (...) {}

  // Rotate Y Pos
  maPushButtons[CmdRotateYPos]->SetText( "Rotate" );
  maPushButtons[CmdRotateYPos]->SetBalloonHelpString( "Rotate" );
  maPushButtons[CmdRotateYPos]->
    SetCommand( mView, "AnimateCameraAzimuthNegative" );
  try {
    IconLoader::SetPushButtonIcon( "RotateYPos", maPushButtons[CmdRotateYPos]);
  } catch (...) {}

  // Rotate Y Neg
  maPushButtons[CmdRotateYNeg]->SetText( "Rotate" );
  maPushButtons[CmdRotateYNeg]->SetBalloonHelpString( "Rotate" );
  maPushButtons[CmdRotateYNeg]->
    SetCommand( mView, "AnimateCameraAzimuthPositive" );
  try {
    IconLoader::SetPushButtonIcon( "RotateYNeg", maPushButtons[CmdRotateYNeg]);
  } catch (...) {}

  // Rotate Z Pos
  maPushButtons[CmdRotateZPos]->SetText( "Rotate" );
  maPushButtons[CmdRotateZPos]->SetBalloonHelpString( "Rotate" );
  maPushButtons[CmdRotateZPos]->
    SetCommand( mView, "AnimateCameraRollNegative" );
  try {
    IconLoader::SetPushButtonIcon( "RotateZPos", maPushButtons[CmdRotateZPos]);
  } catch (...) {}

  // Rotate Z Neg
  maPushButtons[CmdRotateZNeg]->SetText( "Rotate" );
  maPushButtons[CmdRotateZNeg]->SetBalloonHelpString( "Rotate" );
  maPushButtons[CmdRotateZNeg]->
    SetCommand( mView, "AnimateCameraRollPositive" );
  try {
    IconLoader::SetPushButtonIcon( "RotateZNeg", maPushButtons[CmdRotateZNeg]);
  } catch (...) {}

  // Build the menus. File menu.
  // Load Volume.
  this->GetFileMenu()->
    InsertCommand( nFilePos, "L&oad Volume...", this, "LoadVolumeFromDlog" );
  this->GetFileMenu()->SetItemCompoundModeToLeft( nFilePos );
  this->GetFileMenu()->SetItemAccelerator( nFilePos, "Ctrl+O" );
  try {
    IconLoader::SetMenuItemIcon( "LoadVolume", GetFileMenu(), nFilePos );
  } catch (...) {}
  maMenuItems[CmdLoadVolume].menu = GetFileMenu();
  maMenuItems[CmdLoadVolume].nItem = nFilePos;
  nFilePos++;

  this->GetFileMenu()->InsertSeparator( nFilePos++ );

  // Save Volume.
  this->GetFileMenu()->
    InsertCommand( nFilePos, "&Save Volume", this, "SaveVolumeWithConfirm" );
  this->GetFileMenu()->SetItemCompoundModeToLeft( nFilePos );
  this->GetFileMenu()->SetItemAccelerator( nFilePos, "Ctrl+S" );
  try {
    IconLoader::SetMenuItemIcon( "SaveVolume", GetFileMenu(), nFilePos );
  } catch (...) {}
  maMenuItems[CmdSaveVolume].menu = GetFileMenu();
  maMenuItems[CmdSaveVolume].nItem = nFilePos;
  nFilePos++;

  // Save Volume As.
  this->GetFileMenu()->
    InsertCommand( nFilePos, "Save Volume As...", this,"SaveVolumeAsFromDlog");
  this->GetFileMenu()->SetItemCompoundModeToLeft( nFilePos );
  try {
    IconLoader::SetMenuItemIcon( "SaveVolume", GetFileMenu(), nFilePos );
  } catch (...) {}
  maMenuItems[CmdSaveVolumeAs].menu = GetFileMenu();
  maMenuItems[CmdSaveVolumeAs].nItem = nFilePos;
  nFilePos++;

  this->GetFileMenu()->InsertSeparator( nFilePos++ );

  this->InsertRecentFilesMenu( nFilePos++, this );

  this->GetFileMenu()->InsertSeparator( nFilePos++ );

  // Edit menu.
  // Transform Volume.
  this->GetEditMenu()->
    InsertCommand( nEditPos, "Trans&form Volume", this, "TransformVolume" );
  this->GetEditMenu()->SetItemCompoundModeToLeft( nEditPos );
  this->GetEditMenu()->SetItemAccelerator( nEditPos, "Ctrl+F" );
  try {
    IconLoader::SetMenuItemIcon( "TransformVolume", GetEditMenu(), nEditPos );
  } catch (...) {}
  maMenuItems[CmdTransformVolume].menu = GetEditMenu();
  maMenuItems[CmdTransformVolume].nItem = nEditPos;
  nEditPos++;

  this->GetEditMenu()->InsertSeparator( nEditPos++ );

  // Revert.
  this->GetEditMenu()->
    InsertCommand( nEditPos, "&Revert Volume", this, "RevertToSavedTransform");
  this->GetEditMenu()->SetItemCompoundModeToLeft( nEditPos );
  this->GetEditMenu()->SetItemAccelerator( nEditPos, "Ctrl+R" );
  try {
    IconLoader::SetMenuItemIcon( "RevertVolume", GetEditMenu(), nEditPos );
  } catch (...) {}
  maMenuItems[CmdRevertVolume].menu = GetEditMenu();
  maMenuItems[CmdRevertVolume].nItem = nEditPos;
  nEditPos++;

  // View menu.
  // Restore view.
  this->GetViewMenu()->
    InsertCommand( nViewPos, "Restore &View", this, "RestoreView");
  this->GetViewMenu()->SetItemCompoundModeToLeft( nViewPos );
  this->GetViewMenu()->SetItemAccelerator( nViewPos, "Ctrl+V" );
  try {
    IconLoader::SetMenuItemIcon( "RestoreView", GetViewMenu(), nViewPos );
  } catch (...) {}
  maMenuItems[CmdRestoreView].menu = GetViewMenu();
  maMenuItems[CmdRestoreView].nItem = nViewPos;
  nViewPos++;

  this->GetViewMenu()->InsertSeparator( nViewPos++ );

  // Zoom Out.
  this->GetViewMenu()->
    InsertCommand( nViewPos, "Zoom Out", this, "ZoomOut");
  this->GetViewMenu()->SetItemCompoundModeToLeft( nViewPos );
  //  GetViewMenu()->SetItemAccelerator( nViewPos, "Ctrl+Minus" );
  try {
    IconLoader::SetMenuItemIcon( "ZoomOut", GetViewMenu(), nViewPos );
  } catch (...) {}
  maMenuItems[CmdZoomOut].menu = GetViewMenu();
  maMenuItems[CmdZoomOut].nItem = nViewPos;
  nViewPos++;

  // Zoom In.
  this->GetViewMenu()->
    InsertCommand( nViewPos, "Zoom In", this, "ZoomIn");
  this->GetViewMenu()->SetItemCompoundModeToLeft( nViewPos );
  //  GetViewMenu()->SetItemAccelerator( nViewPos, "Ctrl+Plus" );
  try {
    IconLoader::SetMenuItemIcon( "ZoomIn", GetViewMenu(), nViewPos );
  } catch (...) {}
  maMenuItems[CmdZoomIn].menu = GetViewMenu();
  maMenuItems[CmdZoomIn].nItem = nViewPos;
  nViewPos++;


  // Update our menu and buttons.
  this->UpdateCommandStatus();
}

void
vtkKWOrientMRIWindow::LoadVolumeFromDlog () {

  vtkKWLoadSaveDialog* dialog = vtkKWLoadSaveDialog::New();
  dialog->SetApplication( GetApplication() );
  dialog->Create();
  dialog->SetFileTypes( "{MGH {.mgh .mgz}} {Binary {.bshort .bfloat}} {All {*}}" );
  dialog->RetrieveLastPathFromRegistry( "LoadVolume" );
  dialog->SetDefaultExtension( ".mgz" );
  if ( dialog->Invoke() ) {
    dialog->SaveLastPathToRegistry( "LoadVolume" );
    string fnVolume( dialog->GetFileName() );
    this->LoadVolume( fnVolume.c_str() );
  }
}

void
vtkKWOrientMRIWindow::LoadVolume ( const char* ifnVolume ) {

  if( !mView ) throw runtime_error( "mView was NULL" );

  try {
    
    // Try to load the volume.
    vtkFSVolumeSource* volume = vtkFSVolumeSource::New();
    volume->MRIRead( ifnVolume );
    volume->Update();
      
    // Delete existing one and save a reference.
    if ( NULL != mVolume )
      mVolume->Delete();
    mVolume = volume;
    
    // Set it in the view.
    mView->SetCurrentVolume( mVolume );
    
    // Create a basic LUT if we don't have one yet.
    if( !mLUT ) 
      mLUT = vtkLookupTable::New();
    mLUT->SetTableRange( mVolume->GetMinValue(), mVolume->GetMaxValue() );
    mLUT->SetSaturationRange( 0, 0 );
    mLUT->SetHueRange( 0, 0 );
    mLUT->SetValueRange( 0, 1 );
    mLUT->Build();
    for ( int nEntry = 0; nEntry < mLUT->GetIndex(10); nEntry++ )
      mLUT->SetTableValue( nEntry, 0, 0, 0, 0 );

    // Set it in the view.
    mView->SetCurrentVolumeColors( mLUT );

    // Calculate the inverse of our starting viewing transform here so
    // can use it later.
    if ( NULL == mOriginalView )
      mOriginalView = vtkMatrix4x4::New();
    mOriginalView->
      DeepCopy( mView->GetRenderer()->GetActiveCamera()->GetViewTransformMatrix() );
    (*mOriginalView)[0][3] = 0;
    (*mOriginalView)[1][3] = 0;
    (*mOriginalView)[2][3] = 0;
    
    if ( NULL == mOriginalViewI )
      mOriginalViewI = vtkMatrix4x4::New();
    mOriginalViewI->DeepCopy( mOriginalViewI );
    mOriginalViewI->Invert();
    (*mOriginalViewI)[0][3] = 0;
    (*mOriginalViewI)[1][3] = 0;
    (*mOriginalViewI)[2][3] = 0;
    
    // Get and save the original VoxelToRAS.
    if ( NULL == mOriginalVoxelToRASMatrix )
      mOriginalVoxelToRASMatrix = vtkMatrix4x4::New();
    mOriginalVoxelToRASMatrix->DeepCopy( mVolume->GetVoxelToRASMatrix() );
    
    this->SetStatusText( "Volume loaded." );
    this->AddRecentFile( ifnVolume, this, "LoadVolume" );
    
  } catch ( exception& e ) {
    this->GetApplication()->ErrorMessage( e.what() );
  }

  // Update our menu and buttons.
  this->UpdateCommandStatus();
}

void
vtkKWOrientMRIWindow::SaveVolumeWithConfirm () {

  if ( mbDirty ) {
    if ( vtkKWMessageDialog::PopupYesNo
         ( this->GetApplication(), this,
           "Save Volume",
           "Are you sure you want to save changes?" ) ) {
      try {
        this->SaveVolume();
        this->SetStatusText( "Volume saved." );
      } catch ( exception& e ) {
	this->GetApplication()->ErrorMessage( e.what() );
      }
    }
  }

  // Update our menu and buttons.
  this->UpdateCommandStatus();
}

void
vtkKWOrientMRIWindow::SaveVolumeAsFromDlog () {

  vtkKWLoadSaveDialog* dialog = vtkKWLoadSaveDialog::New();
  dialog->SaveDialogOn();
  dialog->SetApplication( GetApplication() );
  dialog->Create();
  dialog->SetFileTypes( "{MGH {.mgh .mgz}} {Binary {.bshort .bfloat}} {All {*}}" );
  dialog->RetrieveLastPathFromRegistry( "LoadVolume" );
  dialog->SetDefaultExtension( ".mgz" );
  if ( dialog->Invoke() ) {
    dialog->SaveLastPathToRegistry( "LoadVolume" );
    string fnVolume( dialog->GetFileName() );
    this->SaveVolume( fnVolume.c_str() );
  }
}


void
vtkKWOrientMRIWindow::SaveVolume () {

  if( !mVolume ) throw runtime_error( "mVolume was NULL" );

  // Tell the volume to write itself.
  mVolume->MRIWrite();

  // Get the new VoxelToRAS as our original.
  if ( NULL == mOriginalVoxelToRASMatrix )
    mOriginalVoxelToRASMatrix = vtkMatrix4x4::New();
  mOriginalVoxelToRASMatrix->DeepCopy( mVolume->GetVoxelToRASMatrix() );

  // No longer dirty.
  mbDirty = false;
}

void
vtkKWOrientMRIWindow::SaveVolume ( const char* ifnVolume ) {

  if( !mVolume ) throw runtime_error( "mVolume was NULL" );

  // Set the volume filename.
  mVolume->MRIWrite( ifnVolume );

  // Get the new VoxelToRAS as our original.
  if ( NULL == mOriginalVoxelToRASMatrix )
    mOriginalVoxelToRASMatrix = vtkMatrix4x4::New();
  mOriginalVoxelToRASMatrix->DeepCopy( mVolume->GetVoxelToRASMatrix() );

  // No longer dirty.
  mbDirty = false;
}

void
vtkKWOrientMRIWindow::RevertToSavedTransform () {

  if( !mView ) throw runtime_error( "mView was NULL" );
  if( !mVolume ) throw runtime_error( "mVolume was NULL" );
  if( !mOriginalVoxelToRASMatrix ) 
    throw runtime_error( "mOriginalVoxelToRASMatrix was NULL" );

  try {
    
    // Set our save original matrix.
    mVolume->SetVoxelToRASMatrix( *mOriginalVoxelToRASMatrix );
      
    // Notify the view that the matrix has changed.
    mView->VolumeToRASTransformChanged();

    // Restore the view.
    mView->RestoreView();

    // No longer dirty.
    mbDirty = false;

  } catch ( exception& e ) {
    this->GetApplication()->ErrorMessage( e.what() );
  }
  
  // Update our menu and buttons.
  this->UpdateCommandStatus();
}


void
vtkKWOrientMRIWindow::TransformVolume () {

  if( !mVolume ) throw runtime_error( "mVolume was NULL" );
  if( !mView ) throw runtime_error( "mView was NULL" );
  if( !mOriginalVoxelToRASMatrix )
    throw runtime_error( "mOriginalVoxelToRASMatrix was NULL" );
  if( !mOriginalView ) throw runtime_error( "mOriginalView was NULL" );
  if( !mOriginalViewI ) throw runtime_error( "mOriginalViewI was NULL" );

  // Get the current VoxelToRAS matrix
  vtkMatrix4x4* currentVoxelToRAS = vtkMatrix4x4::New();
  currentVoxelToRAS->DeepCopy( mVolume->GetVoxelToRASMatrix() );

  // This is the current view transform without the translation
  // component.
  vtkMatrix4x4* viewTransform = vtkMatrix4x4::New();
  viewTransform->
    DeepCopy( mView->GetRenderer()->GetActiveCamera()->GetViewTransformMatrix() );
  (*viewTransform)[0][3] = 0;
  (*viewTransform)[1][3] = 0;
  (*viewTransform)[2][3] = 0;

  // Currently it includes the original camera transform. Since we
  // point the cmaera at the front of the face in RAS space, it
  // introduces a non-identity view transform as a basis, so we undo
  // that here to get the view transforms in relative terms.
  // viewTransform = mOriginalViewI * viewTransform
  vtkMatrix4x4::Multiply4x4( viewTransform, mOriginalViewI, viewTransform );


  // To get voxelToRAS we need the view transform in RAS space, which
  // we get by transforming our view transform to abosulte RAS space,
  // and applying it to the voxelToRAS transform.
  // newVoxelToRAS = currentVoxelToRAS * original-1 * viewTransform * original
  vtkMatrix4x4* newVoxelToRAS = vtkMatrix4x4::New();
  vtkMatrix4x4::Multiply4x4( mOriginalView, viewTransform, newVoxelToRAS );
  vtkMatrix4x4::Multiply4x4( newVoxelToRAS, mOriginalViewI, newVoxelToRAS );
  vtkMatrix4x4::Multiply4x4( newVoxelToRAS, currentVoxelToRAS, newVoxelToRAS );

  // Set the matrix in the volume.
  mVolume->SetVoxelToRASMatrix( *newVoxelToRAS );

  viewTransform->Delete();
  currentVoxelToRAS->Delete();
  newVoxelToRAS->Delete();

  // Recalc our reslice transform.
  mView->VolumeToRASTransformChanged();

  // Restore the view.
  mView->RestoreView();

  // Now dirty.
  mbDirty = true;

  // Update our menu and buttons.
  this->UpdateCommandStatus();
}


void
vtkKWOrientMRIWindow::RestoreView () {

  if( !mView ) throw runtime_error( "mView was NULL" );

  // Tell the view to restore.
  mView->RestoreView();

  // Update our menu and buttons.
  this->UpdateCommandStatus();
}


void
vtkKWOrientMRIWindow::ZoomBy ( float iFactor ) {

  if( !mView ) throw runtime_error( "mView was NULL" );

  // Tell the view to zoom.
  mView->ZoomBy( iFactor );

  // Update our menu and buttons.
  UpdateCommandStatus();
}

void
vtkKWOrientMRIWindow::ZoomIn () {

  this->ZoomBy( 2.0 );
}

void
vtkKWOrientMRIWindow::ZoomOut () {

  this->ZoomBy( 0.5 );
}

void
vtkKWOrientMRIWindow::UpdateCommandStatus () {

  // Determine the enabled state of our commands.
  maCommandEnabled[CmdLoadVolume] = true;
  maCommandEnabled[CmdSaveVolume] = (mView && mbDirty);
  maCommandEnabled[CmdSaveVolumeAs] = (mView && mbDirty);
  maCommandEnabled[CmdTransformVolume] = (mView && mVolume);
  maCommandEnabled[CmdRevertVolume] = (mView && mbDirty);
  maCommandEnabled[CmdRestoreView] = true;
  maCommandEnabled[CmdZoomIn] = true;
  maCommandEnabled[CmdZoomOut] = true;
  maCommandEnabled[CmdRotateXPos] = (mView && mVolume);
  maCommandEnabled[CmdRotateXNeg] = (mView && mVolume);
  maCommandEnabled[CmdRotateYPos] = (mView && mVolume);
  maCommandEnabled[CmdRotateYNeg] = (mView && mVolume);
  maCommandEnabled[CmdRotateZPos] = (mView && mVolume);
  maCommandEnabled[CmdRotateZNeg] = (mView && mVolume);

  // Set the state in the menus and buttons.
  for ( int nCmd = 0; nCmd < kcCommands; nCmd++ ) {
    if ( maCommandEnabled[(Command)nCmd] ) {
      if ( NULL != maPushButtons[(Command)nCmd] )
        maPushButtons[(Command)nCmd]->SetStateToNormal();
      if ( NULL != maMenuItems[(Command)nCmd].menu )
        maMenuItems[(Command)nCmd].menu->
        SetItemStateToNormal( maMenuItems[(Command)nCmd].nItem );
    } else {
      if ( NULL != maPushButtons[(Command)nCmd] )
        maPushButtons[(Command)nCmd]->SetStateToDisabled();
      if ( NULL != maMenuItems[(Command)nCmd].menu )
        maMenuItems[(Command)nCmd].menu->
        SetItemStateToDisabled( maMenuItems[(Command)nCmd].nItem );
    }
  }
}
