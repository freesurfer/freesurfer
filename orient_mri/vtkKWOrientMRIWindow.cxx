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
 *    $Date: 2007/04/19 21:50:40 $
 *    $Revision: 1.10 $
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
#include "vtkFreesurferLookupTable.h"
#include "vtkFSVolumeSource.h"
#include "vtkKWApplication.h"
#include "vtkKWIcon.h"
#include "vtkKWLoadSaveDialog.h"
#include "vtkKWMenu.h"
#include "vtkKWMessageDialog.h"
#include "vtkKWPushButton.h"
#include "vtkKWRadioButton.h"
#include "vtkKWRadioButtonSet.h"
#include "vtkKWToolbar.h"
#include "vtkKWToolbarSet.h"
#include "vtkLookupTable.h"
#include "vtkMatrix4x4.h"
#include "vtkObjectFactory.h"
#include "vtkRenderer.h"

using namespace std;

vtkStandardNewMacro( vtkKWOrientMRIWindow );
vtkCxxRevisionMacro( vtkKWOrientMRIWindow, "$Revision: 1.10 $" );

vtkKWOrientMRIWindow::vtkKWOrientMRIWindow () :
    vtkKWWindow(),
    mView( NULL ),
    mbDirty( false ),
    mVolume( NULL ),
    mGrayScaleColors( NULL ),
    mLUTColors( NULL ),
    mOriginalVoxelToRASMatrix( NULL ),
    mOriginalView( NULL ),
    mOriginalViewI( NULL ) {
}

vtkKWOrientMRIWindow::~vtkKWOrientMRIWindow () {

  if( mView )
    mView->Delete();
  if( mVolume )
    mVolume->Delete();
  if( mGrayScaleColors )
    mGrayScaleColors->Delete();
  if( mLUTColors )
    mLUTColors->Delete();
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

  // Populate the toolbar.

  // Load Volume
  mBtnLoadVolume = vtkKWPushButton::New();
  mBtnLoadVolume->SetParent( toolbar->GetFrame() );
  mBtnLoadVolume->Create();
  toolbar->AddWidget( mBtnLoadVolume );
  mBtnLoadVolume->SetText( "Load Volume" );
  mBtnLoadVolume->SetBalloonHelpString( "Load Volume" );
  mBtnLoadVolume->SetCommand( this, "LoadVolumeFromDlog" );
  try { IconLoader::SetPushButtonIcon( "LoadVolume", mBtnLoadVolume ); }
  catch (...) {}

  // Save Volume
  mBtnSaveVolume = vtkKWPushButton::New();
  mBtnSaveVolume->SetParent( toolbar->GetFrame() );
  mBtnSaveVolume->Create();
  toolbar->AddWidget( mBtnSaveVolume );
  mBtnSaveVolume->SetText( "Save Volume" );
  mBtnSaveVolume->SetBalloonHelpString( "Save Volume" );
  mBtnSaveVolume->SetCommand( this, "SaveVolumeWithConfirm" );
  try { IconLoader::SetPushButtonIcon( "SaveVolume", mBtnSaveVolume ); }
  catch (...) {}

  vtkKWFrame* spacer = vtkKWFrame::New();
  spacer->SetParent( toolbar->GetFrame() );
  spacer->Create();
  spacer->SetWidth( 5 );
  toolbar->AddWidget( spacer );
  spacer->Delete();

  // Transform Volume
  mBtnTransformVolume = vtkKWPushButton::New();
  mBtnTransformVolume->SetParent( toolbar->GetFrame() );
  mBtnTransformVolume->Create();
  toolbar->AddWidget( mBtnTransformVolume );
  mBtnTransformVolume->SetText( "Transform Volume" );
  mBtnTransformVolume->SetBalloonHelpString( "Transform Volume");
  mBtnTransformVolume->SetCommand( this, "TransformVolume" );
  try { IconLoader::SetPushButtonIcon( "TransformVolume", mBtnTransformVolume ); }
  catch (...) {}

  // Revert.
  mBtnRevertVolume = vtkKWPushButton::New();
  mBtnRevertVolume->SetParent( toolbar->GetFrame() );
  mBtnRevertVolume->Create();
  toolbar->AddWidget( mBtnRevertVolume );
  mBtnRevertVolume->SetText( "Revert Volume" );
  mBtnRevertVolume->SetBalloonHelpString( "Revert Volume" );
  mBtnRevertVolume->SetCommand( this, "RevertToSavedTransform" );
  try { IconLoader::SetPushButtonIcon( "RevertVolume", mBtnRevertVolume ); }
  catch (...) {}

  spacer = vtkKWFrame::New();
  spacer->SetParent( toolbar->GetFrame() );
  spacer->Create();
  spacer->SetWidth( 5 );
  toolbar->AddWidget( spacer );
  spacer->Delete();

  // Color table radio button set.
  vtkKWRadioButtonSet* radBtnSetColors = vtkKWRadioButtonSet::New();
  radBtnSetColors->SetParent( toolbar->GetFrame() );
  radBtnSetColors->Create();
  toolbar->AddWidget( radBtnSetColors );
  radBtnSetColors->PackHorizontallyOn();

  mRadBtnUseGrayScaleColors = radBtnSetColors->AddWidget( 1 );
  mRadBtnUseGrayScaleColors->IndicatorVisibilityOff();
  try { IconLoader::SetCheckButtonIcon( "UseGrayScaleColors", mRadBtnUseGrayScaleColors ); }
  catch (...) {}
  mRadBtnUseGrayScaleColors->SelectedStateOn();
  mRadBtnUseGrayScaleColors->SetCommand( this, "UseGrayScaleColors" );

  mRadBtnUseLUTColors = radBtnSetColors->AddWidget( 0 );
  mRadBtnUseLUTColors->SetCommand( this, "UseLUTColors" );
  mRadBtnUseLUTColors->IndicatorVisibilityOff();
  try { IconLoader::SetCheckButtonIcon( "UseLUTColors", mRadBtnUseLUTColors ); }
  catch (...) {}

  spacer = vtkKWFrame::New();
  spacer->SetParent( toolbar->GetFrame() );
  spacer->Create();
  spacer->SetWidth( 5 );
  toolbar->AddWidget( spacer );
  spacer->Delete();

  // Restore View
  mBtnRestoreView = vtkKWPushButton::New();
  mBtnRestoreView->SetParent( toolbar->GetFrame() );
  mBtnRestoreView->Create();
  toolbar->AddWidget( mBtnRestoreView );
  mBtnRestoreView->SetText( "Restore View" );
  mBtnRestoreView->SetBalloonHelpString( "Restore View" );
  mBtnRestoreView->SetCommand( this, "RestoreView" );
  try { IconLoader::SetPushButtonIcon( "RestoreView", mBtnRestoreView ); }
  catch (...) {}

  spacer = vtkKWFrame::New();
  spacer->SetParent( toolbar->GetFrame() );
  spacer->Create();
  spacer->SetWidth( 5 );
  toolbar->AddWidget( spacer );
  spacer->Delete();

  // Zoom Out
  mBtnZoomOut = vtkKWPushButton::New();
  mBtnZoomOut->SetParent( toolbar->GetFrame() );
  mBtnZoomOut->Create();
  toolbar->AddWidget( mBtnZoomOut );
  mBtnZoomOut->SetText( "Zoom Out" );
  mBtnZoomOut->SetBalloonHelpString( "Zoom Out" );
  mBtnZoomOut->SetCommand( this, "ZoomOut" );
  try { IconLoader::SetPushButtonIcon( "ZoomOut", mBtnZoomOut ); }
  catch (...) {}

  // Zoom In
  mBtnZoomIn = vtkKWPushButton::New();
  mBtnZoomIn->SetParent( toolbar->GetFrame() );
  mBtnZoomIn->Create();
  toolbar->AddWidget( mBtnZoomIn );
  mBtnZoomIn->SetText( "Zoom In" );
  mBtnZoomIn->SetBalloonHelpString( "Zoom In" );
  mBtnZoomIn->SetCommand( this, "ZoomIn" );
  try { IconLoader::SetPushButtonIcon( "ZoomIn", mBtnZoomIn ); }
  catch (...) {}

  spacer = vtkKWFrame::New();
  spacer->SetParent( toolbar->GetFrame() );
  spacer->Create();
  spacer->SetWidth( 5 );
  toolbar->AddWidget( spacer );
  spacer->Delete();

  // Rotate X Pos
  mBtnRotateXPos = vtkKWPushButton::New();
  mBtnRotateXPos->SetParent( toolbar->GetFrame() );
  mBtnRotateXPos->Create();
  toolbar->AddWidget( mBtnRotateXPos );
  mBtnRotateXPos->SetText( "Rotate" );
  mBtnRotateXPos->SetBalloonHelpString( "Rotate" );
  mBtnRotateXPos->SetCommand( mView, "AnimateCameraElevatePositive" );
  try { IconLoader::SetPushButtonIcon( "RotateXPos", mBtnRotateXPos); }
  catch (...) {}

  // Rotate X Neg
  mBtnRotateXNeg = vtkKWPushButton::New();
  mBtnRotateXNeg->SetParent( toolbar->GetFrame() );
  mBtnRotateXNeg->Create();
  toolbar->AddWidget( mBtnRotateXNeg );
  mBtnRotateXNeg->SetText( "Rotate" );
  mBtnRotateXNeg->SetBalloonHelpString( "Rotate" );
  mBtnRotateXNeg->SetCommand( mView, "AnimateCameraElevateNegative" );
  try { IconLoader::SetPushButtonIcon( "RotateXNeg", mBtnRotateXNeg); }
  catch (...) {}

  // Rotate Y Pos
  mBtnRotateYPos = vtkKWPushButton::New();
  mBtnRotateYPos->SetParent( toolbar->GetFrame() );
  mBtnRotateYPos->Create();
  toolbar->AddWidget( mBtnRotateYPos );
  mBtnRotateYPos->SetText( "Rotate" );
  mBtnRotateYPos->SetBalloonHelpString( "Rotate" );
  mBtnRotateYPos->SetCommand( mView, "AnimateCameraAzimuthNegative" );
  try { IconLoader::SetPushButtonIcon( "RotateYPos", mBtnRotateYPos); }
  catch (...) {}

  // Rotate Y Neg
  mBtnRotateYNeg = vtkKWPushButton::New();
  mBtnRotateYNeg->SetParent( toolbar->GetFrame() );
  mBtnRotateYNeg->Create();
  toolbar->AddWidget( mBtnRotateYNeg );
  mBtnRotateYNeg->SetText( "Rotate" );
  mBtnRotateYNeg->SetBalloonHelpString( "Rotate" );
  mBtnRotateYNeg->SetCommand( mView, "AnimateCameraAzimuthPositive" );
  try { IconLoader::SetPushButtonIcon( "RotateYNeg", mBtnRotateYNeg); }
  catch (...) {}

  // Rotate Z Pos
  mBtnRotateZPos = vtkKWPushButton::New();
  mBtnRotateZPos->SetParent( toolbar->GetFrame() );
  mBtnRotateZPos->Create();
  toolbar->AddWidget( mBtnRotateZPos );
  mBtnRotateZPos->SetText( "Rotate" );
  mBtnRotateZPos->SetBalloonHelpString( "Rotate" );
  mBtnRotateZPos->SetCommand( mView, "AnimateCameraRollNegative" );
  try { IconLoader::SetPushButtonIcon( "RotateZPos", mBtnRotateZPos); }
  catch (...) {}

  // Rotate Z Neg
  mBtnRotateZNeg = vtkKWPushButton::New();
  mBtnRotateZNeg->SetParent( toolbar->GetFrame() );
  mBtnRotateZNeg->Create();
  toolbar->AddWidget( mBtnRotateZNeg );
  mBtnRotateZNeg->SetText( "Rotate" );
  mBtnRotateZNeg->SetBalloonHelpString( "Rotate" );
  mBtnRotateZNeg->SetCommand( mView, "AnimateCameraRollPositive" );
  try { IconLoader::SetPushButtonIcon( "RotateZNeg", mBtnRotateZNeg); }
  catch (...) {}

  // Start making our menu bars.
  int nFilePos = GetFileMenuInsertPosition();
  int nEditPos = 0;
  int nViewPos = GetViewMenuInsertPosition();

  // File menu.
  // Load Volume.
  this->GetFileMenu()->
    InsertCommand( nFilePos, "L&oad Volume...", this, "LoadVolumeFromDlog" );
  this->GetFileMenu()->SetItemCompoundModeToLeft( nFilePos );
  this->GetFileMenu()->SetItemAccelerator( nFilePos, "Ctrl+O" );
  try { IconLoader::SetMenuItemIcon( "LoadVolume", GetFileMenu(), nFilePos ); }
  catch (...) {}
  mMenuLoadVolume.menu = GetFileMenu();
  mMenuLoadVolume.nItem = nFilePos;
  nFilePos++;

  // Load LUT.
  this->GetFileMenu()->
    InsertCommand( nFilePos, "Load LUT...", this, "LoadLUTFromDlog" );
  this->GetFileMenu()->SetItemCompoundModeToLeft( nFilePos );
  try { IconLoader::SetMenuItemIcon( "LoadLUT", GetFileMenu(), nFilePos ); }
  catch (...) {}
  mMenuLoadLUT.menu = GetFileMenu();
  mMenuLoadLUT.nItem = nFilePos;
  nFilePos++;

  this->GetFileMenu()->InsertSeparator( nFilePos++ );

  // Save Volume.
  this->GetFileMenu()->
    InsertCommand( nFilePos, "&Save Volume", this, "SaveVolumeWithConfirm" );
  this->GetFileMenu()->SetItemCompoundModeToLeft( nFilePos );
  this->GetFileMenu()->SetItemAccelerator( nFilePos, "Ctrl+S" );
  try {IconLoader::SetMenuItemIcon( "SaveVolume", GetFileMenu(), nFilePos );}
  catch (...) {}
  mMenuSaveVolume.menu = GetFileMenu();
  mMenuSaveVolume.nItem = nFilePos;
  nFilePos++;

  // Save Volume As.
  this->GetFileMenu()->
    InsertCommand( nFilePos, "Save Volume As...", this,"SaveVolumeAsFromDlog");
  this->GetFileMenu()->SetItemCompoundModeToLeft( nFilePos );
  try {IconLoader::SetMenuItemIcon( "SaveVolume", GetFileMenu(), nFilePos );}
  catch (...) {}
  mMenuSaveVolumeAs.menu = GetFileMenu();
  mMenuSaveVolumeAs.nItem = nFilePos;
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
  try {IconLoader::SetMenuItemIcon("TransformVolume",GetEditMenu(),nEditPos);}
  catch (...) {}
  mMenuTransformVolume.menu = GetEditMenu();
  mMenuTransformVolume.nItem = nEditPos;
  nEditPos++;

  this->GetEditMenu()->InsertSeparator( nEditPos++ );

  // Revert.
  this->GetEditMenu()->
    InsertCommand( nEditPos, "&Revert Volume", this, "RevertToSavedTransform");
  this->GetEditMenu()->SetItemCompoundModeToLeft( nEditPos );
  this->GetEditMenu()->SetItemAccelerator( nEditPos, "Ctrl+R" );
  try {IconLoader::SetMenuItemIcon( "RevertVolume", GetEditMenu(), nEditPos );}
  catch (...) {}
  mMenuRevertVolume.menu = GetEditMenu();
  mMenuRevertVolume.nItem = nEditPos;
  nEditPos++;

  // View menu.
  // Restore view.
  this->GetViewMenu()->
    InsertCommand( nViewPos, "Restore &View", this, "RestoreView");
  this->GetViewMenu()->SetItemCompoundModeToLeft( nViewPos );
  this->GetViewMenu()->SetItemAccelerator( nViewPos, "Ctrl+V" );
  try {IconLoader::SetMenuItemIcon( "RestoreView", GetViewMenu(), nViewPos );}
  catch (...) {}
  mMenuRestoreView.menu = GetViewMenu();
  mMenuRestoreView.nItem = nViewPos;
  nViewPos++;

  this->GetViewMenu()->InsertSeparator( nViewPos++ );

  // Color scale items.
  this->GetViewMenu()->
    InsertCommand( nViewPos, "Use GrayScale Colors", this, "UseGrayScaleColors");
  this->GetViewMenu()->SetItemCompoundModeToLeft( nViewPos );
  try {IconLoader::SetMenuItemIcon( "UseGrayScaleColors", GetViewMenu(), nViewPos );}
  catch (...) {}
  mMenuUseGrayScaleColors.menu = GetViewMenu();
  mMenuUseGrayScaleColors.nItem = nViewPos;
  nViewPos++;

  this->GetViewMenu()->
    InsertCommand( nViewPos, "Use LUT Colors", this, "UseLUTColors");
  this->GetViewMenu()->SetItemCompoundModeToLeft( nViewPos );
  try {IconLoader::SetMenuItemIcon( "UseLUTColors", GetViewMenu(), nViewPos );}
  catch (...) {}
  mMenuUseLUTColors.menu = GetViewMenu();
  mMenuUseLUTColors.nItem = nViewPos;
  nViewPos++;

  this->GetViewMenu()->InsertSeparator( nViewPos++ );

  // Zoom Out.
  this->GetViewMenu()->
    InsertCommand( nViewPos, "Zoom Out", this, "ZoomOut");
  this->GetViewMenu()->SetItemCompoundModeToLeft( nViewPos );
  try { IconLoader::SetMenuItemIcon( "ZoomOut", GetViewMenu(), nViewPos ); }
  catch (...) {}
  mMenuZoomOut.menu = GetViewMenu();
  mMenuZoomOut.nItem = nViewPos;
  nViewPos++;

  // Zoom In.
  this->GetViewMenu()->
    InsertCommand( nViewPos, "Zoom In", this, "ZoomIn");
  this->GetViewMenu()->SetItemCompoundModeToLeft( nViewPos );
  try { IconLoader::SetMenuItemIcon( "ZoomIn", GetViewMenu(), nViewPos ); }
  catch (...) {}
  mMenuZoomIn.menu = GetViewMenu();
  mMenuZoomIn.nItem = nViewPos;
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
    if( !mGrayScaleColors ) 
      mGrayScaleColors = vtkLookupTable::New();
    mGrayScaleColors->SetTableRange( mVolume->GetMinValue(), mVolume->GetMaxValue() );
    mGrayScaleColors->SetSaturationRange( 0, 0 );
    mGrayScaleColors->SetHueRange( 0, 0 );
    mGrayScaleColors->SetValueRange( 0, 1 );
    mGrayScaleColors->Build();
    for ( int nEntry = 0; nEntry < mGrayScaleColors->GetIndex(10); nEntry++ )
      mGrayScaleColors->SetTableValue( nEntry, 0, 0, 0, 0 );

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

    // Use grayscale for now.
    this->UseGrayScaleColors();
    
    this->SetStatusText( "Volumnne loaded." );
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
vtkKWOrientMRIWindow::LoadLUTFromDlog () {

  vtkKWLoadSaveDialog* dialog = vtkKWLoadSaveDialog::New();
  dialog->SetApplication( GetApplication() );
  dialog->Create();
  dialog->SetFileTypes( "{LUT {.txt *Colors* *LUT*}} {All {*}}" );
  dialog->RetrieveLastPathFromRegistry( "LoadLUT" );
  dialog->SetDefaultExtension( ".txt" );
  if ( dialog->Invoke() ) {
    dialog->SaveLastPathToRegistry( "LoadLUT" );
    string fnLUT( dialog->GetFileName() );
    this->LoadLUT( fnLUT.c_str() );
  }
}

void
vtkKWOrientMRIWindow::LoadLUT ( const char* ifnLUT ) {

  if( !mView ) throw runtime_error( "mView was NULL" );

  try {
    
    // Try to load the table.
    char fnLUT[1024];
    strncpy( fnLUT, ifnLUT, sizeof(fnLUT) );
    COLOR_TABLE* ctab = CTABreadASCII( fnLUT );
    if ( NULL == ctab ) {
      throw new runtime_error( string("Couldn't open color table file ") +
			       ifnLUT );
    }
      
    // Create a table if not already.
    if ( NULL == mLUTColors )
      mLUTColors = vtkFreesurferLookupTable::New();
    mLUTColors->BuildFromCTAB( ctab );
    
    // Free the ctab.
    CTABfree( &ctab );

    // Use this LUT.
    this->UseLUTColors();

    this->SetStatusText( "LUT loaded." );
    this->AddRecentFile( ifnLUT, this, "LoadLUT" );
    
  } catch ( exception& e ) {
    this->GetApplication()->ErrorMessage( e.what() );
  }

  // Update our menu and buttons.
  this->UpdateCommandStatus();
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
vtkKWOrientMRIWindow::UseGrayScaleColors () {

  if( !mView ) throw runtime_error( "mView was NULL" );
  if( !mGrayScaleColors ) throw runtime_error( "mGrayScaleColors was NULL" );
  if( !mRadBtnUseGrayScaleColors ) throw runtime_error( "mRadBtnUseGrayScaleColors was NULL" );

  // Change the colors.
  mView->SetCurrentVolumeColors( mGrayScaleColors );
  mView->Render();
  
  // Turn on our radio buttion.
  mRadBtnUseGrayScaleColors->SelectedStateOn();
}

void
vtkKWOrientMRIWindow::UseLUTColors () {

  if( !mView ) throw runtime_error( "mView was NULL" );
  if( !mLUTColors ) throw runtime_error( "mLUTColors was NULL" );
  if( !mRadBtnUseLUTColors ) throw runtime_error( "mRadBtnUseLUTColors was NULL" );

  // Change the colors.
  mView->SetCurrentVolumeColors( mLUTColors );
  mView->Render();

  // Turn on our radio buttion.
  mRadBtnUseLUTColors->SelectedStateOn();
}

void
vtkKWOrientMRIWindow::UpdateCommandStatus () {

  if( !mView ) throw runtime_error( "mView was NULL" );

  // Determine the enabled state of our commands.
  if( mVolume ) {
    mBtnTransformVolume->SetStateToNormal();
    mMenuTransformVolume.menu->SetItemStateToNormal( mMenuTransformVolume.nItem );
  } else {
    mBtnTransformVolume->SetStateToDisabled();
    mMenuTransformVolume.menu->SetItemStateToDisabled( mMenuTransformVolume.nItem );
  }

  if( mVolume && mbDirty ) {
    mBtnSaveVolume->SetStateToNormal();
    mMenuSaveVolume.menu->SetItemStateToNormal( mMenuSaveVolume.nItem );
    mMenuSaveVolumeAs.menu->SetItemStateToNormal( mMenuSaveVolumeAs.nItem );
    mBtnRevertVolume->SetStateToNormal();
    mMenuRevertVolume.menu->SetItemStateToNormal( mMenuRevertVolume.nItem );
  } else {
    mBtnSaveVolume->SetStateToDisabled();
    mMenuSaveVolume.menu->SetItemStateToDisabled( mMenuSaveVolume.nItem );
    mMenuSaveVolumeAs.menu->SetItemStateToDisabled( mMenuSaveVolumeAs.nItem );
    mBtnRevertVolume->SetStateToDisabled();
    mMenuRevertVolume.menu->SetItemStateToDisabled( mMenuRevertVolume.nItem );
  }

  if( mVolume && mGrayScaleColors ) {
    mRadBtnUseGrayScaleColors->SetStateToNormal();
    mMenuUseGrayScaleColors.menu->SetItemStateToNormal( mMenuUseGrayScaleColors.nItem );
  } else {
    mRadBtnUseGrayScaleColors->SetStateToDisabled();
    mMenuUseGrayScaleColors.menu->SetItemStateToDisabled( mMenuUseGrayScaleColors.nItem );
  }

  if( mVolume && mLUTColors ) {
    mRadBtnUseLUTColors->SetStateToNormal();
    mMenuUseLUTColors.menu->SetItemStateToNormal( mMenuUseLUTColors.nItem );
  } else {
    mRadBtnUseLUTColors->SetStateToDisabled();
    mMenuUseLUTColors.menu->SetItemStateToDisabled( mMenuUseLUTColors.nItem );
  }
}
