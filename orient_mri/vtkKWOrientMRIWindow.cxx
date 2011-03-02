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
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:35 $
 *    $Revision: 1.20 $
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
#include <sstream>
#include <limits>
#include <assert.h>

#include "vtkKWOrientMRIWindow.h"

#include "OrientMRIEvents.h"
#include "IconLoader.h"
#include "vtkAbstractTransform.h"
#include "vtkCallbackCommand.h"
#include "vtkCamera.h"
#include "vtkFreesurferLookupTable.h"
#include "vtkFSVolumeSource.h"
#include "vtkKWApplication.h"
#include "vtkKWEntry.h"
#include "vtkKWEntryWithLabel.h"
#include "vtkKWIcon.h"
#include "vtkKWLoadSaveDialog.h"
#include "vtkKWMenu.h"
#include "vtkKWMessageDialog.h"
#include "vtkKWOrientMRIView2D.h"
#include "vtkKWOrientMRIView3D.h"
#include "vtkKWPushButton.h"
#include "vtkKWRadioButton.h"
#include "vtkKWRadioButtonSet.h"
#include "vtkKWScale.h"
#include "vtkKWSimpleEntryDialog.h"
#include "vtkKWToolbar.h"
#include "vtkKWToolbarSet.h"
#include "vtkLookupTable.h"
#include "vtkMath.h"
#include "vtkMatrix4x4.h"
#include "vtkObjectFactory.h"
#include "vtkRenderer.h"
#include "vtkTransform.h"

using namespace std;

vtkStandardNewMacro( vtkKWOrientMRIWindow );
vtkCxxRevisionMacro( vtkKWOrientMRIWindow, "$Revision: 1.20 $" );

vtkKWOrientMRIWindow::vtkKWOrientMRIWindow () :
  mbDirty( false ),
  mMenuLoadVolume( NULL ),
  mMenuSaveVolume( NULL ),
  mMenuSaveVolumeAs( NULL ),
  mMenuLoadLUT( NULL ),
  mMenuTransformVolume( NULL ),
  mMenuReorthogonalizeVolume( NULL ),
  mMenuRevertVolume( NULL ),
  mMenuRestoreView( NULL ),
  mMenuZoomOut( NULL ),
  mMenuZoomIn( NULL ),
  mMenuUseGrayScaleColors( NULL ),
  mMenuUseLUTColors( NULL ) {

}

vtkKWOrientMRIWindow::~vtkKWOrientMRIWindow () {

}

void
vtkKWOrientMRIWindow::Create () {

  this->SupportHelpOn();

  this->Superclass::Create();

  this->SetPanelLayoutToSecondaryBelowView();
  this->SecondaryPanelVisibilityOff();
  this->MainPanelVisibilityOff();

  // Create our interior views (3 2ds and a 3d).
  for( int nView = 0; nView < 3; nView++ ) {

    // Create the veiw.
    mView2D[nView] = vtkSmartPointer<vtkKWOrientMRIView2D>::New();
    mView2D[nView]->SetParent( this->GetViewFrame() );
    mView2D[nView]->Create();
    mView2D[nView]->SetOrientation( nView );

    // Create scales to control the through plane of the 2D
    // views. Connect it to the view to set the through plane.
    mScaleThroughPlane[nView] = vtkSmartPointer<vtkKWScale>::New();
    mScaleThroughPlane[nView]->SetParent( this->GetViewFrame() );
    mScaleThroughPlane[nView]->Create();
    mScaleThroughPlane[nView]->SetCommand( mView2D[nView], "SetThroughPlane" );
    mScaleThroughPlane[nView]->ValueVisibilityOff();
  }

  // 3D view.
  mView3D = vtkSmartPointer<vtkKWOrientMRIView3D>::New();
  mView3D->SetParent( this->GetViewFrame() );
  mView3D->Create();

  // Set them up in a grid, with the scales underneath the 2D views.
  this->Script( "grid %s -column 0 -row 0 -sticky news -padx 2",
		mView2D[X]->GetWidgetName() );
  this->Script( "grid %s -column 0 -row 1 -sticky news -padx 2",
		mScaleThroughPlane[X]->GetWidgetName() );
  this->Script( "grid %s -column 1 -row 0 -sticky news -padx 2",
		mView2D[Y]->GetWidgetName() );
  this->Script( "grid %s -column 1 -row 1 -sticky news -padx 2",
		mScaleThroughPlane[Y]->GetWidgetName() );
  this->Script( "grid %s -column 0 -row 2 -sticky news -padx 2 -pady 2",
		mView2D[Z]->GetWidgetName() );
  this->Script( "grid %s -column 0 -row 3 -sticky news -padx 2",
		mScaleThroughPlane[Z]->GetWidgetName() );
  this->Script( "grid %s -column 1 -row 2 -rowspan 2 -sticky news -padx 2 -pady 2",
		mView3D->GetWidgetName() );
  this->Script( "grid columnconfigure %s 0 -weight 1", 
		this->GetViewFrame()->GetWidgetName() );
  this->Script( "grid columnconfigure %s 1 -weight 1", 
		this->GetViewFrame()->GetWidgetName() );
  this->Script( "grid rowconfigure %s 0 -weight 1",
		this->GetViewFrame()->GetWidgetName() );
  this->Script( "grid rowconfigure %s 2 -weight 1",
		this->GetViewFrame()->GetWidgetName() );

  // Make a callback and observe our view for UserTransformChanged events.
  vtkSmartPointer<vtkCallbackCommand> callback =
    vtkSmartPointer<vtkCallbackCommand>::New();
  callback->SetCallback( UserTransformChanged );
  callback->SetClientData( this );
  mView3D->AddObserver( OrientMRIEvents::UserTransformChanged, callback );

  // Let the view3D observe us for VolumeToRASTransformChanged
  // events. This is kind of an awkward way to do this.
  callback = vtkSmartPointer<vtkCallbackCommand>::New();
  callback->
    SetCallback( vtkKWOrientMRIView3D::VolumeToRASTransformChangedCallback );
  callback->SetClientData( mView3D );
  this->AddObserver( OrientMRIEvents::VolumeToRASTransformChanged, callback );

  // Make a toolbar.
  vtkSmartPointer<vtkKWToolbar> toolbar = 
    vtkSmartPointer<vtkKWToolbar>::New();
  toolbar->SetName( "Main" );
  toolbar->SetParent( this->GetMainToolbarSet()->GetToolbarsFrame() );
  toolbar->Create();
  this->GetMainToolbarSet()->AddToolbar( toolbar );

  // Populate the toolbar.

  // Load Volume
  mBtnLoadVolume = vtkSmartPointer<vtkKWPushButton>::New();
  mBtnLoadVolume->SetParent( toolbar->GetFrame() );
  mBtnLoadVolume->Create();
  toolbar->AddWidget( mBtnLoadVolume );
  mBtnLoadVolume->SetText( "Load Volume" );
  mBtnLoadVolume->SetBalloonHelpString( "Load Volume" );
  mBtnLoadVolume->SetCommand( this, "LoadVolumeFromDlog" );
  try { IconLoader::SetPushButtonIcon( "LoadVolume", mBtnLoadVolume ); }
  catch (...) {}

  // Save Volume
  mBtnSaveVolume = vtkSmartPointer<vtkKWPushButton>::New();
  mBtnSaveVolume->SetParent( toolbar->GetFrame() );
  mBtnSaveVolume->Create();
  toolbar->AddWidget( mBtnSaveVolume );
  mBtnSaveVolume->SetText( "Save Volume" );
  mBtnSaveVolume->SetBalloonHelpString( "Save Volume" );
  mBtnSaveVolume->SetCommand( this, "SaveVolumeWithConfirm" );
  try { IconLoader::SetPushButtonIcon( "SaveVolume", mBtnSaveVolume ); }
  catch (...) {}

  this->AddSpacerToToolbar( *toolbar );

  // Transform Volume
  mBtnTransformVolume = vtkSmartPointer<vtkKWPushButton>::New();
  mBtnTransformVolume->SetParent( toolbar->GetFrame() );
  mBtnTransformVolume->Create();
  toolbar->AddWidget( mBtnTransformVolume );
  mBtnTransformVolume->SetText( "Transform Volume" );
  mBtnTransformVolume->SetBalloonHelpString( "Transform Volume");
  mBtnTransformVolume->SetCommand( this, "TransformVolumeWithUserTransform" );
  try { IconLoader::SetPushButtonIcon( "TransformVolume", mBtnTransformVolume ); }
  catch (...) {}

  // Revert.
  mBtnRevertVolume = vtkSmartPointer<vtkKWPushButton>::New();
  mBtnRevertVolume->SetParent( toolbar->GetFrame() );
  mBtnRevertVolume->Create();
  toolbar->AddWidget( mBtnRevertVolume );
  mBtnRevertVolume->SetText( "Revert Volume" );
  mBtnRevertVolume->SetBalloonHelpString( "Revert Volume" );
  mBtnRevertVolume->SetCommand( this, "RevertToSavedTransform" );
  try { IconLoader::SetPushButtonIcon( "RevertVolume", mBtnRevertVolume ); }
  catch (...) {}

  this->AddSpacerToToolbar( *toolbar );

  // Color table radio button set.
  vtkSmartPointer<vtkKWRadioButtonSet> radBtnSetColors =
    vtkSmartPointer<vtkKWRadioButtonSet>::New();
  radBtnSetColors->SetParent( toolbar->GetFrame() );
  radBtnSetColors->Create();
  toolbar->AddWidget( radBtnSetColors );
  radBtnSetColors->PackHorizontallyOn();

  mRadBtnUseGrayScaleColors.TakeReference( radBtnSetColors->AddWidget( 1 ) );
  mRadBtnUseGrayScaleColors->IndicatorVisibilityOff();
  try { IconLoader::SetCheckButtonIcon( "UseGrayScaleColors", mRadBtnUseGrayScaleColors ); }
  catch (...) {}
  mRadBtnUseGrayScaleColors->SelectedStateOn();
  mRadBtnUseGrayScaleColors->SetCommand( this, "UseGrayScaleColors" );

  mRadBtnUseLUTColors.TakeReference( radBtnSetColors->AddWidget( 0 ) );
  mRadBtnUseLUTColors->SetCommand( this, "UseLUTColors" );
  mRadBtnUseLUTColors->IndicatorVisibilityOff();
  try { IconLoader::SetCheckButtonIcon( "UseLUTColors", mRadBtnUseLUTColors ); }
  catch (...) {}

  this->AddSpacerToToolbar( *toolbar );

  // Restore View
  mBtnRestoreView = vtkSmartPointer<vtkKWPushButton>::New();
  mBtnRestoreView->SetParent( toolbar->GetFrame() );
  mBtnRestoreView->Create();
  toolbar->AddWidget( mBtnRestoreView );
  mBtnRestoreView->SetText( "Restore View" );
  mBtnRestoreView->SetBalloonHelpString( "Restore View" );
  mBtnRestoreView->SetCommand( this, "RestoreView" );
  try { IconLoader::SetPushButtonIcon( "RestoreView", mBtnRestoreView ); }
  catch (...) {}

  this->AddSpacerToToolbar( *toolbar );

  // Zoom Out
  mBtnZoomOut = vtkSmartPointer<vtkKWPushButton>::New();
  mBtnZoomOut->SetParent( toolbar->GetFrame() );
  mBtnZoomOut->Create();
  toolbar->AddWidget( mBtnZoomOut );
  mBtnZoomOut->SetText( "Zoom Out" );
  mBtnZoomOut->SetBalloonHelpString( "Zoom Out" );
  mBtnZoomOut->SetCommand( this, "ZoomOut" );
  try { IconLoader::SetPushButtonIcon( "ZoomOut", mBtnZoomOut ); }
  catch (...) {}

  // Zoom In
  mBtnZoomIn = vtkSmartPointer<vtkKWPushButton>::New();
  mBtnZoomIn->SetParent( toolbar->GetFrame() );
  mBtnZoomIn->Create();
  toolbar->AddWidget( mBtnZoomIn );
  mBtnZoomIn->SetText( "Zoom In" );
  mBtnZoomIn->SetBalloonHelpString( "Zoom In" );
  mBtnZoomIn->SetCommand( this, "ZoomIn" );
  try { IconLoader::SetPushButtonIcon( "ZoomIn", mBtnZoomIn ); }
  catch (...) {}

  this->AddSpacerToToolbar( *toolbar );

  // Rotate X Pos
  mBtnRotateXPos = vtkSmartPointer<vtkKWPushButton>::New();
  mBtnRotateXPos->SetParent( toolbar->GetFrame() );
  mBtnRotateXPos->Create();
  toolbar->AddWidget( mBtnRotateXPos );
  mBtnRotateXPos->SetText( "Rotate" );
  mBtnRotateXPos->SetBalloonHelpString( "Rotate" );
  mBtnRotateXPos->SetCommand( mView3D, "RotateUserTransform 0 -90" );
  stringstream ssCommand;
  ssCommand << "DoRotateUserTransformDialog 0 -1; "
	    << mBtnRotateXPos->GetTclName() << " SetReliefToRaised";
  mBtnRotateXPos->AddBinding( "<Control-ButtonRelease-1>", this, 
			      ssCommand.str().c_str() );
  try { IconLoader::SetPushButtonIcon( "RotateXPos", mBtnRotateXPos); }
  catch (...) {}

  // Rotate X Neg
  mBtnRotateXNeg = vtkSmartPointer<vtkKWPushButton>::New();
  mBtnRotateXNeg->SetParent( toolbar->GetFrame() );
  mBtnRotateXNeg->Create();
  toolbar->AddWidget( mBtnRotateXNeg );
  mBtnRotateXNeg->SetText( "Rotate" );
  mBtnRotateXNeg->SetBalloonHelpString( "Rotate" );
  mBtnRotateXNeg->SetCommand( mView3D, "RotateUserTransform 0 90" );
  mBtnRotateXNeg->AddBinding( "<Control-ButtonRelease-1>", this,
			      "DoRotateUserTransformDialog 0 1" );
  try { IconLoader::SetPushButtonIcon( "RotateXNeg", mBtnRotateXNeg); }
  catch (...) {}

  // Rotate Y Pos
  mBtnRotateYPos = vtkSmartPointer<vtkKWPushButton>::New();
  mBtnRotateYPos->SetParent( toolbar->GetFrame() );
  mBtnRotateYPos->Create();
  toolbar->AddWidget( mBtnRotateYPos );
  mBtnRotateYPos->SetText( "Rotate" );
  mBtnRotateYPos->SetBalloonHelpString( "Rotate" );
  mBtnRotateYPos->SetCommand( mView3D, "RotateUserTransform 2 90" );
  mBtnRotateYPos->AddBinding( "<Control-ButtonRelease-1>", this,
				   "DoRotateUserTransformDialog 2 1" );
  try { IconLoader::SetPushButtonIcon( "RotateYPos", mBtnRotateYPos); }
  catch (...) {}

  // Rotate Y Neg
  mBtnRotateYNeg = vtkSmartPointer<vtkKWPushButton>::New();
  mBtnRotateYNeg->SetParent( toolbar->GetFrame() );
  mBtnRotateYNeg->Create();
  toolbar->AddWidget( mBtnRotateYNeg );
  mBtnRotateYNeg->SetText( "Rotate" );
  mBtnRotateYNeg->SetBalloonHelpString( "Rotate" );
  mBtnRotateYNeg->SetCommand( mView3D, "RotateUserTransform 2 -90" );
  mBtnRotateYNeg->AddBinding( "<Control-ButtonRelease-1>", this,
			      "DoRotateUserTransformDialog 2 -1" );
  try { IconLoader::SetPushButtonIcon( "RotateYNeg", mBtnRotateYNeg); }
  catch (...) {}

  // Rotate Z Pos
  mBtnRotateZPos = vtkSmartPointer<vtkKWPushButton>::New();
  mBtnRotateZPos->SetParent( toolbar->GetFrame() );
  mBtnRotateZPos->Create();
  toolbar->AddWidget( mBtnRotateZPos );
  mBtnRotateZPos->SetText( "Rotate" );
  mBtnRotateZPos->SetBalloonHelpString( "Rotate" );
  mBtnRotateZPos->SetCommand( mView3D, "RotateUserTransform 1 -90" );
  mBtnRotateZPos->AddBinding( "<Control-ButtonRelease-1>", this,
			      "DoRotateUserTransformDialog 1 -1" );
  try { IconLoader::SetPushButtonIcon( "RotateZPos", mBtnRotateZPos); }
  catch (...) {}

  // Rotate Z Neg
  mBtnRotateZNeg = vtkSmartPointer<vtkKWPushButton>::New();
  mBtnRotateZNeg->SetParent( toolbar->GetFrame() );
  mBtnRotateZNeg->Create();
  toolbar->AddWidget( mBtnRotateZNeg );
  mBtnRotateZNeg->SetText( "Rotate" );
  mBtnRotateZNeg->SetBalloonHelpString( "Rotate" );
  mBtnRotateZNeg->SetCommand( mView3D, "RotateUserTransform 1 90" );
  mBtnRotateZNeg->AddBinding( "<Control-ButtonRelease-1>", this,
			      "DoRotateUserTransformDialog 1 1" );
  try { IconLoader::SetPushButtonIcon( "RotateZNeg", mBtnRotateZNeg); }
  catch (...) {}

  // Start making our menu bars.
  int nItem = this->GetFileMenuInsertPosition();

  // File menu.
  // Load Volume.
  mMenuLoadVolume = new MenuItem();
  mMenuLoadVolume->
    MakeCommand( this->GetFileMenu(), nItem++,
		 "L&oad Volume...", this, "LoadVolumeFromDlog",
		 "Ctrl+O", "LoadVolume" );

  // Load LUT.
  mMenuLoadLUT = new MenuItem();
  mMenuLoadLUT->
    MakeCommand( this->GetFileMenu(), nItem++,
		 "Load LUT...", this, "LoadLUTFromDlog",
		 "", "LoadVolume" );

  this->GetFileMenu()->InsertSeparator( nItem++ );

  // Save Volume.
  mMenuSaveVolume = new MenuItem();
  mMenuSaveVolume->
    MakeCommand( this->GetFileMenu(), nItem++,
		 "&Save Volume...", this, "SaveVolumeWithConfirm",
		 "Ctrl+S", "SaveVolume" );

  // Save Volume As.
  mMenuSaveVolumeAs = new MenuItem();
  mMenuSaveVolumeAs->
    MakeCommand( this->GetFileMenu(), nItem++,
		 "Save Volume As...", this, "SaveVolumeAsFromDlog",
		 "", "SaveVolume" );

  this->GetFileMenu()->InsertSeparator( nItem++ );

  this->InsertRecentFilesMenu( nItem++, this );

  this->GetFileMenu()->InsertSeparator( nItem++ );

  // Edit menu.
  nItem = 0;

  // Transform Volume.
  mMenuTransformVolume = new MenuItem();
  mMenuTransformVolume->
    MakeCommand( this->GetEditMenu(), nItem++,
		 "Trans&form Volume", this, "TransformVolumeWithUserTransform",
		 "Ctrl+F", "TransformVolume" );

  mMenuReorthogonalizeVolume = new MenuItem();
  mMenuReorthogonalizeVolume->
    MakeCommand( this->GetEditMenu(), nItem++,
		 "Reorthogonalize Volume", this, "ReorthogonalizeVolume",
		 NULL, NULL );

  this->GetEditMenu()->InsertSeparator( nItem++ );

  // Revert.
  mMenuRevertVolume = new MenuItem();
  mMenuRevertVolume->
    MakeCommand( this->GetEditMenu(), nItem++,
		 "&Revert Volume", this, "RevertToSavedTransform",
		 "Ctrl+R", "RevertVolume" );

  // View menu.
  nItem = this->GetViewMenuInsertPosition();

  // Restore view.
  mMenuRestoreView = new MenuItem();
  mMenuRestoreView->
    MakeCommand( this->GetViewMenu(), nItem++,
		 "Restore &View", this, "RestoreView",
		 "Ctrl+V", "RestoreView" );

  this->GetViewMenu()->InsertSeparator( nItem++ );

  // Color scale items.
  mMenuUseGrayScaleColors = new MenuItem();
  mMenuUseGrayScaleColors->
    MakeCommand( this->GetViewMenu(), nItem++,
		 "Use Grayscale Colors", this, "UseGrayScaleColors",
		 "", "UseGrayScaleColors" );

  mMenuUseLUTColors = new MenuItem();
  mMenuUseLUTColors->
    MakeCommand( this->GetViewMenu(), nItem++,
		 "Use LUT Colors", this, "UseLUTColors",
		 "", "UseLUTColors" );

  this->GetViewMenu()->InsertSeparator( nItem++ );

  // Zoom Out.
  mMenuZoomOut = new MenuItem();
  mMenuZoomOut->
    MakeCommand( this->GetViewMenu(), nItem++,
		 "Zoom Out", this, "ZoomOut",
		 "", "ZoomOut" );

  // Zoom In.
  mMenuZoomIn = new MenuItem();
  mMenuZoomIn->
    MakeCommand( this->GetViewMenu(), nItem++,
		 "Zoom In", this, "ZoomIn",
		 "", "ZoomIn" );

  // Update our menu and buttons.
  this->UpdateCommandStatus();
}

void
vtkKWOrientMRIWindow::LoadVolumeFromDlog () {

  // The smart pointer version is commented out because of a bug that
  // occurs when the dialog is deleted.
  //  vtkSmartPointer<vtkKWLoadSaveDialog> dialog = 
  //    vtkSmartPointer<vtkKWLoadSaveDialog>::New();
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

  assert( mView2D[X].GetPointer() );
  assert( mView2D[Y].GetPointer() );
  assert( mView2D[Z].GetPointer() );
  assert( mView3D.GetPointer() );

  try {
    
    // Try to load the volume.
    vtkSmartPointer<vtkFSVolumeSource> volume = 
      vtkSmartPointer<vtkFSVolumeSource>::New();
    //volume->ActualSpacingOn();
    volume->MRIRead( ifnVolume );
    volume->Update();
      
    // Save the new volume.
    mVolume = volume;

    // Initialize our transforms.
    mVolumeToRASTransform = vtkSmartPointer<vtkImageReslice>::New();
#if 0   
    vtkSmartPointer<vtkImageReslice> pixelSizer =
      vtkSmartPointer<vtkImageReslice>::New();
    pixelSizer->SetInputConnection( mVolume->GetOutputPort() );
    pixelSizer->SetOutputDimensionality( 3 );

    vtkSmartPointer<vtkTransform> pixelSizerTransform =
      vtkSmartPointer<vtkTransform>::New();
    pixelSizerTransform->Scale( volume->GetPixelSizeX(),
				volume->GetPixelSizeY(),
				volume->GetPixelSizeZ() );
    pixelSizer->SetResliceTransform( pixelSizerTransform );



    mVolumeToRASTransform->SetInputConnection( pixelSizer->GetOutputPort() );
#else
    mVolumeToRASTransform->SetInputConnection( mVolume->GetOutputPort() );
#endif

    mVolumeToRASTransform->SetOutputDimensionality( 3 );
    mVolumeToRASTransform->AutoCropOutputOff();

    mUserTransform = vtkSmartPointer<vtkImageReslice>::New();
    mUserTransform->
      SetInputConnection( mVolumeToRASTransform->GetOutputPort() );
    mUserTransform->SetOutputDimensionality( 3 );
    mUserTransform->AutoCropOutputOff();

    // Update these objects before passing them into our view.
    mVolumeToRASTransform->Update();
    mUserTransform->Update();

    // Set the resliced (transformed to RAS) data in the view. THe 2D
    // views get the output of the user transform, so they are always
    // updated with the combination of the RAS and the user
    // transform. The 3D view only gets the RAS transform, as it
    // defines the user transform.
    mView2D[X]->SetFSVolumeAndImage( *mVolume, *mUserTransform->GetOutput() );
    mView2D[Y]->SetFSVolumeAndImage( *mVolume, *mUserTransform->GetOutput() );
    mView2D[Z]->SetFSVolumeAndImage( *mVolume, *mUserTransform->GetOutput() );
    mView3D->SetFSVolumeAndImage( *mVolume, 
				  *mVolumeToRASTransform->GetOutput() );

    // Create a basic LUT.
    mGrayScaleColors = vtkSmartPointer<vtkLookupTable>::New();
    mGrayScaleColors->SetTableRange( mVolume->GetMinValue(), mVolume->GetMaxValue() );
    mGrayScaleColors->SetSaturationRange( 0, 0 );
    mGrayScaleColors->SetHueRange( 0, 0 );
    mGrayScaleColors->SetValueRange( 0, 1 );
    mGrayScaleColors->Build();
    for ( int nEntry = 0; nEntry < mGrayScaleColors->GetIndex(10); nEntry++ )
      mGrayScaleColors->SetTableValue( nEntry, 0, 0, 0, 0 );

    // Get and save the original VoxelToRAS.
    mOriginalVoxelToRASMatrix = vtkSmartPointer<vtkMatrix4x4>::New();
    mOriginalVoxelToRASMatrix->DeepCopy( mVolume->GetVoxelToRASMatrix() );

    // Use grayscale for now.
    this->UseGrayScaleColors();

    // Build the initial volume transform.
    this->RebuildVolumeToRASTransform();
    this->UpdateUserTransform();

    // Render everything.
    mView2D[X]->Render();
    mView2D[Y]->Render();
    mView2D[Z]->Render();
    mView3D->Render();
    
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

  // The smart pointer version is commented out because of a bug that
  // occurs when the dialog is deleted.
  //  vtkSmartPointer<vtkKWLoadSaveDialog> dialog = 
  //    vtkSmartPointer<vtkKWLoadSaveDialog>::New();
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

  assert( mVolume.GetPointer() );

  // Tell the volume to write itself.
  mVolume->MRIWrite();

  // Get the new VoxelToRAS as our original.
  mOriginalVoxelToRASMatrix = vtkSmartPointer<vtkMatrix4x4>::New();
  mOriginalVoxelToRASMatrix->DeepCopy( mVolume->GetVoxelToRASMatrix() );

  // No longer dirty.
  mbDirty = false;
}

void
vtkKWOrientMRIWindow::SaveVolume ( const char* ifnVolume ) {

  assert( mVolume.GetPointer() );

  // Set the volume filename.
  mVolume->MRIWrite( ifnVolume );

  // Get the new VoxelToRAS as our original.
  mOriginalVoxelToRASMatrix = vtkSmartPointer<vtkMatrix4x4>::New();
  mOriginalVoxelToRASMatrix->DeepCopy( mVolume->GetVoxelToRASMatrix() );

  // No longer dirty.
  mbDirty = false;
}

void
vtkKWOrientMRIWindow::LoadLUTFromDlog () {

  // The smart pointer version is commented out because of a bug that
  // occurs when the dialog is deleted.
  //  vtkSmartPointer<vtkKWLoadSaveDialog> dialog = 
  //    vtkSmartPointer<vtkKWLoadSaveDialog>::New();
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

  assert( mView2D[X].GetPointer() );
  assert( mView2D[Y].GetPointer() );
  assert( mView2D[Z].GetPointer() );
  assert( mView3D.GetPointer() );

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
    mLUTColors = vtkSmartPointer<vtkFreesurferLookupTable>::New();
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

  assert( mVolume.GetPointer() );
  assert( mOriginalVoxelToRASMatrix.GetPointer() );

  try {

    // Set our save original matrix.
    mVolume->SetVoxelToRASMatrix( *mOriginalVoxelToRASMatrix );
      
    // Rebuild our transform.
    this->RebuildVolumeToRASTransform();

    // Restore the view.
    this->RestoreView();

    // No longer dirty.
    mbDirty = false;

  } catch ( exception& e ) {
    this->GetApplication()->ErrorMessage( e.what() );
  }
  
  // Update our menu and buttons.
  this->UpdateCommandStatus();
}


void
vtkKWOrientMRIWindow::TransformVolumeWithUserTransform () {

  assert( mVolume.GetPointer() );
  assert( mOriginalVoxelToRASMatrix.GetPointer() );

  // orient_mri only seeks to change the rotation of a volume, not its
  // location in space. We want to keep the same CRAS. So what we do
  // is take the existing volume in RAS space, translate to the
  // center, perform the rotation, and then translate it back, so that
  // the volume was rotated around the existing center.

  // Get the current VoxelToRAS matrix
  vtkSmartPointer<vtkMatrix4x4> currentVoxelToRASM = 
    vtkSmartPointer<vtkMatrix4x4>::New();
  currentVoxelToRASM->DeepCopy( mVolume->GetVoxelToRASMatrix() );

  // We need to transform this volume so that the center is at the RAS
  // origin. Get the RAS bounds. and calculate the amount we need to
  // translate to do this.
  float RASBounds[6];
  mVolume->GetRASBounds( RASBounds );
  double trans[3];
  trans[0] = (RASBounds[1]-RASBounds[0])/2.0 + RASBounds[0];
  trans[1] = (RASBounds[3]-RASBounds[2])/2.0 + RASBounds[2];
  trans[2] = (RASBounds[5]-RASBounds[4])/2.0 + RASBounds[4];

  // Get the user transform. Get it as a vtkTransform as well because
  // we want to extract the orientation.
  vtkSmartPointer<vtkMatrix4x4> userTransformM =
    vtkSmartPointer<vtkMatrix4x4>::New();
  vtkSmartPointer<vtkTransform> userTransform =
    vtkSmartPointer<vtkTransform>::New();
  userTransformM->DeepCopy( &this->GetUserTransform() );
  userTransformM->Invert();
  userTransform->SetMatrix( userTransformM );

  // Extract the orientation from thse user transform.
  double orientation[3];
  userTransform->GetOrientation( orientation );

  // Make a new transform. Set it to post multiply because we're going
  // to start with an initial transform and then stack stuff on after
  // that.
  vtkSmartPointer<vtkTransform> transform = 
    vtkSmartPointer<vtkTransform>::New();
  transform->PostMultiply();

  // Initialize it with the current VoxelToRAS. 
  transform->SetMatrix( currentVoxelToRASM );

  // First, undo the translation to get it to the center.
  transform->Translate( -trans[0],-trans[1], -trans[2] );
   
  // Now rotate it with the orientation from the user transform.
  transform->RotateX( orientation[0] );
  transform->RotateY( orientation[1] );
  transform->RotateZ( orientation[2] );

  // Now translate it back to where it was before.
  transform->Translate( trans[0], trans[1], trans[2] );
  
  // Set the matrix in the volume.
  mVolume->SetVoxelToRASMatrix( *transform->GetMatrix() );
  mVolume->Update();

#if 0
  cerr << "current " << endl << *currentVoxelToRASM << endl
       << "user " << endl << *userTransformM << endl
       << "composed " << endl << *transform->GetMatrix() << endl;
#endif
  
  // Recalc our reslice transform.
  this->RebuildVolumeToRASTransform();

  // Update our menu and buttons.
  this->UpdateCommandStatus();
}

void
vtkKWOrientMRIWindow::ReorthogonalizeVolume () {

  double const* xStart = NULL;
  double const* xEnd = NULL;
  double const* yStart = NULL;
  double const* yEnd = NULL;
  double const* zStart = NULL;
  double const* zEnd = NULL;

  vtkKWOrientMRIView2D::GetReorthoPoints( xStart, xEnd, 
					  yStart, yEnd,
					  zStart, zEnd );

  // Make sure we have good points; they must all be different.
  if( (fabs(xStart[0] - xEnd[0]) < numeric_limits<double>::epsilon() &&
       fabs(xStart[1] - xEnd[2]) < numeric_limits<double>::epsilon() &&
       fabs(xStart[2] - xEnd[1]) < numeric_limits<double>::epsilon()) ||
      (fabs(yStart[0] - yEnd[0]) < numeric_limits<double>::epsilon() &&
       fabs(yStart[1] - yEnd[2]) < numeric_limits<double>::epsilon() &&
       fabs(yStart[2] - yEnd[1]) < numeric_limits<double>::epsilon()) ||
      (fabs(zStart[0] - zEnd[0]) < numeric_limits<double>::epsilon() &&
       fabs(zStart[1] - zEnd[2]) < numeric_limits<double>::epsilon() &&
       fabs(zStart[2] - zEnd[1]) < numeric_limits<double>::epsilon()) ) {
    this->GetApplication()->ErrorMessage( "You must draw all three orthogonalization arrows before reorthogonalizing the volume." );
    return;
  }

  // First make vectors out of our points and normalize them.
  double x[3], y[3], z[3];
  for( int n = 0; n < 3; n++ ) {
    x[n] = xEnd[n] - xStart[n];
    y[n] = yEnd[n] - yStart[n];
    z[n] = zEnd[n] - zStart[n];
  }
  vtkMath::Normalize( x );
  vtkMath::Normalize( y );
  vtkMath::Normalize( z );

  // The rest of this is from some matlab code I got from Bruce.
  // d = x'*y ;
  double d = vtkMath::Dot( x, y );

  // y = y - d*x ;
  for( int n = 0; n < 3; n++ )
    y[n] = y[n] - d*x[n];

  vtkMath::Normalize( y );

  // t = cross(x,y);
  double t[3];
  vtkMath::Cross( x, y, t );
  
  // d = t'*z ;
  d = vtkMath::Dot( t, z );

  // if (d < 0)
  //   z = -1*t ;
  // else
  //   z = t;
  // end
  if( d < 0 )
    for( int n = 0; n < 3; n++ )
      z[n] = -t[n];
  else
    for( int n = 0; n < 3; n++ )
      z[n] = t[n];

  vtkMath::Normalize( z );

  // Now the 3x3 matrix made but putting x in the first row, y in the
  // second, and z in the third, is our transform.
  vtkSmartPointer<vtkMatrix4x4> matrix = 
    vtkSmartPointer<vtkMatrix4x4>::New();

  // 0,0 0,1 0,2 0,3
  // 1,0 1,1 1,2 1,3
  // 2,0 2,1 2,2 2,3
  // 3,0 3,1 3,2 3,3
  matrix->SetElement( 0, 0, x[0] );
  matrix->SetElement( 0, 1, x[1] );
  matrix->SetElement( 0, 2, x[2] );
  matrix->SetElement( 0, 3, 0 );
  matrix->SetElement( 1, 0, y[0] );
  matrix->SetElement( 1, 1, y[1] );
  matrix->SetElement( 1, 2, y[2] );
  matrix->SetElement( 1, 3, 0 );
  matrix->SetElement( 2, 0, z[0] );
  matrix->SetElement( 2, 1, z[1] );
  matrix->SetElement( 2, 2, z[2] );
  matrix->SetElement( 2, 3, 0 );
  matrix->SetElement( 3, 0, 0 );
  matrix->SetElement( 3, 1, 0 );
  matrix->SetElement( 3, 2, 0 );
  matrix->SetElement( 3, 3, 1 );

  // Note that the following code looks similar to the code in
  // TransformVolumeWithUserTransform, except in that function, we
  // extract the rotation elements from the user transform and rotate
  // the voxelToRAS by them, where as here, we actually concatenate
  // the volume created above, including skewing elements. This is
  // done because we actually want to allow the volume to be skewed.

  // Make a new transform. Set it to post multiply because we're going
  // to start with an initial transform and then stack stuff on after
  // that.
  vtkSmartPointer<vtkTransform> transform = 
    vtkSmartPointer<vtkTransform>::New();
  transform->PostMultiply();

  // Get the current VoxelToRAS matrix
  vtkSmartPointer<vtkMatrix4x4> currentVoxelToRASM = 
    vtkSmartPointer<vtkMatrix4x4>::New();
  currentVoxelToRASM->DeepCopy( mVolume->GetVoxelToRASMatrix() );

  // Initialize it with the current VoxelToRAS. 
  transform->SetMatrix( currentVoxelToRASM );

  // Concatenate the matrix we got.
  transform->Concatenate( matrix );
  
  // Set the matrix in the volume.
  mVolume->SetVoxelToRASMatrix( *transform->GetMatrix() );
  mVolume->Update();

  // Recalc our reslice transform.
  this->RebuildVolumeToRASTransform();

  // Reset the reortho points.
  vtkKWOrientMRIView2D::ResetOrthoPoints();
}

void
vtkKWOrientMRIWindow::RestoreView () {

  assert( mView3D.GetPointer() );

  // Tell the views to restore.
  mView3D->RestoreView();

  mView2D[X]->Render();
  mView2D[Y]->Render();
  mView2D[Z]->Render();

  // Update our menu and buttons.
  this->UpdateCommandStatus();
}


void
vtkKWOrientMRIWindow::ZoomBy ( float iFactor ) {

  assert( mView3D.GetPointer() );

  // Tell the view to zoom.
  mView3D->ZoomBy( iFactor );

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

  assert( mView3D.GetPointer() );
  assert( mGrayScaleColors.GetPointer() );
  assert( mRadBtnUseGrayScaleColors.GetPointer() );

  // Change the colors.
  mView3D->SetImageColors( *mGrayScaleColors );
  mView3D->Render();
  
  // Turn on our radio buttion.
  mRadBtnUseGrayScaleColors->SelectedStateOn();
}

void
vtkKWOrientMRIWindow::UseLUTColors () {

  assert( mView3D.GetPointer() );
  assert( mLUTColors.GetPointer() );
  assert( mRadBtnUseLUTColors.GetPointer() );

  // Change the colors.
  mView3D->SetImageColors( *mLUTColors );
  mView3D->Render();

  // Turn on our radio buttion.
  mRadBtnUseLUTColors->SelectedStateOn();
}

void
vtkKWOrientMRIWindow::UpdateUserTransform () {

  assert( mUserTransform.GetPointer() );
  assert( mView2D[X].GetPointer() );
  assert( mView2D[Y].GetPointer() );
  assert( mView2D[Z].GetPointer() );

  // Get the user transform and copy it into mUserTransform.
  vtkSmartPointer<vtkTransform> userTransform = 
    vtkSmartPointer<vtkTransform>::New();
  userTransform->SetMatrix( &this->GetUserTransform() );
  mUserTransform->SetResliceTransform( userTransform );
  mUserTransform->Update();
  
  // Redraw our 2D views.
  for( int nView = 0; nView < 3; nView++ )
    mView2D[nView]->Render();
}

void
vtkKWOrientMRIWindow::UpdateCommandStatus () {

  assert( mView3D.GetPointer() );

  // Determine the enabled state of our commands.
  if( mVolume.GetPointer() ) {
    mBtnTransformVolume->SetStateToNormal();
    mMenuTransformVolume->SetStateToNormal();

    mMenuReorthogonalizeVolume->SetStateToNormal();
  } else {
    mBtnTransformVolume->SetStateToDisabled();
    mMenuTransformVolume->SetStateToDisabled();

    mMenuReorthogonalizeVolume->SetStateToDisabled();
  }

  if( mVolume.GetPointer() && mbDirty ) {
    mBtnSaveVolume->SetStateToNormal();
    mMenuSaveVolume->SetStateToNormal();
    mMenuSaveVolumeAs->SetStateToNormal();
    mBtnRevertVolume->SetStateToNormal();
    mMenuRevertVolume->SetStateToNormal();
  } else {
    mBtnSaveVolume->SetStateToDisabled();
    mMenuSaveVolume->SetStateToDisabled();
    mMenuSaveVolumeAs->SetStateToDisabled();
    mBtnRevertVolume->SetStateToDisabled();
    mMenuRevertVolume->SetStateToDisabled();
  }

  if( mVolume.GetPointer() && mGrayScaleColors.GetPointer() ) {
    mRadBtnUseGrayScaleColors->SetStateToNormal();
    mMenuUseGrayScaleColors->SetStateToNormal();
  } else {
    mRadBtnUseGrayScaleColors->SetStateToDisabled();
    mMenuUseGrayScaleColors->SetStateToDisabled();
  }

  if( mVolume.GetPointer() && mLUTColors.GetPointer() ) {
    mRadBtnUseLUTColors->SetStateToNormal();
    mMenuUseLUTColors->SetStateToNormal();
  } else {
    mRadBtnUseLUTColors->SetStateToDisabled();
    mMenuUseLUTColors->SetStateToDisabled();
  }
}

void
vtkKWOrientMRIWindow::RebuildVolumeToRASTransform () {

  assert( mVolume.GetPointer() );
  assert( mVolumeToRASTransform.GetPointer() );
  assert( mScaleThroughPlane[X].GetPointer() );
  assert( mScaleThroughPlane[Y].GetPointer() );
  assert( mScaleThroughPlane[Z].GetPointer() );
  
  // This rotates the volume to the proper orientation. From
  // ImageReslice: "applying a transform to the resampling grid (which
  // lies in the output coordinate system) is equivalent to applying
  // the inverse of that transform to the input volume."
  double* rtv = mVolume->GetRASToVoxelMatrix();

  // 0,0 0,1 0,2 0,3      rtv[0]  rtv[1]  rtv[2]  0
  // 1,0 1,1 1,2 1,3  =>  rtv[4]  rtv[5]  rtv[6]  0
  // 2,0 2,1 2,2 2,3      rtv[8]  rtv[9]  rtv[10] 0
  // 3,0 3,1 3,2 3,3        0       0       0     1
  vtkSmartPointer<vtkMatrix4x4> matrix = 
    vtkSmartPointer<vtkMatrix4x4>::New();
  matrix->SetElement( 0, 0, rtv[0] );
  matrix->SetElement( 0, 1, rtv[1] );
  matrix->SetElement( 0, 2, rtv[2] );
  matrix->SetElement( 0, 3, 0 );
  matrix->SetElement( 1, 0, rtv[4] );
  matrix->SetElement( 1, 1, rtv[5] );
  matrix->SetElement( 1, 2, rtv[6] );
  matrix->SetElement( 1, 3, 0 );
  matrix->SetElement( 2, 0, rtv[8] );
  matrix->SetElement( 2, 1, rtv[9] );
  matrix->SetElement( 2, 2, rtv[10] );
  matrix->SetElement( 2, 3, 0 );
  matrix->SetElement( 3, 0, 0 );
  matrix->SetElement( 3, 1, 0 );
  matrix->SetElement( 3, 2, 0 );
  matrix->SetElement( 3, 3, 1 );

  vtkSmartPointer<vtkTransform> transform = 
    vtkSmartPointer<vtkTransform>::New();
  transform->SetMatrix( matrix );

  // Set the transform in our reslicer.
  mVolumeToRASTransform->SetResliceTransform( transform );
  mVolumeToRASTransform->Update();
  
  // Get the RAS bounds from the volume and set the ranges of the
  // scales.
  float RASBounds[6];
  mVolume->GetRASBounds( RASBounds );
  mScaleThroughPlane[X]->SetRange( RASBounds[0], RASBounds[1] );
  mScaleThroughPlane[X]->SetValue( RASBounds[0] +
				   (RASBounds[1]-RASBounds[0])/2.0 );
  mScaleThroughPlane[Y]->SetRange( RASBounds[2], RASBounds[3] );
  mScaleThroughPlane[Y]->SetValue( RASBounds[2] +
				   (RASBounds[3]-RASBounds[2])/2.0 );
  mScaleThroughPlane[Z]->SetRange( RASBounds[4], RASBounds[5] );
  mScaleThroughPlane[Z]->SetValue( RASBounds[4] + 
				   (RASBounds[5]-RASBounds[4])/2.0 );
  

#if 0
  int* extent;
  double* bounds;
  cerr << "RAS bounds " << RASBounds[0] << " " << RASBounds[1] << " "
       << RASBounds[2] << " " << RASBounds[3] << " "
       << RASBounds[4] << " " << RASBounds[5] << " " << endl;
  extent = mVolumeToRASTransform->GetOutput()->GetWholeExtent();
  cerr << "in window, vtr output whole extent "
       << extent[0] << " " << extent[1] << " "
       << extent[2] << " " << extent[3] << " "
       << extent[4] << " " << extent[5] << " " << endl;
  extent = mVolumeToRASTransform->GetOutput()->GetExtent();
  cerr << "in window, vtr output extent "
       << extent[0] << " " << extent[1] << " "
       << extent[2] << " " << extent[3] << " "
       << extent[4] << " " << extent[5] << " " << endl;
  extent = mVolumeToRASTransform->GetOutputExtent();
  cerr << "in window, vtr OutputExtent "
       << extent[0] << " " << extent[1] << " "
       << extent[2] << " " << extent[3] << " "
       << extent[4] << " " << extent[5] << " " << endl;
  bounds = mVolumeToRASTransform->GetOutput()->GetBounds();
  cerr << "in window, vtr output bounds "
       << bounds[0] << " " << bounds[1] << " "
       << bounds[2] << " " << bounds[3] << " "
       << bounds[4] << " " << bounds[5] << " " << endl;
  extent = mUserTransform->GetOutput()->GetExtent();
  cerr << "in window, user extent " << extent[0] << " " << extent[1] << " "
      << extent[2] << " " << extent[3] << " "
      << extent[4] << " " << extent[5] << " " << endl;

  bounds = mUserTransform->GetOutput()->GetBounds();
  cerr << "in window, user bounds " << bounds[0] << " " << bounds[1] << " "
      << bounds[2] << " " << bounds[3] << " "
      << bounds[4] << " " << bounds[5] << " " << endl;
#endif
  
  // Broadcast our event.
  this->InvokeEvent( OrientMRIEvents::VolumeToRASTransformChanged );
  
  // Restore the view.
  this->RestoreView();
  
  // Now dirty.
  mbDirty = true;
  
  // Update our menu and buttons.
  this->UpdateCommandStatus();
 
}

vtkMatrix4x4&
vtkKWOrientMRIWindow::GetUserTransform () {

  assert( mView3D.GetPointer() );

  // Just get it from the 3D view.
  return mView3D->GetUserTransform();
}


void
vtkKWOrientMRIWindow::AddSpacerToToolbar ( vtkKWToolbar& iToolbar, int iWidth ) {
  
  // Create a frame of the specified width inside the toolbar.
  vtkSmartPointer<vtkKWFrame> spacer = 
    vtkSmartPointer<vtkKWFrame>::New();
  spacer->SetParent( &iToolbar );
  spacer->Create();
  spacer->SetWidth( iWidth );
  iToolbar.AddWidget( spacer );

}

void
vtkKWOrientMRIWindow::UserTransformChanged ( vtkObject* iCaller, 
					     unsigned long iEventId,
					     void* iClientData,
					     void* iCallData ) {

  // We're listening for UserTransformChanged events. When we get one,
  // call the window's UpdateUserTransform function.
  if( OrientMRIEvents::UserTransformChanged == iEventId ) {

    // Get our window pointer from the client data.
    assert( iClientData );
    vtkKWOrientMRIWindow* window = 
      static_cast<vtkKWOrientMRIWindow*>( iClientData );
    
    window->UpdateUserTransform();
  }

}

void
vtkKWOrientMRIWindow::DoRotateUserTransformDialog ( int iAxis, 
						    int iMultiplier ) {
  
  assert( mView3D );

  // Default amount to rotate by.
  static double degrees = 90;

  // Make our dialog.
  vtkSmartPointer<vtkKWSimpleEntryDialog> dialog =
    vtkSmartPointer<vtkKWSimpleEntryDialog>::New();
  dialog->SetApplication( this->GetApplication() );
  dialog->Create();
  dialog->SetOptions( vtkKWMessageDialog::OkDefault );
  dialog->SetTitle( "Rotate" );

  // Configure the entry so that it only accepts doubles. Fill it in
  // with the amount of degrees.
  vtkSmartPointer<vtkKWEntryWithLabel> labeledEntry = dialog->GetEntry();
  labeledEntry->SetLabelText( "Degrees to rotate: " );
  labeledEntry->GetWidget()->SetRestrictValueToDouble();
  labeledEntry->GetWidget()->SetValueAsDouble( degrees );
  labeledEntry->GetWidget()->SetWidth( 5 );
  labeledEntry->GetWidget()->Focus();

  if( dialog->Invoke() ) {
    // Get the new value and rotate the transform.
    degrees = labeledEntry->GetWidget()->GetValueAsDouble();
    mView3D->RotateUserTransform( iAxis, degrees * (double)iMultiplier );
  }
  
}

vtkKWOrientMRIWindow::MenuItem::MenuItem ()  :
    mMenu( NULL ),
    mnItem( -1 ) {}


void
vtkKWOrientMRIWindow::MenuItem::MakeCommand ( vtkKWMenu* iMenu,
					 int inItem,
					 const char* isText,
					 vtkObject* iCommandObject,
					 const char* isCommand,
					 const char* isAccelerator,
					 const char* isIconKey ) {
  
  if( NULL == iMenu )
    throw( "MakeCommand: iMenu was NULL" );
  if( NULL == isText )
    throw( "MakeCommand: isText was NULL" );
  if( NULL == isCommand )
    throw( "MakeCommand: isCommand was NULL" );

  // Try inserting the command. If it returns -1 it failed.
  int nIndex =
    iMenu->InsertCommand( inItem, isText, iCommandObject, isCommand );

  if( -1 == nIndex )
    throw runtime_error( "Couldn't create menu item" );

  // It's in. Save the info.
  mMenu = iMenu;
  mnItem = nIndex;

  // Set the accelerator if available.
  if( isAccelerator )
    mMenu->SetItemAccelerator( mnItem, isAccelerator );

  // Try to load the icon. If it fails, use the KW built-in question
  // mark icon.
  if( isIconKey ) {
    try {
      mMenu->SetItemCompoundModeToLeft( mnItem );
      IconLoader::SetMenuItemIcon( isIconKey, mMenu, mnItem );
    } catch (exception& e) {
      stringstream ssError;
      ssError << "Error loading icon: " << e.what();
      mMenu->GetApplication()->ErrorMessage( ssError.str().c_str() );
      mMenu->SetItemImageToPredefinedIcon( mnItem, vtkKWIcon::IconQuestion );
    }
  }

}

void
vtkKWOrientMRIWindow::MenuItem::MakeCheckButton ( vtkKWMenu* iMenu,
					     int inItem,
					     const char* isText,
					     vtkObject* iCommandObject,
					     const char* isCommand,
					     const char* isAccelerator,
					     const char* isIconKey ) {
  
  if( NULL == iMenu )
    throw( "MakeCheckButton: iMenu was NULL" );
  if( NULL == isText )
    throw( "MakeCheckButton: isText was NULL" );
  if( NULL == isCommand )
    throw( "MakeCheckButton: isCommand was NULL" );

  // Try inserting the check button. If it returns -1 it failed.
  int nIndex =
    iMenu->InsertCheckButton( inItem, isText, iCommandObject, isCommand );

  if( -1 == nIndex )
    throw runtime_error( "Couldn't create menu item" );

  // It's in. Save the info.
  mMenu = iMenu;
  mnItem = nIndex;

  // Set the accelerator if available.
  if( isAccelerator )
    mMenu->SetItemAccelerator( mnItem, isAccelerator );

  // Try to load the icon. If it fails, use the KW built-in question
  // mark icon.
  if( isIconKey ) {
    try {
      mMenu->SetItemCompoundModeToLeft( mnItem );
      IconLoader::SetMenuItemIcon( isIconKey, mMenu, mnItem );
    } catch (exception& e) {
      stringstream ssError;
      ssError << "Error loading icon: " << e.what();
      mMenu->GetApplication()->ErrorMessage( ssError.str().c_str() );
      mMenu->SetItemImageToPredefinedIcon( mnItem, vtkKWIcon::IconQuestion );
    }
  }

}

void
vtkKWOrientMRIWindow::MenuItem::SetStateToDisabled () {

  if( NULL == mMenu.GetPointer() )
    throw runtime_error( "MenuItem::SetStateToDisabled: mMenu was NULL" );
  if( -1 == mnItem )
    throw runtime_error( "MenuItem::SetStateToDisabled: mnItem was -1" );

  mMenu->SetItemStateToDisabled( mnItem );
}

void
vtkKWOrientMRIWindow::MenuItem::SetStateToNormal () {

  if( NULL == mMenu.GetPointer() )
    throw runtime_error( "MenuItem::SetStateToNormal: mMenu was NULL" );
  if( -1 == mnItem )
    throw runtime_error( "MenuItem::SetStateToNormal: mnItem was -1" );

  mMenu->SetItemStateToNormal( mnItem );
}

void
vtkKWOrientMRIWindow::MenuItem::SetSelectedState ( int ibOn ) {

  if( NULL == mMenu.GetPointer() )
    throw runtime_error( "MenuItem::SetSelectedState: mMenu was NULL" );
  if( -1 == mnItem )
    throw runtime_error( "MenuItem::SetSelectedState: mnItem was -1" );

  mMenu->SetItemSelectedState( mnItem, ibOn );
}

int
vtkKWOrientMRIWindow::MenuItem::GetSelectedState () {

  if( NULL == mMenu.GetPointer() )
    throw runtime_error( "MenuItem::GetSelectedState: mMenu was NULL" );
  if( -1 == mnItem )
    throw runtime_error( "MenuItem::GetSelectedState: mnItem was -1" );

  return mMenu->GetItemSelectedState( mnItem );
}

int
vtkKWOrientMRIWindow::MenuItem::GetIndex () const {

  return mnItem;
}
