/**
 * @file  vtkKWScubaWindow.cxx
 * @brief The main window
 *
 * Manages the 'current' tool, view, and layer. Holds one or more
 * views and lays them out. Top level document loading and creation of
 * layers. Runs the info area and gets updates from mouseovered
 * objects. Runs the toolbars. Manages the UI panels for settings.
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

#include <assert.h>
#include <string>
#include "vtkKWScubaWindow.h"
#include "IconLoader.h"
#include "vtkKWApplication.h"
#include "vtkKWFrameWithLabel.h"
#include "vtkKWIcon.h"
#include "vtkKWLoadSaveDialog.h"
#include "vtkKWMenu.h"
#include "vtkKWMenuButton.h"
#include "vtkKWMenuButtonWithSpinButtons.h"
#include "vtkKWMenuButtonWithSpinButtonsWithLabel.h"
#include "vtkKWMessageDialog.h"
#include "vtkKWMultiColumnList.h"
#include "vtkKWPushButton.h"
#include "vtkKWRadioButton.h"
#include "vtkKWRadioButtonSet.h"
#include "vtkKWScaleWithEntry.h"
#include "vtkKWScubaApplicationSettingsInterface.h"
#include "vtkKWScubaLayerCollectionDTI.h"
#include "vtkKWScubaLayerCollectionMRI.h"
#include "vtkKWScubaLayerCollectionMRIS.h"
#include "vtkKWScubaLayerCollectionODF.h"
#include "vtkKWScubaLayerCollectionPath.h"
#include "vtkKWScubaToolEdit2DMRI.h"
#include "vtkKWScubaToolNavigate.h"
#include "vtkKWSplitFrame.h"
#include "vtkKWToolbar.h"
#include "vtkKWToolbarSet.h"
#include "vtkKWUserInterfacePanel.h"
#include "vtkObjectFactory.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkScubaInteractorStyle.h"

using namespace std;

vtkStandardNewMacro( vtkKWScubaWindow );
vtkCxxRevisionMacro( vtkKWScubaWindow, "$Revision: 1.9 $" );

int const vtkKWScubaWindow::kToolbarSpacerWidth = 5;
const string vtkKWScubaWindow::sDefaultVolumeFileExtension = ".mgz";
const char* vtkKWScubaWindow::sRegistryKey = "WindowSettings";
const char* vtkKWScubaWindow::sAutoSizeInfoAreaKey = "AutoSizeInfoArea";


vtkKWScubaWindow::vtkKWScubaWindow () :
    Listener( "vtkKWScubaWindow" ),
    mMenuLoadVolume( NULL ),
    mMenuLoadSurface( NULL ),
    mMenuLoadDTI( NULL ),
    mMenuLoadPath( NULL ),
    mMenuLoadODF( NULL ),
    mMenuSaveVolume( NULL ),
    mMenuZoomOut( NULL ),
    mMenuZoomIn( NULL ),
    mCurrentViewLayout( NoViewLayout ),
    mbAutoSizeInfoArea( 1 ) {

}

vtkKWScubaWindow::~vtkKWScubaWindow () {
}

void
vtkKWScubaWindow::CreateWidget () {

  // Restore our preferences. This needs to be done before we call the
  // superclass Create because that will instantiate our
  // vtkKWScubaApplicationSettingsInterface, which will ask us for
  // values.
  if( this->GetApplication()->
      HasRegistryValue( 2, sRegistryKey, sAutoSizeInfoAreaKey ) )
    mbAutoSizeInfoArea = this->GetApplication()->
      GetIntRegistryValue( 2, sRegistryKey, sAutoSizeInfoAreaKey );

  // Our default view layout.
  this->SetPanelLayoutToSecondaryBelowMainAndView();

  // Create the window. The layout must be set before this.
  this->Superclass::CreateWidget();

  // Get our geometry. This needs to be done after the superclass
  // Create because we can't set geometry before that.
  this->RestoreWindowGeometryFromRegistry();

  // Turn on our panels.
  this->SecondaryPanelVisibilityOn();
  this->MainPanelVisibilityOn();

  // Start out with a small secondary frame. We'll expand it when we
  // get items to put into it.
  this->GetSecondarySplitFrame()->SetFrame1Size( 0 );

  //
  // Create the secondary panel. This is our 'info area' and contains
  // two tables that we populate with data from the view and
  // layers. Note that we don't use the secondary UserInterface
  // Manager, we just pack these things right into the frame. This is
  // because we don't need any kind of paged interface here.
  mCursorInfoTable = vtkSmartPointer<vtkKWMultiColumnList>::New();
  mCursorInfoTable->SetParent( this->GetSecondaryPanelFrame() );
  mCursorInfoTable->Create();
  mCursorInfoTable->ColumnLabelsVisibilityOff();
  mCursorInfoTable->SetWidth( 1 );
  mCursorInfoTable->AddColumn( "" );
  mCursorInfoTable->ColumnStretchableOff( 0 );
  mCursorInfoTable->AddColumn( "" );
  mCursorInfoTable->ColumnStretchableOn( 1 );

  mMouseOverInfoTable = vtkSmartPointer<vtkKWMultiColumnList>::New();
  mMouseOverInfoTable->SetParent( this->GetSecondaryPanelFrame() );
  mMouseOverInfoTable->Create();
  mMouseOverInfoTable->ColumnLabelsVisibilityOff();
  mMouseOverInfoTable->SetWidth( 1 );
  mMouseOverInfoTable->AddColumn( "" );
  mMouseOverInfoTable->ColumnStretchableOff( 0 );
  mMouseOverInfoTable->AddColumn( "" );
  mMouseOverInfoTable->ColumnStretchableOn( 1 );

  this->Script( "pack %s %s -side left -expand yes -fill both",
                mCursorInfoTable->GetWidgetName(),
                mMouseOverInfoTable->GetWidgetName() );

  //
  // Create the main panel. We use the interface manager here which
  // implements basically a tabbed notebook. Then we create pages to
  // go in the manager, and get the frame for them. We'll pack all the
  // labeled frames inside those. One page for Tools (gets tools
  // frame) and one page for Display (gets views and layers frames).
  //

  vtkKWUserInterfacePanel* panel = vtkKWUserInterfacePanel::New();
  panel->SetUserInterfaceManager( this->GetMainUserInterfaceManager() );
  panel->Create();
  panel->AddPage( "Tool", "Tool page", NULL );
  panel->AddPage( "Display", "Display page", NULL );

  //
  // Start in the tool frame.
  vtkKWWidget* panelFrame = panel->GetPageWidget( "Tool" );

  //
  // The tool frame.
  vtkSmartPointer<vtkKWFrameWithLabel> toolFrame = 
    vtkSmartPointer<vtkKWFrameWithLabel>::New();
  toolFrame->SetParent( panelFrame );
  toolFrame->Create();
  toolFrame->SetLabelText( "Tool" );
  this->Script( "pack %s -side top -fill x -anchor nw",
		toolFrame->GetWidgetName() );

  // Menu for tools.
  vtkSmartPointer<vtkKWMenuButtonWithSpinButtonsWithLabel> labeledMenu =
    vtkSmartPointer<vtkKWMenuButtonWithSpinButtonsWithLabel>::New();
  labeledMenu->SetParent( toolFrame->GetFrame() );
  labeledMenu->Create();
  labeledMenu->SetLabelText( "Tool: " );

  mMenuTool = labeledMenu->GetWidget()->GetWidget();

  // A subframe for the tool.
  mToolUIFrame = vtkSmartPointer<vtkKWFrame>::New();
  mToolUIFrame->SetParent( toolFrame->GetFrame() );
  mToolUIFrame->Create();

  this->Script( "pack %s %s -side top -fill x -anchor nw",
                labeledMenu->GetWidgetName(),
                mToolUIFrame->GetWidgetName() );

  //
  // Now start using using the display page.
  panelFrame = panel->GetPageWidget( "Display" );

  //
  // Create the View frame. Make our view-selecting menu for it,
  // and then another frame underneath that. The views will populate
  // that second frame. Leave it blank for now.
  vtkSmartPointer<vtkKWFrameWithLabel> viewFrame =
    vtkSmartPointer<vtkKWFrameWithLabel>::New();
  viewFrame->SetParent( panelFrame );
  viewFrame->Create();
  viewFrame->SetLabelText( "View" );
  this->Script( "pack %s -side top -fill x -anchor nw",
		viewFrame->GetWidgetName() );

  labeledMenu = 
    vtkSmartPointer<vtkKWMenuButtonWithSpinButtonsWithLabel>::New();
  labeledMenu->SetParent( viewFrame->GetFrame() );
  labeledMenu->Create();
  labeledMenu->SetLabelText( "View: " );

  mMenuView = labeledMenu->GetWidget()->GetWidget();

  mViewUIFrame = vtkSmartPointer<vtkKWFrame>::New();
  mViewUIFrame->SetParent( viewFrame->GetFrame() );
  mViewUIFrame->Create();

  this->Script( "pack %s %s -side top -fill x -anchor nw -pady 5",
                labeledMenu->GetWidgetName(),
                mViewUIFrame->GetWidgetName() );


  //
  // Create the Layer frame. Make our layer-selecting menu for it,
  // and then another frame underneath that. The layers will populate
  // that second frame. Leave it blank for now.
  vtkSmartPointer<vtkKWFrameWithLabel> layerFrame =
    vtkSmartPointer<vtkKWFrameWithLabel>::New();
  layerFrame->SetParent( panelFrame );
  layerFrame->Create();
  layerFrame->SetLabelText( "Layer" );
  this->Script( "pack %s -side top -fill x -anchor nw",
		layerFrame->GetWidgetName() );

  labeledMenu = 
    vtkSmartPointer<vtkKWMenuButtonWithSpinButtonsWithLabel>::New();
  labeledMenu->SetParent( layerFrame->GetFrame() );
  labeledMenu->Create();
  labeledMenu->SetLabelText( "Layer: " );

  mMenuLayer = labeledMenu->GetWidget()->GetWidget();

  mLayerUIFrame = vtkSmartPointer<vtkKWFrame>::New();
  mLayerUIFrame->SetParent( layerFrame->GetFrame() );
  mLayerUIFrame->Create();

  this->Script( "pack %s %s -side top -fill x -anchor nw -pady 5",
                labeledMenu->GetWidgetName(),
                mLayerUIFrame->GetWidgetName() );

  // ====================================================================

  //
  // Make our tool toolbar.
  mToolbarTools = vtkSmartPointer<vtkKWToolbar>::New();
  mToolbarTools->SetName( "Tools" );
  mToolbarTools->SetParent( this->GetMainToolbarSet()->GetToolbarsFrame() );
  mToolbarTools->Create();
  this->GetMainToolbarSet()->AddToolbar( mToolbarTools );
  
  //
  // Make our window command toolbar.
  mToolbarWindow = vtkSmartPointer<vtkKWToolbar>::New();
  mToolbarWindow->SetName( "Main" );
  mToolbarWindow->SetParent( this->GetMainToolbarSet()->GetToolbarsFrame() );
  mToolbarWindow->Create();
  this->GetMainToolbarSet()->AddToolbar( mToolbarWindow );
  
  //
  // Make the view toolbar, we'll fill it in later.
  mToolbarView = vtkSmartPointer<vtkKWToolbar>::New();
  mToolbarView->SetName( "View" );
  mToolbarView->SetParent( this->GetMainToolbarSet()->GetToolbarsFrame() );
  mToolbarView->Create();
  this->GetMainToolbarSet()->AddToolbar( mToolbarView );

  // Make our buttons. For each one, we create a button, set the
  // parent, give it some basic characterists, an icon, and set its
  // command. We finish by adding it to the toolbar, which will pack
  // it for us.

  // Load Volume
  mBtnLoadVolume = vtkSmartPointer<vtkKWPushButton>::New();
  mBtnLoadVolume->SetParent( mToolbarWindow->GetFrame() );
  mBtnLoadVolume->Create();
  mBtnLoadVolume->SetText( "Load Volume" );
  mBtnLoadVolume->SetBalloonHelpString( "Load Volume" );
  mBtnLoadVolume->SetImageToPredefinedIcon( vtkKWIcon::IconFileOpen );
  mBtnLoadVolume->SetCommand( this, "LoadVolumeFromDlog" );
  try { IconLoader::SetPushButtonIcon( "LoadVolume", mBtnLoadVolume ); }
  catch( exception& e ) { cerr << "Error loading icon: " << e.what() << endl; }
  mToolbarWindow->AddWidget( mBtnLoadVolume );
    
  // Save Volume
  mBtnSaveVolume = vtkSmartPointer<vtkKWPushButton>::New();
  mBtnSaveVolume->SetParent( mToolbarWindow->GetFrame() );
  mBtnSaveVolume->Create();
  mBtnSaveVolume->SetText( "Save Volume" );
  mBtnSaveVolume->SetBalloonHelpString( "Save Volume" );
  mBtnSaveVolume->SetImageToPredefinedIcon( vtkKWIcon::IconFloppy );
  mBtnSaveVolume->SetCommand( this, "SaveVolumeWithConfirm" );
  try { IconLoader::SetPushButtonIcon( "SaveVolume", mBtnSaveVolume ); }
  catch( exception& e ) { cerr << "Error loading icon: " << e.what() << endl; }
  mToolbarWindow->AddWidget( mBtnSaveVolume );

  // View Layout radio button set.
  mRadBtnSetViewLayout = vtkSmartPointer<vtkKWRadioButtonSet>::New();
  mRadBtnSetViewLayout->SetParent( mToolbarWindow->GetFrame() );
  mRadBtnSetViewLayout->Create();
  mRadBtnSetViewLayout->PackHorizontallyOn();
  mToolbarWindow->AddWidget( mRadBtnSetViewLayout );

  vtkSmartPointer<vtkKWRadioButton> radBtnViewLayoutSingle;
  radBtnViewLayoutSingle.
    TakeReference( mRadBtnSetViewLayout->AddWidget( (int)Single ) );
  radBtnViewLayoutSingle->SetCommand( this, "SetViewLayoutToSingle" );
  radBtnViewLayoutSingle->IndicatorVisibilityOff();
  try { IconLoader::SetCheckButtonIcon("ViewSingle",radBtnViewLayoutSingle); }
  catch( exception& e ) { cerr << "Error loading icon: " << e.what() << endl; }
  if( Single == mCurrentViewLayout )
    radBtnViewLayoutSingle->SelectedStateOn();
  
  vtkSmartPointer<vtkKWRadioButton> radBtnViewLayoutTwoByTwo;
  radBtnViewLayoutTwoByTwo.
    TakeReference( mRadBtnSetViewLayout->AddWidget( (int)TwoByTwo ) );
  radBtnViewLayoutTwoByTwo->SetCommand( this, "SetViewLayoutToTwoByTwo" );
  radBtnViewLayoutTwoByTwo->IndicatorVisibilityOff();
  try { IconLoader::SetCheckButtonIcon("ViewTwoByTwo",radBtnViewLayoutTwoByTwo); }
  catch( exception& e ) { cerr << "Error loading icon: " << e.what() << endl; }
  if( TwoByTwo == mCurrentViewLayout )
    radBtnViewLayoutTwoByTwo->SelectedStateOn();

  //
  // Build the menus. File menu.

  int nItem = this->GetFileMenuInsertPosition();

  // Load Volume.
  mMenuLoadVolume = new MenuItem();
  mMenuLoadVolume->
    MakeCommand( this->GetFileMenu(), nItem++,
		 "L&oad Volume...", this, "LoadVolumeFromDlog",
		 "Ctrl+O", "LoadVolume" );

  // Load Surface.
  mMenuLoadSurface = new MenuItem();
  mMenuLoadSurface->
    MakeCommand( this->GetFileMenu(), nItem++,
		 "Load Surface...", this, "LoadSurfaceFromDlog",
		 "", "LoadSurface" );

  // Load DTI.
  mMenuLoadDTI = new MenuItem();
  mMenuLoadDTI->
    MakeCommand( this->GetFileMenu(), nItem++,
		 "Load DTI...", this, "LoadDTIFromDlog",
		 "", "LoadDTI" );
  
  // Load Path
  mMenuLoadPath = new MenuItem();
  mMenuLoadPath->
    MakeCommand( this->GetFileMenu(), nItem++,
		 "Load Path...", this, "LoadPathFromDlog",
		 "", "LoadPath" );
  
  // Load ODF
  mMenuLoadODF = new MenuItem();
  mMenuLoadODF->
    MakeCommand( this->GetFileMenu(), nItem++,
		 "Load ODF...", this, "LoadODFFromDlog",
		 "", "LoadODF" );

  this->GetFileMenu()->InsertSeparator( nItem++ );

  // Save Volume.
  mMenuSaveVolume = new MenuItem();
  mMenuSaveVolume->
    MakeCommand( this->GetFileMenu(), nItem++,
		 "&Save Volume...", this, "SaveVolumeWithConfirm",
		 "Ctrl+S", "SaveVolume" );

  this->GetFileMenu()->InsertSeparator( nItem++ );

  this->InsertRecentFilesMenu( nItem++, this );

  this->GetFileMenu()->InsertSeparator( nItem++ );

  // View menu.
  nItem = this->GetViewMenuInsertPosition();

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

  //
  // Load initial tools. This will just create them and then we will
  // get them later with the IDTracker methods when we build our
  // menu. Don't Delete them tho, as there are no other references to
  // them.
  vtkSmartPointer<vtkKWScubaTool> navigateTool = 
    vtkSmartPointer<vtkKWScubaToolNavigate>::New();
  navigateTool->SetApplication( this->GetApplication() );

  vtkSmartPointer<vtkKWScubaTool> tool = 
    vtkSmartPointer<vtkKWScubaToolEdit2DMRI>::New();
  tool->SetApplication( this->GetApplication() );

  // Make a radio button set for our toolbar.
  mRadBtnSetTool = vtkSmartPointer<vtkKWRadioButtonSet>::New();
  mRadBtnSetTool->SetParent( mToolbarWindow->GetFrame() );
  mRadBtnSetTool->Create();
  
  // Build the menu and make one too the current one for the initial
  // selection.
  this->SetCurrentTool( *navigateTool );
  this->UpdateToolMenu();

  
  //
  // Make our initial view setup.
  this->SetViewLayoutToSingle();
  this->InitializeViewSettingsForLayout();

  //
  // Update our menu and buttons.
  this->UpdateCommandStatus();
}

void
vtkKWScubaWindow::PrepareForDelete () {

  // Save our geometry.
  this->SaveWindowGeometryToRegistry();
}

vtkKWApplicationSettingsInterface*
vtkKWScubaWindow::GetApplicationSettingsInterface() {

  // If we don't have one already...
  if( NULL == ApplicationSettingsInterface ) {

    // Create our custom interface.
    vtkSmartPointer<vtkKWScubaApplicationSettingsInterface> settings = 
      vtkSmartPointer<vtkKWScubaApplicationSettingsInterface>::New();
    settings->SetWindow( this );
    settings->SetScubaWindow( this );
    settings->SetUserInterfaceManager
      ( this->GetApplicationSettingsUserInterfaceManager() );

    ApplicationSettingsInterface = settings;
  }
  
  return ApplicationSettingsInterface;
}

void
vtkKWScubaWindow::LoadVolumeFromDlog () {

  // Create a Load dialog and set it up.
  // The smart pointer version is commented out because of a bug that
  // occurs when the dialog is deleted.
  //  vtkSmartPointer<vtkKWLoadSaveDialog> dialog = 
  //    vtkSmartPointer<vtkKWLoadSaveDialog>::New();
  vtkKWLoadSaveDialog* dialog = vtkKWLoadSaveDialog::New();
  dialog->SetApplication( this->GetApplication() );
  dialog->Create();
  dialog->SetTitle( "Load a volume" );
  dialog->SetFileTypes( "{MGH {.mgh .mgz}} "
                        "{Binary {.bshort .bfloat}} {All {*}}" );
  dialog->RetrieveLastPathFromRegistry( "LoadVolume" );
  dialog->SetDefaultExtension( sDefaultVolumeFileExtension.c_str() );

  // Show the dialog, and when it returns, Invoke() will be true if
  // they clicked OK and gave us a filename.
  if ( dialog->Invoke() ) {
    dialog->SaveLastPathToRegistry( "LoadVolume" );
    string fnVolume( dialog->GetFileName() );
    this->LoadVolume( fnVolume.c_str() );
  }
}

void
vtkKWScubaWindow::LoadSurfaceFromDlog () {

  // Create a Load dialog and set it up.
  // The smart pointer version is commented out because of a bug that
  // occurs when the dialog is deleted.
  //  vtkSmartPointer<vtkKWLoadSaveDialog> dialog = 
  //    vtkSmartPointer<vtkKWLoadSaveDialog>::New();
  vtkKWLoadSaveDialog* dialog = vtkKWLoadSaveDialog::New();
  dialog->SetApplication( this->GetApplication() );
  dialog->Create();
  dialog->SetTitle( "Load a surface" );
  dialog->SetFileTypes( "{\"Known surfaces\" {*.inflated* *.orig* *sphere* *.pial* *.smooth* *.white}} "
                        "{\"All surfaces\" {lh.* rh.*}} {All {*}}" );
  dialog->RetrieveLastPathFromRegistry( "LoadSurface" );

  // Show the dialog, and when it returns, Invoke() will be true if
  // they clicked OK and gave us a filename.
  if ( dialog->Invoke() ) {
    dialog->SaveLastPathToRegistry( "LoadSurface" );
    string fnSurface( dialog->GetFileName() );
    this->LoadSurface( fnSurface.c_str() );
  }
}

void
vtkKWScubaWindow::LoadDTIFromDlog () {

  // Create a Load dialog and set it up.
  // The smart pointer version is commented out because of a bug that
  // occurs when the dialog is deleted.
  //  vtkSmartPointer<vtkKWLoadSaveDialog> dialog = 
  //    vtkSmartPointer<vtkKWLoadSaveDialog>::New();
  vtkKWLoadSaveDialog* dialog = vtkKWLoadSaveDialog::New();
  dialog->SetApplication( this->GetApplication() );
  dialog->Create();
  dialog->SetTitle( "Load a DTI" );
  dialog->SetFileTypes( "{MGH {.mgh .mgz}} "
                        "{Binary {.bshort .bfloat}} {All {*}}" );
  dialog->RetrieveLastPathFromRegistry( "LoadDTI" );
  dialog->SetDefaultExtension( sDefaultVolumeFileExtension.c_str() );

  // Show the dialog, and when it returns, Invoke() will be true if
  // they clicked OK and gave us a filename.
  if ( dialog->Invoke() ) {
    dialog->SaveLastPathToRegistry( "LoadDTI" );
    string fnDTI( dialog->GetFileName() );
    this->LoadDTI( fnDTI.c_str() );
  }
}

void
vtkKWScubaWindow::LoadPathFromDlog () {

  // Create a Load dialog and set it up.
  // The smart pointer version is commented out because of a bug that
  // occurs when the dialog is deleted.
  //  vtkSmartPointer<vtkKWLoadSaveDialog> dialog = 
  //    vtkSmartPointer<vtkKWLoadSaveDialog>::New();
  vtkKWLoadSaveDialog* dialog = vtkKWLoadSaveDialog::New();
  dialog->SetApplication( this->GetApplication() );
  dialog->Create();
  dialog->SetTitle( "Load a path" );
  dialog->SetFileTypes( "{MGH {.mgh .mgz}} "
                        "{Binary {.bshort .bfloat}} {All {*}}" );
  dialog->RetrieveLastPathFromRegistry( "LoadPath" );
  dialog->SetDefaultExtension( sDefaultVolumeFileExtension.c_str() );

  // Show the dialog, and when it returns, Invoke() will be true if
  // they clicked OK and gave us a filename.
  if ( dialog->Invoke() ) {
    dialog->SaveLastPathToRegistry( "LoadPath" );
    string fnPath( dialog->GetFileName() );
    this->LoadPath( fnPath.c_str() );
  }
  
}

void
vtkKWScubaWindow::LoadODFFromDlog () {

  // Create a Load dialog and set it up.
  // The smart pointer version is commented out because of a bug that
  // occurs when the dialog is deleted.
  //  vtkSmartPointer<vtkKWLoadSaveDialog> dialog = 
  //    vtkSmartPointer<vtkKWLoadSaveDialog>::New();
  vtkKWLoadSaveDialog* dialog = vtkKWLoadSaveDialog::New();
  dialog->SetApplication( this->GetApplication() );
  dialog->Create();
  dialog->SetTitle( "Load ODF" );
  dialog->SetFileTypes( "{MGH {.mgh .mgz}} "
                        "{Binary {.bshort .bfloat}} {All {*}}" );
  dialog->RetrieveLastPathFromRegistry( "LoadODF" );
  dialog->SetDefaultExtension( sDefaultVolumeFileExtension.c_str() );

  // Show the dialog, and when it returns, Invoke() will be true if
  // they clicked OK and gave us a filename.
  if ( dialog->Invoke() ) {
    dialog->SaveLastPathToRegistry( "LoadODF" );
    string fnPath( dialog->GetFileName() );
    this->LoadODF( fnPath.c_str() );
  }
  
}

void
vtkKWScubaWindow::SaveVolumeWithConfirm () {

  // If the view has dirty data...
  if ( mCurrentView.GetPointer() && mCurrentView->IsDataDirty() ) {

    // Bring up a yes/no dialog, and if they click Yes...
    if ( vtkKWMessageDialog::PopupYesNo
         ( this->GetApplication(), this,
           "Save Volume",
           "Are you sure you want to save changes?" ) ) {

      // Try to save the view.
      try {
        // SAVE COMMAND HERE
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
vtkKWScubaWindow::LoadVolume ( const char* ifnVolume ) {
  
  try {

    // Make a MRI layer collection and set it up.
    vtkSmartPointer<vtkKWScubaLayerCollectionMRI> col = 
      vtkSmartPointer<vtkKWScubaLayerCollectionMRI>::New();
    col->SetApplication( this->GetApplication() );

    // Set the file name in the collectoin.
    col->SetVolumeFileName( ifnVolume );
    
    // Select this layer.
    this->SetCurrentLayerCollection( *col );
    
    // Status message.
    this->SetStatusText( "Volume loaded." );

    this->AddLayerCollectionToViews( col );
    
    // Add this file to the Recent menu.
    this->AddRecentFile( ifnVolume, this, "LoadVolume" );
    
    // Update our layer menu.
    this->UpdateLayerMenu();
    
  } catch ( exception& e ) {
    this->GetApplication()->ErrorMessage( e.what() );
  }
  
  // Update our menu and buttons.
  this->UpdateCommandStatus();
}

void
vtkKWScubaWindow::LoadSurface ( const char* ifnSurface ) {

  try {
    
    // Make a MRIS layer collectio and set it up.
    vtkSmartPointer<vtkKWScubaLayerCollectionMRIS> col =
      vtkSmartPointer<vtkKWScubaLayerCollectionMRIS>::New();
    col->SetApplication( this->GetApplication() );
    
    // Set the file name in the collection.
    col->SetSurfaceFileName( ifnSurface );
    
    // Select this layer.
    this->SetCurrentLayerCollection( *col );
    
    // Status message.
    this->SetStatusText( "Surface loaded." );

    this->AddLayerCollectionToViews( col );

    // Add this file to the Recent menu.
    this->AddRecentFile( ifnSurface, this, "LoadSurface" );
    
    // Update our layer menu.
    this->UpdateLayerMenu();
    
  } catch ( exception& e ) {
    this->GetApplication()->ErrorMessage( e.what() );
  }

  // Update our menu and buttons.
  this->UpdateCommandStatus();
}

void
vtkKWScubaWindow::LoadDTI ( const char* ifnDTI ) {
  
  try {

    // Make a DTI layer collection and set it up.
    vtkSmartPointer<vtkKWScubaLayerCollectionDTI> col =
      vtkSmartPointer<vtkKWScubaLayerCollectionDTI>::New();
    col->SetApplication( this->GetApplication() );

    // Set the file name in the collectoin.
    col->SetFAVolumeFileName( ifnDTI );
    
    // Select this layer.
    this->SetCurrentLayerCollection( *col );
    
    // Status message.
    this->SetStatusText( "DTI loaded." );

    this->AddLayerCollectionToViews( col );
    
    // Add this file to the Recent menu.
    this->AddRecentFile( ifnDTI, this, "LoadDTI" );
    
    // Update our layer menu.
    this->UpdateLayerMenu();
    
  } catch ( exception& e ) {
    this->GetApplication()->ErrorMessage( e.what() );
  }
  
  // Update our menu and buttons.
  this->UpdateCommandStatus();
}

void
vtkKWScubaWindow::LoadPath ( const char* ifnPath ) {
  try {

    // Make a path layer collection and set it up.
    vtkSmartPointer<vtkKWScubaLayerCollectionPath> col =
      vtkSmartPointer<vtkKWScubaLayerCollectionPath>::New();
    col->SetApplication( this->GetApplication() );

    // Set the file name in the collection.
    col->SetVolumeFileName( ifnPath );
    
    // Select this layer.
    this->SetCurrentLayerCollection( *col );
    
    // Status message.
    this->SetStatusText( "Path loaded." );

    this->AddLayerCollectionToViews( col );
    
    // Add this file to the Recent menu.
    this->AddRecentFile( ifnPath, this, "LoadPath" );
    
    // Update our layer menu.
    this->UpdateLayerMenu();
    
  } catch ( exception& e ) {
    this->GetApplication()->ErrorMessage( e.what() );
  }
  
  // Update our menu and buttons.
  this->UpdateCommandStatus();
}

void
vtkKWScubaWindow::LoadODF ( const char* ifnODF ) {
  try {

    // Make an ODF layer collection and set it up.
    vtkSmartPointer<vtkKWScubaLayerCollectionODF> col =
      vtkSmartPointer<vtkKWScubaLayerCollectionODF>::New();
    col->SetApplication( this->GetApplication() );

    // Set the file name in the collection.
    col->SetODFVolumeFileName( ifnODF );
    
    // Select this layer.
    this->SetCurrentLayerCollection( *col );
    
    // Status message.
    this->SetStatusText( "ODF loaded." );

    this->AddLayerCollectionToViews( col );
    
    // Add this file to the Recent menu.
    this->AddRecentFile( ifnODF, this, "LoadODF" );
    
    // Update our layer menu.
    this->UpdateLayerMenu();
    
  } catch ( exception& e ) {
    this->GetApplication()->ErrorMessage( e.what() );
  }
  
  // Update our menu and buttons.
  this->UpdateCommandStatus();
}


void
vtkKWScubaWindow::ZoomBy ( float iFactor ) {

  // Go through all our views.
  map<int,vtkSmartPointer<vtkKWScubaView> >::iterator tView;
  for( tView = maView.begin(); tView != maView.end(); ++tView ) {
    vtkKWScubaView* view = tView->second;
    view->Set2DZoomLevel( view->Get2DZoomLevel() * iFactor );
  }

  // Update our menu and buttons.
  this->UpdateCommandStatus();
}

void
vtkKWScubaWindow::ZoomIn () {
  this->ZoomBy( 2.0 );
}

void
vtkKWScubaWindow::ZoomOut () {
  this->ZoomBy( 0.5 );
}

void
vtkKWScubaWindow::SetViewLayout ( ViewLayout iLayout ) {

  // If we're not already in this view layout...
  if( mCurrentViewLayout != iLayout ) {

    mCurrentViewLayout = iLayout;

    // Unpack all our vurrent views.
    this->GetViewPanelFrame()->UnpackChildren();

    // Switch on the kind of layout.
    switch( mCurrentViewLayout ) {
    case Single: {
      
      // Just make one view the size of the whole frame.
      vtkKWScubaView* view = this->GetOrMakeNthView( 0 );
      this->Script( "pack %s -expand yes -fill both -anchor c",
		    view->GetWidgetName() );

      // Set its label.
      view->SetLabel( "Main" );

      // Select this view.
      this->SetCurrentView( *view );
      
    } 
      break;
    case TwoByTwo: {

      // Get four views.
      vtkKWScubaView* topLeft = this->GetOrMakeNthView( 0 );
      vtkKWScubaView* topRight = this->GetOrMakeNthView( 1 );
      vtkKWScubaView* bottomLeft = this->GetOrMakeNthView( 2 );
      vtkKWScubaView* bottomRight = this->GetOrMakeNthView( 3 );

      // Pack them in a 2x2 grid.
      this->Script( "grid %s -column 0 -row 0 -sticky news",
		    topLeft->GetWidgetName() );
      this->Script( "grid %s -column 1 -row 0 -sticky news",
		    topRight->GetWidgetName() );
      this->Script( "grid %s -column 0 -row 1 -sticky news",
		    bottomLeft->GetWidgetName() );
      this->Script( "grid %s -column 1 -row 1 -sticky news",
		    bottomRight->GetWidgetName() );

      // Give all columns and rows equal weight.
      vtkKWWidget* parent = topLeft->GetParent();
      this->Script( "grid rowconfigure %s 0 -weight 1", 
 		    parent->GetWidgetName() );
      this->Script( "grid rowconfigure %s 1 -weight 1", 
 		    parent->GetWidgetName() );
      this->Script( "grid columnconfigure %s 0 -weight 1", 
 		    parent->GetWidgetName() );
      this->Script( "grid columnconfigure %s 1 -weight 1", 
 		    parent->GetWidgetName() );
      
      // Give them descriptive labels.
      topLeft->SetLabel( "Top Left" );
      topRight->SetLabel( "Top Right" );
      bottomLeft->SetLabel( "Bottom Left" );
      bottomRight->SetLabel( "Bottom Right" );

      // Select the top left view.
      this->SetCurrentView( *topLeft );
    }
      break;
    default:
      throw runtime_error( "Invalid view layout" );
    }

    // Set the current view settings (like RAS positions and
    // orientation) for this layout.
    this->InitializeViewSettingsForLayout();
  
    // Update the view menu with our new views and names.
    this->UpdateViewMenu();
  }
}

vtkKWScubaWindow::ViewLayout
vtkKWScubaWindow::GetViewLayout () const {
  return mCurrentViewLayout;
}

void
vtkKWScubaWindow::SetViewLayoutToSingle () {

  this->SetViewLayout( Single );
}

void
vtkKWScubaWindow::SetViewLayoutToTwoByTwo () {

  this->SetViewLayout( TwoByTwo );
}


vtkKWScubaTool*
vtkKWScubaWindow::GetCurrentTool () {
  return mCurrentTool;
}

void
vtkKWScubaWindow::SetCurrentTool ( vtkKWScubaTool& iTool ) {

  if( mCurrentTool.GetPointer() == &iTool )
    return;
  
  // If we have a current tool, tell it to unpopulate the UI frame.
  if ( mCurrentTool.GetPointer() ) {
    mCurrentTool->DepopulateControlPage();
  }

  // Unpack the tool UI frame and tell the new tool to populate it.
  mToolUIFrame->UnpackChildren();
  mToolUIFrame->SetHeight( 1 );
  iTool.PopulateControlPage( mToolUIFrame );

  // Set the text in the menu.
  mMenuTool->SetValue( iTool.GetLabel() );

  // Save the current tool.
  mCurrentTool = &iTool;
}

void
vtkKWScubaWindow::SetCurrentToolFromMenu () {

  assert( mMenuTool.GetPointer() );

  // Try to get an entry. If it's in our map...
  int nEntry = mMenuTool->GetMenu()->GetIndexOfItem( mMenuTool->GetValue() );
  if ( mToolMenuIndexToPointerMap.end() !=
       mToolMenuIndexToPointerMap.find( nEntry ) ) {
    
      // Set the current tool from the pointer our map.
    this->SetCurrentTool( *mToolMenuIndexToPointerMap[nEntry] );
  } else {
    throw runtime_error( "Tool table doens't have that entry" );
  }
}

void
vtkKWScubaWindow::SetCurrentView ( vtkKWScubaView& iView ) {

  if( mCurrentView.GetPointer() == &iView )
    return;

  // If we have a current view, tell it to unpopulate the UI frame.
  if ( mCurrentView.GetPointer() ) {
    mCurrentView->DepopulateControlPage();
    mCurrentView->DepopulateToolbar();
  }

  // Unpack the view UI frame and tell the new view to populate it.
  mViewUIFrame->UnpackChildren();
  mViewUIFrame->SetHeight( 1 );
  iView.PopulateControlPage( mViewUIFrame );

  // Show the view's toolbar.
  mToolbarView->RemoveAllWidgets();
  iView.PopulateToolbar( mToolbarView );

  // Set the text in the menu.
  mMenuView->SetValue( iView.GetLabel() );

  // Save the current view.
  mCurrentView = &iView;
}

void
vtkKWScubaWindow::SetCurrentViewFromMenu () {
 
  assert( mMenuView.GetPointer() );
 
  // Try to get an entry. If it's in our map...
  int nEntry = mMenuView->GetMenu()->GetIndexOfItem( mMenuView->GetValue() );
  if ( mViewMenuIndexToPointerMap.end() !=
       mViewMenuIndexToPointerMap.find( nEntry ) ) {
    
    // Set the current view from the pointer our map.
    this->SetCurrentView( *mViewMenuIndexToPointerMap[nEntry] );
  } else {
    throw runtime_error( "View table doens't have that entry" );
  }
}


void
vtkKWScubaWindow::SetCurrentLayerCollection ( vtkKWScubaLayerCollection& iCol ) {

  // If we have a current layer collection, tell it to unpopulate the UI frame.
  if ( mCurrentLayerCollection.GetPointer() )
    mCurrentLayerCollection->DepopulateControlPage();

  // Unpack the layer UI frame and tell the new layer to populate it.
  mLayerUIFrame->UnpackChildren();
  mLayerUIFrame->SetHeight( 1 );
  iCol.PopulateControlPage( mLayerUIFrame );

  // Set the text in the menu box.
  mMenuLayer->SetValue( iCol.GetLabel() );

  // Save the current layer.
  mCurrentLayerCollection = &iCol;
}

void
vtkKWScubaWindow::SetCurrentLayerCollectionFromMenu () {

  assert( mMenuLayer.GetPointer() );

  // Try to get an entry. If it's in our map...
  int nEntry = mMenuLayer->GetMenu()->GetIndexOfItem( mMenuLayer->GetValue() );
  if ( mLayerMenuIndexToPointerMap.end() !=
       mLayerMenuIndexToPointerMap.find( nEntry ) ) {
    
    // Set the current layer from the pointer in our map.
    this->SetCurrentLayerCollection( *mLayerMenuIndexToPointerMap[nEntry] );
  } else {
    throw runtime_error( "Layer table doens't have that entry" );
  }
}

void
vtkKWScubaWindow::EventDone () {

  // This is called every time an event is done. If we have a view, if
  // it's info has changed, update our info area. This makes sure we
  // always update when the mosue moves or any other event happens.
  if ( mCurrentView.GetPointer() ) {
    if ( mCurrentView->IsInfoChanged() ) {
      this->UpdateInfoArea();
    }
  }
}

void
vtkKWScubaWindow::DoListenToMessage ( string const isMessage, 
				      void* const iData ) {

  // If the layer label changed, update our layer menu.
  if ( isMessage == "LayerLabelChanged" ) {
    this->UpdateLayerMenu();
  }
}

void
vtkKWScubaWindow::UpdateToolMenu () {

  // Clear out the menu and the map.
  mMenuTool->GetMenu()->DeleteAllItems();
  mToolMenuIndexToPointerMap.clear();

  // Get a list of tools.
  list<vtkKWScubaTool*> lTools;
  vtkKWScubaTool::GetPointerList( lTools );
  list<vtkKWScubaTool*>::iterator tTool;

  // For each one...
  int nEntry = 0;
  for ( tTool = lTools.begin(); tTool != lTools.end(); ++tTool ) {
    vtkKWScubaTool* tool = *tTool;

    // Add this entry to the menu, and associate this entry with
    // the pointer to the layer.
    mMenuTool->GetMenu()->AddRadioButton( tool->GetLabel(), 
					  this, "SetCurrentToolFromMenu" );
    mToolMenuIndexToPointerMap[nEntry++] = tool;
  }
}

int
vtkKWScubaWindow::GetAutoSizeInfoArea () {

  return mbAutoSizeInfoArea;
}

void
vtkKWScubaWindow::SetAutoSizeInfoArea ( int ibSize ) {

  if( ibSize != mbAutoSizeInfoArea ) {
    mbAutoSizeInfoArea = ibSize;

    // Write our pref.
    this->GetApplication()->
      SetRegistryValue( 2, sRegistryKey, sAutoSizeInfoAreaKey, 
			"%d", mbAutoSizeInfoArea );

    this->UpdateInfoArea();
  }
}


void
vtkKWScubaWindow::UpdateViewMenu () {

  assert( mMenuView.GetPointer() );

  // Clear out the menu and the map.
  mMenuView->GetMenu()->DeleteAllItems();
  mViewMenuIndexToPointerMap.clear();

  // Get a list of views.
  list<vtkKWScubaView*> lViews;
  vtkKWScubaView::GetPointerList( lViews );
  list<vtkKWScubaView*>::iterator tView;

  // For each one...
  int nEntry = 0;
  for ( tView = lViews.begin(); tView != lViews.end(); ++tView ) {
    vtkKWScubaView* view = *tView;

    // Add this entry to the menu, and associate this entry with
    // the pointer to the view.
    mMenuView->GetMenu()->AddRadioButton( view->GetLabel(),
					  this, "SetCurrentViewFromMenu" );
    mViewMenuIndexToPointerMap[nEntry++] = view;
  }

  // Set the value of the currently selected view.
  if( NULL != mCurrentView )
    mMenuView->SetValue( mCurrentView->GetLabel() );
}

void
vtkKWScubaWindow::UpdateLayerMenu () {

  assert( mMenuLayer.GetPointer() );

  // Clear out the menu and the map.
  mMenuLayer->GetMenu()->DeleteAllItems();
  mLayerMenuIndexToPointerMap.clear();

  // Get a list of collections.
  list<vtkKWScubaLayerCollection*> lCols;
  vtkKWScubaLayerCollection::GetPointerList( lCols );
  list<vtkKWScubaLayerCollection*>::iterator tCol;

  // For each one...
  int nEntry = 0;
  for ( tCol = lCols.begin(); tCol != lCols.end(); ++tCol ) {
    vtkKWScubaLayerCollection* col = *tCol;

    // Add this entry to the menu, and associate this entry with
    // the pointer to the collection.
    mMenuLayer->GetMenu()->AddRadioButton( col->GetLabel(),
				this, "SetCurrentLayerCollectionFromMenu" );
    mLayerMenuIndexToPointerMap[nEntry++] = col;
  }
}


void
vtkKWScubaWindow::UpdateCommandStatus () {

  assert( mMenuLoadVolume );
  assert( mBtnLoadVolume.GetPointer() );
  assert( mMenuZoomOut );
  assert( mMenuZoomIn );
  assert( mMenuSaveVolume );
  assert( mBtnSaveVolume.GetPointer() );

  // Enable/disable our menus and toolbar buttons accordingly.
  mMenuLoadVolume->SetStateToNormal();
  mBtnLoadVolume->SetStateToNormal();
  mMenuZoomOut->SetStateToNormal();
  mMenuZoomIn->SetStateToNormal();

  if ( mCurrentView && mCurrentView->IsDataDirty() ) {
    mMenuSaveVolume->SetStateToNormal();
    mBtnSaveVolume->SetStateToNormal();
  } else {
    mMenuSaveVolume->SetStateToDisabled();
    mBtnSaveVolume->SetStateToDisabled();
  }

}

void
vtkKWScubaWindow::UpdateInfoArea () {

  assert( mMouseOverInfoTable.GetPointer() );
  assert( mCursorInfoTable.GetPointer() );

  if ( mCurrentView ) {

    // Cursor items. Get the list of cursor items from the view,
    // then for each one, make an entry in our table.
    mCursorInfoTable->InsertCellText( 0, 0, "Cursor" );
    list<ScubaInfoItem> lCursorItems;
    mCurrentView->GetCursorInfoItems( lCursorItems );
    int nItem = 1;
    list<ScubaInfoItem>::iterator tItem;
    for ( tItem = lCursorItems.begin();
	  tItem != lCursorItems.end(); ++tItem ) {
      ScubaInfoItem& item = *tItem;
      mCursorInfoTable->InsertCellText( nItem, 0, item.GetLabel() );
      mCursorInfoTable->InsertCellText( nItem, 1, item.GetValue() );
      nItem++;
    }
    
    // Mouseover items. Get the list of mouseover items from the
    // view, then for each one, make an entry in our table.
    list<ScubaInfoItem> lMouseOverItems;
    mCurrentView->GetMouseOverInfoItems( lMouseOverItems );
    mMouseOverInfoTable->InsertCellText( 0, 0, "Mouse Over" );
    nItem = 1;
    for ( tItem = lMouseOverItems.begin();
	  tItem != lMouseOverItems.end(); ++tItem ) {
      ScubaInfoItem& item = *tItem;
      mMouseOverInfoTable->InsertCellText( nItem, 0, item.GetLabel() );
      mMouseOverInfoTable->InsertCellText( nItem, 1, item.GetValue() );
      nItem++;
    }
    
    // Keep track of how many items we had last time. If we have fewer...
    static int cLastItems = 0;
    if ( nItem < cLastItems ) {
      
      // Clear the rows we no longer need. Note that we count from
      // nItem to nLastDeepestRow, but keep clearing nItem, since
      // our number of rows change when we delete one. So we delete
      // the next row as many times as we need to.
      for ( int nRow = nItem; nRow <= cLastItems; nRow++ ) {
	mMouseOverInfoTable->DeleteRow( nItem );
	mCursorInfoTable->DeleteRow( nItem );
      }
    }
    // Remember how many items we had.
    cLastItems = nItem;
    
    // Set our panel height to something decent. There doesn't seem
    // to be any way to get the actual row height from the table, so
    // this is hopefully a good estimate.
    if( mbAutoSizeInfoArea ) 
      this->GetSecondarySplitFrame()->
	SetFrame1Size( mCursorInfoTable->GetHeight() * 10 );
  }
  
  // Tell the view that we have updated our info.
  mCurrentView->InfoUpdated();
}

vtkKWScubaView* 
vtkKWScubaWindow::GetOrMakeNthView ( int inView ) {
  
  map<int,vtkSmartPointer<vtkKWScubaView> >::iterator tView;
  tView = maView.find( inView );
  if( tView != maView.end() )
    if( tView->second )
      return tView->second;

  vtkSmartPointer<vtkKWScubaView> view = 
    vtkSmartPointer<vtkKWScubaView>::New();
  view->SetParent( this->GetViewPanelFrame() );
  view->Create();

  maView[inView] = view;

  // We listen to the view. This means the view can call SendBroadcast
  // and we'll get it in DoListenToMessage.
  view->AddListener( this );

  // Create our custom interactor style (listener) for this
  // view. Associate it with the view and our tool. Tell the view's
  // window to use this style.
  vtkSmartPointer<vtkScubaInteractorStyle> style = 
    vtkSmartPointer<vtkScubaInteractorStyle>::New();
  style->SetWindow( this );
  view->GetRenderWindow()->GetInteractor()->SetInteractorStyle( style );

  // We'll get the layers from the current view and set them in this
  // view.
  if( NULL != mCurrentView.GetPointer() ) {
    int cSlots = mCurrentView->GetHighestFilledLayerSlot();
    for( int nSlot = 0; nSlot <= cSlots; nSlot++ ) {
      // Try to get a collection at this slot.
      vtkKWScubaLayerCollection* col = 
	mCurrentView->GetCollectionAtSlot( nSlot );
      // Set the collection if we got one.
      if( col )
        view->SetLayerCollectionAtSlot( nSlot, col );
    }
    
    // Notify views that the layers have changed.
    view->LayerListChanged();
    
    // Reset our camera to fit all the data now.
    view->ResetAllCameras();
  }

  return view;
}

void
vtkKWScubaWindow::InitializeViewSettingsForLayout () {

  switch( mCurrentViewLayout ) {
  case Single: {

    // Single view, put it in 2D mode and on the 0 RAS point.
    vtkKWScubaView* view = this->GetOrMakeNthView( 0 );

    view->SetDisplayMode( vtkKWScubaView::TwoDee );
    view->Set2DRASZ( 0 );
    view->Set2DInPlane( 0 );

    this->SetCurrentView( *view );

  }
    break;
  case TwoByTwo: {

    // Four views.
    vtkKWScubaView* topLeft = this->GetOrMakeNthView( 0 );
    vtkKWScubaView* topRight = this->GetOrMakeNthView( 1 );
    vtkKWScubaView* bottomLeft = this->GetOrMakeNthView( 2 );
    vtkKWScubaView* bottomRight = this->GetOrMakeNthView( 3 );

    // Put the first three in 2d mode with X, Y, and Z orientations,
    // and the last one in 3d mode.
    topLeft->SetDisplayMode( vtkKWScubaView::TwoDee );
    topLeft->Set2DRASZ( 0 );
    topLeft->Set2DInPlane( 0 );
    topRight->SetDisplayMode( vtkKWScubaView::TwoDee );
    topRight->Set2DRASZ( 0 );
    topRight->Set2DInPlane( 1 );
    bottomLeft->SetDisplayMode( vtkKWScubaView::TwoDee );
    bottomLeft->Set2DRASZ( 0 );
    bottomLeft->Set2DInPlane( 2 );
    bottomRight->SetDisplayMode( vtkKWScubaView::ThreeDee );
    bottomRight->Set3DRASX( 0 );
    bottomRight->Set3DRASY( 0 );
    bottomRight->Set3DRASZ( 0 );

    this->SetCurrentView( *topLeft );
  }
    break;
  default:
    throw runtime_error( "Invalid view layout" );
  }
}

void
vtkKWScubaWindow::AddLayerCollectionToViews ( vtkKWScubaLayerCollection* col ) {
  
  // Go through all our views.
  map<int,vtkSmartPointer<vtkKWScubaView> >::iterator tView;
  for( tView = maView.begin(); tView != maView.end(); ++tView ) {
    vtkKWScubaView* view = tView->second;

    // Find a good slot number. If this returns -1, we'll add a new slot.
    int nSlot = view->GetFirstUnusedLayerSlot();
    if( -1 == nSlot ) nSlot = 0;
  
    // Set the layer in the view at a slot.
    view->SetLayerCollectionAtSlot( nSlot, col );
    
    // Notify views that the layers have changed.
    view->LayerListChanged();

    // Reset our camera to fit all the data now.
    view->ResetAllCameras();
  }
  
}

void
vtkKWScubaWindow::AddSpacerToToolbar ( vtkKWToolbar* iToolbar, int iWidth ) {
  
  // Create a frame of the specified width inside the toolbar.
  vtkSmartPointer<vtkKWFrame> spacer = 
    vtkSmartPointer<vtkKWFrame>::New();
  spacer->SetParent( iToolbar );
  spacer->Create();
  spacer->SetWidth( iWidth );
  iToolbar->AddWidget( spacer );

}

vtkKWScubaWindow::MenuItem::MenuItem ()  :
    mMenu( NULL ),
    mnItem( -1 ) {}


void
vtkKWScubaWindow::MenuItem::MakeCommand ( vtkKWMenu* iMenu,
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
vtkKWScubaWindow::MenuItem::MakeCheckButton ( vtkKWMenu* iMenu,
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
vtkKWScubaWindow::MenuItem::SetStateToDisabled () {

  if( NULL == mMenu.GetPointer() )
    throw runtime_error( "MenuItem::SetStateToDisabled: mMenu was NULL" );
  if( -1 == mnItem )
    throw runtime_error( "MenuItem::SetStateToDisabled: mnItem was -1" );

  mMenu->SetItemStateToDisabled( mnItem );
}

void
vtkKWScubaWindow::MenuItem::SetStateToNormal () {

  if( NULL == mMenu.GetPointer() )
    throw runtime_error( "MenuItem::SetStateToNormal: mMenu was NULL" );
  if( -1 == mnItem )
    throw runtime_error( "MenuItem::SetStateToNormal: mnItem was -1" );

  mMenu->SetItemStateToNormal( mnItem );
}

void
vtkKWScubaWindow::MenuItem::SetSelectedState ( int ibOn ) {

  if( NULL == mMenu.GetPointer() )
    throw runtime_error( "MenuItem::SetSelectedState: mMenu was NULL" );
  if( -1 == mnItem )
    throw runtime_error( "MenuItem::SetSelectedState: mnItem was -1" );

  mMenu->SetItemSelectedState( mnItem, ibOn );
}

int
vtkKWScubaWindow::MenuItem::GetSelectedState () {

  if( NULL == mMenu.GetPointer() )
    throw runtime_error( "MenuItem::GetSelectedState: mMenu was NULL" );
  if( -1 == mnItem )
    throw runtime_error( "MenuItem::GetSelectedState: mnItem was -1" );

  return mMenu->GetItemSelectedState( mnItem );
}

int
vtkKWScubaWindow::MenuItem::GetIndex () const {

  return mnItem;
}

