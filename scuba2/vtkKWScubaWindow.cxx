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
 *    $Author: kteich $
 *    $Date: 2007/04/06 22:23:06 $
 *    $Revision: 1.1 $
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


#include <string>
#include "IconLoader.h"
#include "vtkKWApplication.h"
#include "vtkKWFrameWithLabel.h"
#include "vtkKWIcon.h"
#include "vtkKWLoadSaveDialog.h"
#include "vtkKWScubaLayerCollectionDTI.h"
#include "vtkKWScubaLayerCollectionMRI.h"
#include "vtkKWScubaLayerCollectionMRIS.h"
#include "vtkKWScubaLayerCollectionPath.h"
#include "vtkKWScubaToolEdit2DMRI.h"
#include "vtkKWScubaToolNavigate.h"
#include "vtkKWScubaWindow.h"
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
#include "vtkKWToolbar.h"
#include "vtkKWToolbarSet.h"
#include "vtkKWUserInterfacePanel.h"
#include "vtkScubaInteractorStyle.h"
#include "vtkObjectFactory.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"

using namespace std;

vtkStandardNewMacro( vtkKWScubaWindow );
vtkCxxRevisionMacro( vtkKWScubaWindow, "$Revision: 1.1 $" );

const string vtkKWScubaWindow::DEFAULT_VOLUME_FILE_EXTENSION = ".mgz";

vtkKWScubaWindow::vtkKWScubaWindow () :
    Listener( "vtkKWScubaWindow" ),
    mToolUIFrame( NULL ),
    mViewUIFrame( NULL ),
    mLayerUIFrame( NULL ),
    mMenuTool( NULL ),
    mMenuView( NULL ),
    mMenuLayer( NULL ),
    mBtnLoadVolume( NULL ),
    mBtnSaveVolume( NULL ),
    mCursorInfoTable( NULL ),
    mMouseOverInfoTable( NULL ),
    mCurrentTool( NULL ),
    mCurrentView( NULL ),
    mCurrentLayerCollection( NULL ),
    mCurrentViewLayout( NoViewLayout ) {

}

vtkKWScubaWindow::~vtkKWScubaWindow () {

  if ( mCursorInfoTable ) mCursorInfoTable->Delete();
  if ( mMouseOverInfoTable ) mMouseOverInfoTable->Delete();
  if ( mToolUIFrame ) mToolUIFrame->Delete();
  if ( mViewUIFrame ) mLayerUIFrame->Delete();
  if ( mLayerUIFrame ) mLayerUIFrame->Delete();
  if ( mBtnLoadVolume ) mBtnLoadVolume->Delete();
  if ( mBtnSaveVolume ) mBtnSaveVolume->Delete();

  map<int,vtkKWScubaView*>::iterator tView;
  for( tView = maView.begin(); tView != maView.end(); ++tView ) {
    vtkKWScubaView* view = tView->second;
    if( view )
      view->Delete();
  }
}

void
vtkKWScubaWindow::Create () {

  // Our default view layout.
  SetPanelLayoutToSecondaryBelowMainAndView();

  // Create the window. The layout must be set before this.
  this->Superclass::Create();

  // Turn on our panels.
  SecondaryPanelVisibilityOn();
  MainPanelVisibilityOn();

  //
  // Create the secondary panel. This is our 'info area' and contains
  // two tables that we populate with data from the view and
  // layers. Note that we don't use the secondary UserInterface
  // Manager, we just pack these things right into the frame. This is
  // because we don't need any kind of paged interface here.
  mCursorInfoTable = vtkKWMultiColumnList::New();
  mCursorInfoTable->SetParent( this->GetSecondaryPanelFrame() );
  mCursorInfoTable->Create();
  mCursorInfoTable->ColumnLabelsVisibilityOff();
  mCursorInfoTable->SetWidth( 1 );
  mCursorInfoTable->AddColumn( "" );
  mCursorInfoTable->ColumnStretchableOff( 0 );
  mCursorInfoTable->AddColumn( "" );
  mCursorInfoTable->ColumnStretchableOn( 1 );

  mMouseOverInfoTable = vtkKWMultiColumnList::New();
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
  vtkKWFrameWithLabel* toolFrame = vtkKWFrameWithLabel::New();
  toolFrame->SetParent( panelFrame );
  toolFrame->Create();
  toolFrame->SetLabelText( "Tool" );
  this->Script( "pack %s -side top -fill x -anchor nw",
		toolFrame->GetWidgetName() );

  // Menu for tools.
  vtkKWMenuButtonWithSpinButtonsWithLabel* labeledMenu =
    vtkKWMenuButtonWithSpinButtonsWithLabel::New();
  labeledMenu->SetParent( toolFrame->GetFrame() );
  labeledMenu->Create();
  labeledMenu->SetLabelText( "Tool: " );

  mMenuTool = labeledMenu->GetWidget()->GetWidget();

  // A subframe for the tool.
  mToolUIFrame = vtkKWFrame::New();
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
  vtkKWFrameWithLabel* viewFrame = vtkKWFrameWithLabel::New();
  viewFrame->SetParent( panelFrame );
  viewFrame->Create();
  viewFrame->SetLabelText( "View" );
  this->Script( "pack %s -side top -fill x -anchor nw",
		viewFrame->GetWidgetName() );

  labeledMenu = vtkKWMenuButtonWithSpinButtonsWithLabel::New();
  labeledMenu->SetParent( viewFrame->GetFrame() );
  labeledMenu->Create();
  labeledMenu->SetLabelText( "View: " );

  mMenuView = labeledMenu->GetWidget()->GetWidget();

  mViewUIFrame = vtkKWFrame::New();
  mViewUIFrame->SetParent( viewFrame->GetFrame() );
  mViewUIFrame->Create();

  this->Script( "pack %s %s -side top -fill x -anchor nw -pady 5",
                labeledMenu->GetWidgetName(),
                mViewUIFrame->GetWidgetName() );


  //
  // Create the Layer frame. Make our layer-selecting menu for it,
  // and then another frame underneath that. The layers will populate
  // that second frame. Leave it blank for now.
  vtkKWFrameWithLabel* layerFrame = vtkKWFrameWithLabel::New();
  layerFrame->SetParent( panelFrame );
  layerFrame->Create();
  layerFrame->SetLabelText( "Layer" );
  this->Script( "pack %s -side top -fill x -anchor nw",
		layerFrame->GetWidgetName() );

  labeledMenu = vtkKWMenuButtonWithSpinButtonsWithLabel::New();
  labeledMenu->SetParent( layerFrame->GetFrame() );
  labeledMenu->Create();
  labeledMenu->SetLabelText( "Layer: " );

  mMenuLayer = labeledMenu->GetWidget()->GetWidget();

  mLayerUIFrame = vtkKWFrame::New();
  mLayerUIFrame->SetParent( layerFrame->GetFrame() );
  mLayerUIFrame->Create();

  this->Script( "pack %s %s -side top -fill x -anchor nw -pady 5",
                labeledMenu->GetWidgetName(),
                mLayerUIFrame->GetWidgetName() );

  // ====================================================================

  //
  // Make our tool toolbar.
  mToolbarTools = vtkKWToolbar::New();
  mToolbarTools->SetName( "Tools" );
  mToolbarTools->SetParent( this->GetMainToolbarSet()->GetToolbarsFrame() );
  mToolbarTools->Create();
  this->GetMainToolbarSet()->AddToolbar( mToolbarTools );
  
  //
  // Make our window command toolbar.
  mToolbarWindow = vtkKWToolbar::New();
  mToolbarWindow->SetName( "Main" );
  mToolbarWindow->SetParent( this->GetMainToolbarSet()->GetToolbarsFrame() );
  mToolbarWindow->Create();
  this->GetMainToolbarSet()->AddToolbar( mToolbarWindow );
  
  //
  // Make the view toolbar, we'll fill it in later.
  mToolbarView = vtkKWToolbar::New();
  mToolbarView->SetName( "View" );
  mToolbarView->SetParent( this->GetMainToolbarSet()->GetToolbarsFrame() );
  mToolbarView->Create();
  this->GetMainToolbarSet()->AddToolbar( mToolbarView );

  // Make our buttons. For each one, we create a button, set the
  // parent, give it some basic characterists, an icon, and set its
  // command. We finish by adding it to the toolbar, which will pack
  // it for us.

  // Load Volume
  mBtnLoadVolume = vtkKWPushButton::New();
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
  mBtnSaveVolume = vtkKWPushButton::New();
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
  vtkKWRadioButtonSet* mRadBtnSetViewLayout = vtkKWRadioButtonSet::New();
  mRadBtnSetViewLayout->SetParent( mToolbarWindow->GetFrame() );
  mRadBtnSetViewLayout->Create();
  mRadBtnSetViewLayout->PackHorizontallyOn();
  mToolbarWindow->AddWidget( mRadBtnSetViewLayout );

  vtkKWRadioButton* mRadBtnViewLayoutSingle = 
    mRadBtnSetViewLayout->AddWidget( (int)Single );
  mRadBtnViewLayoutSingle->SetCommand( this, "SetViewLayoutToSingle" );
  mRadBtnViewLayoutSingle->IndicatorVisibilityOff();
  try { IconLoader::SetCheckButtonIcon("ViewSingle",mRadBtnViewLayoutSingle); }
  catch( exception& e ) { cerr << "Error loading icon: " << e.what() << endl; }
  if( Single == mCurrentViewLayout )
    mRadBtnViewLayoutSingle->SelectedStateOn();
  
  vtkKWRadioButton* mRadBtnViewLayoutTwoByTwo = 
    mRadBtnSetViewLayout->AddWidget( (int)TwoByTwo );
  mRadBtnViewLayoutTwoByTwo->SetCommand( this, "SetViewLayoutToTwoByTwo" );
  mRadBtnViewLayoutTwoByTwo->IndicatorVisibilityOff();
  try { IconLoader::SetCheckButtonIcon("ViewTwoByTwo",mRadBtnViewLayoutTwoByTwo); }
  catch( exception& e ) { cerr << "Error loading icon: " << e.what() << endl; }
  if( TwoByTwo == mCurrentViewLayout )
    mRadBtnViewLayoutTwoByTwo->SelectedStateOn();

  //
  // Build the menus. File menu.

  int nFilePos = GetFileMenuInsertPosition();

  // Load Volume.
  GetFileMenu()->
  InsertCommand( nFilePos, "L&oad Volume...", this, "LoadVolumeFromDlog" );
  GetFileMenu()->
  SetItemImageToPredefinedIcon( nFilePos, vtkKWIcon::IconFileOpen );
  GetFileMenu()->SetItemCompoundModeToLeft( nFilePos );
  GetFileMenu()->SetItemAccelerator( nFilePos, "Ctrl+O" );
  mMenuLoadVolume.menu = GetFileMenu();
  mMenuLoadVolume.nItem = nFilePos;
  try { IconLoader::SetMenuItemIcon( "LoadVolume", GetFileMenu(), nFilePos ); }
  catch( exception& e ) { cerr << "Error loading icon: " << e.what() << endl; }
  nFilePos++;

  // Load Surface.
  GetFileMenu()->
  InsertCommand( nFilePos, "Load Surface...", this, "LoadSurfaceFromDlog" );
  GetFileMenu()->
  SetItemImageToPredefinedIcon( nFilePos, vtkKWIcon::IconFileOpen );
  GetFileMenu()->SetItemCompoundModeToLeft( nFilePos );
  mMenuLoadSurface.menu = GetFileMenu();
  mMenuLoadSurface.nItem = nFilePos;
  try { IconLoader::SetMenuItemIcon("LoadSurface", GetFileMenu(), nFilePos ); }
  catch( exception& e ) { cerr << "Error loading icon: " << e.what() << endl; }
  nFilePos++;

  // Load DTI.
  GetFileMenu()->
  InsertCommand( nFilePos, "Load DTI...", this, "LoadDTIFromDlog" );
  GetFileMenu()->
  SetItemImageToPredefinedIcon( nFilePos, vtkKWIcon::IconFileOpen );
  GetFileMenu()->SetItemCompoundModeToLeft( nFilePos );
  mMenuLoadDTI.menu = GetFileMenu();
  mMenuLoadDTI.nItem = nFilePos;
  try { IconLoader::SetMenuItemIcon( "LoadDTI", GetFileMenu(), nFilePos ); }
  catch( exception& e ) { cerr << "Error loading icon: " << e.what() << endl; }
  nFilePos++;
  
  // Load Path
  GetFileMenu()->
  InsertCommand( nFilePos, "Load Path...", this, "LoadPathFromDlog" );
  GetFileMenu()->
  SetItemImageToPredefinedIcon( nFilePos, vtkKWIcon::IconFileOpen );
  GetFileMenu()->SetItemCompoundModeToLeft( nFilePos );
  mMenuLoadPath.menu = GetFileMenu();
  mMenuLoadPath.nItem = nFilePos;
  try { IconLoader::SetMenuItemIcon( "LoadPath", GetFileMenu(), nFilePos ); }
  catch( exception& e ) { cerr << "Error loading icon: " << e.what() << endl; }
  nFilePos++;

  GetFileMenu()->InsertSeparator( nFilePos++ );

  // Save Volume.
  GetFileMenu()->
  InsertCommand( nFilePos, "&Save Volume", this, "SaveVolumeWithConfirm" );
  GetFileMenu()->
  SetItemImageToPredefinedIcon( nFilePos, vtkKWIcon::IconFloppy );
  GetFileMenu()->SetItemCompoundModeToLeft( nFilePos );
  GetFileMenu()->SetItemAccelerator( nFilePos, "Ctrl+S" );
  mMenuSaveVolume.menu = GetFileMenu();
  mMenuSaveVolume.nItem = nFilePos;
  try { IconLoader::SetMenuItemIcon( "SaveVolume", GetFileMenu(), nFilePos ); }
  catch( exception& e ) { cerr << "Error loading icon: " << e.what() << endl; }
  nFilePos++;

  GetFileMenu()->InsertSeparator( nFilePos++ );

  InsertRecentFilesMenu( nFilePos++, this );

  GetFileMenu()->InsertSeparator( nFilePos++ );

  // View menu.
  int nViewPos = GetViewMenuInsertPosition();

  // Zoom Out.
  GetViewMenu()->InsertCommand( nViewPos, "Zoom Out", this, "ZoomOut");
  GetViewMenu()->SetItemCompoundModeToLeft( nViewPos );
  //  GetViewMenu()->SetItemAccelerator( nViewPos, "Ctrl+Minus" );
  mMenuZoomOut.menu = GetViewMenu();
  mMenuZoomOut.nItem = nViewPos;
  try { IconLoader::SetMenuItemIcon( "ZoomOut", GetViewMenu(), nViewPos ); }
  catch( exception& e ) { cerr << "Error loading icon: " << e.what() << endl; }
  nViewPos++;

  // Zoom In.
  GetViewMenu()->InsertCommand( nViewPos, "Zoom In", this, "ZoomIn");
  GetViewMenu()->SetItemCompoundModeToLeft( nViewPos );
  //  GetViewMenu()->SetItemAccelerator( nViewPos, "Ctrl+Plus" );
  mMenuZoomIn.menu = GetViewMenu();
  mMenuZoomIn.nItem = nViewPos;
  try { IconLoader::SetMenuItemIcon( "ZoomIn", GetViewMenu(), nViewPos ); }
  catch( exception& e ) { cerr << "Error loading icon: " << e.what() << endl; }
  nViewPos++;

  //
  // Load initial tools. This will just create them and then we will
  // get them later with the IDTracker methods when we build our
  // menu. Don't Delete them tho, as there are no other references to
  // them.
  vtkKWScubaTool* navigateTool = vtkKWScubaToolNavigate::New();
  navigateTool->SetApplication( this->GetApplication() );

  vtkKWScubaTool* tool = vtkKWScubaToolEdit2DMRI::New();
  tool->SetApplication( this->GetApplication() );

  // Make a radio button set for our toolbar.
  mRadBtnSetTool = vtkKWRadioButtonSet::New();
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
vtkKWScubaWindow::LoadVolumeFromDlog () {

  // Create a Load dialog and set it up.
  vtkKWLoadSaveDialog* dialog = vtkKWLoadSaveDialog::New();
  dialog->SetApplication( this->GetApplication() );
  dialog->Create();
  dialog->SetTitle( "Load a volume" );
  dialog->SetFileTypes( "{MGH {.mgh .mgz}} "
                        "{Binary {.bshort .bfloat}} {All {*}}" );
  dialog->RetrieveLastPathFromRegistry( "LoadVolume" );
  dialog->SetDefaultExtension( DEFAULT_VOLUME_FILE_EXTENSION.c_str() );

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
  vtkKWLoadSaveDialog* dialog = vtkKWLoadSaveDialog::New();
  dialog->SetApplication( this->GetApplication() );
  dialog->Create();
  dialog->SetTitle( "Load a DTI" );
  dialog->SetFileTypes( "{MGH {.mgh .mgz}} "
                        "{Binary {.bshort .bfloat}} {All {*}}" );
  dialog->RetrieveLastPathFromRegistry( "LoadDTI" );
  dialog->SetDefaultExtension( DEFAULT_VOLUME_FILE_EXTENSION.c_str() );

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
  vtkKWLoadSaveDialog* dialog = vtkKWLoadSaveDialog::New();
  dialog->SetApplication( this->GetApplication() );
  dialog->Create();
  dialog->SetTitle( "Load a path" );
  dialog->SetFileTypes( "{MGH {.mgh .mgz}} "
                        "{Binary {.bshort .bfloat}} {All {*}}" );
  dialog->RetrieveLastPathFromRegistry( "LoadPath" );
  dialog->SetDefaultExtension( DEFAULT_VOLUME_FILE_EXTENSION.c_str() );

  // Show the dialog, and when it returns, Invoke() will be true if
  // they clicked OK and gave us a filename.
  if ( dialog->Invoke() ) {
    dialog->SaveLastPathToRegistry( "LoadPath" );
    string fnPath( dialog->GetFileName() );
    this->LoadPath( fnPath.c_str() );
  }
  
}


void
vtkKWScubaWindow::SaveVolumeWithConfirm () {

  // If the view has dirty data...
  if ( mCurrentView && mCurrentView->IsDataDirty() ) {

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
    vtkKWScubaLayerCollectionMRI* col = vtkKWScubaLayerCollectionMRI::New();
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
    vtkKWScubaLayerCollectionMRIS* col =vtkKWScubaLayerCollectionMRIS::New();
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
    vtkKWScubaLayerCollectionDTI* col = vtkKWScubaLayerCollectionDTI::New();
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
    vtkKWScubaLayerCollectionPath* col = vtkKWScubaLayerCollectionPath::New();
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
vtkKWScubaWindow::ZoomBy ( float iFactor ) {

  // Go through all our views.
  map<int,vtkKWScubaView*>::iterator tView;
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

  if( mCurrentViewLayout != iLayout ) {

    mCurrentViewLayout = iLayout;

    this->GetViewPanelFrame()->UnpackChildren();

    if( Single == mCurrentViewLayout ) {
      
      vtkKWScubaView* view = this->GetNthView( 0 );
      this->Script( "pack %s -expand yes -fill both -anchor c",
		    view->GetWidgetName() );
      view->SetLabel( "Main" );

      this->SetCurrentView( *view );

      //      view->Delete();

    } else if( TwoByTwo == mCurrentViewLayout ) {

      vtkKWScubaView* topLeft = this->GetNthView( 0 );
      vtkKWScubaView* topRight = this->GetNthView( 1 );
      vtkKWScubaView* bottomLeft = this->GetNthView( 2 );
      vtkKWScubaView* bottomRight = this->GetNthView( 3 );

      this->Script( "grid %s -column 0 -row 0 -sticky news",
		    topLeft->GetWidgetName() );
      this->Script( "grid %s -column 1 -row 0 -sticky news",
		    topRight->GetWidgetName() );
      this->Script( "grid %s -column 0 -row 1 -sticky news",
		    bottomLeft->GetWidgetName() );
      this->Script( "grid %s -column 1 -row 1 -sticky news",
		    bottomRight->GetWidgetName() );

       vtkKWWidget* parent = topLeft->GetParent();
       this->Script( "grid rowconfigure %s 0 -weight 1", 
 		    parent->GetWidgetName() );
       this->Script( "grid rowconfigure %s 1 -weight 1", 
 		    parent->GetWidgetName() );
       this->Script( "grid columnconfigure %s 0 -weight 1", 
 		    parent->GetWidgetName() );
       this->Script( "grid columnconfigure %s 1 -weight 1", 
 		    parent->GetWidgetName() );

      topLeft->SetLabel( "Top Left" );
      topRight->SetLabel( "Top Right" );
      bottomLeft->SetLabel( "Bottom Left" );
      bottomRight->SetLabel( "Bottom Right" );

      this->SetCurrentView( *topLeft );

//       topLeft->Delete();
//       topRight->Delete();
//       bottomLeft->Delete();
//       bottomRight->Delete();

    }

    this->InitializeViewSettingsForLayout();
  
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

  if( mCurrentTool == &iTool )
    return;
  
  // If we have a current view, tell it to unpopulate the UI frame.
  if ( mCurrentTool ) {
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

#include "vtkRenderer.h"

void
vtkKWScubaWindow::SetCurrentView ( vtkKWScubaView& iView ) {

  if( mCurrentView == &iView )
    return;

  // If we have a current view, tell it to unpopulate the UI frame.
  if ( mCurrentView ) {
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
  if ( mCurrentLayerCollection )
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
  if ( mCurrentView ) {
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
    UpdateLayerMenu();
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

void
vtkKWScubaWindow::UpdateViewMenu () {

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

  // Enable/disable our menus and toolbar buttons accordingly.
  mMenuLoadVolume.menu->SetItemStateToNormal( mMenuLoadVolume.nItem );
  mBtnLoadVolume->SetStateToNormal();
  mMenuZoomOut.menu->SetItemStateToNormal( mMenuZoomOut.nItem );
  mMenuZoomIn.menu->SetItemStateToNormal( mMenuZoomIn.nItem );

  if ( mCurrentView && mCurrentView->IsDataDirty() ) {
    mMenuSaveVolume.menu->SetItemStateToNormal( mMenuSaveVolume.nItem );
    mBtnSaveVolume->SetStateToNormal();
  } else {
    mMenuSaveVolume.menu->SetItemStateToDisabled( mMenuSaveVolume.nItem );
    mBtnSaveVolume->SetStateToDisabled();
  }

}

void
vtkKWScubaWindow::UpdateInfoArea () {

  if ( mCurrentView ) {

    // Make sure we have the tables...
    if ( mMouseOverInfoTable && mCursorInfoTable ) {

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

    }

    // Tell the view that we have updated our info.
    mCurrentView->InfoUpdated();
  }

}

vtkKWScubaView* 
vtkKWScubaWindow::GetNthView ( int inView ) {

  map<int,vtkKWScubaView*>::iterator tView;
  tView = maView.find( inView );
  if( tView != maView.end() )
    if( tView->second )
      return tView->second;

  vtkKWScubaView* view = vtkKWScubaView::New();
  view->SetParent( this->GetViewPanelFrame() );
  view->Create();

  maView[inView] = view;

  // We listen to the view. This means the view can call SendBroadcast
  // and we'll get it in DoListenToMessage.
  view->AddListener( this );

  // Create our custom interactor style (listener) for this
  // view. Associate it with the view and our tool. Tell the view's
  // window to use this style.
  vtkScubaInteractorStyle* style = vtkScubaInteractorStyle::New();
  style->SetWindow( this );
  view->GetRenderWindow()->GetInteractor()->SetInteractorStyle( style );

  // We'll get the layers from the current view and set them in this
  // view.
  if( NULL != mCurrentView ) {
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

  if( Single == mCurrentViewLayout ) {
    
    vtkKWScubaView* view = this->GetNthView( 0 );

    view->SetDisplayMode( vtkKWScubaView::TwoDee );
    view->Set2DRASZ( 0 );
    view->Set2DInPlane( 0 );

    this->SetCurrentView( *view );

  } else if( TwoByTwo == mCurrentViewLayout ) {

    vtkKWScubaView* topLeft = this->GetNthView( 0 );
    vtkKWScubaView* topRight = this->GetNthView( 1 );
    vtkKWScubaView* bottomLeft = this->GetNthView( 2 );
    vtkKWScubaView* bottomRight = this->GetNthView( 3 );

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
}

void
vtkKWScubaWindow::AddLayerCollectionToViews ( vtkKWScubaLayerCollection* col ) {
  
  // Go through all our views.
  map<int,vtkKWScubaView*>::iterator tView;
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
