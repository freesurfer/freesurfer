/**
 * @file  vtkKWQdecWindow.cxx
 * @brief Main QDEC logic
 *
 * Loads in all types of data and manages them. Main logic for running
 * analysis and handling results. Composes display objects and sends
 * it to the main view. Holds settings for color scales and currently
 * displayed objects.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2008/01/10 23:51:48 $
 *    $Revision: 1.3 $
 *
 * Copyright (C) 2007,
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

#include <fstream>
#include <vector>
#include <stdexcept>
#include <sstream>

#include "vtkKWQdecWindow.h"
#include "IconLoader.h"
#include "QdecEvents.h"
#include "QdecUtilities.h"
#include "vtkCallbackCommand.h"
#include "vtkCellArray.h"
#include "vtkColorTransferFunction.h"
#include "vtkDataArray.h"
#include "vtkFreesurferLookupTable.h"
#include "vtkFSSurfaceSource.h"
#include "vtkFSSurfaceLabelSource.h"
#include "vtkFloatArray.h"
#include "vtkKWApplication.h"
#include "vtkKWBltGraph.h"
#include "vtkKWCheckButton.h"
#include "vtkKWCheckButtonWithLabel.h"
#include "vtkKWEntry.h"
#include "vtkKWEntryWithLabel.h"
#include "vtkKWEvent.h"
#include "vtkKWFileBrowserDialog.h"
#include "vtkKWFrame.h"
#include "vtkKWFrameWithLabel.h"
#include "vtkKWFrameWithScrollbar.h"
#include "vtkKWHistogram.h"
#include "vtkKWIcon.h"
#include "vtkKWLabel.h"
#include "vtkKWListBox.h"
#include "vtkKWListBoxWithScrollbars.h"
#include "vtkKWListBoxWithScrollbarsWithLabel.h"
#include "vtkKWLoadSaveDialog.h"
#include "vtkKWLoadSaveButton.h"
#include "vtkKWMenu.h"
#include "vtkKWMenuButton.h"
#include "vtkKWMessageDialog.h"
#include "vtkKWMultiColumnList.h"
#include "vtkKWMultiColumnListWithScrollbars.h"
#include "vtkKWNotebook.h"
#include "vtkKWProgressGauge.h"
#include "vtkKWPushButton.h"
#include "vtkKWRadioButton.h"
#include "vtkKWRadioButtonSet.h"
#include "vtkKWRGBATransferFunctionEditor.h"
#include "vtkKWScale.h"
#include "vtkKWScaleWithLabel.h"
#include "vtkKWSimpleEntryDialog.h"
#include "vtkKWSplitFrame.h"
#include "vtkKWToolbar.h"
#include "vtkKWToolbarSet.h"
#include "vtkKWUserInterfacePanel.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkRenderLargeImage.h"
#include "vtkRGBATransferFunction.h"
#include "vtkSelectPolyData.h"
#include "vtkSmartPointer.h"
#include "vtkTIFFWriter.h"
#include "vtkWindowToImageFilter.h"
#include "QdecProject.h"

extern "C" {
#include "colortab.h"
#include "mrisutils.h"
#include "fsenv.h"
}

using namespace std;

vtkStandardNewMacro( vtkKWQdecWindow );
vtkCxxRevisionMacro( vtkKWQdecWindow, "$Revision: 1.3 $" );

const char* vtkKWQdecWindow::ksSubjectsPanelName = "Subjects";
const char* vtkKWQdecWindow::ksDesignPanelName = "Design";
const char* vtkKWQdecWindow::ksDisplayPanelName = "Display";

vtkKWQdecWindow::vtkKWQdecWindow () :
  vtkKWWindow(),
  mbUseHistogramEditor( true ),
  mMenuLoadDataTable( NULL ),
  mMenuLoadProjectFile( NULL ),
  mMenuLoadLabel( NULL ),
  mMenuSaveProjectFile( NULL ),
  mMenuSaveTIFF( NULL ),
  mMenuSaveGDFPostscript( NULL ),
  mMenuSaveLabel( NULL ),
  mMenuMapLabel( NULL ),
  mMenuClearCurvature( NULL ),
  mMenuClearSurfaceScalars( NULL ),
  mMenuRestoreView( NULL ),
  mMenuZoomOut( NULL ),
  mMenuZoomIn( NULL ),
  mMenuShowCursor( NULL ),
  mMenuAddSelectionToROI( NULL ),
  mMenuRemoveSelectionFromROI( NULL ),
  mMenuClearROI( NULL ),
  mMenuSmoothCurvatureScalars( NULL ),
  mMenuSmoothSurfaceScalars( NULL ),
  mMenuGraphAverageROI( NULL ),
  mbFrameSurfaceInited( false ),
  mbFrameCurvatureInited( false ),
  mbFrameOverlayInited( false ),
  mbFrameSurfaceScalarsInited( false ),
  mPlotContinuousFactorSelection( -1 ),
  mQdecProject( NULL ),
  mcVertices( -1 ),
  msCurrentSurfaceSource( "" ),
  mGDFID( -1 ),
  mbGDFLoaded( false ),
  mnCurrentSurfaceScalars( -1 ),
  maAnnotationIndicies( NULL ),
  mAnnotationTable( NULL ),
  mbShowCurvature( true ),
  mSurfaceScalarsColorMin( 2.0 ),
  mSurfaceScalarsColorMid( 4.0 ),
  mSurfaceScalarsColorMax( 5.0 ),
  mSurfaceScalarsColorOffset( 0.0 ),
  mbSurfaceScalarsColorReverse( false ),
  mbSurfaceScalarsColorShowPositive( true ),
  mbSurfaceScalarsColorShowNegative( true ),
  mbDrawCurvatureGreenRedIfNoScalars( false ),
  mSurfaceScalarsColorsFDRRate( 0.05 ),
  msOverlayDescription( "" ),
  mbViewInitialized( false ) {

  maDiscreteFactorSelection[0] = -1;
  maDiscreteFactorSelection[1] = -1;
  maContinuousFactorSelection[0] = -1;
  maContinuousFactorSelection[1] = -1;
}

vtkKWQdecWindow::~vtkKWQdecWindow () {
  
  delete mMenuLoadDataTable;
  delete mMenuLoadProjectFile;
  delete mMenuLoadLabel;
  delete mMenuSaveProjectFile;
  delete mMenuSaveTIFF;
  delete mMenuSaveGDFPostscript;
  delete mMenuSaveLabel;
  delete mMenuMapLabel;
  delete mMenuClearCurvature;
  delete mMenuClearSurfaceScalars;
  delete mMenuRestoreView;
  delete mMenuZoomOut;
  delete mMenuZoomIn;
  delete mMenuShowCursor;
  delete mMenuAddSelectionToROI;
  delete mMenuRemoveSelectionFromROI;
  delete mMenuClearROI;
  delete mMenuSmoothCurvatureScalars;
  delete mMenuSmoothSurfaceScalars;
  delete mMenuGraphAverageROI;
  delete mQdecProject;
}

void
vtkKWQdecWindow::SetUseHistogramEditor ( bool ibUse ) {

  if( this->IsCreated() )
    vtkErrorMacro( << "SetUseHistogramEditor() should be called before Create()" );

  mbUseHistogramEditor = ibUse;
}

void
vtkKWQdecWindow::CreateWidget () {

  this->SupportHelpOn();

  this->Superclass::CreateWidget();

  this->SetTitle( "qdec" );

  this->SetPanelLayoutToSecondaryBelowView();
  this->SecondaryPanelVisibilityOff();

  // Create the primary data storage and working class
  this->mQdecProject = new QdecProject();

  // Create our interior view.
  mView = vtkSmartPointer<vtkKWQdecView>::New();
  mView->SetParent( this->GetViewFrame() );
  mView->Create();

  // Create our interior graph.
  mGraph = vtkSmartPointer<vtkKWBltGraph>::New();
  mGraph->SetParent( this->GetViewFrame() );
  mGraph->Create();
  mGraph->SetLegendVisibleToOff();

  // Start with the graph packed.
  this->Script( "pack %s -expand yes -fill both -anchor c",
                mGraph->GetWidgetName() );

  // ---------------------------------------------------------------------
  // Make a toolbar.
  vtkSmartPointer<vtkKWToolbar> toolbar = 
    vtkSmartPointer<vtkKWToolbar>::New();
  toolbar->SetName( "Main" );
  toolbar->SetParent( GetMainToolbarSet()->GetToolbarsFrame() );
  toolbar->Create();
  this->GetMainToolbarSet()->AddToolbar( toolbar );

  // Buttons for the toolbar.
  // Load data table
  mBtnLoadDataTable.TakeReference( 
    this->MakeToolbarButton( toolbar,
                             "Load a .dat file into the Design tab",
                             this, "LoadDataTableFromDlog", "LoadDataTable" )
    );
    
  // Load project file
  mBtnLoadProjectFile.TakeReference(
    this->MakeToolbarButton( toolbar,
                             "Load a project file",
                             this, "LoadProjectFileFromDlog",
                             "LoadAnalyzedData" )
    );

  // Load label
  mBtnLoadLabel.TakeReference(
    this->MakeToolbarButton( toolbar,
                             "Load a label as an ROI",
                             this, "LoadLabelFromDlog",
                             "LoadLabel" )
    );

  this->AddSpacerToToolbar( toolbar );

  // Save TIFF
  mBtnSaveTIFF.TakeReference(
    this->MakeToolbarButton( toolbar,
                             "Save a picture of the view as a TIFF",
                             this, "SaveTIFFImageFromDlog", "SaveTIFF" )
    );

  // Save Label
  mBtnSaveLabel.TakeReference(
    this->MakeToolbarButton( toolbar,
                             "Save the current label",
                             this, "SaveLabelFromDlog", "SaveLabel" )
    );

  this->AddSpacerToToolbar( toolbar );

  // Restore View
  mBtnRestoreView.TakeReference(
    this->MakeToolbarButton( toolbar,
                             "Restore the camera to a lateral view",
                             this, "RestoreView", "RestoreView" )
    );

  this->AddSpacerToToolbar( toolbar );

  // Zoom out and in
  mBtnZoomOut.TakeReference(
    this->MakeToolbarButton( toolbar, "Zoom out a step",
                             this, "ZoomOut", "ZoomOut" )
    );
  mBtnZoomIn.TakeReference(
    this->MakeToolbarButton( toolbar, "Zoom in a step",
                             this, "ZoomIn", "ZoomIn" )
    );

  this->AddSpacerToToolbar( toolbar );

  // Rotate buttons.
  mBtnCameraElevateNegative.TakeReference(
    this->MakeToolbarButton( toolbar, "Rotate the surface (pitch forward)",
                             mView, "AnimateCameraElevateNegative",
                             "AnimateCameraElevateNegative" )
    );
  mBtnCameraElevatePositive.TakeReference(
    this->MakeToolbarButton( toolbar, "Rotate the surface (pitch backward)",
                             mView, "AnimateCameraElevatePositive",
                             "AnimateCameraElevatePositive" )
    );
  mBtnCameraAzimuthNegative.TakeReference(
    this->MakeToolbarButton( toolbar, "Rotate the surface (yaw ccw)",
                             mView, "AnimateCameraAzimuthNegative",
                             "AnimateCameraAzimuthNegative" )
    );
  mBtnCameraAzimuthPositive.TakeReference(
    this->MakeToolbarButton( toolbar, "Rotate the surface (yaw cw)",
                             mView, "AnimateCameraAzimuthPositive",
                             "AnimateCameraAzimuthPositive" )
    );
  mBtnCameraRollPositive.TakeReference(
    this->MakeToolbarButton( toolbar, "Rotate the surface (roll ccw)",
                             mView, "AnimateCameraRollPositive",
                             "AnimateCameraRollPositive" )
    );
  mBtnCameraRollNegative.TakeReference(
    this->MakeToolbarButton( toolbar, "Rotate the surface (roll cw)",
                             mView, "AnimateCameraRollNegative",
                             "AnimateCameraRollNegative" )
    );

  this->AddSpacerToToolbar( toolbar );

  // Show cursor
  mBtnShowCursor.TakeReference(
    this->MakeToolbarCheckButton( toolbar, "Show the cursor",
                                  this, "ShowCursor", "ShowCursor" )
    );
  if( mView )
    mBtnShowCursor->SetSelectedState( mView->GetShowCursor() );

  this->AddSpacerToToolbar( toolbar );

  // Goto vertex entry.
  vtkSmartPointer<vtkKWEntryWithLabel> labeledEntryVertex =
    vtkSmartPointer<vtkKWEntryWithLabel>::New();
  labeledEntryVertex->SetParent( toolbar->GetFrame() );
  labeledEntryVertex->Create();
  labeledEntryVertex->SetLabelText( "Jump to vertex number: " );
  toolbar->AddWidget( labeledEntryVertex );

  mEntrySelectVertex = labeledEntryVertex->GetWidget();
  mEntrySelectVertex->SetWidth( 6 );
  mEntrySelectVertex->SetRestrictValueToInteger();
  mEntrySelectVertex->SetValueAsInt( 0 );
  mEntrySelectVertex->SetCommandTrigger( vtkKWEntry::TriggerOnReturnKey );
  // I don't think this mapping is useful here, and it triggers it
  // when we switch window spaces, so undo it here.
  mEntrySelectVertex->
    RemoveBinding( "<Unmap>", mEntrySelectVertex, "ValueCallback" );

  this->AddSpacerToToolbar( toolbar );

  // Selection buttons.
  mBtnAddSelectionToROI.TakeReference(
    this->MakeToolbarButton( toolbar, "Add the selection to the ROI",
                             this, "AddSelectionToROI",
                             "AddSelectionToROI" )
    );

  mBtnRemoveSelectionFromROI.TakeReference(
    this->MakeToolbarButton( toolbar, "Remove the selection from the ROI",
                             this, "RemoveSelectionFromROI",
                             "RemoveSelectionFromROI" )
    );

  // ---------------------------------------------------------------------
  // Build the menus. File menu.
  // Load Data table..
  //  this->GetFileMenu()->SetFont( "helvetica 12" );
  int nItem = this->GetFileMenuInsertPosition();
  mMenuLoadDataTable = new MenuItem();
  mMenuLoadDataTable->
    MakeCommand( this->GetFileMenu(), nItem++,
                 "L&oad Data Table...", this, "LoadDataTableFromDlog",
                 "Ctrl+O", "LoadDataTable" );

  // Load Project File..
  mMenuLoadProjectFile = new MenuItem();
  mMenuLoadProjectFile->
    MakeCommand( this->GetFileMenu(), nItem++,
                 "Load Pro&ject File...", this, "LoadProjectFileFromDlog",
                 "Ctrl+J", "LoadAnalyzedData" );

  // Load Label
  mMenuLoadLabel = new MenuItem();
  mMenuLoadLabel->
    MakeCommand( this->GetFileMenu(), nItem++,
                 "Load &Label...", this, "LoadLabelFromDlog",
                 "Ctrl+L", "LoadLabel" );

  // These menu items are for loading pieces of data individuall, and
  // can be enabled for debugging.
#if 0
  this->GetFileMenu()->InsertCommand( nItem++, "Load Surface...",
                                      this, "LoadSurfaceFromDlog" );
  this->GetFileMenu()->InsertCommand( nItem++, "Load GDF...",
                                      this, "LoadGDFFromDlog" );
  this->GetFileMenu()->InsertCommand( nItem++, "Load Scalars...",
                                      this, "LoadSurfaceScalarsFromDlog" );
  this->GetFileMenu()->InsertCommand( nItem++, "Load Curvature...",
                                      this, "LoadCurvatureFromDlog" );
  this->GetFileMenu()->InsertCommand( nItem++, "Load Annotation...",
                                      this, "LoadAnnotationFromDlog" );
#endif

  this->GetFileMenu()->InsertSeparator( nItem++ );

  // Save Project File.
  mMenuSaveProjectFile = new MenuItem();
  mMenuSaveProjectFile->
    MakeCommand( this->GetFileMenu(), nItem++,
                 "Save Project File...", this, "SaveProjectFileFromDlog",
                 NULL, NULL );

  // Save TIFF
  mMenuSaveTIFF = new MenuItem();
  mMenuSaveTIFF->
    MakeCommand( this->GetFileMenu(), nItem++,
                 "Save TIFF...", this, "SaveTIFFImageFromDlog",
                 NULL, "SaveTIFF" );

  // Save Postcript of GDF
  mMenuSaveGDFPostscript = new MenuItem();
  mMenuSaveGDFPostscript->
    MakeCommand( this->GetFileMenu(), nItem++,
                 "Save GDF Plot as Postcript...", this,
                 "SaveGDFPostscriptFromDlog", NULL, NULL );

  // Save Label
  mMenuSaveLabel = new MenuItem();
  mMenuSaveLabel->
    MakeCommand( this->GetFileMenu(), nItem++,
                 "Save Label...", this, "SaveLabelFromDlog",
                 NULL, "SaveLabel" );

  // Map label
  mMenuMapLabel = new MenuItem();
  mMenuMapLabel->
    MakeCommand( this->GetFileMenu(), nItem++,
                 "Map Label to Subjects...", this, "MapLabelFromDlog",
                 NULL, NULL );
  
  // Insert a separator, then the recent files, then another
  // separator. The quit items will be after that.
  this->GetFileMenu()->InsertSeparator( nItem++ );

  this->InsertRecentFilesMenu( nItem++, this );

  this->GetFileMenu()->InsertSeparator( nItem++ );

  // Edit menu.
  nItem = 0;

  // Clear Curvatures.
  mMenuClearCurvature = new MenuItem();
  mMenuClearCurvature->
    MakeCommand( this->GetEditMenu(), nItem++,
                 "Clear Curvature", this, "ClearCurvature", NULL, NULL );

  // Clear scalars.
  mMenuClearSurfaceScalars = new MenuItem();
  mMenuClearSurfaceScalars->
    MakeCommand( this->GetEditMenu(), nItem++,
                 "Clear Surface Scalars", this, "ClearSurfaceScalars",
                 NULL, NULL );

  this->GetEditMenu()->InsertSeparator( nItem++ );

  // Add selection to ROI
  mMenuAddSelectionToROI = new MenuItem();
  mMenuAddSelectionToROI->
    MakeCommand( this->GetEditMenu(), nItem++,
                 "Add Selection to ROI", this, "AddSelectionToROI",
                 NULL, NULL );
  
  // Remove selection remove from ROI
  mMenuRemoveSelectionFromROI = new MenuItem();
  mMenuRemoveSelectionFromROI->
    MakeCommand( this->GetEditMenu(), nItem++,
                 "Remove Selection from ROI", this, "RemoveSelectionFromROI",
                 NULL, NULL );
  
  // Clear ROI
  mMenuClearROI = new MenuItem();
  mMenuClearROI->
    MakeCommand( this->GetEditMenu(), nItem++,
                 "Clear ROI", this, "ClearROI",
                 NULL, NULL );
  
  this->GetEditMenu()->InsertSeparator( nItem++ );

  // Smooth curvature.
  mMenuSmoothCurvatureScalars = new MenuItem();
  mMenuSmoothCurvatureScalars->
    MakeCommand( this->GetEditMenu(), nItem++, "Smooth Curvature Scalars", 
                 this, "SmoothCurvatureScalarsFromDlog", NULL, NULL );

  // Smooth surface scalars.
  mMenuSmoothSurfaceScalars = new MenuItem();
  mMenuSmoothSurfaceScalars->
    MakeCommand( this->GetEditMenu(), nItem++, "Smooth Surface Scalars", 
                 this, "SmoothSurfaceScalarsFromDlog", NULL, NULL );

  // View menu.
  // Restore view.
  nItem = this->GetViewMenuInsertPosition();
  mMenuRestoreView = new MenuItem();
  mMenuRestoreView->
    MakeCommand( this->GetViewMenu(), nItem++,
                 "Restore &View", this, "RestoreView", 
                 "Ctrl+V", "RestoreView" );

  this->GetViewMenu()->InsertSeparator( nItem++ );

  // Zoom Out.
  mMenuZoomOut = new MenuItem();
  mMenuZoomOut->
    MakeCommand( this->GetViewMenu(), nItem++,
                 "Zoom Out", this, "ZoomOut", NULL, "ZoomOut" );

  // Zoom In.
  mMenuZoomIn = new MenuItem();
  mMenuZoomIn->
    MakeCommand( this->GetViewMenu(), nItem++,
                 "Zoom In", this, "ZoomIn", NULL, "ZoomIn" );

  // Show Cursor
  mMenuShowCursor = new MenuItem();
  mMenuShowCursor->
    MakeCheckButton( this->GetViewMenu(), nItem++,
                     "Show Cursor",
                     this, "SetShowCursorFromMenu", NULL, "ShowCursor" );

  this->GetViewMenu()->InsertSeparator( nItem++ );

  // Graph average ROI.
  mMenuGraphAverageROI = new MenuItem();
  mMenuGraphAverageROI->
    MakeCommand( this->GetViewMenu(), nItem++,
                 "Graph Average ROI",
                 this, "GraphAverageROIInGDF", NULL, NULL );


  // ---------------------------------------------------------------------
  //
  // Create the main panel. We use the interface manager here which
  // implements basically a tabbed notebook. Then we create pages to
  // go in the manager, and get the frame for them. We'll pack all the
  // labeled frames inside those.
  //

  this->GetMainSplitFrame()->SetFrame1Size( 325 );

  mPanel = vtkSmartPointer<vtkKWUserInterfacePanel>::New();
  mPanel->SetUserInterfaceManager( this->GetMainUserInterfaceManager() );
  mPanel->Create();
  mPanel->AddPage(
    ksSubjectsPanelName, "Explore data from the subjects data table", NULL );
  mPanel->AddPage(
    ksDesignPanelName, "Design the GLM query and begin the analysis", NULL );
  mPanel->AddPage(
    ksDisplayPanelName, "Configure the view of the results", NULL );

  //
  // Use the Subjects pane and put some sample form stuff in it
  vtkKWWidget* panelFrame = mPanel->GetPageWidget( "Subjects" );

  // Make the top level frame in this page a scrolling frame.
  vtkSmartPointer<vtkKWFrameWithScrollbar> scrolledFrame = 
    vtkSmartPointer<vtkKWFrameWithScrollbar>::New();
  scrolledFrame->SetParent( panelFrame );
  scrolledFrame->Create();
  scrolledFrame->VerticalScrollbarVisibilityOn();
  scrolledFrame->HorizontalScrollbarVisibilityOff();
  this->Script( "pack %s -expand yes -fill both",
                scrolledFrame->GetWidgetName() );

  // Get the inner scrolled frame as the parent for everything else.
  panelFrame = scrolledFrame->GetFrame();

  // We're packing this frame in a grid, keep track of our rows.
  int nRow = 0;

  vtkSmartPointer<vtkKWEntryWithLabel> labeledEntry =
    vtkSmartPointer<vtkKWEntryWithLabel>::New();
  labeledEntry->SetParent( panelFrame );
  labeledEntry->SetLabelText( "SUBJECTS_DIR:" );
  labeledEntry->Create();
  labeledEntry->SetLabelPositionToTop();
  mEntrySubjectsDir = labeledEntry->GetWidget();
  mEntrySubjectsDir->SetCommandTrigger ( vtkKWEntry::TriggerOnReturnKey );
  mEntrySubjectsDir->RemoveBinding( "<Unmap>", 
                                    mEntrySubjectsDir, 
                                    "ValueCallback" );
  mEntrySubjectsDir->SetCommand ( this, "SetSubjectsDir" );
  if( this->mQdecProject )
  {
    string sSubjectsDir = this->mQdecProject->GetSubjectsDir();
    mEntrySubjectsDir->SetValue( sSubjectsDir.c_str() );
  }
  mEntrySubjectsDir->SetReadOnly( 0 );
  this->Script( "grid %s -column 0 -columnspan 2 -row %d -sticky new",
                labeledEntry->GetWidgetName(), nRow );
  nRow++;

  labeledEntry = vtkSmartPointer<vtkKWEntryWithLabel>::New();
  labeledEntry->SetParent( panelFrame );
  labeledEntry->SetLabelText( "Data Table:" );
  labeledEntry->Create();
  labeledEntry->SetLabelPositionToTop();
  mEntryDataTable = labeledEntry->GetWidget();
  mEntryDataTable->SetValue( "<not loaded>" );
  mEntryDataTable->SetReadOnly( 1 );
  this->Script( "grid %s -column 0 -columnspan 2 -row %d -sticky new",
                labeledEntry->GetWidgetName(), nRow );
  nRow++;

  labeledEntry = vtkSmartPointer<vtkKWEntryWithLabel>::New();
  labeledEntry->SetParent( panelFrame );
  labeledEntry->SetLabelText( "Number of Subjects:" );
  labeledEntry->Create();
  mEntryNumberOfSubjects = labeledEntry->GetWidget();
  mEntryNumberOfSubjects->SetValue( "0" );
  mEntryNumberOfSubjects->SetReadOnly( 1 );
  this->Script( "grid %s -column 0 -columnspan 2 -row %d -sticky new",
                labeledEntry->GetWidgetName(), nRow );
  nRow++;

  // Create the subject data scatter-plot exploration frame.
  vtkSmartPointer<vtkKWFrameWithLabel> exploreFrame = 
    vtkSmartPointer<vtkKWFrameWithLabel>::New();
  exploreFrame->SetParent( panelFrame );
  exploreFrame->Create();
  exploreFrame->SetLabelText( "Scatter Plot" );
  this->Script( "grid %s -column 0 -columnspan 2 -row %d -sticky new",
                exploreFrame->GetWidgetName(), nRow );
  nRow++;

  vtkSmartPointer<vtkKWListBoxWithScrollbarsWithLabel> listBox =
    vtkSmartPointer<vtkKWListBoxWithScrollbarsWithLabel>::New();
  listBox->SetParent( exploreFrame->GetFrame() );
  listBox->SetLabelText( "Continuous Factors:" );
  listBox->Create();
  listBox->SetLabelPositionToTop();
  mListPlotContinuousFactors = listBox->GetWidget()->GetWidget();
  mListPlotContinuousFactors->ExportSelectionOff();
  this->Script( "pack %s -side top -fill x -expand yes",
                listBox->GetWidgetName() );
  mListPlotContinuousFactors->
    SetSelectionCommand( this, "PlotContinuousFactorsListBoxCallback" );

  labeledEntry = vtkSmartPointer<vtkKWEntryWithLabel>::New();
  labeledEntry->SetParent( panelFrame );
  labeledEntry->SetLabelText( "Common-space Subject:" );
  labeledEntry->Create();
  mEntryAverageSubject = labeledEntry->GetWidget();
  mEntryAverageSubject->SetValue( "fsaverage" );
  mEntryAverageSubject->SetCommandTrigger ( vtkKWEntry::TriggerOnReturnKey );
  mEntryAverageSubject->RemoveBinding( "<Unmap>", 
                                       mEntryAverageSubject, 
                                       "ValueCallback" );
  mEntryAverageSubject->SetCommand ( this, "SetAverageSubject" );
  this->Script( "grid %s -column 0 -columnspan 2 -row %d -sticky new",
                labeledEntry->GetWidgetName(), nRow );
  nRow++;

  this->Script( "grid columnconfigure %s 0 -weight 0",
                panelFrame->GetWidgetName() );
  this->Script( "grid columnconfigure %s 1 -weight 1",
                panelFrame->GetWidgetName() );
  this->Script( "grid rowconfigure %s %d -weight 1",
                panelFrame->GetWidgetName(), nRow );
  for( int nRowConfigure = 0; nRowConfigure <= nRow; nRowConfigure++ )
    this->Script( "grid rowconfigure %s %d -pad 4",
                  panelFrame->GetWidgetName(), nRowConfigure );

  // ---------------------------------------------------------------------
  //
  // Design pane gets the control for configuring the input to
  // mri_glmfit.
  panelFrame = mPanel->GetPageWidget( "Design" );

  // Make the top level frame in this page a scrolling frame.
  scrolledFrame = vtkSmartPointer<vtkKWFrameWithScrollbar>::New();
  scrolledFrame->SetParent( panelFrame );
  scrolledFrame->Create();
  scrolledFrame->VerticalScrollbarVisibilityOn();
  scrolledFrame->HorizontalScrollbarVisibilityOff();
  this->Script( "pack %s -expand yes -fill both",
                scrolledFrame->GetWidgetName() );

  // Get the inner scrolled frame as the parent for everything else.
  panelFrame = scrolledFrame->GetFrame();

  labeledEntry = vtkSmartPointer<vtkKWEntryWithLabel>::New();
  labeledEntry->SetParent( panelFrame );
  labeledEntry->SetLabelText( "Design Name:" );
  labeledEntry->Create();
  mEntryDesignName = labeledEntry->GetWidget();
  mEntryDesignName->SetValue( "Untitled" );
  mEntryDesignName->SetCommandTrigger ( vtkKWEntry::TriggerOnReturnKey );
  mEntryDesignName->RemoveBinding( "<Unmap>", 
                                   mEntryDesignName, 
                                   "ValueCallback" );
  mEntryDesignName->SetCommand ( this, "SetDesignName" );
  this->Script( "pack %s -fill x", labeledEntry->GetWidgetName() );

  //
  // Create the Factors frame.
  vtkSmartPointer<vtkKWFrameWithLabel> factorsFrame =
    vtkSmartPointer<vtkKWFrameWithLabel>::New();
  factorsFrame->SetParent( panelFrame );
  factorsFrame->Create();
  factorsFrame->SetLabelText( "Factors" );
  this->Script( "pack %s -fill x", factorsFrame->GetWidgetName() );

  // Discrete and continuous list boxes are inside the Factors frame.
  listBox = vtkSmartPointer<vtkKWListBoxWithScrollbarsWithLabel>::New();
  listBox->SetParent( factorsFrame->GetFrame() );
  listBox->SetLabelText( "Discrete (choose up to two):" );
  listBox->Create();
  listBox->SetLabelPositionToTop();
  mListDiscreteFactors = listBox->GetWidget()->GetWidget();
  mListDiscreteFactors->ExportSelectionOff();
  mListDiscreteFactors->SetHeight( 3 );
  this->Script( "pack %s -fill x -expand y", listBox->GetWidgetName() );
  mListDiscreteFactors->SetSelectionModeToMultiple();
  mListDiscreteFactors->
    SetSelectionCommand( this, "DiscreteFactorsListBoxCallback" );

  listBox = vtkSmartPointer<vtkKWListBoxWithScrollbarsWithLabel>::New();
  listBox->SetParent( factorsFrame->GetFrame() );
  listBox->SetLabelText( "Continuous (choose up to two):" );
  listBox->Create();
  listBox->SetLabelPositionToTop();
  mListContinuousFactors = listBox->GetWidget()->GetWidget();
  mListContinuousFactors->ExportSelectionOff();
  this->Script( "pack %s -fill x -expand y", listBox->GetWidgetName() );
  mListContinuousFactors->SetSelectionModeToMultiple();
  mListContinuousFactors->
    SetSelectionCommand( this, "ContinuousFactorsListBoxCallback" );

  //
  // Create the Measures frame, below the Factors frame.
  vtkSmartPointer<vtkKWFrameWithLabel> measuresFrame =
    vtkSmartPointer<vtkKWFrameWithLabel>::New();
  measuresFrame->SetParent( panelFrame );
  measuresFrame->Create();
  measuresFrame->SetLabelText( "Measures" );
  this->Script( "pack %s -fill x", measuresFrame->GetWidgetName() );

  //
  // Create the Surface-based frame, inside the Measures frame.
  vtkSmartPointer<vtkKWFrameWithLabel> surfaceMeasuresFrame = 
    vtkSmartPointer<vtkKWFrameWithLabel>::New();
  surfaceMeasuresFrame->SetParent( measuresFrame->GetFrame() );
  surfaceMeasuresFrame->Create();
  surfaceMeasuresFrame->SetLabelText( "Surface-based" );
  this->Script( "pack %s -fill x", surfaceMeasuresFrame->GetWidgetName() );

  // The widgets inside the Surface-based Measures frame. We'll pack
  // these in a grid so we can get the actual menus as wide as
  // possible and aligned nicely.
  vtkSmartPointer<vtkKWFrameWithLabel> morphMeasuresFrame =
    vtkSmartPointer<vtkKWFrameWithLabel>::New();
  morphMeasuresFrame->SetParent( surfaceMeasuresFrame->GetFrame() );
  morphMeasuresFrame->Create();
  morphMeasuresFrame->SetLabelText( "Morphometric" );
  this->Script( "pack %s -fill x", morphMeasuresFrame->GetWidgetName() );

  nRow = 0;
  vtkSmartPointer<vtkKWLabel> label = 
    vtkSmartPointer<vtkKWLabel>::New();
  label->SetParent( morphMeasuresFrame->GetFrame() );
  label->Create();
  label->SetText( "Measure: " );
  label->SetJustificationToRight();
  this->Script( "grid %s -column 0 -row %d -sticky ne",
                label->GetWidgetName(), nRow );

  mMenuMeasure = vtkSmartPointer<vtkKWMenuButton>::New();
  mMenuMeasure->SetParent( morphMeasuresFrame->GetFrame() );
  mMenuMeasure->Create();
  mMenuMeasure->GetMenu()->AddRadioButton( "thickness" );
  mMenuMeasure->GetMenu()->AddRadioButton( "area" );
  mMenuMeasure->GetMenu()->AddRadioButton( "area.pial" );
  mMenuMeasure->GetMenu()->AddRadioButton( "volume" );
  mMenuMeasure->GetMenu()->AddRadioButton( "sulc" );
  mMenuMeasure->GetMenu()->AddRadioButton( "curv" );
  mMenuMeasure->GetMenu()->AddRadioButton( "white.K" );
  mMenuMeasure->GetMenu()->AddRadioButton( "white.H" );
  mMenuMeasure->GetMenu()->AddRadioButton( "jacobian_white" );
  mMenuMeasure->SetValue( "thickness" );
  this->Script( "grid %s -column 1 -row %d -sticky nw",
                mMenuMeasure->GetWidgetName(), nRow );

  nRow++;
  label = vtkSmartPointer<vtkKWLabel>::New();
  label->SetParent( morphMeasuresFrame->GetFrame() );
  label->Create();
  label->SetText( "Hemisphere: " );
  label->SetJustificationToRight();
  this->Script( "grid %s -column 0 -row %d -sticky ne",
                label->GetWidgetName(), nRow );

  mMenuHemisphere = vtkSmartPointer<vtkKWMenuButton>::New();
  mMenuHemisphere->SetParent( morphMeasuresFrame->GetFrame() );
  mMenuHemisphere->Create();
  mMenuHemisphere->GetMenu()->AddRadioButton( "lh" );
  mMenuHemisphere->GetMenu()->AddRadioButton( "rh" );
  mMenuHemisphere->SetValue( "lh" );
  this->Script( "grid %s -column 1 -row %d -sticky nw",
                mMenuHemisphere->GetWidgetName(), nRow );

  nRow++;
  label = vtkSmartPointer<vtkKWLabel>::New();
  label->SetParent( morphMeasuresFrame->GetFrame() );
  label->Create();
  label->SetText( "Smoothing (FWHM): " );
  label->SetJustificationToRight();
  this->Script( "grid %s -column 0 -row %d -sticky ne",
                label->GetWidgetName(), nRow );

  mMenuSmoothness = vtkSmartPointer<vtkKWMenuButton>::New();
  mMenuSmoothness->SetParent( morphMeasuresFrame->GetFrame() );
  mMenuSmoothness->Create();
  mMenuSmoothness->GetMenu()->AddRadioButton( "0" );
  mMenuSmoothness->GetMenu()->AddRadioButton( "5" );
  mMenuSmoothness->GetMenu()->AddRadioButton( "10" );
  mMenuSmoothness->GetMenu()->AddRadioButton( "15" );
  mMenuSmoothness->GetMenu()->AddRadioButton( "20" );
  mMenuSmoothness->GetMenu()->AddRadioButton( "25" );
  mMenuSmoothness->SetValue( "10" );
  this->Script( "grid %s -column 1 -row %d -sticky nw",
                mMenuSmoothness->GetWidgetName(), nRow );

  // Weight our grid properly.
  this->Script( "grid columnconfigure %s 0 -weight 0",
                morphMeasuresFrame->GetFrame()->GetWidgetName() );
  this->Script( "grid columnconfigure %s 1 -weight 1",
                morphMeasuresFrame->GetFrame()->GetWidgetName() );
  this->Script( "grid rowconfigure %s %d -weight 1",
                morphMeasuresFrame->GetFrame()->GetWidgetName(), nRow );
  for( int nRowConfigure = 0; nRowConfigure <= nRow; nRowConfigure++ )
    this->Script( "grid rowconfigure %s %d -pad 4",
                  surfaceMeasuresFrame->GetFrame()->GetWidgetName(),
                  nRowConfigure );

  //
  // Create the Functional based frame, inside the Measures frame.
  vtkSmartPointer<vtkKWFrameWithLabel> functionalMeasuresFrame = 
    vtkSmartPointer<vtkKWFrameWithLabel>::New();
  functionalMeasuresFrame->SetParent( surfaceMeasuresFrame->GetFrame() );
  functionalMeasuresFrame->Create();
  functionalMeasuresFrame->SetLabelText( "Functional" );
  functionalMeasuresFrame->CollapseFrame();
  this->Script( "pack %s -fill x", functionalMeasuresFrame->GetWidgetName() );

  //
  // Create the Volume-based frame, inside the Measures frame.
  vtkSmartPointer<vtkKWFrameWithLabel> volumeMeasuresFrame =
    vtkSmartPointer<vtkKWFrameWithLabel>::New();
  volumeMeasuresFrame->SetParent( measuresFrame->GetFrame() );
  volumeMeasuresFrame->Create();
  volumeMeasuresFrame->SetLabelText( "Volume-based" );
  volumeMeasuresFrame->CollapseFrame();
  this->Script( "pack %s -fill x", volumeMeasuresFrame->GetWidgetName() );

  //
  // Now for the 'go' button (to start analysis)
  vtkSmartPointer<vtkKWPushButton> button =
    vtkSmartPointer<vtkKWPushButton>::New();
  button->SetParent( panelFrame );
  button->SetText( "Analyze" );
  button->Create();
  button->SetWidth( 20 );
  button->SetCommand( this, "AnalyzeDesign" );
  this->Script( "pack %s -pady 20", button->GetWidgetName() );


  // ---------------------------------------------------------------------
  //
  // Display pane gets the settings for the surface display.
  panelFrame = mPanel->GetPageWidget( "Display" );


  // Make the top level frame in this page a scrolling frame.
  scrolledFrame = vtkSmartPointer<vtkKWFrameWithScrollbar>::New();
  scrolledFrame->SetParent( panelFrame );
  scrolledFrame->Create();
  scrolledFrame->VerticalScrollbarVisibilityOn();
  scrolledFrame->HorizontalScrollbarVisibilityOff();
  this->Script( "pack %s -expand yes -fill both",
                scrolledFrame->GetWidgetName() );

  // Get the inner scrolled frame as the parent for everything else.
  panelFrame = scrolledFrame->GetFrame();

  // We create four placeholder frames here, which we'll populate
  // later in the UpdateDisplayPage() function later on if we get data
  // for them.

  //
  // Surface frame. This is a placeholder frame in which we'll stuff a
  // label frame and a few radio buttons if we load more than one
  // surface. Otherwise it just remains empty.
  mFrameSurface = vtkSmartPointer<vtkKWFrame>::New();
  mFrameSurface->SetParent( panelFrame );
  mFrameSurface->Create();
  mbFrameSurfaceInited = false;

  // Curvature frame. This is a placeholder frame in which we'll stuff
  // a checkbox for controlling curvature visibility if we load
  // curvature. Leave it empty for now.
  mFrameCurvature = vtkSmartPointer<vtkKWFrame>::New();
  mFrameCurvature->SetParent( panelFrame );
  mFrameCurvature->Create();
  mbFrameCurvatureInited = false;

  //
  // Overlay frame. This is a placeholder frame in which we'll stuff a
  // slider for controlling the overlay opacity if we load
  // one. Otherwise it just remains empty.
  mFrameOverlay = vtkSmartPointer<vtkKWFrame>::New();
  mFrameOverlay->SetParent( panelFrame );
  mFrameOverlay->Create();
  mbFrameOverlayInited = false;

  // Scalars frame. This is a placeholder frame which will contain the
  // scalars controls, including the listbox, checkboxes, and the
  // color table editor.
  mFrameSurfaceScalars = vtkSmartPointer<vtkKWFrame>::New();
  mFrameSurfaceScalars->SetParent( panelFrame );
  mFrameSurfaceScalars->Create();
  mbFrameSurfaceScalarsInited = false;

  this->Script( "pack %s %s %s %s -side top -fill x -pady 5",
                mFrameSurface->GetWidgetName(),
                mFrameCurvature->GetWidgetName(),
                mFrameOverlay->GetWidgetName(),
                mFrameSurfaceScalars->GetWidgetName() );

  // ---------------------------------------------------------------------
  //
  // Make a callback and observe our view for UserSelectedVertex events.
  vtkSmartPointer<vtkCallbackCommand> callback =
    vtkSmartPointer<vtkCallbackCommand>::New();
  callback->SetCallback( UserSelectedVertexCallback );
  callback->SetClientData( this );
  mView->AddObserver( QdecEvents::UserSelectedVertex, callback );

  // Make a callback for our notebook.
  callback = vtkSmartPointer<vtkCallbackCommand>::New();
  callback->SetClientData( this );
  callback->SetCallback( NotebookRaisePageCallback );
  this->GetMainNotebook()->AddObserver( vtkKWEvent::NotebookRaisePageEvent, 
                                        callback );

  // Add callbacks for the graph.
  callback = vtkSmartPointer<vtkCallbackCommand>::New();
  callback->SetClientData( this );
  callback->SetCallback( ContinuousPlotGraphMouseoverEnterElementCallback );
  mGraph->AddObserver( vtkKWBltGraph::MouseoverEnterElementEvent, callback );

  callback = vtkSmartPointer<vtkCallbackCommand>::New();
  callback->SetClientData( this );
  callback->SetCallback( ContinuousPlotGraphMouseoverExitElementCallback );
  mGraph->AddObserver( vtkKWBltGraph::MouseoverExitElementEvent, callback );

  callback = vtkSmartPointer<vtkCallbackCommand>::New();
  callback->SetClientData( this );
  callback->SetCallback( ContinuousPlotGraphContextualMenuOpeningCallback );
  mGraph->AddObserver( vtkKWBltGraph::ContextualMenuOpening, callback );

}

void
vtkKWQdecWindow::FinishCreating () {

  // If we set up our menu item visibilities in CreateWidget(), they
  // are overridden, and I can't find a good way to get around that
  // (PostCreate doesn't work beucase the superclass's CreateWidget
  // calls PostCreate before CreateWidget is done...???). So the app
  // calls this to let us finish setting up our menu items.

  // Update our menu and buttons.
  this->UpdateCommandStatus();
}

// ------------------------------------
// Implement ProgressUpdateGUI virtuals
void 
vtkKWQdecWindow::BeginActionWithProgress( const char* isTitle ) {
  this->SetStatusText( isTitle );
}
void 
vtkKWQdecWindow::UpdateProgressMessage( const char* isMessage ) {
  this->SetStatusText( isMessage );
}
void
vtkKWQdecWindow::UpdateProgressPercent( float iPercent ) {
  this->GetProgressGauge()->SetValue( iPercent );
}
void
vtkKWQdecWindow::EndActionWithProgress() {
}
// ------------------------------------



void
vtkKWQdecWindow::LoadDataTableFromDlog () {

  vtkSmartPointer<vtkKWLoadSaveDialog> dialog =
    vtkSmartPointer<vtkKWLoadSaveDialog>::New();
  dialog->SetApplication( this->GetApplication() );
  dialog->Create();
  dialog->SetTitle( "Load Data Table" );
  dialog->SetFileTypes( "{{Data file} {*.dat}} {{All} {*}}" );
  if( this->mQdecProject )
  {
    string sWorkingDir = this->mQdecProject->GetWorkingDir();
    dialog->SetLastPath( sWorkingDir.c_str() );
    string sSubjectsDir = this->mQdecProject->GetSubjectsDir();
    mEntrySubjectsDir->SetValue( sSubjectsDir.c_str() );
  }
  if( dialog->Invoke() ) {
    string fnDataTable( dialog->GetFileName() );
    this->LoadDataTable( fnDataTable.c_str() );
  }
}

void
vtkKWQdecWindow::LoadProjectFileFromDlog () {

  vtkSmartPointer<vtkKWLoadSaveDialog> dialog =
    vtkSmartPointer<vtkKWLoadSaveDialog>::New();
  dialog->SetApplication( this->GetApplication() );
  dialog->Create();
  dialog->SetTitle( "Load Project File" );
  dialog->SetFileTypes( "{{Qdec Project file} {*.qdec}} {{All} {*}}" );

  // Set the default dir to the subjects dir / qdec.
  string sQdecDir = this->mQdecProject->GetSubjectsDir() + "/qdec";
  dialog->SetLastPath( sQdecDir.c_str() );

  if( dialog->Invoke() ) {
    string fnProject( dialog->GetFileName() );
    this->LoadProjectFile( fnProject.c_str() );
  }

  //  dialog->Delete();
}

void
vtkKWQdecWindow::LoadLabelFromDlog () {

  vtkSmartPointer<vtkKWLoadSaveDialog> dialog =
    vtkSmartPointer<vtkKWLoadSaveDialog>::New();
  dialog->SetApplication( this->GetApplication() );
  dialog->Create();
  dialog->SetTitle( "Load Label" );
  dialog->SetFileTypes( "{{Label file} {*.label}} {{All} {*}}" );
  dialog->RetrieveLastPathFromRegistry( "LoadLabel" );
  if( dialog->Invoke() ) {
    dialog->SaveLastPathToRegistry( "LoadLabel" );
    string fnLabel( dialog->GetFileName() );
    this->LoadLabel( fnLabel.c_str() );
  }
}

void
vtkKWQdecWindow::LoadSurfaceFromDlog () {

  vtkSmartPointer<vtkKWLoadSaveDialog> dialog =
    vtkSmartPointer<vtkKWLoadSaveDialog>::New();
  dialog->SetApplication( this->GetApplication() );
  dialog->Create();
  dialog->SetTitle( "Load Surface" );
  dialog->SetFileTypes( "{{Surface} {lh.* rh.*}} {{All} {*}}" );
  dialog->RetrieveLastPathFromRegistry( "LoadSurface" );
  if( dialog->Invoke() ) {
    dialog->SaveLastPathToRegistry( "LoadSurface" );
    string fn( dialog->GetFileName() );
    this->LoadSurface( fn.c_str() );
  }
}

void
vtkKWQdecWindow::LoadGDFFromDlog () {

  vtkSmartPointer<vtkKWLoadSaveDialog> dialog =
    vtkSmartPointer<vtkKWLoadSaveDialog>::New();
  dialog->SetApplication( this->GetApplication() );
  dialog->Create();
  dialog->SetTitle( "Load GDF" );
  dialog->SetFileTypes( "{{FSGDF} {*.fsgd}} {{All} {*}}" );
  dialog->RetrieveLastPathFromRegistry( "LoadGDF" );
  if( dialog->Invoke() ) {
    dialog->SaveLastPathToRegistry( "LoadGDF" );
    string fn( dialog->GetFileName() );
    this->LoadGDFFile( fn.c_str() );
  }
}

void
vtkKWQdecWindow::LoadSurfaceScalarsFromDlog () {

  vtkSmartPointer<vtkKWLoadSaveDialog> dialog =
    vtkSmartPointer<vtkKWLoadSaveDialog>::New();
  dialog->SetApplication( this->GetApplication() );
  dialog->Create();
  dialog->SetTitle( "Load Scalars" );
  dialog->SetFileTypes( 
    "{{Volume encoded scalars} {*.mgh *.mgz}} {{All} {*}}" );
  dialog->RetrieveLastPathFromRegistry( "LoadSurfaceScalars" );
  if( dialog->Invoke() ) {
    dialog->SaveLastPathToRegistry( "LoadSurfaceScalars" );
    string fn( dialog->GetFileName() );
    this->LoadSurfaceScalars( fn.c_str() );
  }
}

void
vtkKWQdecWindow::LoadCurvatureFromDlog () {

  vtkSmartPointer<vtkKWLoadSaveDialog> dialog =
    vtkSmartPointer<vtkKWLoadSaveDialog>::New();
  dialog->SetApplication( this->GetApplication() );
  dialog->Create();
  dialog->SetTitle( "Load Curvature" );
  dialog->SetFileTypes( "{{Curvature} {*.curv}} {{All} {*}}" );
  dialog->RetrieveLastPathFromRegistry( "LoadCurvature" );
  if( dialog->Invoke() ) {
    dialog->SaveLastPathToRegistry( "LoadCurvature" );
    string fn( dialog->GetFileName() );
    this->LoadSurfaceCurvatureScalars( fn.c_str() );
  }
}

void
vtkKWQdecWindow::LoadAnnotationFromDlog () {

  vtkSmartPointer<vtkKWLoadSaveDialog> dialog =
    vtkSmartPointer<vtkKWLoadSaveDialog>::New();
  dialog->SetApplication( this->GetApplication() );
  dialog->Create();
  dialog->SetTitle( "Load Annotation" );
  dialog->SetFileTypes( "{{Annotation} {*.annot}} {{All} {*}}" );
  dialog->RetrieveLastPathFromRegistry( "LoadAnnotation" );
  if( dialog->Invoke() ) {
    dialog->SaveLastPathToRegistry( "LoadAnnotation" );
    string fn( dialog->GetFileName() );
    this->LoadAnnotation( fn.c_str() );
  }
}

void
vtkKWQdecWindow::SaveProjectFileFromDlog () {

  assert( mQdecProject );
  assert( mQdecProject->GetGlmDesign() );

  vtkSmartPointer<vtkKWLoadSaveDialog> dialog =
    vtkSmartPointer<vtkKWLoadSaveDialog>::New();
  dialog->SetApplication( this->GetApplication() );
  dialog->SaveDialogOn();
  dialog->Create();
  dialog->SetTitle( "Save Project File" );
  dialog->SetFileTypes( "{{Qdec Project file} {*.qdec}} {{All} {*}}" );
  dialog->SetDefaultExtension( ".qdec" );

  // Make a file name based on the analysis name.
  string fnDefault = mQdecProject->GetGlmDesign()->GetName() + ".qdec";
  dialog->SetInitialFileName( fnDefault.c_str() );

  // Set the default dir to the subjects dir / qdec.
  string sQdecDir = this->mQdecProject->GetSubjectsDir() + "/qdec";
  dialog->SetLastPath( sQdecDir.c_str() );

  if( dialog->Invoke() ) {
    string fnProject( dialog->GetFileName() );
    this->SaveProjectFile( fnProject.c_str() );
  }
}

void
vtkKWQdecWindow::SaveTIFFImageFromDlog () {

  vtkSmartPointer<vtkKWLoadSaveDialog> dialog =
    vtkSmartPointer<vtkKWLoadSaveDialog>::New();
  dialog->SetApplication( this->GetApplication() );
  dialog->SaveDialogOn();
  dialog->Create();
  dialog->SetTitle( "Save TIFF Image" );
  dialog->SetFileTypes( "{{All} {*}}" );
  dialog->RetrieveLastPathFromRegistry( "SaveTIFFImage" );

  // Turn on the visibility frame so we can add the widgets to specify
  // the magnification rate.
  dialog->PreviewFrameVisibilityOn();

  vtkSmartPointer<vtkKWFrameWithLabel> labeledFrame = 
    vtkSmartPointer<vtkKWFrameWithLabel>::New();
  labeledFrame->SetParent( dialog->GetPreviewFrame() );
  labeledFrame->Create();
  labeledFrame->SetLabelText( "Magnification" );

  vtkSmartPointer<vtkKWLabel> sizeLabel =
    vtkSmartPointer<vtkKWLabel>::New();
  sizeLabel->SetParent( labeledFrame->GetFrame() );
  sizeLabel->Create();
  stringstream ssCurrentSize;
  ssCurrentSize << "Current size: " << mView->GetWidth() 
                << ", " << mView->GetHeight() << endl;
  sizeLabel->SetText( ssCurrentSize.str().c_str() );

  vtkSmartPointer<vtkKWEntryWithLabel> labeledEntry = 
    vtkSmartPointer<vtkKWEntryWithLabel>::New();
  labeledEntry->SetParent( labeledFrame->GetFrame() );
  labeledEntry->Create();
  labeledEntry->SetLabelText( "Magnification level: " );
  labeledEntry->GetWidget()->SetRestrictValueToInteger();
  labeledEntry->GetWidget()->SetValueAsInt( 1 );

  this->Script( "pack %s %s -side top -anchor w",
                sizeLabel->GetWidgetName(),
                labeledEntry->GetWidgetName() );

  this->Script( "pack %s -side top -fill x",
                labeledFrame->GetWidgetName() );

  if( dialog->Invoke() ) {
    dialog->SaveLastPathToRegistry( "SaveTIFFImage" );
    string fnTIFF( dialog->GetFileName() );
    this->SaveTIFFImage( fnTIFF.c_str(), 
                         labeledEntry->GetWidget()->GetValueAsInt() );
  }
}

void
vtkKWQdecWindow::SaveGDFPostscriptFromDlog () {

  assert( mbGDFLoaded );
  assert( mGDFID >= 0  );
  assert( this->GetApplication()->
          EvaluateBooleanExpression( "FsgdfPlot_IsWindowShowing %d", mGDFID ));

  vtkSmartPointer<vtkKWLoadSaveDialog> dialog =
    vtkSmartPointer<vtkKWLoadSaveDialog>::New();
  dialog->SetApplication( this->GetApplication() );
  dialog->SaveDialogOn();
  dialog->Create();
  dialog->SetTitle( "Save GDF Plot" );

  dialog->SetFileTypes( "{{All} {*}}" );
  dialog->RetrieveLastPathFromRegistry( "SaveGDF" );
  if( dialog->Invoke() ) {
    dialog->SaveLastPathToRegistry( "SaveGDF" );
    string fnPostscript( dialog->GetFileName() );
    this->Script( "FsgdfPlot_SaveToPostscript %d %s", 
                  mGDFID, fnPostscript.c_str() );
  }
}

void
vtkKWQdecWindow::SaveLabelFromDlog () {

  assert( mROISource );

  vtkSmartPointer<vtkKWLoadSaveDialog> dialog =
    vtkSmartPointer<vtkKWLoadSaveDialog>::New();
  dialog->SetApplication( this->GetApplication() );
  dialog->SaveDialogOn();
  dialog->Create();
  dialog->SetTitle( "Save Label" );

  string fnLabel( mROISource->GetLabelFileName() );
  string::size_type lastSlash = fnLabel.rfind( "/" );
  fnLabel = fnLabel.substr( lastSlash+1, string::npos );

  dialog->SetInitialFileName( fnLabel.c_str() );
  dialog->SetFileTypes( "{{All} {*}}" );
  dialog->RetrieveLastPathFromRegistry( "LoadLabel" );
  if( dialog->Invoke() ) {
    dialog->SaveLastPathToRegistry( "LoadLabel" );
    string fnLabel( dialog->GetFileName() );
    this->SaveLabel( fnLabel.c_str() );
  }
}

void
vtkKWQdecWindow::MapLabelFromDlog () {

  assert( mROISource );

  vtkKWSimpleEntryDialog* dialog = vtkKWSimpleEntryDialog::New();
  dialog->SetApplication( this->GetApplication() );
  dialog->SetStyleToOkCancel();
  dialog->Create();
  dialog->SetTitle( "Map Label" );
  dialog->SetText( "This will save the current label and map it to each "
                   "indiviudal subject's space, saving it in each one's "
                   "label directory." );
  dialog->SetOKButtonText( "Save" );

  dialog->GetEntry()->SetLabelText( "Label name: " );

  if( dialog->Invoke() ) {
    string fnLabel = dialog->GetEntry()->GetWidget()->GetValue();
    this->MapLabelToSubjects( fnLabel.c_str() );
  }
}

void
vtkKWQdecWindow::SmoothCurvatureScalarsFromDlog () {

  assert( mCurvatureScalars );

  vtkKWSimpleEntryDialog* dialog = vtkKWSimpleEntryDialog::New();
  dialog->SetApplication( this->GetApplication() );
  dialog->SetStyleToOkCancel();
  dialog->Create();
  dialog->SetTitle( "Smooth Curvature" );
  dialog->SetText( "This will smooth the curvature values using the surface "
                   "topography by the number of steps given." );
  dialog->SetOKButtonText( "Smooth" );

  dialog->GetEntry()->SetLabelText( "Number of steps: " );
  dialog->GetEntry()->GetWidget()->SetRestrictValueToInteger();
  dialog->GetEntry()->GetWidget()->SetValueAsInt( 1 );

  if( dialog->Invoke() ) {
    int cSteps = dialog->GetEntry()->GetWidget()->GetValueAsInt();
    this->SmoothCurvatureScalars( cSteps );
  }

}

void
vtkKWQdecWindow::SmoothSurfaceScalarsFromDlog () {

  assert( maSurfaceScalars.find(mnCurrentSurfaceScalars) !=
          maSurfaceScalars.end() );
  assert( maSurfaceScalars[mnCurrentSurfaceScalars].mValues );

  vtkKWSimpleEntryDialog* dialog = vtkKWSimpleEntryDialog::New();
  dialog->SetApplication( this->GetApplication() );
  dialog->SetStyleToOkCancel();
  dialog->Create();
  dialog->SetTitle( "Smooth Scalars" );
  dialog->SetText( "This will smooth the current scalars using the surface "
                   "topography by the number of steps given." );
  dialog->SetOKButtonText( "Smooth" );

  dialog->GetEntry()->SetLabelText( "Number of steps: " );
  dialog->GetEntry()->GetWidget()->SetRestrictValueToInteger();
  dialog->GetEntry()->GetWidget()->SetValueAsInt( 1 );

  if( dialog->Invoke() ) {
    int cSteps = dialog->GetEntry()->GetWidget()->GetValueAsInt();
    this->SmoothSurfaceScalars( cSteps );
  }
}

void
vtkKWQdecWindow::LoadDataTable ( const char* ifnDataTable ) {

  if( mView ) {
    try {

      if( this->mQdecProject->LoadDataTable( ifnDataTable ) )
        throw runtime_error( "Error loading the data table." );
      else
        this->mQdecProject->DumpDataTable( stdout );

      // data table may have set SUBJECTS_DIR, so update our mEntry display
      string sSubjectsDir = this->mQdecProject->GetSubjectsDir();
      mEntrySubjectsDir->SetValue( sSubjectsDir.c_str() );

      // We need to update our tabs.
      this->UpdateSubjectsPage();
      this->UpdateDesignPage();

      this->SetStatusText( "Data table loaded." );
      this->AddRecentFile( ifnDataTable, this, "LoadDataTable" );

    } catch ( exception& e ) {
      this->GetApplication()->ErrorMessage( e.what() );
    }
  }

  // Update our menu and buttons.
  this->UpdateCommandStatus();
}

void
vtkKWQdecWindow::LoadProjectFile ( const char* ifnProject ) {

  if( mView ) {
    try {
      if( ifnProject ) { 

        this->SetStatusText( "Loading project file..." );

        if( this->mQdecProject->LoadProjectFile( ifnProject ) )
          throw runtime_error( "Error loading the project file." );

        // the project file always has a subjects dir, so make sure the GUI
        // Subjects Display panel has that updated info
        string sSubjectsDir = this->mQdecProject->GetSubjectsDir();
        this->mEntrySubjectsDir->SetValue( sSubjectsDir.c_str() );

        QdecGlmFitResults* results = mQdecProject->GetGlmFitResults();
        assert( results );

        // We need to update our tabs.
        this->UpdateDesignPage();
        this->UpdateSubjectsPage();

        // Load in the data.
        this->LoadAnalyzedData( results );
	
        this->SetStatusText( "Project file loaded." );
        this->AddRecentFile( ifnProject, this, "LoadProjectFile" );

      }
    } catch ( exception& e ) {
      this->GetApplication()->ErrorMessage( e.what() );
      this->SetStatusText( "Error loading project file" );
    }
  }
}

void
vtkKWQdecWindow::SaveProjectFile ( const char* ifnProject ) {

  if( mView ) {
    try {
      if( ifnProject ) { 

        this->SetStatusText( "Saving project file..." );

        if( this->mQdecProject->SaveProjectFile( ifnProject ) )
          throw runtime_error( "Error saving the project file." );
	
        this->SetStatusText( "Project file saved." );
        this->AddRecentFile( ifnProject, this, "SaveProjectFile" );
	
      }
    } catch ( exception& e ) {
      this->GetApplication()->ErrorMessage( e.what() );
      this->SetStatusText( "Error saving project file" );
    }
  }
}

void
vtkKWQdecWindow::LoadAnalyzedData (  QdecGlmFitResults* iGlmResults ) {

  if( mView ) {
    try {

      if( NULL == iGlmResults ) {
        throw runtime_error( "No results have been generated!" );
      }

      // We'll run a progress bar that just tracks steps.
      const float cSteps = 8.0;
      float stepIncrement = 100.0 / cSteps;
      int nStep = 1;

      // Get SUBJECTS_DIR and average subj
      string sSubjectsDir = iGlmResults->GetGlmDesign()->GetSubjectsDir();
      string sAverageSubj = iGlmResults->GetGlmDesign()->GetAverageSubject();
      string sHemisphere  = iGlmResults->GetGlmDesign()->GetHemi();

      // Load the white surface at
      // SUBJECTS_DIR/fsaverage/surf/$hemi.white if it's available.
      string fnWhite = sSubjectsDir + "/" +  sAverageSubj +
        "/surf/" + sHemisphere + ".white";
      if( QdecUtilities::IsFileReadable( fnWhite ) ) {
        this->SetStatusText( "Loading white surface..." );
        this->GetProgressGauge()->SetValue( (float)nStep++ * stepIncrement );
        this->LoadSurface( fnWhite.c_str(), "white" );
      }

      // Load the pial surface at
      // SUBJECTS_DIR/fsaverage/surf/$hemi.pial if it's available.
      string fnPial = sSubjectsDir + "/" +  sAverageSubj +
        "/surf/" + sHemisphere + ".pial";
      if( QdecUtilities::IsFileReadable( fnPial ) ) {
        this->SetStatusText( "Loading pial surface..." );
        this->GetProgressGauge()->SetValue( (float)nStep++ * stepIncrement );
        this->LoadSurface( fnPial.c_str(), "pial" );
      }

      // Make sure we have an fsaverage subject. The surface is
      // SUBJECTS_DIR/fsaverage/surf/$hemi.inflated
      this->SetStatusText( "Loading inflated surface..." );
      this->GetProgressGauge()->SetValue( (float)nStep++ * stepIncrement );
      string fnInflated = sSubjectsDir + "/" + sAverageSubj +
        "/surf/" + sHemisphere + ".inflated";
      QdecUtilities::AssertFileIsReadable( fnInflated );
      this->LoadSurface( fnInflated.c_str(), "inflated" );

      // We also want to load the curvature in
      // SUBJECTS_DIR/fsaverage/surf/$hemi.curv.
      this->SetStatusText( "Loading curvature..." );
      this->GetProgressGauge()->SetValue( (float)nStep++ * stepIncrement );
      string fnCurvature =
        sSubjectsDir + "/" +  sAverageSubj +
        "/surf/" + sHemisphere + ".curv";
      QdecUtilities::AssertFileIsReadable( fnCurvature );
      this->LoadSurfaceCurvatureScalars( fnCurvature.c_str() );

      // Try to load the annot file if it's available.
      this->SetStatusText( "Trying to load annotation..." );
      this->GetProgressGauge()->SetValue( (float)nStep++ * stepIncrement );
      string fnAnnot =
        sSubjectsDir + "/" +  sAverageSubj +
        "/label/" + sHemisphere + ".aparc.annot";
      if( QdecUtilities::IsFileReadable( fnAnnot ) ) {
        this->LoadAnnotation( fnAnnot.c_str() );
      }

      // There are three scalars to load. The first is the "Contrasts"
      // and it is fnWorkingDir/<contrast>/sig.mgh. We also need to get the
      // contrast question and make it the label.
      this->SetStatusText( "Loading contrasts..." );
      this->GetProgressGauge()->SetValue( (float)nStep++ * stepIncrement );
      vector<string> fnContrastSigFiles = iGlmResults->GetContrastSigFiles();
      vector<string> sContrastQuestions = iGlmResults->GetContrastQuestions();
      assert( fnContrastSigFiles.size() == sContrastQuestions.size() );
      for( unsigned int i=0;
           i < iGlmResults->GetContrastSigFiles().size(); i++) {
        QdecUtilities::AssertFileIsReadable( fnContrastSigFiles[i] );
        int nScalar = this->LoadSurfaceScalars( fnContrastSigFiles[i].c_str(),
                                                sContrastQuestions[i].c_str());
        if( nScalar < 0 )
          throw runtime_error( string("Couldn't load scalars from ") + 
                               fnContrastSigFiles[i] );
      }

      // The second is the "Residual Error StdDev" in fnWorkingDir/rstd.mgh.
      this->SetStatusText( "Loading error stddev..." );
      this->GetProgressGauge()->SetValue( (float)nStep++ * stepIncrement );
      string fnRSTD = iGlmResults->GetResidualErrorStdDevFile();
      QdecUtilities::AssertFileIsReadable( fnRSTD );
      int nScalar = this->LoadSurfaceScalars( fnRSTD.c_str(),
                                              "Residual Error StdDev" );
      if( nScalar < 0 )
        throw runtime_error( string("Couldn't load scalars from ") + fnRSTD );

      // The third is the "Regression Coefficients" in fnWorkingDir/beta.mgh.
      this->SetStatusText( "Loading regression coefficients..." );
      this->GetProgressGauge()->SetValue( (float)nStep++ * stepIncrement );
      string fnBeta = iGlmResults->GetRegressionCoefficientsFile();
      QdecUtilities::AssertFileIsReadable( fnBeta );
      nScalar = this->LoadSurfaceScalars( fnBeta.c_str(),
                                          "Regression Coefficients" );
      if( nScalar < 0 )
        throw runtime_error( string("Couldn't load scalars from ") + fnBeta );

      // The GDF to load is fnWorkingDir/y.fsgd.
      this->SetStatusText( "Loading GDF..." );
      this->GetProgressGauge()->SetValue( (float)nStep++ * stepIncrement );
      string fnGDF = iGlmResults->GetFsgdFile();
      QdecUtilities::AssertFileIsReadable( fnGDF );
      this->LoadGDFFile( fnGDF.c_str() );

      this->SetStatusText( "Analyzed data loaded." );
      this->GetProgressGauge()->SetValue( 100 );
      this->AddRecentFile(iGlmResults->GetGlmDesign()->GetWorkingDir().c_str(),
                          this, "LoadAnalyzedData" );

    } catch ( exception& e ) {
      this->GetApplication()->ErrorMessage( e.what() );
      this->SetStatusText( "Error loading analyzed data" );
    }
  }

  // Update our menu and buttons.
  this->UpdateCommandStatus();

  // Update the display page.
  this->UpdateDisplayPage();

  printf("\n============================================================\n"
         "Completed loading of analyzed data.\n");
}

void
vtkKWQdecWindow::LoadSurface ( const char* ifnSurface, const char* isLabel ) {

  if( mView ) {
    try {

      // Try to load the surface object.
      vtkSmartPointer<vtkFSSurfaceSource> surface =
        vtkSmartPointer<vtkFSSurfaceSource>::New();
      surface->MRISRead( ifnSurface );
      surface->Update();

      // If we already have a number of vertices, check it. Otherwise,
      // save this as the number of vertices.
      if( -1 != mcVertices ) {
        if( mcVertices != surface->GetNumberOfVertices() ) {
          throw runtime_error(
            "Multiple surfaces must have the same number of vertices" );
        }
      } else {
        mcVertices = surface->GetNumberOfVertices();
      }

      // Get a good label.
      string sLabel;
      if( NULL == isLabel ) {
        sLabel = "user";
      } else {
        sLabel = isLabel;
      }

      // If we already have a surface with this label, delete it.
      if( maSurfaceSource.find( sLabel ) != maSurfaceSource.end() )
        maSurfaceSource.erase( sLabel );

      // Save a pointer to it.
      maSurfaceSource[sLabel] = surface;

      // Select this one.
      this->SetCurrentSurface( sLabel.c_str() );

      // First time we do this, reset the view to initialize it.
      if( !mbViewInitialized ) {
        mView->ResetView();
        mbViewInitialized = true;
      }

      this->SetStatusText( "Surface loaded." );
      this->AddRecentFile( ifnSurface, this, "LoadSurface" );

    } catch ( exception& e ) {
      this->GetApplication()->ErrorMessage( e.what() );
    }
  }

  // Update our menu and buttons.
  this->UpdateCommandStatus();
}

void
vtkKWQdecWindow::LoadGDFFile ( const char* ifnGDFFile ) {

  if( mView ) {
    try {

      // Try to load the GDF.
      const char* sResult = this->Script( "FsgdfPlot_Read %s", ifnGDFFile );

      int id = strtol( sResult, (char**)NULL, 10 );
      if( ERANGE == errno ) {
        cerr << "Couldn't load GDF: " << sResult;
        return;
      }

      mGDFID = id;
      mbGDFLoaded = true;

      // Tell the view.
      mView->SetGDFID( mGDFID );

      this->SetStatusText( "GDF file loaded." );
      this->AddRecentFile( ifnGDFFile, this, "LoadGDFFile" );

    } catch ( exception& e ) {
      this->GetApplication()->ErrorMessage( e.what() );
    }
  }

  // Update our menu and buttons.
  this->UpdateCommandStatus();
}

// Load a surface scalar file and return the entry index associated
// with it. Optionally takes a label name to associate with this label file.
// Optionally takes a frame number if reading a multi-frame file. 
// This does not change the active scalars file.
int  
vtkKWQdecWindow::LoadSurfaceScalars ( const char* ifnScalars, 
                                      const char* isLabel,
                                      int inFrame ) {
  int rnEntry = -1;

  try {

    if( maSurfaceSource.size() < 1 )
      throw runtime_error( "Must load a surface before loading scalars." );

    // Try to load the scalars.
    vtkSmartPointer<vtkFloatArray> scalars;
    scalars.TakeReference( 
      this->NewScalarsFromSurfaceScalarsFile( ifnScalars, inFrame ) );

    // Make an entry into our table of scalars;
    int nEntry = 0;
    while ( maSurfaceScalars.find(nEntry) != maSurfaceScalars.end() )
      nEntry++;
    maSurfaceScalars[nEntry].mValues = scalars;
    maSurfaceScalars[nEntry].mnEntry = nEntry;
    maSurfaceScalars[nEntry].mfnSource = ifnScalars;
    if( isLabel ) {
      maSurfaceScalars[nEntry].msLabel = isLabel;
    } else {
      // if a label was not supplied, just use the filename
      string fnScalars = ifnScalars;
      string::size_type lastSlash = fnScalars.rfind( "/" );
      maSurfaceScalars[nEntry].msLabel =
        fnScalars.substr( lastSlash+1, string::npos ).c_str();
    }

    // Return the index of the scalars we loaded.
    rnEntry = nEntry;

    this->SetStatusText( "Scalars loaded." );
    this->AddRecentFile( ifnScalars, this, "LoadSurfaceScalars" );

  } catch ( exception& e ) {
    this->GetApplication()->ErrorMessage( e.what() );
  }

  // Update our menu and buttons.
  this->UpdateCommandStatus();
  this->UpdateDisplayPage();

  return rnEntry;
}

void
vtkKWQdecWindow::LoadSurfaceCurvatureScalars ( const char* ifnScalars ) {

  try {

    if( maSurfaceSource.size() < 1 )
      throw runtime_error( "Must load a surface before loading scalars." );

    // Try to load the scalars.
    vtkSmartPointer<vtkFloatArray> scalars;
    scalars.TakeReference(
      this->NewScalarsFromSurfaceScalarsFile( ifnScalars ) );

    // Save our curvature.
    mCurvatureScalars = scalars;

    this->ComposeSurfaceScalarsAndShow();

    this->SetStatusText( "Curvature loaded." );
    this->AddRecentFile( ifnScalars, this, "LoadSurfaceCurvatureScalars" );

  } catch ( exception& e ) {
    this->GetApplication()->ErrorMessage( e.what() );
  }

  // Update our menu and buttons.
  this->UpdateCommandStatus();
  this->UpdateDisplayPage();
}

void
vtkKWQdecWindow::LoadAnnotation ( const char* ifnAnnotation ) {

  assert( ifnAnnotation );
  if( maSurfaceSource.size() < 1 )
    throw runtime_error( "Must load a surface before loading annotation." );

  int* indices = NULL;
  COLOR_TABLE* table = NULL;

  try {

    int cValues = mcVertices;

    // Try to read the annotation.
    int eRead =
      MRISreadAnnotationIntoArray( ifnAnnotation, cValues, &indices );
    if( 0 != eRead )
      throw runtime_error ("Could not read annotation file");

    // Try to read a color table too.
    eRead = MRISreadCTABFromAnnotationIfPresent( ifnAnnotation, &table );
    if( 0 != eRead )
      throw runtime_error ("Could not read color table file");

    // If we got both, save the new ones.
    if( maAnnotationIndicies )
      free( maAnnotationIndicies );
    maAnnotationIndicies = indices;
    indices = NULL;

    if( mAnnotationTable )
      CTABfree( &mAnnotationTable );
    mAnnotationTable = table;
    table = NULL;

    // Set us as a lookup object in the view.
    mView->SetVertexAnnotationLookup( this );

    // Also load this as an overlay with an initial opacity of 0.
    this->LoadSurfaceOverlayScalars( ifnAnnotation );
    this->SetSurfaceOverlayOpacity( 0.0 );

  } catch (exception& e) {
    this->GetApplication()->ErrorMessage( e.what() );
  }

  if( table )
    free( table );
  if( table )
    CTABfree( &table );

  // Update our menu and buttons.
  this->UpdateCommandStatus();
  this->UpdateDisplayPage();

}

void
vtkKWQdecWindow::LoadSurfaceOverlayScalars ( const char* ifnScalars,
                                             const char* ifnColors ) {

  assert( ifnScalars );
  if( maSurfaceSource.size() < 1 )
    throw runtime_error( "Must load a surface before loading annotation." );

  int* lookupTable = NULL;
  COLOR_TABLE* table = NULL;

  try {

    int cValues = mcVertices;

    // Try to read the annotation.
    int eRead =
      MRISreadAnnotationIntoArray( ifnScalars, cValues, &lookupTable );
    if( 0 != eRead )
      throw runtime_error ("Could not read annotation file");

    // Try to read a color table too.
    eRead = MRISreadCTABFromAnnotationIfPresent( ifnScalars, &table );
    if( 0 != eRead ) {

      // If we didn't get a table, they have to have passed in a file
      // name for one.
      if( ifnColors == NULL ) {
        throw runtime_error( "No color table present in overlay file; "
                             "need to specify an external file." );
      }

      //  Try to load the table.
      char* fnColors = strdup( ifnColors );
      table = CTABreadASCII( fnColors );
      free( fnColors );
      if( NULL == table )
        throw runtime_error( string("Couldn't read color table ") + ifnColors);
    }

    // Init a float array.
    vtkSmartPointer<vtkFloatArray> scalars = 
      vtkSmartPointer<vtkFloatArray>::New();

    // Allocate our scalars.
    scalars->Allocate( cValues );
    scalars->SetNumberOfComponents( 1 );

    // Copy our array into the scalars.
    for( int nValue = 0; nValue < cValues; nValue ++ ) {
      int nEntry = 0;
      CTABfindAnnotation( table, lookupTable[nValue], &nEntry );
      scalars->InsertNextValue( static_cast<float>(nEntry) );
    }

    // Convert the CTAB to a vtkLookupTable.
    vtkSmartPointer<vtkFreesurferLookupTable> colors =
      vtkSmartPointer<vtkFreesurferLookupTable>::New();
    colors->BuildFromCTAB( table );

    // Save pointers.
    mOverlayScalars = scalars;
    mOverlayColors = colors;

    // Set it in the view.
    mView->
      SetSurfaceOverlayScalarsAndColors( mOverlayScalars, mOverlayColors );

    // If the suffix is .annot, use "annotation" as a label, otherwise
    // use its name as the label.
    string fnScalars = ifnScalars;
    if( fnScalars.find( ".annot" ) != string::npos ) {
      msOverlayDescription = "Annotation";
    } else {
      string::size_type lastSlash = fnScalars.rfind( "/" );
      msOverlayDescription =
        fnScalars.substr( lastSlash+1, string::npos ).c_str();
    }

  } catch ( exception& e ) {
    this->GetApplication()->ErrorMessage( e.what() );
  }

  if( lookupTable )
    free( lookupTable );
  if( table )
    CTABfree( &table );

  // Update our menu and buttons.
  this->UpdateCommandStatus();
  this->UpdateDisplayPage();

}

void
vtkKWQdecWindow::LoadLabel ( const char* ifnLabel ) {

  assert( ifnLabel );
  if( maSurfaceSource.size() < 1 )
    throw runtime_error( "Must load a surface before loading labels." );

  try {

    // Get the first surface to use to read the label.
    MRIS* mris = maSurfaceSource.begin()->second->GetMRIS();

    vtkSmartPointer<vtkFSSurfaceLabelSource> labelSource =
      vtkSmartPointer<vtkFSSurfaceLabelSource>::New();
    labelSource->SetLabelFileName( ifnLabel );
    labelSource->SetMris( mris );
    labelSource->Update();

    // Assign the new one.
    mROISource = labelSource;

    // Send the poly data to the view.
    mView->SetROI ( mROISource->GetOutput() );

  } catch (exception& e) {
    this->GetApplication()->ErrorMessage( e.what() );
  }

  this->UpdateCommandStatus();
}

void
vtkKWQdecWindow::SaveLabel ( const char* ifnLabel ) {

  assert( ifnLabel );
  assert( mROISource );

  try {

    // Try to write the label.
    mROISource->SetLabelFileName( ifnLabel );
    mROISource->WriteLabelFile();

    this->SetStatusText( "Label saved." );

  } catch (exception& e) {
    this->GetApplication()->ErrorMessage( e.what() );
  }

  this->UpdateCommandStatus();
}

void
vtkKWQdecWindow::MapLabelToSubjects ( const char* ifnLabel ) {

  assert( ifnLabel );
  assert( mROISource );
  assert( mQdecProject );

  try {

    // Try to write the label to the average subject.
    mROISource->SetLabelFileName( ifnLabel );
    mROISource->WriteLabelFile();

    // Tell the project to run this for all our subjects.
    mQdecProject->GenerateMappedLabelForAllSubjects( ifnLabel, this );

    this->SetStatusText( "Labels mapped and saved." );

  } catch (exception& e) {
    this->GetApplication()->ErrorMessage( e.what() );
  }

  this->UpdateCommandStatus();
}

void
vtkKWQdecWindow::SaveTIFFImage ( const char* ifnTIFF, 
                                 int iMagnificationLevel ) {

  assert( iMagnificationLevel >= 1 );

  if( mView ) {
    try {

      // Turn the legend and annotation off if we're doing any
      // magnification.
      if( iMagnificationLevel != 1 ) {
        mView->SetShowLegend( 0 );
        mView->SetShowAnnotation( 0 );
      }
      
      // Create our pipeline objects. vtkWindowToImageFilter just
      // takes a normal screenshot, while vtkRenderLargeImage actually
      // does an offscreen render.
      vtkSmartPointer<vtkWindowToImageFilter> shot;
      vtkSmartPointer<vtkRenderLargeImage> imageMag;

      // If we have a magnification level, make the imageMag,
      // otherwise make the shot.
      if( iMagnificationLevel != 1 ) {
        imageMag = vtkSmartPointer<vtkRenderLargeImage>::New();
        imageMag->SetInput( (vtkRenderer*)mView->GetRenderer() );
        imageMag->SetMagnification( iMagnificationLevel );
      } else {
        shot = vtkSmartPointer<vtkWindowToImageFilter>::New();
        shot->SetInput( (vtkWindow*)mView->GetRenderWindow() );
      }
      
      // Make the writer object and connect whichever filter object we
      // just made to it.
      vtkSmartPointer<vtkTIFFWriter> writer = 
        vtkSmartPointer<vtkTIFFWriter>::New();
      if( imageMag.GetPointer() ) 
        writer->SetInput( imageMag->GetOutput() );
      else 
        writer->SetInput( shot->GetOutput() );

      // Write the file.
      writer->SetFileName( ifnTIFF );
      writer->Write();

      // Reshow the legend and annotation.
      if( iMagnificationLevel != 1 ) {
        mView->SetShowLegend( 1 );
        mView->SetShowAnnotation( 1 );
      }

      this->SetStatusText( "TIFF saved." );

    } catch ( exception& e ) {
      this->GetApplication()->ErrorMessage( e.what() );
    }
  }
}


void
vtkKWQdecWindow::RestoreView () {

  if( mView.GetPointer() )
    mView->RestoreView();

  // Update our menu and buttons.
  this->UpdateCommandStatus();
}


void
vtkKWQdecWindow::SetCurrentSurfaceScalars ( int inEntry ) {

  if( mView.GetPointer() ) {
    try {

      // Make sure we have these scalars.
      if( maSurfaceScalars.find(inEntry) == maSurfaceScalars.end() )
        return;

      // Save the new current index.
      mnCurrentSurfaceScalars = inEntry;

#if 0
      // Quick histogram.
      if( mnCurrentSurfaceScalars >= 0 ) {
        vtkFloatArray* a =
          maSurfaceScalars[mnCurrentSurfaceScalars].mValues;
        int const cBins = 10;
        int* aBins = (int*) calloc( cBins, sizeof(int) );
        double binRange =
          (a->GetRange()[1] - a->GetRange()[0] + 0.00001) / (double)(cBins);
        for( int n = 0; n < a->GetNumberOfTuples(); n++ )
          aBins[(int)floor((a->GetTuple1(n)-a->GetRange()[0])/binRange)]++;
        cerr << "Histogram " << a->GetNumberOfTuples()
             << " values over [" << a->GetRange()[0]
             << "," << a->GetRange()[1] << "]" << endl;
        for( double n = 0; n < cBins; n += 1.0 )
          cerr << "[" << a->GetRange()[0] + (n*binRange) << ","
               << a->GetRange()[0] + ((n+1.0)*binRange)
               << "): " << aBins[(int)n] << endl;
        free( aBins );
      };
#endif

      // Call this to display these scalars.
      this->ComposeSurfaceScalarsAndShow();

      this->UpdateSurfaceScalarsColorsEditor();

      // Update our menu.
      this->UpdateCommandStatus();

    } catch ( exception& e ) {
      this->GetApplication()->ErrorMessage( e.what() );
    }
  }

}

void
vtkKWQdecWindow::ClearSurfaceScalars () {

  // Delete all the scalar arrays.
  maSurfaceScalars.clear();

  // No selected scalar.
  mnCurrentSurfaceScalars = -1;

  // Redraw.
  this->ComposeSurfaceScalarsAndShow();

  // Update the colors editor.
  this->UpdateSurfaceScalarsColorsEditor();

  // Update our menu and page.
  this->UpdateDisplayPage();
  this->UpdateCommandStatus();
}

void
vtkKWQdecWindow::ClearCurvature () {

  if( NULL == mCurvatureScalars )
    return;

  // Delete the curvature scalars.
  mCurvatureScalars = NULL;

  // Redraw.
  this->ComposeSurfaceScalarsAndShow();

  // Update our menu.
  this->UpdateCommandStatus();
}

void
vtkKWQdecWindow::ZoomBy ( float iFactor ) {

  if( mView.GetPointer() )
    mView->ZoomBy( iFactor );

  // Update our menu and buttons.
  this->UpdateCommandStatus();
}

void
vtkKWQdecWindow::ZoomIn () {

  this->ZoomBy( 2.0 );
}

void
vtkKWQdecWindow::ZoomOut () {

  this->ZoomBy( 0.5 );
}

void
vtkKWQdecWindow::ShowCursor ( int ibShow ) {

  assert( mView.GetPointer() );

  mView->SetShowCursor( ibShow );

  this->UpdateCommandStatus();
}

void
vtkKWQdecWindow::SetShowCursorFromMenu () {

  assert( mMenuShowCursor );

  this->ShowCursor( mMenuShowCursor->GetSelectedState() );
}

void
vtkKWQdecWindow::SetCurrentSurface ( const char* isLabel ) {

  assert( isLabel );

  // If we are setting to NULL, that means clear the current surface.
  if( NULL == isLabel ) {

    // No current surface.
    msCurrentSurfaceSource = "";

    // Tell the view to clear its surface.
    mView->SetSurface( NULL );

  } else {
    
    // Make sure it exists.
    if( maSurfaceSource.find( string(isLabel) ) == maSurfaceSource.end() )
      throw runtime_error( "SetCurrentSurface: isLabel surface was not found");
    
    // Save this is the current label.
    msCurrentSurfaceSource = isLabel;
    
    // Set it in the view.
    mView->SetSurface( maSurfaceSource[msCurrentSurfaceSource] );
    
    // Update the source for the ROI if we have one.
    if( NULL != mROISource ) {
      mROISource->
        SetMris( maSurfaceSource[msCurrentSurfaceSource]->GetMRIS() );
      mView->SetROI( mROISource->GetOutput() );
    }
  }
}

void
vtkKWQdecWindow::SetSurfaceOverlayOpacity ( double iOpacity ) {

  assert( mView.GetPointer() );

  if( iOpacity < 0.0 || iOpacity > 1.0 )
    throw runtime_error( "Invalid opacity" );

  // Set our scale if we have it.
  if( mScaleOverlay ) {
    mScaleOverlay->SetValue( iOpacity );
  }
  
  // Set it in the view.
  mView->SetSurfaceOverlayOpacity( iOpacity );
}

void
vtkKWQdecWindow::SetCurrentSurfaceScalarsFromTableSelection () {

  assert( mTableSurfaceScalars.GetPointer() );

  int nEntry = mTableSurfaceScalars->GetIndexOfFirstSelectedRow();

  this->SetCurrentSurfaceScalars( nEntry );
}

void
vtkKWQdecWindow::DiscreteFactorsListBoxCallback () {

  assert( mListDiscreteFactors.GetPointer() );

  this->ManageFactorListBoxSelections( mListDiscreteFactors,
                                       maDiscreteFactorSelection );
}

void
vtkKWQdecWindow::ContinuousFactorsListBoxCallback () {

  assert( mListContinuousFactors.GetPointer() );

  this->ManageFactorListBoxSelections( mListContinuousFactors,
                                       maContinuousFactorSelection );
}

void
vtkKWQdecWindow::PlotContinuousFactorsListBoxCallback () {
  assert( mListPlotContinuousFactors.GetPointer() );

  mPlotContinuousFactorSelection =
    mListPlotContinuousFactors->GetSelectionIndex();

  this->UpdateContinuousFactorPlot();
}

void
vtkKWQdecWindow::ManageFactorListBoxSelections ( vtkKWListBox* iListBox,
                                                 int iaSelections[2] ) {

  // Make sure that we only have two items selected, max. If more than
  // two, unselect the first one.
  int cSelected = 0;
  for( int nItem = 0; nItem < iListBox->GetNumberOfItems(); nItem++ )
    if( iListBox->GetSelectState( nItem ) )
      cSelected++;

  if( cSelected > 2 ) {

    // Unselect the first one.
    int nUnselect = iaSelections[0];
    iListBox->SetSelectState( nUnselect, false );

    // Bump the second selection into the first slot.
    iaSelections[0] = iaSelections[1];
    iaSelections[1] = -1;
  }

  // Find the selected items. If an item is selected but not in our
  // array of selections, make the first or second selection,
  // depending on which ones we have. If it's not selected, make sure
  // it is no longer in our array. If we had to clear the first
  // selection, bump the second one into the first slot.
  for( int nItem = 0; nItem < iListBox->GetNumberOfItems(); nItem++ ) {

    if( iListBox->GetSelectState( nItem ) ) {

      if( -1 == iaSelections[0] )
        iaSelections[0] = nItem;
      // Make sure it's not already marked in the first selection.
      else if( nItem != iaSelections[0] )
        iaSelections[1] = nItem;

    } else {

      if( iaSelections[0] == nItem ) {
        iaSelections[0] = iaSelections[1];
        iaSelections[1] = -1;
      } else if( iaSelections[1] == nItem ) {
        iaSelections[1] = -1;
      }

    }
  }
}

void
vtkKWQdecWindow::UpdateSubjectsPage () {

  assert( mEntryDataTable.GetPointer() );
  assert( mEntryNumberOfSubjects.GetPointer() );
  assert( mQdecProject );
  assert( mQdecProject->GetDataTable() );
  assert( mListPlotContinuousFactors.GetPointer() );

  QdecDataTable* dTable = this->mQdecProject->GetDataTable();

  mEntryDataTable->SetValue( dTable->GetFileName().c_str() );
  mEntryNumberOfSubjects->SetValueAsInt( dTable->GetSubjectIDs().size() );
  this->mQdecProject->SetAverageSubject( mEntryAverageSubject->GetValue() );
  this->mQdecProject->SetSubjectsDir( mEntrySubjectsDir->GetValue() );

  mListPlotContinuousFactors->DeleteAll();

  vector< string > factors =
    this->mQdecProject->GetDataTable()->GetContinuousFactors();
  for(unsigned int i=0; i < factors.size(); i++) {
    mListPlotContinuousFactors->Append( factors[i].c_str() );
  }
}

void
vtkKWQdecWindow::UpdateDesignPage () {

  assert( mListDiscreteFactors.GetPointer() );
  assert( mListContinuousFactors.GetPointer() );
  assert( mQdecProject ) ;
  assert( mQdecProject->GetDataTable() );
  assert( mQdecProject->GetGlmDesign() );

  // Fill our our entries with the values from the design.
  QdecGlmDesign* design =  mQdecProject->GetGlmDesign();
  mEntryDesignName->SetValue( design->GetName().c_str() );
  mMenuMeasure->SetValue( design->GetMeasure().c_str() );
  mMenuHemisphere->SetValue( design->GetHemi().c_str() );
  stringstream ssSmoothness;
  ssSmoothness << design->GetSmoothness();
  mMenuSmoothness->SetValue( ssSmoothness.str().c_str() );

  // Clear the factor lists.
  mListDiscreteFactors->DeleteAll();
  mListContinuousFactors->DeleteAll();

  // Get the current factors.
  vector<QdecFactor*> const& lDiscreteFactors =
    design->GetDiscreteFactors();
  vector<QdecFactor*> const& lContinuousFactors = 
    design->GetContinuousFactors();
  int nDiscreteSelection = 0;
  int nContinuousSelection = 0;

  vector< string > factors =
    this->mQdecProject->GetDataTable()->GetDiscreteFactors();
  for(unsigned int i=0; i < factors.size(); i++) {
    mListDiscreteFactors->Append( factors[i].c_str() );

    // If this factor is one of our chosen discrete factors from the
    // design, enter its index in our selections, and select this
    // item.
    if( (lDiscreteFactors.size() > 0 && 
         lDiscreteFactors[0]->GetFactorName() == factors[i]) ||
        (lDiscreteFactors.size() > 1 &&
         lDiscreteFactors[1]->GetFactorName() == factors[i]) ) {
      maDiscreteFactorSelection[nDiscreteSelection++] = i;
      mListDiscreteFactors->SetSelectState( i, 1 );
    }
  }

  factors = this->mQdecProject->GetDataTable()->GetContinuousFactors();
  for(unsigned int i=0; i < factors.size(); i++) {
    mListContinuousFactors->Append( factors[i].c_str() );

    if( (lContinuousFactors.size() > 0 && 
         lContinuousFactors[0]->GetFactorName() == factors[i] ) ||
        (lContinuousFactors.size() > 1 && 
         lContinuousFactors[1]->GetFactorName() == factors[i] ) ) {
      maContinuousFactorSelection[nContinuousSelection++] = i;
      mListContinuousFactors->SetSelectState( i, 1 );
    }
  }
}

void
vtkKWQdecWindow::UpdateDisplayPage () {

  assert( mFrameSurface.GetPointer() );
  assert( mFrameCurvature.GetPointer() );
  assert( mFrameOverlay.GetPointer() );
  assert( mFrameSurfaceScalars.GetPointer() );

  // We want to check if we have to populate any of the empty frames
  // we asserted above. For each one, if we have data for it and it's
  // uninited, init it with the specific widgets it needs. If we no
  // longer have data for it and it's initied, unpack it.

  // If we have more than one surface loaded, we need to populate the
  // surface frame with a label frame and some radio buttons for
  // choosing the surface.
  if( maSurfaceSource.size() > 1 && !mbFrameSurfaceInited ) {

    vtkSmartPointer<vtkKWFrameWithLabel> surfaceFrame =
      vtkSmartPointer<vtkKWFrameWithLabel>::New();
    surfaceFrame->SetParent( mFrameSurface );
    surfaceFrame->Create();
    surfaceFrame->SetLabelText( "Surfaces" );
    this->Script( "pack %s -side top -fill both",
                  surfaceFrame->GetWidgetName() );

    mRadBtnSetSurface = vtkSmartPointer<vtkKWRadioButtonSet>::New();
    mRadBtnSetSurface->SetParent( surfaceFrame->GetFrame() );
    mRadBtnSetSurface->Create();
    mRadBtnSetSurface->PackHorizontallyOn();
    this->Script( "pack %s -side top -fill both",
                  mRadBtnSetSurface->GetWidgetName() );

    int nButton = 0;
    map<string,vtkSmartPointer<vtkFSSurfaceSource> >::iterator tSurface;
    for( tSurface = maSurfaceSource.begin();
         tSurface != maSurfaceSource.end(); ++tSurface ) {

      string sSurface = tSurface->first;

      // Make a button for this surface.
      vtkSmartPointer<vtkKWRadioButton> radBtn;
      radBtn.TakeReference( mRadBtnSetSurface->AddWidget( nButton++ ) );
      radBtn->SetText( sSurface.c_str() );

      // Make the command for this button.
      string sCmd = string("SetCurrentSurface ") + sSurface;
      radBtn->SetCommand( this, sCmd.c_str() );

      // If this is our currently displayed surface, select the button.
      if( sSurface == msCurrentSurfaceSource ) radBtn->SelectedStateOn();
    }

    mbFrameSurfaceInited = true;
  }

  // Make sure our currently displayed surface's button is selected
  if( maSurfaceSource.size() > 1 && mbFrameSurfaceInited ) {
    int nButton = 0;
    vtkSmartPointer<vtkKWRadioButton> radBtn = NULL;
    do {
      radBtn = mRadBtnSetSurface->GetWidget ( nButton++ );
      if ( radBtn ) {
        string sSurface = radBtn->GetText();
        // If this is our currently displayed surface, select the button.
        if( sSurface == msCurrentSurfaceSource ) radBtn->SelectedStateOn();
      }
    }
    while ( NULL != radBtn);
  }

  // If we don't have more than one surface and the frame is inited,
  // unpack it and shrink it to hide it.
  if( maSurfaceSource.size() <= 1 && mbFrameSurfaceInited ) {
    mFrameSurface->UnpackChildren();
    mFrameSurface->SetHeight( 1 );
    mRadBtnSetSurface = NULL;
    mbFrameSurfaceInited = false;
  }


  // If we have curvature, create a labeled frame for it and put a
  // checkbox in it controlling the visibility.
  if( mCurvatureScalars.GetPointer() && !mbFrameCurvatureInited ) {

    vtkSmartPointer<vtkKWFrameWithLabel> curvatureFrame =
      vtkSmartPointer<vtkKWFrameWithLabel>::New();
    curvatureFrame->SetParent( mFrameCurvature );
    curvatureFrame->Create();
    curvatureFrame->SetLabelText( "Curvature" );
    this->Script( "pack %s -side top -fill both",
                  curvatureFrame->GetWidgetName() );

    vtkSmartPointer<vtkKWCheckButtonWithLabel> labeledCheck =
      vtkSmartPointer<vtkKWCheckButtonWithLabel>::New();
    labeledCheck->SetParent( curvatureFrame->GetFrame() );
    labeledCheck->Create();
    labeledCheck->SetLabelText( "Show Curvature" );
    labeledCheck->SetLabelPositionToRight();
    this->Script( "pack %s -side top",
                  labeledCheck->GetWidgetName() );

    mCheckShowCurvature = labeledCheck->GetWidget();
    mCheckShowCurvature->SetSelectedState( mbShowCurvature );
    mCheckShowCurvature->SetCommand( this, "SetShowCurvature" );

    labeledCheck = vtkSmartPointer<vtkKWCheckButtonWithLabel>::New();
    labeledCheck->SetParent( curvatureFrame->GetFrame() );
    labeledCheck->Create();
    labeledCheck->
      SetLabelText( "Draw Curvature in green/red\nif no scalars loaded" );
    labeledCheck->SetLabelPositionToRight();
    this->Script( "pack %s -side top",
                  labeledCheck->GetWidgetName() );

    mCheckDrawCurvatureGreenRed = labeledCheck->GetWidget();
    mCheckDrawCurvatureGreenRed->
      SetSelectedState( mbDrawCurvatureGreenRedIfNoScalars );
    mCheckDrawCurvatureGreenRed->
      SetCommand( this, "SetDrawCurvatureGreenRed" );

    mbFrameCurvatureInited = true;
  }
  if( NULL == mCurvatureScalars.GetPointer() && mbFrameCurvatureInited ) {
    mFrameCurvature->UnpackChildren();
    mFrameCurvature->SetHeight( 1 );
    mCheckShowCurvature = NULL;
    mCheckDrawCurvatureGreenRed = NULL;
    mbFrameCurvatureInited = false;
  }


  // If we have an overlay, create a scale for controlling the
  // opacity. We'll put it in a labeled frame and set the label to the
  // description we got from how we loaded the overlay; from the
  // annotation or from a file directly.
  if( mOverlayScalars.GetPointer() && 
      mOverlayColors.GetPointer() &&
      !mbFrameOverlayInited ) {

    vtkSmartPointer<vtkKWFrameWithLabel> overlayFrame =
      vtkSmartPointer<vtkKWFrameWithLabel>::New();
    overlayFrame->SetParent( mFrameOverlay );
    overlayFrame->Create();
    overlayFrame->SetLabelText( msOverlayDescription.c_str() );
    this->Script( "pack %s -side top -fill both",
                  overlayFrame->GetWidgetName() );

    vtkSmartPointer<vtkKWScaleWithLabel> labeledScale =
      vtkSmartPointer<vtkKWScaleWithLabel>::New();
    labeledScale->SetParent( overlayFrame->GetFrame() );
    labeledScale->Create();
    labeledScale->SetLabelText( "Opacity: " );
    this->Script( "pack %s -side top -fill both",
                  labeledScale->GetWidgetName() );

    mScaleOverlay = labeledScale->GetWidget();
    mScaleOverlay->SetRange( 0, 1 );
    mScaleOverlay->SetResolution( 0.1 );
    mScaleOverlay->ValueVisibilityOff();
    mScaleOverlay->SetValue( mView->GetSurfaceOverlayOpacity() );
    mScaleOverlay->SetEndCommand( mView, "SetSurfaceOverlayOpacity");

    mbFrameOverlayInited = true;
  }
  if( NULL == mOverlayScalars.GetPointer() && mbFrameOverlayInited ) {
    mFrameOverlay->UnpackChildren();
    mFrameOverlay->SetHeight( 1 );
    mScaleOverlay = NULL;
    mbFrameOverlayInited = false;
  }


  // If we have any scalars, populate the surface scalars frame.
  if( maSurfaceScalars.size() > 0 && !mbFrameSurfaceScalarsInited ) {

    vtkSmartPointer<vtkKWFrameWithLabel> scalarsFrame =
      vtkSmartPointer<vtkKWFrameWithLabel>::New();
    scalarsFrame->SetParent( mFrameSurfaceScalars );
    scalarsFrame->Create();
    scalarsFrame->SetLabelText( "Scalars" );
    this->Script( "pack %s -side top -fill both",
                  scalarsFrame->GetWidgetName() );

    // Make an multi column list; we'll populate with a list of our
    // scalars here.
    vtkSmartPointer<vtkKWMultiColumnListWithScrollbars> scrolls =
      vtkSmartPointer<vtkKWMultiColumnListWithScrollbars>::New();
    scrolls->SetParent( scalarsFrame->GetFrame() );
    scrolls->Create();
    scrolls->HorizontalScrollbarVisibilityOff();
    this->Script( "pack %s -side top -fill both",
                  scrolls->GetWidgetName() );

    // Save a pointer to the table inside. Set it up with callbacks
    // and nice colors.
    mTableSurfaceScalars = scrolls->GetWidget();
    mTableSurfaceScalars->AddColumn( "Description" );
    mTableSurfaceScalars->SetSelectionModeToSingle();
    mTableSurfaceScalars->SetSelectionChangedCommand
      ( this,
        "SetCurrentSurfaceScalarsFromTableSelection" );
    mTableSurfaceScalars->SetSelectionBackgroundColor( 0.5, 0.8, 0.5 );
    mTableSurfaceScalars->SetSelectionForegroundColor( 0.0, 0.0, 0.0 );
    mTableSurfaceScalars->SetColumnFormatCommandToEmptyOutput( 0 );
    mTableSurfaceScalars->
      SetPotentialCellColorsChangedCommand
      ( mTableSurfaceScalars,
        "ScheduleRefreshColorsOfAllCellsWithWindowCommand" );

    // Checkbox for reverse values.
    vtkSmartPointer<vtkKWCheckButtonWithLabel> labeledCheckReverse =
      vtkSmartPointer<vtkKWCheckButtonWithLabel>::New();
    labeledCheckReverse->SetParent( scalarsFrame->GetFrame() );
    labeledCheckReverse->Create();
    labeledCheckReverse->SetLabelText( "Reverse Values" );
    labeledCheckReverse->SetLabelPositionToRight();
    this->Script( "pack %s -side top",
                  labeledCheckReverse->GetWidgetName() );

    mCheckSurfaceScalarsColorReverse = labeledCheckReverse->GetWidget();
    mCheckSurfaceScalarsColorReverse->
      SetSelectedState( mbSurfaceScalarsColorReverse );
    mCheckSurfaceScalarsColorReverse->
      SetCommand( this, "SetSurfaceScalarsColorReverse" );


    // Checkboxes for positve and negative values (inside an inner
    // frame).
    vtkSmartPointer<vtkKWFrame> framePosNeg = 
      vtkSmartPointer<vtkKWFrame>::New();
    framePosNeg->SetParent( scalarsFrame->GetFrame() );
    framePosNeg->Create();
    this->Script( "pack %s -side top -fill both",
                  framePosNeg->GetWidgetName() );

    vtkSmartPointer<vtkKWLabel> labelShowValues = 
      vtkSmartPointer<vtkKWLabel>::New();
    labelShowValues->SetParent( framePosNeg );
    labelShowValues->Create();
    labelShowValues->SetText( "Show Values: " );

    vtkSmartPointer<vtkKWCheckButtonWithLabel> labeledCheckPos =
      vtkSmartPointer<vtkKWCheckButtonWithLabel>::New();
    labeledCheckPos->SetParent( framePosNeg );
    labeledCheckPos->Create();
    labeledCheckPos->SetLabelText( "Positive" );
    labeledCheckPos->SetLabelPositionToRight();

    mCheckSurfaceScalarsColorShowPositive = labeledCheckPos->GetWidget();
    mCheckSurfaceScalarsColorShowPositive->
      SetSelectedState(mbSurfaceScalarsColorShowPositive);
    mCheckSurfaceScalarsColorShowPositive->
      SetCommand( this, "SetSurfaceScalarsColorShowPositive" );

    vtkSmartPointer<vtkKWCheckButtonWithLabel> labeledCheckNeg =
      vtkSmartPointer<vtkKWCheckButtonWithLabel>::New();
    labeledCheckNeg->SetParent( framePosNeg );
    labeledCheckNeg->Create();
    labeledCheckNeg->SetLabelText( "Negative" );
    labeledCheckNeg->SetLabelPositionToRight();

    mCheckSurfaceScalarsColorShowNegative = labeledCheckNeg->GetWidget();
    mCheckSurfaceScalarsColorShowNegative->
      SetSelectedState(mbSurfaceScalarsColorShowNegative);
    mCheckSurfaceScalarsColorShowNegative->
      SetCommand( this, "SetSurfaceScalarsColorShowNegative" );
    
    // Pack the inner frame, then the frame.
    this->Script( "pack %s %s %s -side left -fill both",
                  labelShowValues->GetWidgetName(),
                  labeledCheckPos->GetWidgetName(),
                  labeledCheckNeg->GetWidgetName() );
    this->Script( "pack %s -side top -fill both",
                  framePosNeg->GetWidgetName() );

    // We switch here on the flag to use the histogram editor. If we
    // use the editor, make it now and pack all the entries in the
    // user frame. If not, just create a frame for the entries.

    // Here begins the histogram. It takes a bunch of code. We'll also
    // create it with a user frame, which we'll populate after this.
    vtkSmartPointer<vtkKWFrame> editorUserFrame;
    if( mbUseHistogramEditor ) {
      
      mEditorSurfaceScalarColors = 
        vtkSmartPointer<vtkKWRGBATransferFunctionEditor>::New();
      mEditorSurfaceScalarColors->SetParent( scalarsFrame->GetFrame() );
      
      mEditorSurfaceScalarColors->ParameterRangeVisibilityOff();
      mEditorSurfaceScalarColors->ValueRangeVisibilityOff();
      mEditorSurfaceScalarColors->ColorSpaceOptionMenuVisibilityOff();
      mEditorSurfaceScalarColors->ColorRampVisibilityOff();
      mEditorSurfaceScalarColors->PointEntriesVisibilityOff();
      mEditorSurfaceScalarColors->ValueEntriesVisibilityOff();
      mEditorSurfaceScalarColors->PointIndexVisibilityOff();
      mEditorSurfaceScalarColors->SharpnessEntryVisibilityOff();
      mEditorSurfaceScalarColors->MidPointEntryVisibilityOff();
      mEditorSurfaceScalarColors->FunctionLineVisibilityOff();
      mEditorSurfaceScalarColors->UserFrameVisibilityOn();
      
      mEditorSurfaceScalarColors->Create();
      
      mEditorSurfaceScalarColors->ExpandCanvasWidthOn();
      //  mEditorSurfaceScalarColors->SetCanvasWidth( 450 );
      mEditorSurfaceScalarColors->SetCanvasHeight( 150 );
      mEditorSurfaceScalarColors->SetLabelText("Colors");
      mEditorSurfaceScalarColors->SetRangeLabelPositionToTop();
      mEditorSurfaceScalarColors->SetPointPositionInValueRangeToTop();
      mEditorSurfaceScalarColors->SetPointStyleToCursorDown();
      mEditorSurfaceScalarColors->PointGuidelineVisibilityOn();
      mEditorSurfaceScalarColors->SelectedPointIndexVisibilityOn();
      mEditorSurfaceScalarColors->SetLabelPositionToTop();
      mEditorSurfaceScalarColors->ParameterTicksVisibilityOn();
      mEditorSurfaceScalarColors->ComputeValueTicksFromHistogramOn();
      mEditorSurfaceScalarColors->SetParameterTicksFormat("%-#6.0f");
      mEditorSurfaceScalarColors->
        SetFunctionChangedCommand( this, "SurfaceScalarColorsEditorChanged" );
      
      // This is the histogram that will go in the editor. We make it with
      // some dummy data so that it will take up space in the beginning.
      mHistogramSurfaceScalarColors = vtkSmartPointer<vtkKWHistogram>::New();
      vtkSmartPointer<vtkFloatArray> dummy =
        vtkSmartPointer<vtkFloatArray>::New();
      dummy->Allocate( 2 );
      dummy->SetNumberOfComponents( 1 );
      dummy->InsertNextValue( -mSurfaceScalarsColorMax );
      dummy->InsertNextValue( mSurfaceScalarsColorMax );
      mHistogramSurfaceScalarColors->BuildHistogram( dummy, 0 );
      mEditorSurfaceScalarColors->SetHistogram
        ( mHistogramSurfaceScalarColors );
      // Pack the editor.
      this->Script( "pack %s -side top -fill both",
                    mEditorSurfaceScalarColors->GetWidgetName() );

      // Now we'll get the user frame and populate it.
      editorUserFrame = mEditorSurfaceScalarColors->GetUserFrame();

    } else {

      editorUserFrame = vtkSmartPointer<vtkKWFrame>::New();
      editorUserFrame->SetParent( scalarsFrame->GetFrame() );
      editorUserFrame->Create();

      // Pack the editor.
      this->Script( "pack %s -side top -fill both",
                    editorUserFrame->GetWidgetName() );
    }

    vtkSmartPointer<vtkKWFrame> frameEntries = 
      vtkSmartPointer<vtkKWFrame>::New();
    frameEntries->SetParent( editorUserFrame );
    frameEntries->Create();

    // First do three entries with our min/mid/max values. These can
    // be used to set the values directly, and will also show the
    // values in the editor.
    vtkSmartPointer<vtkKWEntryWithLabel> labeledEntryMin =
      vtkSmartPointer<vtkKWEntryWithLabel>::New();
    labeledEntryMin->SetParent( frameEntries );
    labeledEntryMin->SetLabelText( "Min: " );
    labeledEntryMin->Create();
    mEntrySurfaceScalarsColorMin = labeledEntryMin->GetWidget();
    mEntrySurfaceScalarsColorMin->SetWidth( 6 );
    mEntrySurfaceScalarsColorMin->SetRestrictValueToDouble();
    mEntrySurfaceScalarsColorMin->
      SetCommand( this, "SetSurfaceScalarsColorMin");
    mEntrySurfaceScalarsColorMin->SetValueAsDouble( mSurfaceScalarsColorMin );

    vtkSmartPointer<vtkKWEntryWithLabel> labeledEntryMid =
      vtkSmartPointer<vtkKWEntryWithLabel>::New();
    labeledEntryMid->SetParent( frameEntries );
    labeledEntryMid->SetLabelText( "Mid: " );
    labeledEntryMid->Create();
    mEntrySurfaceScalarsColorMid = labeledEntryMid->GetWidget();
    mEntrySurfaceScalarsColorMid->SetWidth( 6 );
    mEntrySurfaceScalarsColorMid->SetRestrictValueToDouble();
    mEntrySurfaceScalarsColorMid->
      SetCommand( this, "SetSurfaceScalarsColorMid");
    mEntrySurfaceScalarsColorMid->SetValueAsDouble( mSurfaceScalarsColorMid );

    vtkSmartPointer<vtkKWEntryWithLabel> labeledEntryMax = 
      vtkSmartPointer<vtkKWEntryWithLabel>::New();
    labeledEntryMax->SetParent( frameEntries );
    labeledEntryMax->SetLabelText( "Max: " );
    labeledEntryMax->Create();
    mEntrySurfaceScalarsColorMax = labeledEntryMax->GetWidget();
    mEntrySurfaceScalarsColorMax->SetWidth( 6 );
    mEntrySurfaceScalarsColorMax->SetRestrictValueToDouble();
    mEntrySurfaceScalarsColorMax->
      SetCommand( this, "SetSurfaceScalarsColorMax");
    mEntrySurfaceScalarsColorMax->SetValueAsDouble( mSurfaceScalarsColorMax );
    this->Script( "pack %s %s %s -side left -padx 5 -fill x",
                  labeledEntryMin->GetWidgetName(),
                  labeledEntryMid->GetWidgetName(),
                  labeledEntryMax->GetWidgetName() );

    // Now a frame for the button to set the color scale using FDR,
    // and its rate entry.
    vtkSmartPointer<vtkKWFrame> frameFDR =
      vtkSmartPointer<vtkKWFrame>::New();
    frameFDR->SetParent( editorUserFrame );
    frameFDR->Create();

    vtkSmartPointer<vtkKWPushButton> buttonFDR = 
      vtkSmartPointer<vtkKWPushButton>::New();
    buttonFDR->SetParent( frameFDR );
    buttonFDR->Create();
    buttonFDR->SetText( "Set Using FDR" );
    buttonFDR->SetCommand( this, "SetSurfaceScalarsColorsUsingFDR" );

    vtkSmartPointer<vtkKWEntryWithLabel> labeledEntryFDR =
      vtkSmartPointer<vtkKWEntryWithLabel>::New();
    labeledEntryFDR->SetParent( frameFDR );
    labeledEntryFDR->SetLabelText( "Rate: " );
    labeledEntryFDR->Create();
    labeledEntryFDR->GetWidget()->SetWidth( 6 );
    labeledEntryFDR->GetWidget()->SetRestrictValueToDouble();
    labeledEntryFDR->GetWidget()->
      SetCommand( this, "SetSurfaceScalarsColorsFDRRate" );
    labeledEntryFDR->GetWidget()->
      SetValueAsDouble( mSurfaceScalarsColorsFDRRate );
    this->Script( "pack %s -side left -expand y -fill x -padx 5",
                  buttonFDR->GetWidgetName() );
    this->Script( "pack %s -side left -fill x -padx 5",
                  labeledEntryFDR->GetWidgetName() );
    
    // Single entry for the offet.
    vtkSmartPointer<vtkKWEntryWithLabel> labeledEntryOffset =
      vtkSmartPointer<vtkKWEntryWithLabel>::New();
    labeledEntryOffset->SetParent( editorUserFrame );
    labeledEntryOffset->SetLabelText( "Value offset: " );
    labeledEntryOffset->Create();
    mEntrySurfaceScalarsColorOffset = labeledEntryOffset->GetWidget();
    mEntrySurfaceScalarsColorOffset->SetWidth( 6 );
    mEntrySurfaceScalarsColorOffset->SetRestrictValueToDouble();
    mEntrySurfaceScalarsColorOffset->
      SetCommand( this, "SetSurfaceScalarsColorOffset");
    mEntrySurfaceScalarsColorOffset->
      SetValueAsDouble( mSurfaceScalarsColorOffset );

    // Pack the three items in the user frame.
    this->Script( "pack %s %s %s -side top -fill x -pady 2",
                  frameEntries->GetWidgetName(), 
                  frameFDR->GetWidgetName(),
                  labeledEntryOffset->GetWidgetName() );

    // Update the editors with the current min/mid/max values and a
    // dummy histogram.
    this->UpdateSurfaceScalarsColorsEditor();

    mbFrameSurfaceScalarsInited = true;
  }
  if( maSurfaceScalars.size() == 0 && mbFrameSurfaceScalarsInited ) {
    mFrameSurfaceScalars->UnpackChildren();
    mFrameSurfaceScalars->SetHeight( 1 );
    mTableSurfaceScalars = NULL;
    mCheckSurfaceScalarsColorShowPositive = NULL;
    mCheckSurfaceScalarsColorShowNegative = NULL;
    mEditorSurfaceScalarColors = NULL;
    mHistogramSurfaceScalarColors = NULL;
    mEntrySurfaceScalarsColorMin = NULL;
    mEntrySurfaceScalarsColorMid = NULL;
    mEntrySurfaceScalarsColorMax = NULL;
    mEntrySurfaceScalarsColorOffset = NULL;
    mbFrameSurfaceScalarsInited = false;
  }

  // If we have our table of scalars, refresh the table of scalars
  // with the current scalars.
  if( mTableSurfaceScalars.GetPointer() && mbFrameSurfaceScalarsInited ) {

    mTableSurfaceScalars->DeleteAllRows();
    map<int,SurfaceScalar>::iterator tScalar;
    for( tScalar = maSurfaceScalars.begin();
         tScalar != maSurfaceScalars.end();
         ++tScalar ) {

      // Set the cell text.
      mTableSurfaceScalars->
        InsertCellText( (*tScalar).second.mnEntry, 0,
                        (*tScalar).second.msLabel.c_str() );

      // When this cell is actually created, call CreateScalarTableEntry.
      mTableSurfaceScalars->SetCellWindowCommand
        ( (*tScalar).second.mnEntry, 0,
          this, "CreateScalarTableEntry" );

    }
  }
}

void
vtkKWQdecWindow::UpdateSurfaceScalarsColorsEditor () {

  if( NULL == mEditorSurfaceScalarColors.GetPointer() )
    return;

  assert( mHistogramSurfaceScalarColors.GetPointer() );

  // Create a color table and add our poitns. This isn't our actual
  // color table that we'll use; we just use this to show the settable
  // values in the editor the values. When the values change, we
  // update our internal min/mid/max variables, and compose the
  // scalars.
  vtkSmartPointer<vtkRGBATransferFunction> colors =
    vtkSmartPointer<vtkRGBATransferFunction>::New();
  mEditorSurfaceScalarColors->SetRGBATransferFunction( colors );
  mnNegativeMaxValue =
    colors->AddRGBAPoint( -mSurfaceScalarsColorMax, 0, 1, 1, 1 );
  mnNegativeMidValue =
    colors->AddRGBAPoint( -mSurfaceScalarsColorMid, 0, 0, 1, 1 );
  mnNegativeMinValue =
    colors->AddRGBAPoint( -mSurfaceScalarsColorMin, 0, 0, 1, 1 );
  mnPositiveMinValue =
    colors->AddRGBAPoint(  mSurfaceScalarsColorMin, 1, 0, 0, 1 );
  mnPositiveMidValue =
    colors->AddRGBAPoint(  mSurfaceScalarsColorMid, 1, 0, 0, 1 );
  mnPositiveMaxValue =
    colors->AddRGBAPoint(  mSurfaceScalarsColorMax, 1, 1, 0, 1 );
  colors->Build();

  // Set up the point symmetry in the colors.
  mEditorSurfaceScalarColors->SetPointCountMinimum( 6 );
  mEditorSurfaceScalarColors->SetPointCountMaximum( 6 );
  mEditorSurfaceScalarColors->
    SetPointSymmetry( mnNegativeMaxValue, mnPositiveMaxValue );
  mEditorSurfaceScalarColors->
    SetPointSymmetry( mnNegativeMidValue, mnPositiveMidValue );
  mEditorSurfaceScalarColors->
    SetPointSymmetry( mnNegativeMinValue, mnPositiveMinValue );

  // If we have scalars right now, make a histogram of them, set it in
  // the editor, and make sure the range for the editor is correct.
  if( maSurfaceScalars.find( mnCurrentSurfaceScalars ) !=
      maSurfaceScalars.end() ) {

    vtkFloatArray* surfaceScalars =
      maSurfaceScalars[mnCurrentSurfaceScalars].mValues;

    // Make the histogram.
    mHistogramSurfaceScalarColors->BuildHistogram( surfaceScalars, 0 );
    mEditorSurfaceScalarColors->
      SetHistogram( mHistogramSurfaceScalarColors );

    // Set the range to show the full abs of the range, so that
    // we'll see the -max as well as the max. Make sure that this
    // also includes the max color.
    double range[2];
    surfaceScalars->GetRange( range );
    double highestAbs = MAX( MAX( fabs(range[0]), fabs(range[1]) ),
                             mSurfaceScalarsColorMax );
    mEditorSurfaceScalarColors->
      SetWholeParameterRange( -highestAbs, highestAbs );
    mEditorSurfaceScalarColors->
      SetVisibleParameterRangeToWholeParameterRange();
    
  } else {

    // No scalars, so just make some dummy data.
    vtkSmartPointer<vtkFloatArray> dummy = 
      vtkSmartPointer<vtkFloatArray>::New();
    dummy->Allocate( 2 );
    dummy->SetNumberOfComponents( 1 );
    dummy->InsertNextValue( -mSurfaceScalarsColorMax );
    dummy->InsertNextValue( mSurfaceScalarsColorMax );
    mHistogramSurfaceScalarColors->BuildHistogram( dummy, 0 );

    // Set the dummy histogram.
    mEditorSurfaceScalarColors->SetHistogram( mHistogramSurfaceScalarColors );

    // Adjust the range.
    mEditorSurfaceScalarColors->
      SetVisibleParameterRangeToWholeParameterRange();
  }
}

void
vtkKWQdecWindow::CreateScalarTableEntry ( const char* iTable,
                                          int iRow, int iColumn,
                                          const char* iWidget ) {
  
  assert( mTableSurfaceScalars.GetPointer() );
  assert( iTable );
  assert( iWidget );

  // In this function, we are given a widget name, and we can create a
  // widget for that widget name. This lets us customize the contents
  // of a cell. We'll create a label widget so that we can use line
  // breaks to the contents from scrolling.

  // Create a label. Set the widget name to the one we got.
  vtkSmartPointer<vtkKWLabel> label =
    vtkSmartPointer<vtkKWLabel>::New();
  label->SetParent( mTableSurfaceScalars );
  label->SetWidgetName( iWidget );
  label->Create();

  // Match the colors with the table default.
  double r = 0, g = 0, b = 0;
  mTableSurfaceScalars->GetCellBackgroundColor( iRow, iColumn, &r, &g, &b );
  label->SetBackgroundColor( r, g, b );
  
  mTableSurfaceScalars->GetCellForegroundColor( iRow, iColumn, &r, &g, &b );
  label->SetForegroundColor( r, g, b );

  // Let the table add bindings to this widget.
  mTableSurfaceScalars->AddBindingsToWidgetName( iWidget );

  // Set the string. The row is the scalar index.
  string sLabel = mTableSurfaceScalars->GetCellText( iRow, iColumn );
  label->SetText( sLabel.c_str() );

  // Set wrapping and font.
  label->SetWrapLength( "200p" );
  label->SetFont( "helvetica 9" );
}

void
vtkKWQdecWindow::SurfaceScalarColorsEditorChanged () {

  assert( mEditorSurfaceScalarColors.GetPointer() );

  // Get the value from the table and see if it is different from the
  // value we have. If so, set it. Note that the value at
  // mnPositiveMinValue will be the same as the value at
  // -mnNegativeMinValue beause the editor keeps them symmetric.
  double value;
  mEditorSurfaceScalarColors->
    GetFunctionPointParameter( mnPositiveMinValue, &value );
  if( mSurfaceScalarsColorMin != value )
    this->SetSurfaceScalarsColorMin( value );
  
  mEditorSurfaceScalarColors->
    GetFunctionPointParameter( mnPositiveMidValue, &value );
  if( mSurfaceScalarsColorMid != value )
    this->SetSurfaceScalarsColorMid( value );
  
  mEditorSurfaceScalarColors->
    GetFunctionPointParameter( mnPositiveMaxValue, &value );
  if( mSurfaceScalarsColorMax != value )
    this->SetSurfaceScalarsColorMax( value );
}

void
vtkKWQdecWindow::SetSurfaceScalarsColorMin ( double iMin ) {

  assert( mEntrySurfaceScalarsColorMin.GetPointer() );

  // If not already equal...
  if( fabs( iMin - mSurfaceScalarsColorMin ) >
      numeric_limits<double>::epsilon() ) {

    // Save our new value.
    mSurfaceScalarsColorMin = iMin;

    // Set the value in the entry if not equal.
    if( fabs( mEntrySurfaceScalarsColorMin->GetValueAsDouble() -
              mSurfaceScalarsColorMin ) > numeric_limits<double>::epsilon() )
      mEntrySurfaceScalarsColorMin->SetValueAsDouble
        ( mSurfaceScalarsColorMin );

    // Update the editor.
    this->UpdateSurfaceScalarsColorsEditor();

    // Draw with the new values.
    this->ComposeSurfaceScalarsAndShow();
  }
}

void
vtkKWQdecWindow::SetSurfaceScalarsColorMid ( double iMid ) {

  assert( mEntrySurfaceScalarsColorMin.GetPointer() );

  // If not already equal...
  if( fabs( iMid - mSurfaceScalarsColorMid ) >
      numeric_limits<double>::epsilon() ) {

    // Save our new value.
    mSurfaceScalarsColorMid = iMid;

    // Set the value in the entry if not equal.
    if( fabs( mEntrySurfaceScalarsColorMid->GetValueAsDouble() -
              mSurfaceScalarsColorMid ) > numeric_limits<double>::epsilon() )
      mEntrySurfaceScalarsColorMid->SetValueAsDouble
        ( mSurfaceScalarsColorMid );

    // Update the editor.
    this->UpdateSurfaceScalarsColorsEditor();

    // Draw with the new values.
    this->ComposeSurfaceScalarsAndShow();
  }
}

void
vtkKWQdecWindow::SetSurfaceScalarsColorMax ( double iMax ) {

  assert( mEntrySurfaceScalarsColorMax.GetPointer() );

  // If not already equal...
  if( fabs( iMax - mSurfaceScalarsColorMax ) >
      numeric_limits<double>::epsilon() ) {

    // Save our new value.
    mSurfaceScalarsColorMax = iMax;

    // Set the value in the entry if not equal.
    if( fabs( mEntrySurfaceScalarsColorMax->GetValueAsDouble() -
              mSurfaceScalarsColorMax ) > numeric_limits<double>::epsilon() )
      mEntrySurfaceScalarsColorMax->SetValueAsDouble
        ( mSurfaceScalarsColorMax );

    // Update the editor.
    this->UpdateSurfaceScalarsColorsEditor();

    // Draw with the new values.
    this->ComposeSurfaceScalarsAndShow();
  }
}

void
vtkKWQdecWindow::SetSurfaceScalarsColors ( double iMin, double iMid,
                                           double iMax ) {

  assert( mEntrySurfaceScalarsColorMin.GetPointer() );
  assert( mEntrySurfaceScalarsColorMid.GetPointer() );
  assert( mEntrySurfaceScalarsColorMax.GetPointer() );

  // Does the same thing as the above colors, except all at once to
  // avoid three redraws.
  bool bChanged = false;

  // If not already equal...
  if( fabs( iMin - mSurfaceScalarsColorMin ) >
      numeric_limits<double>::epsilon() ) {

    // Save our new value.
    mSurfaceScalarsColorMin = iMin;

    // Set the value in the entry if not equal.
    if( fabs( mEntrySurfaceScalarsColorMin->GetValueAsDouble() -
              mSurfaceScalarsColorMin ) > numeric_limits<double>::epsilon() )
      mEntrySurfaceScalarsColorMin->SetValueAsDouble
        ( mSurfaceScalarsColorMin );

    bChanged = true;
  }


  // If not already equal...
  if( fabs( iMid - mSurfaceScalarsColorMid ) >
      numeric_limits<double>::epsilon() ) {

    // Save our new value.
    mSurfaceScalarsColorMid = iMid;

    // Set the value in the entry if not equal.
    if( fabs( mEntrySurfaceScalarsColorMid->GetValueAsDouble() -
              mSurfaceScalarsColorMid ) > numeric_limits<double>::epsilon() )
      mEntrySurfaceScalarsColorMid->SetValueAsDouble
        ( mSurfaceScalarsColorMid );

    bChanged = true;
  }


  // If not already equal...
  if( fabs( iMax - mSurfaceScalarsColorMax ) >
      numeric_limits<double>::epsilon() ) {

    // Save our new value.
    mSurfaceScalarsColorMax = iMax;

    // Set the value in the entry if not equal.
    if( fabs( mEntrySurfaceScalarsColorMax->GetValueAsDouble() -
              mSurfaceScalarsColorMax ) > numeric_limits<double>::epsilon() )
      mEntrySurfaceScalarsColorMax->SetValueAsDouble
        ( mSurfaceScalarsColorMax );

    bChanged = true;
  }


  if( bChanged ) {

    // Update the editor.
    this->UpdateSurfaceScalarsColorsEditor();

    // Draw with the new values.
    this->ComposeSurfaceScalarsAndShow();
  }
}

void
vtkKWQdecWindow::SetSurfaceScalarsColorOffset ( double iOffset ) {

  assert( mEntrySurfaceScalarsColorOffset.GetPointer() );

  // If not already equal...
  if( fabs( iOffset - mSurfaceScalarsColorOffset ) >
      numeric_limits<double>::epsilon() ) {

    // Save our new value.
    mSurfaceScalarsColorOffset = iOffset;

    // Set the value in the entry if not equal.
    if( fabs( mEntrySurfaceScalarsColorOffset->GetValueAsDouble() -
              mSurfaceScalarsColorOffset ) > 
        numeric_limits<double>::epsilon() )
      mEntrySurfaceScalarsColorOffset->SetValueAsDouble
        ( mSurfaceScalarsColorOffset );

    // Draw with the new values.
    this->ComposeSurfaceScalarsAndShow();
  }
}

void
vtkKWQdecWindow::SetSurfaceScalarsColorsUsingFDR () {

  this->SetStatusText( "Setting threshold using FDR..." );
  try {

    if( maSurfaceSource.find( msCurrentSurfaceSource ) ==
        maSurfaceSource.end() ||
        maSurfaceScalars.find( mnCurrentSurfaceScalars ) ==
        maSurfaceScalars.end() ) {
      throw runtime_error( "Must have a surface loaded and scalars selected.");
    }

    // If they are only looking at positive or negative values, set the
    // sign. If they are reversing, reverse the sign.
    int sign = 0;
    if( !mbSurfaceScalarsColorShowPositive )
      sign = -1;
    if( !mbSurfaceScalarsColorShowNegative )
      sign = 1;
    if( mbSurfaceScalarsColorReverse )
      sign = -sign;

    // The FDR function will look at values in the val field of the
    // surface. Get the MRIS and write our current scalar values into
    // the val field.
    MRIS* mris = maSurfaceSource[msCurrentSurfaceSource]->GetMRIS();
    vtkFloatArray* scalars = maSurfaceScalars[mnCurrentSurfaceScalars].mValues;
    for( int nVertex = 0; nVertex < mris->nvertices; nVertex++ ) {
      mris->vertices[nVertex].val = scalars->GetTuple1( nVertex );
      mris->vertices[nVertex].ripflag = 0;
    }

    // Call the FDR function. We pass our surface, the rate, the sign we
    // got before, FALSE for only_marked (not supported), and get the
    // threshold back.
    double threshold;
    int err = MRISfdr2vwth( mris, mSurfaceScalarsColorsFDRRate, sign,
                            true,  // log10flag
                            false, // maskflag
                            &threshold);
    if( err )
      throw runtime_error( "Error from MRISfdr2vwth" );
    else
      this->SetStatusText( "Completed FDR threshold set" );

    // Set our min to the threshold, and calculate good values for mid
    // and max.
    this->SetSurfaceScalarsColors( threshold,
                                   threshold + 1.5,
                                   threshold + 2.25 );

  } catch (exception& e) {
    stringstream ssError;
    ssError << "Error in calculating FDR: " << e.what();
    this->GetApplication()->ErrorMessage( ssError.str().c_str() );
    this->SetStatusText( "Error during FDR calculation" );
  }
}

void
vtkKWQdecWindow::SetSurfaceScalarsColorsFDRRate ( const char* isValue ) {

  stringstream ssValue( isValue );
  ssValue >> mSurfaceScalarsColorsFDRRate;
}

void
vtkKWQdecWindow::UpdateContinuousFactorPlot () {

  assert( mGraph.GetPointer() );
  assert( mQdecProject );

  // Remove all the elements in the graph.
  mGraph->DeleteAllElements();

  if( mPlotContinuousFactorSelection != -1 ) {
      
    // Start a list of points for our data to graph.
    vector<double> lPoints;
    
    // Get the name of the factor.
    string sFactor( mListPlotContinuousFactors->GetItem
                    (mPlotContinuousFactorSelection) );
    
    // For each subject...
    vector<QdecSubject*>::iterator tSubject;
    vector<QdecSubject*> lSubjects = 
      mQdecProject->GetDataTable()->GetSubjects();
    int nIndex = 0;
    for( tSubject = lSubjects.begin();
         tSubject != lSubjects.end();
         ++tSubject, nIndex++ ) {
      
      // Get the value.
      float value = (*tSubject)->GetContinuousFactor( sFactor.c_str() );
      
      // Add this index,value pair to the list of values to graph.
      lPoints.clear();
      lPoints.push_back( nIndex );
      lPoints.push_back( value );

      // We'll use a blue plus if it's not excluded and an red x if
      // so.
      QdecGlmDesign* design = mQdecProject->GetGlmDesign();
      assert( design );
      string sSymbol;
      float color[3] = {0, 0, 0};
      if( design->GetExcludeSubjectID( (*tSubject)->GetId().c_str() ) ) {
        sSymbol = "cross";
        color[0] = 1.0; // red
      } else {
        sSymbol = "plus";
        color[2] = 1.0; // blue
      }

      // Add the data to the graph.
      mGraph->AddElement( (*tSubject)->GetId().c_str(),
                          lPoints,
                          sSymbol.c_str(),
                          0,
                          color[0], color[1], color[2] );
    }
  }

  // Draw the graph.
  mGraph->Draw();

}

void
vtkKWQdecWindow::SelectSurfaceVertex ( int inVertex ) {

  assert( mView.GetPointer() );
  assert( mEntrySelectVertex.GetPointer() );

  try {

    // Try to go to this vertex.
    mView->SelectSurfaceVertex( inVertex );

  } catch (exception& e) {

    // If we fail, give an error and try to suggest a better index.
    this->GetApplication()->ErrorMessage( e.what() );
    if( inVertex < 0 )
      mEntrySelectVertex->SetValueAsInt( 0 );
    else if( inVertex >= mcVertices )
      mEntrySelectVertex->SetValueAsInt( mcVertices-1 );
  }

}

void
vtkKWQdecWindow::AddSelectionToROI () {

  assert( mView.GetPointer() );
  if( maSurfaceSource.size() < 1 )
    throw runtime_error( "Must load a surface before "
                         "working on the selection." );

  // Get the list of selected points. This is a closed loop.
  vtkPoints* points = mView->GetSurfaceSelectionPoints();
  assert( points );
  if( points->GetNumberOfPoints() < 4 )
    return;

  // Get the vertices in this path.
  vtkSmartPointer<vtkIntArray> vertices = 
    vtkSmartPointer<vtkIntArray>::New();
  this->FindVerticesInsidePath( points, vertices );

  // If we don't have a label source, create one now.
  if( !mROISource.GetPointer() ) {
    mROISource = vtkSmartPointer<vtkFSSurfaceLabelSource>::New();
    string sLabelName = mMenuHemisphere->GetValue();
    sLabelName += ".untitled.label";
    mROISource->SetLabelFileName( sLabelName.c_str() );
    vtkFSSurfaceSource* surface = maSurfaceSource[msCurrentSurfaceSource];
    mROISource->SetMris( surface->GetMRIS() );
    mROISource->InitializeEmptyLabel();
    mView->SetROI( mROISource->GetOutput() );
  }

  // Add this list of vertices to the label.
  mROISource->AddVerticesToLabel( vertices->GetNumberOfTuples(),
                                  vertices->GetPointer(0) );

  // Redraw the ROI.
  mROISource->Update();

  mView->Render();

  this->UpdateCommandStatus();

}

void
vtkKWQdecWindow::RemoveSelectionFromROI () {

  assert( mView.GetPointer() );
  if( maSurfaceSource.size() < 1 )
    throw runtime_error( "Must load a surface "
                         "before working on the selection." );

  // Get the list of selected points. This is a closed loop.
  vtkPoints* points = mView->GetSurfaceSelectionPoints();
  assert( points );
  if( points->GetNumberOfPoints() < 4 )
    return;

  // Get the vertices in this path.
  vtkSmartPointer<vtkIntArray> vertices = 
    vtkSmartPointer<vtkIntArray>::New();
  this->FindVerticesInsidePath( points, vertices );

  // If we don't have a label source, create one now.
  if( !mROISource.GetPointer() ) {
    mROISource = vtkSmartPointer<vtkFSSurfaceLabelSource>::New();
    string sLabelName = mMenuHemisphere->GetValue();
    sLabelName += ".untitled.label";
    mROISource->SetLabelFileName( sLabelName.c_str() );
    vtkFSSurfaceSource* surface = maSurfaceSource.begin()->second;
    mROISource->SetMris( surface->GetMRIS() );
    mROISource->InitializeEmptyLabel();
    mView->SetROI( mROISource->GetOutput() );
  }

  // Add this list of vertices to the label.
  mROISource->RemoveVerticesFromLabel( vertices->GetNumberOfTuples(),
                                       vertices->GetPointer(0) );

  // Redraw the ROI.
  mROISource->Update();

  mView->Render();
}

void
vtkKWQdecWindow::ClearROI () {
  
  if( NULL == mROISource.GetPointer() )
    return;

  assert( mView.GetPointer() );

  // Empty the label.
  mROISource->InitializeEmptyLabel();

  // Redraw the ROI.
  mROISource->Update();

  mView->Render();
}

void
vtkKWQdecWindow::GraphAverageROIInGDF () {

  assert( mROISource.GetPointer() );
  assert( mbGDFLoaded );

  // Create a vertex vector.
  vector<int> lVertices;

  // Get it filled out.
  mROISource->GetLabeledVertices( lVertices );

  // Make sure we got something.
  if( 0 == lVertices.size() ) {
    vtkSmartPointer<vtkKWMessageDialog> dialog =
      vtkSmartPointer<vtkKWMessageDialog>::New();

    dialog->SetStyleToMessage();
    dialog->SetOptions( vtkKWMessageDialog::WarningIcon );
    dialog->SetApplication( this->GetApplication() );
    dialog->Create();
    dialog->SetText( "The ROI has no points in it; "
                     "there is nothing to average." );
    dialog->Invoke();
    return;
  }

  // Start the point list.
  this->Script( "FsgdfPlot_BeginPointList %d", mGDFID );

  // For each vertex, add it to the point list.
  vector<int>::iterator tVertex;
  for( tVertex = lVertices.begin(); tVertex != lVertices.end(); ++tVertex ) {
    int nVertex = *tVertex;
    this->Script( "FsgdfPlot_AddPoint %d %d 0 0", mGDFID, nVertex );
  }

  // End the point list and draw the graph.
  this->Script( "FsgdfPlot_EndPointList %d", mGDFID );

  // Set the info string.
  this->Script( "FsgdfPlot_SetInfo %d \"Average of ROI\"", mGDFID );
}

void
vtkKWQdecWindow::SmoothCurvatureScalars ( int icSteps ) {

  assert( maSurfaceSource.size() > 0 );
  assert( mCurvatureScalars.GetPointer() );

  // Get the first surface to use to do the smoothing, and smooth the
  // values using that surface topography.
  maSurfaceSource.begin()->second->
    SmoothValuesOnSurface( *mCurvatureScalars, icSteps );

  // Redraw.
  this->ComposeSurfaceScalarsAndShow();
}

void
vtkKWQdecWindow::SmoothSurfaceScalars ( int icSteps ) {
  
  assert( maSurfaceSource.size() > 0 );
  assert( maSurfaceScalars.find(mnCurrentSurfaceScalars) !=
          maSurfaceScalars.end() );
  assert( maSurfaceScalars[mnCurrentSurfaceScalars].mValues.GetPointer() );

  // Get the first surface to use to do the smoothing, and smooth the
  // values using that surface topography.
  maSurfaceSource.begin()->second->
    SmoothValuesOnSurface( *maSurfaceScalars[mnCurrentSurfaceScalars].mValues, 
                           icSteps );

  // Redraw.
  this->ComposeSurfaceScalarsAndShow();
}

void
vtkKWQdecWindow::NotebookPageRaised ( const char* isTitle ) {

  assert( isTitle );
  assert( mView.GetPointer() );
  assert( mGraph.GetPointer() );

  if( 0 == strcmp( isTitle, ksSubjectsPanelName ) ) {

    this->GetViewFrame()->UnpackChildren();
    this->Script( "pack %s -expand yes -fill both -anchor c",
                  mGraph->GetWidgetName() );
  
  } else if ( 0 == strcmp( isTitle, ksDesignPanelName ) ) {

    this->GetViewFrame()->UnpackChildren();

  } else if ( 0 == strcmp( isTitle, ksDisplayPanelName ) ) {

    this->GetViewFrame()->UnpackChildren();
    this->Script( "pack %s -expand yes -fill both -anchor c",
                  mView->GetWidgetName() );
  }
}

void
vtkKWQdecWindow::NotebookRaisePageCallback ( vtkObject* iCaller,
                                             unsigned long iEventId,
                                             void* iClientData,
                                             void* iCallData ) {

  assert( vtkKWEvent::NotebookRaisePageEvent == iEventId );
  assert( iCallData );
  assert( iClientData );

  // Extract the client and call data and call the window's
  // ContinuousPlotGraphMouseoverEnterElement function.
  try {

    char** args = static_cast<char**>( iCallData );
    char* sTitle = args[0];

    vtkKWQdecWindow* window = 
      static_cast<vtkKWQdecWindow*>( iClientData );
    
    if( window )
      window->NotebookPageRaised( sTitle );
  }
  catch(...) {
    cerr << "Invalid call or client data in NotebookRaisePageCallback" << endl;
  }
}

void 
vtkKWQdecWindow::ContinuousPlotGraphMouseoverEnterElement 
( const char* isElement ) {

  // Just show it in our status bar.
  this->SetStatusText( isElement );
}

void 
vtkKWQdecWindow::ContinuousPlotGraphMouseoverExitElement () {

  // Clear our status bar.
  this->SetStatusText( "" );
}

void
vtkKWQdecWindow::ContinuousPlotGraphSetUpContextualMenu (const char* isElement,
                                                         vtkKWMenu* iMenu ) {

  assert( isElement );
  assert( iMenu );
  assert( mQdecProject );
  assert( mMenuHemisphere.GetPointer() );

  // Add an inactive item at the top with the element (subject) name.
  iMenu->AddCommand( isElement,  NULL, "" ); 
  iMenu->SetItemStateToDisabled( 1 );
  iMenu->AddSeparator();

  // Check this element. If it shows up in the design's exclusion
  // list, add a command to unexclude it. Otherwise, add a command to
  // exclude it.
  QdecGlmDesign* design = mQdecProject->GetGlmDesign();
  assert( design );
  if( design->GetExcludeSubjectID( isElement ) ) {
    stringstream ssCommand;
    ssCommand << "SetExcludeSubjectID " << isElement << " 0";
    iMenu->AddCommand( "Don't exclude", this, ssCommand.str().c_str() );
  } else {
    stringstream ssCommand;
    ssCommand << "SetExcludeSubjectID " << isElement << " 1";
    iMenu->AddCommand( "Exclude", this, ssCommand.str().c_str() );
  }

  // Check to see if we can find our freesurfer directory.
  FSENV* envVars = FSENVgetenv();
  if( NULL != envVars ) {
    
    // Add commands to load this subject's data in tksurfer and tkmedit.
    iMenu->AddSeparator();
    stringstream ssTkmeditCommand;
    // catch { exec tkmedit $s norm.mgz $hemi.white -segmentation aseg.mgz }
    ssTkmeditCommand << "set err [catch { exec "
                     << envVars->FREESURFER_HOME << "/bin/tkmedit "
                     << isElement << " norm.mgz " 
                     << mMenuHemisphere->GetValue() << ".white "
                     << "-segmentation aseg.mgz"
                     << " }] ; if { $err != 0 } { "
                     << this->GetApplication()->GetTclName() << " "
                     << "ErrorMessage \"Couldn't start tkmedit.\" }";
    iMenu->AddCommand( "Open in tkmedit ", NULL, 
                       ssTkmeditCommand.str().c_str() );
    stringstream ssTksurferCommand;
    // catch { exec tksurfer $s $hemi inflated -annotation aparc.annot }
    ssTksurferCommand << "catch { exec "
                      << envVars->FREESURFER_HOME << "/bin/tksurfer "
                      << isElement << " " 
                      << mMenuHemisphere->GetValue() << " inflated "
                      << "-annotation aparc.annot"
                      << " }";
    iMenu->AddCommand( "Open in tksurfer ", NULL, 
                       ssTksurferCommand.str().c_str() );

  } else {

    // Informative disabled menu commands.
    iMenu->AddCommand( "tkmedit not found", NULL, "" );
    iMenu->SetItemStateToDisabled( 5 );
    iMenu->AddCommand( "tksurfer not found", NULL, "" );
    iMenu->SetItemStateToDisabled( 6 );
  }

}

void
vtkKWQdecWindow::ContinuousPlotGraphMouseoverEnterElementCallback
( vtkObject* iCaller, unsigned long iEventId, void* iClientData,
  void* iCallData ) {

  assert( vtkKWBltGraph::MouseoverEnterElementEvent == iEventId );
  assert( iCallData );
  assert( iClientData );

  // Extract the client and call data and call the window's
  // ContinuousPlotGraphMouseoverEnterElement function.
  try {

    vtkKWBltGraph::SelectedElementAndPoint* foundElement =
      static_cast<vtkKWBltGraph::SelectedElementAndPoint*>( iCallData );

    vtkKWQdecWindow* window = 
      static_cast<vtkKWQdecWindow*>( iClientData );
    
    if( window )
      window->ContinuousPlotGraphMouseoverEnterElement( foundElement->msLabel);
  }
  catch(...) {
    cerr << "Invalid call or client data in "
         << "ContinuousPlotGraphMouseoverEnterElementCallback" << endl;
  }
}

void
vtkKWQdecWindow::ContinuousPlotGraphMouseoverExitElementCallback
( vtkObject* iCaller, unsigned long iEventId, void* iClientData,
  void* iCallData ) {

  assert( vtkKWBltGraph::MouseoverExitElementEvent == iEventId );
  assert( iClientData );

  // Extract the client data and call the window's
  // ContinuousPlotGraphMouseoverExitElement function.
  try {

    vtkKWQdecWindow* window = 
      static_cast<vtkKWQdecWindow*>( iClientData );
    
    if( window )
      window->ContinuousPlotGraphMouseoverExitElement();
  }
  catch(...) {
    cerr << "Invalid call or client data in "
         << "ContinuousPlotGraphMouseoverExitElementCallback" << endl;
  }
}

void
vtkKWQdecWindow::ContinuousPlotGraphContextualMenuOpeningCallback
( vtkObject* iCaller, unsigned long iEventId, void* iClientData,
  void* iCallData ) {

  assert( vtkKWBltGraph::ContextualMenuOpening == iEventId );
  assert( iClientData );
  assert( iCallData );

  // Extract the client and call data and call the window's
  // ContinuousPlotGraphSetUpContextualMenu function.
  try {

    vtkKWBltGraph::ContextualMenuElement* clickedElement =
      static_cast<vtkKWBltGraph::ContextualMenuElement*>( iCallData );

    vtkKWQdecWindow* window = 
      static_cast<vtkKWQdecWindow*>( iClientData );
    
    if( window )
      window->ContinuousPlotGraphSetUpContextualMenu(clickedElement->msElement,
                                                     clickedElement->mMenu);
  }
  catch(...) {
    cerr << "Invalid call or client data in "
         << "ContinuousPlotGraphContextualMenuOpeningCallback" << endl;
  }
}

void
vtkKWQdecWindow::SetExcludeSubjectID ( const char* isElement, 
                                       int ibExclude ) {

  assert( isElement );
  assert( mQdecProject );

  // Get the design from the project.
  QdecGlmDesign* design = mQdecProject->GetGlmDesign();
  assert( design );
  
  // Set the exclude flag in the design.
  design->SetExcludeSubjectID( isElement, ibExclude );

  // Redraw our graph.
  this->UpdateContinuousFactorPlot();
}

char const*
vtkKWQdecWindow::GetAnnotationForVertex ( int inVertex ) {

  // Bail if we have nothing to lookup?
  if( NULL == mAnnotationTable ||
      NULL == maAnnotationIndicies )
    return NULL;

  // Check for invalid vertex.
  if( inVertex < 0 || inVertex >= mcVertices )
    return NULL;

  // Get the table index and then look it up in the table.
  int nEntry = 0;
  CTABfindAnnotation( mAnnotationTable, 
                      maAnnotationIndicies[inVertex], &nEntry );
  
  // Make sure this entry is valid.
  int bValid = false;
  CTABisEntryValid( mAnnotationTable, nEntry, &bValid );
  if( bValid ) {

    // Copy it into temp story and return a pointer to it.
    static char sLabel[1024];
    CTABcopyName( mAnnotationTable, nEntry, sLabel, sizeof(sLabel) );
    return sLabel;
  }

  // Invalid entry.
  return NULL;
}

void
vtkKWQdecWindow::FindVerticesInsidePath ( vtkPoints* iPath,
                                          vtkIntArray* ioVertices ) {

  assert( iPath );
  assert( ioVertices );
  assert( maSurfaceSource.size() > 0 );

  // Get our surface.
  vtkFSSurfaceSource* surface = maSurfaceSource.begin()->second;

  // Make a selector and pass in the list of points and the surface.
  vtkSmartPointer<vtkSelectPolyData> selector =
    vtkSmartPointer<vtkSelectPolyData>::New();
  selector->SetLoop( iPath );
  selector->SetInputConnection( surface->GetOutputPort() );

  // Select the smallest region.
  selector->SetSelectionModeToSmallestRegion();
  selector->GenerateSelectionScalarsOn();

  // Run the selector.
  selector->Update();

  // Go through the scalars from the output of the selector. If the
  // value is 0 it's on line, and <0 is inside the selection. Add
  // those vertex numbers to our list.
  vtkDataArray* scalars = selector->GetOutput()->GetPointData()->GetScalars();
  for( int nVertex = 0; nVertex < scalars->GetNumberOfTuples(); nVertex++ ) {

    double value = scalars->GetTuple1( (vtkIdType)nVertex );
    if( value <= 0.0 )
      ioVertices->InsertNextValue( nVertex );
  }
}

vtkFloatArray*
vtkKWQdecWindow::NewScalarsFromSurfaceScalarsFile ( const char* ifnScalars,
                                                    int inIndex ) {

  if( maSurfaceSource.size() < 1 )
    throw runtime_error( "Must load a surface before loading scalars." );

  // Try to load the scalars.

  // Init a float array.
  vtkFloatArray* scalars = vtkFloatArray::New();

  int cValues = mcVertices;

  // Set the frame to load.
  MRISsetReadFrame( inIndex );

  // Try to read the scalars.
  float* values = NULL;
  int eRead = MRISreadValuesIntoArray( ifnScalars, cValues, &values );
  if( 0 != eRead ) {
    if( values ) free( values );
    throw runtime_error ("Could not read scalar file");
  }

  // Allocate our scalars.
  scalars->Allocate( cValues );
  scalars->SetNumberOfComponents( 1 );

  // Copy our array into the scalars.
  for( int nValue = 0; nValue < cValues; nValue ++ )
    scalars->InsertNextValue( values[nValue] );

  // MRISreadValuesIntoArray allocated the array, so we free it.
  free( values );

  return scalars;
}

void
vtkKWQdecWindow::ComposeSurfaceScalarsAndShow () {

  if( mView.GetPointer() ) {

    // If no scalars to show, set NULL everything and return.
    if( NULL == mCurvatureScalars.GetPointer() &&
        maSurfaceScalars.find(mnCurrentSurfaceScalars) ==
        maSurfaceScalars.end() ) {

      mView->SetSurfaceScalars( NULL );
      mView->SetSurfaceScalarsColors( NULL );
      mView->SetAnnotationMessage( "" );
      mView->SetSurfaceLookupScalars( NULL );
      mView->SetSurfaceLegendColors( NULL );

      // Render the view with lack of scalars.
      mView->Render();

      return;
    }

    // We need to generate a new set of scalars to display based on
    // our current settings. If we're showing the curvature (and
    // it's loaded), we'll use curvature values where the scalars is
    // < the min. If we have to reverse or not show scalars values,
    // we'll do that here too.
    vtkSmartPointer<vtkFloatArray> composedScalars = 
      vtkSmartPointer<vtkFloatArray>::New();

    vtkFloatArray* surfaceScalars = NULL;
    if( maSurfaceScalars.find(mnCurrentSurfaceScalars) !=
        maSurfaceScalars.end() )
      surfaceScalars = maSurfaceScalars[mnCurrentSurfaceScalars].mValues;

    // Get the number of elements and initialize the composed scalar
    // array.
    int cValues = 0;
    if( mCurvatureScalars.GetPointer() )
      cValues = mCurvatureScalars->GetNumberOfTuples();
    else if( surfaceScalars )
      cValues = surfaceScalars->GetNumberOfTuples();

    composedScalars->Allocate( cValues );
    composedScalars->SetNumberOfComponents( 1 );

    // NOTE: epsilon is actually a bit too small here, as adding an
    // entry to the table and then another entry EPS bigger will
    // sometimes make the table think it's the same value.
    const double EPS = numeric_limits<double>::epsilon() * 2.0;

    // For each value, check the surface scalar value. If it's < min,
    // use the background value. If we're reversing, reverse the
    // scalar value. If we're not showing one side, use the
    // background value. If we are showing curvature (and have it),
    // the background value is our curvature value.
    for( int nValue = 0; nValue < cValues; nValue++ ) {

      float background = 0;
      if( mCurvatureScalars.GetPointer() && mbShowCurvature )
        background = mCurvatureScalars->GetTuple1( nValue );

      float surfaceScalar = 0;
      if( surfaceScalars )
        surfaceScalar = surfaceScalars->GetTuple1( nValue ) -
          mSurfaceScalarsColorOffset; // subtract the offset

      if( mbSurfaceScalarsColorReverse )
        surfaceScalar = -surfaceScalar;
      if( surfaceScalar > 0 && !mbSurfaceScalarsColorShowPositive )
        surfaceScalar = 0;

      if( surfaceScalar < 0 && !mbSurfaceScalarsColorShowNegative )
        surfaceScalar = 0;

      // Insert the appropriate color into the composed array.
      if( surfaceScalar < mSurfaceScalarsColorMin &&
          surfaceScalar > -mSurfaceScalarsColorMin ) {
        // Here, it's possible for the background (curvature) value to
        // be higher than the min scalar value. If so, we'll get a
        // background color with a scalar color, and that will be
        // bad. So clamp it to the -min/min. NOTE this should be able
        // to be adjusted with EPS but it was never big enough to make
        // a difference, so use 0.1 for now.
        if( background < -mSurfaceScalarsColorMin )
          background = -mSurfaceScalarsColorMin + 0.1;
        if( background > mSurfaceScalarsColorMin )
          background = mSurfaceScalarsColorMin - 0.1;
        composedScalars->InsertNextValue( background );
      }
      else
        composedScalars->InsertNextValue( surfaceScalar );
    }

    // Now we need to create a color table. If there is no curvature,
    // this is easy; it's a red point at min, a red point at mid, and
    // yellow point at max, and blue and auqa points on the negative
    // side. If we have curvature as well, we need to have points
    // going from dark gray at curvatureMin to light gray at
    // -curvatureMin. Everything between min and curvatureMin and -min
    // and -curvatureMin will be gray. However, if our min <
    // curvatureMin, then we need to set the curvature min do abut the
    // min.
    double curvatureMin = 0;
    if( mCurvatureScalars.GetPointer() && mbShowCurvature ) {
      double range[2];
      mCurvatureScalars->GetRange( range );
      double highestAbs = MAX( fabs(range[0]), fabs(range[1]) );
      curvatureMin = highestAbs;
    }

    bool bUseGray = true;
    if( mSurfaceScalarsColorMin <= curvatureMin ) {
      curvatureMin = mSurfaceScalarsColorMin - EPS;
      bUseGray = false;
    }

    // This is the negative scalar range.
    vtkSmartPointer<vtkColorTransferFunction> composedColors =
      vtkSmartPointer<vtkColorTransferFunction>::New();
    composedColors->AddRGBPoint( -mSurfaceScalarsColorMax, 0, 1, 1 );
    composedColors->AddRGBPoint( -mSurfaceScalarsColorMid, 0, 0, 1 );
    composedColors->AddRGBPoint( -mSurfaceScalarsColorMin, 0, 0, 1 );

    // This pads the range from the -scalar min to the -curv min with
    // flat gray.
    if( bUseGray && mSurfaceScalarsColorMin != 0 ) {
      composedColors->AddRGBPoint
        ( -mSurfaceScalarsColorMin + EPS, 0.5, 0.5, 0.5 );
      if( mCurvatureScalars.GetPointer() && mbShowCurvature ){
        composedColors->AddRGBPoint( -curvatureMin - EPS, 0.5, 0.5, 0.5 );
      }
    }

    // This is our curvature color range -- binary gray or red/green.
    if( mCurvatureScalars.GetPointer() &&
        mbShowCurvature &&
        mSurfaceScalarsColorMin != 0 ) {

      // Do red/green if we have our flag set and there are no
      // scalars, otherwise do gray.
      if( mbDrawCurvatureGreenRedIfNoScalars && 
          NULL == surfaceScalars ) {
        composedColors->AddRGBPoint( -curvatureMin,     0.0, 1.0, 0.0 );
        composedColors->AddRGBPoint( -curvatureMin/2.0, 0.0, 0.8, 0.0 );
        composedColors->AddRGBPoint(  curvatureMin/2.0, 0.8, 0.0, 0.0 );
        composedColors->AddRGBPoint(  curvatureMin,     1.0, 0.0, 0.0 );
      } else {
        composedColors->AddRGBPoint( -curvatureMin, 0.6, 0.6, 0.6 );
        composedColors->AddRGBPoint(  0,            0.6, 0.6, 0.6 );
        composedColors->AddRGBPoint(  EPS,          0.4, 0.4, 0.4 );
        composedColors->AddRGBPoint(  curvatureMin, 0.4, 0.4, 0.4 );
      }
    }

    // This pads the range from the curv min to the scalar min.
    if( bUseGray && mSurfaceScalarsColorMin != 0 ) {
      if( mCurvatureScalars.GetPointer() && mbShowCurvature )
        composedColors->AddRGBPoint( curvatureMin + EPS, 0.5, 0.5, 0.5 );
      composedColors->AddRGBPoint
        ( mSurfaceScalarsColorMin - EPS, 0.5, 0.5, 0.5 );
    }

    // Positive scalar range.
    composedColors->AddRGBPoint( mSurfaceScalarsColorMin, 1, 0, 0 );
    composedColors->AddRGBPoint( mSurfaceScalarsColorMid, 1, 0, 0 );
    composedColors->AddRGBPoint( mSurfaceScalarsColorMax, 1, 1, 0 );

    composedColors->Build();

    // Set the composed scalars and colors in the view to draw on the
    // surface.
    mView->SetSurfaceScalars( composedScalars );
    mView->SetSurfaceScalarsColors( composedColors );

    // If we have surface scalars...
    if( surfaceScalars ) {

      // Show the label in the annoation.
      mView->SetAnnotationMessage
        ( maSurfaceScalars[mnCurrentSurfaceScalars].msLabel.c_str() );

      // Pass it to the view as the lookup scalars. This way we'll
      // only print surface scalar values and not curvature values when
      // clicking on a point.
      mView->SetSurfaceLookupScalars( surfaceScalars );

      // Make a quick color table with just the scalar colors and set
      // that in the view's legend colors.
      vtkSmartPointer<vtkColorTransferFunction> surfaceScalarsColors =
        vtkSmartPointer<vtkColorTransferFunction>::New();
      surfaceScalarsColors->AddRGBPoint( -mSurfaceScalarsColorMax, 0, 1, 1 );
      surfaceScalarsColors->AddRGBPoint( -mSurfaceScalarsColorMid, 0, 0, 1 );
      surfaceScalarsColors->AddRGBPoint( -mSurfaceScalarsColorMin, 0, 0, 1 );
      surfaceScalarsColors->AddRGBPoint
        ( -mSurfaceScalarsColorMin + EPS, 0.5, 0.5, 0.5 );
      surfaceScalarsColors->AddRGBPoint
        ( mSurfaceScalarsColorMin - EPS, 0.5, 0.5, 0.5 );
      surfaceScalarsColors->AddRGBPoint( mSurfaceScalarsColorMin, 1, 0, 0 );
      surfaceScalarsColors->AddRGBPoint( mSurfaceScalarsColorMid, 1, 0, 0 );
      surfaceScalarsColors->AddRGBPoint( mSurfaceScalarsColorMax, 1, 1, 0 );
      surfaceScalarsColors->Build();
      mView->SetSurfaceLegendColors( surfaceScalarsColors );

    } else {

      mView->SetAnnotationMessage( "" );
      mView->SetSurfaceLookupScalars( NULL );
    }

    // Render the view with the new scalars.
    mView->Render();
  }
}

void
vtkKWQdecWindow::UserSelectedVertexCallback ( vtkObject* iCaller, 
                                              unsigned long iEventId,
                                              void* iClientData,
                                              void* iCallData ) {

  assert( QdecEvents::UserSelectedVertex == iEventId );
  assert( iClientData );
  assert( iCallData );

  // Extract the client and call data and call the window's
  // UserSelectedVertex function.
  try {
    vtkKWQdecWindow* window = 
      static_cast<vtkKWQdecWindow*>( iClientData );
    
    vtkKWQdecView::SurfaceVertexInformation* info = 
      static_cast<vtkKWQdecView::SurfaceVertexInformation*>( iCallData );
    
    if( window )
      window->UserSelectedVertex( *info );
  }
  catch(...) {
    cerr << "Invalid client data in UserTransformChanged callback" << endl;
  }
}

void
vtkKWQdecWindow::UserSelectedVertex 
( vtkKWQdecView::SurfaceVertexInformation& iInfo ) {

  assert( mMenuSaveGDFPostscript );

  // If the vertex that was selected is valid and our GDF is loaded
  // and the window is showing, update the status of our Save GDF menu
  // command.
  if( iInfo.mnVertexIndex >= 0 &&
      mbGDFLoaded && 
      this->GetApplication()->
      EvaluateBooleanExpression( "FsgdfPlot_IsWindowShowing %d", mGDFID )) {
    mMenuSaveGDFPostscript->SetStateToNormal();
  } else {
    mMenuSaveGDFPostscript->SetStateToDisabled();
  }
}

void
vtkKWQdecWindow::UpdateCommandStatus () {

  if( mPanel.GetPointer() ) {
    vtkKWWidget* panelSubjects = mPanel->GetPageWidget( "Subjects" );
    if( panelSubjects )
      if( this->mQdecProject->GetDataTable() ) {
        panelSubjects->EnabledOn();
        panelSubjects->UpdateEnableState();
      } else {
        panelSubjects->EnabledOff();
        panelSubjects->UpdateEnableState();
      }

    vtkKWWidget* panelDesign = mPanel->GetPageWidget( "Design" );
    if( panelDesign )
      if( this->mQdecProject->GetDataTable() ) {
        panelDesign->EnabledOn();
        panelDesign->UpdateEnableState();
      } else {
        panelDesign->EnabledOff();
        panelDesign->UpdateEnableState();
      }

    vtkKWWidget* panelDisplay = mPanel->GetPageWidget( "Display" );
    if( panelDisplay )
      if( maSurfaceSource.size() > 0 ) {
        panelDisplay->EnabledOn();
        panelDisplay->UpdateEnableState();
      } else {
        panelDisplay->EnabledOff();
        panelDisplay->UpdateEnableState();
      }
  }

  assert( mMenuSaveProjectFile );
  if( mQdecProject && mQdecProject->GetGlmFitResults() )
    mMenuSaveProjectFile->SetStateToNormal();
  else
    mMenuSaveProjectFile->SetStateToDisabled();

  assert( mMenuSaveTIFF );
  assert( mBtnSaveTIFF.GetPointer() );
  if( maSurfaceSource.size() > 0 ) {
    mBtnSaveTIFF->SetStateToNormal();
    mMenuSaveTIFF->SetStateToNormal();
  } else {
    mBtnSaveTIFF->SetStateToDisabled();
    mMenuSaveTIFF->SetStateToDisabled();
  }

  assert( mMenuSaveGDFPostscript );
  if( mbGDFLoaded && 
      this->GetApplication()->EvaluateBooleanExpression
      ( "FsgdfPlot_IsWindowShowing %d", mGDFID )) {
    mMenuSaveGDFPostscript->SetStateToNormal();
  } else {
    mMenuSaveGDFPostscript->SetStateToDisabled();
  }

  assert( mMenuClearCurvature );

  if( mCurvatureScalars.GetPointer() ) {
    mMenuClearCurvature->SetStateToNormal();
    mMenuSmoothCurvatureScalars->SetStateToNormal();
  } else {
    mMenuClearCurvature->SetStateToDisabled();
    mMenuSmoothCurvatureScalars->SetStateToDisabled();
  }

  assert( mMenuClearSurfaceScalars );

  if( maSurfaceScalars.size() > 0 ) {
    mMenuClearSurfaceScalars->SetStateToNormal();
  } else {
    mMenuClearSurfaceScalars->SetStateToDisabled();
  }

  if( maSurfaceScalars.size() > 0 &&
      maSurfaceScalars.find(mnCurrentSurfaceScalars) !=
      maSurfaceScalars.end() &&
      NULL != maSurfaceScalars[mnCurrentSurfaceScalars].mValues ) {
    mMenuSmoothSurfaceScalars->SetStateToNormal();
  } else {
    mMenuSmoothSurfaceScalars->SetStateToDisabled();
  }

  assert( mBtnShowCursor.GetPointer() );
  assert( mMenuShowCursor );
  assert( mEntrySelectVertex.GetPointer() );

  if( maSurfaceSource.size() > 0 ) {
    mBtnShowCursor->SetStateToNormal();
    mBtnShowCursor->SetSelectedState( mView->GetShowCursor() );
    mMenuShowCursor->SetStateToNormal();
    mMenuShowCursor->SetSelectedState( mView->GetShowCursor() );
    mEntrySelectVertex->SetStateToNormal();
    mEntrySelectVertex->SetCommand( this, "SelectSurfaceVertex" );
  } else {
    mBtnShowCursor->SetStateToDisabled();
    mMenuShowCursor->SetStateToDisabled();
    mEntrySelectVertex->SetStateToDisabled();
    mEntrySelectVertex->SetCommand( NULL, NULL );
  }

  assert( mBtnSaveLabel.GetPointer() );
  assert( mMenuSaveLabel );
  assert( mBtnRemoveSelectionFromROI.GetPointer() );
  assert( mBtnAddSelectionToROI.GetPointer() );
  assert( mMenuRemoveSelectionFromROI );
  assert( mMenuAddSelectionToROI );
  assert( mMenuClearROI );

  if( mROISource.GetPointer() ) {
    mBtnSaveLabel->SetStateToNormal();
    mMenuSaveLabel->SetStateToNormal();
    mBtnRemoveSelectionFromROI->SetStateToNormal();
    mMenuRemoveSelectionFromROI->SetStateToNormal();
    mMenuClearROI->SetStateToNormal();
  } else {
    mBtnSaveLabel->SetStateToDisabled();
    mMenuSaveLabel->SetStateToDisabled();
    mBtnRemoveSelectionFromROI->SetStateToDisabled();
    mMenuRemoveSelectionFromROI->SetStateToDisabled();
    mMenuClearROI->SetStateToDisabled();
  }
  if( maSurfaceSource.size() > 0 ) {
    mBtnAddSelectionToROI->SetStateToNormal();
    mMenuAddSelectionToROI->SetStateToNormal();
  } else {
    mBtnAddSelectionToROI->SetStateToDisabled();
    mMenuAddSelectionToROI->SetStateToDisabled();
  }

  assert( mMenuMapLabel );

  if( mROISource.GetPointer() && 
      maSurfaceSource.size() > 0 &&
      mQdecProject && mQdecProject->GetSubjectIDs().size() > 0 ) {
    mMenuMapLabel->SetStateToNormal();
  } else {
    mMenuMapLabel->SetStateToDisabled();
  }

  if( mbGDFLoaded && mROISource.GetPointer() ) {
    mMenuGraphAverageROI->SetStateToNormal();
  } else {
    mMenuGraphAverageROI->SetStateToDisabled();
  }
}

void
vtkKWQdecWindow::AnalyzeDesign () {

  assert( mListDiscreteFactors.GetPointer() );
  assert( mListContinuousFactors.GetPointer() );
  assert( mMenuSmoothness.GetPointer() );
  assert( mMenuMeasure.GetPointer() );
  assert( mMenuHemisphere.GetPointer() );
  assert( mEntryDesignName.GetPointer() );
  assert( mEntryAverageSubject.GetPointer() );

  try {

    // in case the user changed the average subject and the design name,
    // then update them
    assert( this->mQdecProject );
    this->mQdecProject->SetAverageSubject( mEntryAverageSubject->GetValue() );
    string fnWorkingDir = this->mQdecProject->GetDefaultWorkingDir();
    fnWorkingDir += "/";
    fnWorkingDir += mEntryDesignName->GetValue();
    this->mQdecProject->SetWorkingDir( fnWorkingDir.c_str() );

    // Gather-up our selected design parameters, from user menu selections
    const char* name = strdup( mEntryDesignName->GetValue() );
    const char* measure = strdup( mMenuMeasure->GetValue() );
    int smoothness = atoi( mMenuSmoothness->GetValue() );
    const char* hemi = strdup( mMenuHemisphere->GetValue() );
    const char sNone[10] = "none";
    const char* df1;
    if( -1 != maDiscreteFactorSelection[0] ) {
      df1 =
        strdup( mListDiscreteFactors->GetItem(maDiscreteFactorSelection[0]) );
    } else {
      df1 = strdup( sNone );
    }
    const char* df2;
    if( -1 != maDiscreteFactorSelection[1] ) {
      df2 =
        strdup( mListDiscreteFactors->GetItem(maDiscreteFactorSelection[1]) );
    } else {
      df2 = strdup( sNone );
    }
    const char* cf1;
    if( -1 != maContinuousFactorSelection[0] ) {
      cf1 = strdup( mListContinuousFactors->GetItem
                    (maContinuousFactorSelection[0]) );
    } else {
      cf1 = strdup( sNone );
    }
    const char* cf2;
    if( -1 != maContinuousFactorSelection[1] ) {
      cf2 = strdup( mListContinuousFactors->GetItem
                    (maContinuousFactorSelection[1]) );
    } else {
      cf2 = strdup( sNone );
    }

    // Now create the design input files (checking for validity before running)
    if( this->mQdecProject->CreateGlmDesign
        (name,
         df1, // first discrete factor
         df2, // second discrete factor
         cf1, // first continuous factor
         cf2, // second continuous factor
         measure,
         hemi,
         smoothness,
         this  /* ProgressUpdateGUI */ )
      ) {
      this->SetStatusText( "Error during design creation" );
    } else if( this->mQdecProject->RunGlmFit() ) { // Now actually run the fit
      this->SetStatusText( "Error during design analysis" );
    } else {
      this->SetStatusText( "Design analysis complete." );
      this->GetProgressGauge()->SetValue( 100 );

      // Clear the surface scalars and curvature data we have now.
      this->ClearCurvature();
      this->ClearSurfaceScalars();

      //
      // Load the analyzed data into the view.
      //
      this->LoadAnalyzedData( this->mQdecProject->GetGlmFitResults() );

      // Clear the ROI if there is one.
      this->ClearROI();
      
      // Raise the display page so they can see the data.
      mPanel->RaisePage( ksDisplayPanelName );
      // (Why doesn't it do this automatically?)
      this->NotebookPageRaised( ksDisplayPanelName ); 
    }
  } catch (exception& e) {
    stringstream ssError;
    ssError << "Error in Analyze: " << e.what();
    this->GetApplication()->ErrorMessage( ssError.str().c_str() );
    this->SetStatusText( "Error during design analysis" );
  }
}

void
vtkKWQdecWindow::SetSubjectsDir ( const char* isSubjectsDir ) {

  if( this->mQdecProject )
  {
    this->mQdecProject->SetSubjectsDir( isSubjectsDir );
    string sSubjectsDir = this->mQdecProject->GetSubjectsDir();
    this->mEntrySubjectsDir->SetValue( sSubjectsDir.c_str() );
  }
}

void
vtkKWQdecWindow::SetAverageSubject ( const char* isAverageSubject ) {

  if( this->mQdecProject )
  {
    this->mQdecProject->SetAverageSubject( isAverageSubject );
    string sAverageSubject = this->mQdecProject->GetAverageSubject();
    this->mEntryAverageSubject->SetValue( sAverageSubject.c_str() );
  }
}

void
vtkKWQdecWindow::SetDesignName ( const char* isDesignName ) {

  this->mEntryDesignName->SetValue( isDesignName );
}

void
vtkKWQdecWindow::SetShowCurvature ( int ibShow ) {

  mbShowCurvature = ibShow;
  this->ComposeSurfaceScalarsAndShow();
}

void
vtkKWQdecWindow::SetDrawCurvatureGreenRed ( int ibDraw ) {

  mbDrawCurvatureGreenRedIfNoScalars = ibDraw;
  this->ComposeSurfaceScalarsAndShow();
}

void
vtkKWQdecWindow::SetSurfaceScalarsColorMin ( const char* isMin ) {

  double value;
  stringstream ssMin( isMin );
  ssMin >> value;

  this->SetSurfaceScalarsColorMin( value );
}

void
vtkKWQdecWindow::SetSurfaceScalarsColorMid ( const char* isMid ) {

  double value;
  stringstream ssMid( isMid );
  ssMid >> value;

  this->SetSurfaceScalarsColorMid( value );
}

void
vtkKWQdecWindow::SetSurfaceScalarsColorMax ( const char* isMax ) {

  double value;
  stringstream ssMax( isMax );
  ssMax >> value;

  this->SetSurfaceScalarsColorMax( value );
}

void
vtkKWQdecWindow::SetSurfaceScalarsColorReverse ( int ibReverse ) {

  mbSurfaceScalarsColorReverse = ibReverse;
  this->ComposeSurfaceScalarsAndShow();
}

void
vtkKWQdecWindow::SetSurfaceScalarsColorShowPositive ( int ibShow ) {

  mbSurfaceScalarsColorShowPositive = ibShow;
  this->ComposeSurfaceScalarsAndShow();
}

void
vtkKWQdecWindow::SetSurfaceScalarsColorShowNegative ( int ibShow ) {

  mbSurfaceScalarsColorShowNegative = ibShow;
  this->ComposeSurfaceScalarsAndShow();
}


vtkKWPushButton*
vtkKWQdecWindow::MakeToolbarButton ( vtkKWToolbar* iToolbar,
                                     const char* isBalloonHelpText,
                                     vtkObject* iCommandObject,
                                     const char* isCommand,
                                     const char* isIconKey ) {

  if( NULL == iToolbar )
    throw runtime_error( "MakeToolbarButton: iToolbar was NULL" );

  // Create the button in the toolbar.
  vtkKWPushButton* button = vtkKWPushButton::New();
  button->SetParent( iToolbar->GetFrame() );
  button->Create();

  // Set balloon text if available.
  if( isBalloonHelpText )
    button->SetBalloonHelpString( isBalloonHelpText );

  // Set the command if available.
  if( iCommandObject && isCommand  )
    button->SetCommand( iCommandObject, isCommand );

  // Try to load the icon. If it fails, use the KW built-in question
  // mark icon.
  if( isIconKey ) {
    try {
      IconLoader::SetPushButtonIcon( isIconKey, button );
    } catch (exception& e) {
      stringstream ssError;
      ssError << "Error loading icon: " << e.what();
      this->GetApplication()->ErrorMessage( ssError.str().c_str() );
      button->SetImageToPredefinedIcon( vtkKWIcon::IconQuestion );
    }
  }

  // Add it to the toolbar.
  iToolbar->AddWidget( button );

  // Return the button.
  return button;
}

vtkKWCheckButton*
vtkKWQdecWindow::MakeToolbarCheckButton ( vtkKWToolbar* iToolbar,
                                          const char* isBalloonHelpText,
                                          vtkObject* iCommandObject,
                                          const char* isCommand,
                                          const char* isIconKey ) {

  if( NULL == iToolbar )
    throw runtime_error( "MakeCheckToolbarButton: iToolbar was NULL" );

  // Create the button in the toolbar.
  vtkKWCheckButton* button = vtkKWCheckButton::New();
  button->SetParent( iToolbar->GetFrame() );
  button->Create();

  // Turn the little box off.
  button->IndicatorVisibilityOff();

  // Set balloon text if available.
  if( isBalloonHelpText )
    button->SetBalloonHelpString( isBalloonHelpText );

  // Set the command if available.
  if( iCommandObject && isCommand )
    button->SetCommand( iCommandObject, isCommand );

  // Try to load the icon. If it fails, use the KW built-in question
  // mark icon.
  if( isIconKey ) {
    try {
      IconLoader::SetCheckButtonIcon( isIconKey, button );
    } catch (exception& e) {
      stringstream ssError;
      ssError << "Error loading icon: " << e.what();
      this->GetApplication()->ErrorMessage( ssError.str().c_str() );
      button->SetImageToPredefinedIcon( vtkKWIcon::IconQuestion );
    }
  }

  // Add it to the toolbar.
  iToolbar->AddWidget( button );

  // Return the button.
  return button;
}

void
vtkKWQdecWindow::AddSpacerToToolbar ( vtkKWToolbar* iToolbar, int iWidth ) {

  // Create a frame of the specified width inside the toolbar.
  vtkSmartPointer<vtkKWFrame> spacer = 
    vtkSmartPointer<vtkKWFrame>::New();
  spacer->SetParent( iToolbar );
  spacer->Create();
  spacer->SetWidth( iWidth );
  iToolbar->AddWidget( spacer );

}

vtkKWQdecWindow::MenuItem::MenuItem ()  :
  mMenu( NULL ),
  mnItem( -1 ) {}


void
vtkKWQdecWindow::MenuItem::MakeCommand ( vtkKWMenu* iMenu,
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
vtkKWQdecWindow::MenuItem::MakeCheckButton ( vtkKWMenu* iMenu,
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
vtkKWQdecWindow::MenuItem::SetStateToDisabled () {

  if( NULL == mMenu.GetPointer() )
    throw runtime_error( "MenuItem::SetStateToDisabled: mMenu was NULL" );
  if( -1 == mnItem )
    throw runtime_error( "MenuItem::SetStateToDisabled: mnItem was -1" );

  mMenu->SetItemStateToDisabled( mnItem );
}

void
vtkKWQdecWindow::MenuItem::SetStateToNormal () {

  if( NULL == mMenu.GetPointer() )
    throw runtime_error( "MenuItem::SetStateToNormal: mMenu was NULL" );
  if( -1 == mnItem )
    throw runtime_error( "MenuItem::SetStateToNormal: mnItem was -1" );

  mMenu->SetItemStateToNormal( mnItem );
}

void
vtkKWQdecWindow::MenuItem::SetSelectedState ( int ibOn ) {

  if( NULL == mMenu.GetPointer() )
    throw runtime_error( "MenuItem::SetSelectedState: mMenu was NULL" );
  if( -1 == mnItem )
    throw runtime_error( "MenuItem::SetSelectedState: mnItem was -1" );

  mMenu->SetItemSelectedState( mnItem, ibOn );
}

int
vtkKWQdecWindow::MenuItem::GetSelectedState () {

  if( NULL == mMenu.GetPointer() )
    throw runtime_error( "MenuItem::GetSelectedState: mMenu was NULL" );
  if( -1 == mnItem )
    throw runtime_error( "MenuItem::GetSelectedState: mnItem was -1" );

  return mMenu->GetItemSelectedState( mnItem );
}

int
vtkKWQdecWindow::MenuItem::GetIndex () const {

  return mnItem;
}

