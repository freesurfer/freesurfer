/**
 * @brief Main QDEC logic
 *
 * Loads in all types of data and manages them. Main logic for running
 * analysis and handling results. Composes display objects and sends
 * it to the main view. Holds settings for color scales and currently
 * displayed objects.
 */
/*
 * Original Author: Kevin Teich
 *
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
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
#include <fstream>
#include <vector>
#include <stdexcept>
#include <sstream>

#include <stdio.h>
#include <stdlib.h>

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

#include "mrisutils.h"
#include "fsenv.h"

using namespace std;

vtkStandardNewMacro( vtkKWQdecWindow );

const char* vtkKWQdecWindow::ksSubjectsPanelName = "Subjects";
const char* vtkKWQdecWindow::ksDesignPanelName = "Design";
const char* vtkKWQdecWindow::ksContrastPanelName = "Contrast";
const char* vtkKWQdecWindow::ksDisplayPanelName = "Display";

const int vtkKWQdecWindow::kMaxDiscreteFactors = 2;
const int vtkKWQdecWindow::kMaxContinuousFactors = 1;

// if USE_MID is 1, then enable usage of the 'mid' threshold value (entry
// of it, and inclusion in the histogram box).  in practice, this 'mid'
// value is useless, and the code typically automatically sets it to some
// epsilon greater than the 'min' value.
#define USE_MID 0 

vtkKWQdecWindow::vtkKWQdecWindow () :
  vtkKWWindow(),
  mCurrentNotebookPanelName ( NULL ),
  mbUseHistogramEditor( true ),
  mMenuLoadDataTable( NULL ),
  mMenuLoadProjectFile( NULL ),
  mMenuLoadLabel( NULL ),
  mMenuLoadAnnotation( NULL ),
  mMenuSaveDataTable( NULL ),
  mMenuSaveProjectFile( NULL ),
  mMenuSaveScatterPlotPostscript( NULL ),
  mMenuSaveTIFF( NULL ),
  mMenuQuickSnapsTIFF( NULL ),
  mMenuSaveGDFPostscript( NULL ),
  mMenuSaveLabel( NULL ),
  mMenuMapLabel( NULL ),
  mMenuClearCurvature( NULL ),
  mMenuClearSurfaceScalars( NULL ),
  mMenuRestoreView( NULL ),
  mMenuZoomOut( NULL ),
  mMenuZoomIn( NULL ),
  mMenuShowCursor( NULL ),
  mMenuShowCurvature( NULL ),
  mMenuAddSelectionToROI( NULL ),
  mMenuRemoveSelectionFromROI( NULL ),
  mMenuClearROI( NULL ),
  mMenuSmoothCurvatureScalars( NULL ),
  mMenuSmoothSurfaceScalars( NULL ),
  mMenuGraphAverageROI( NULL ),
  mMenuButtonSimulationThresh( NULL ),
  mMenuButtonSimulationSign( NULL ),
  mbFrameSurfaceInited( false ),
  mbFrameOverlayInited( false ),
  mbFrameSurfaceScalarsInited( false ),
  mScatterPlotSelection( -1 ),
  mScatterPlotLegend( "" ),
  mVertexPlot( NULL ),
  mQdecProject( NULL ),
  mcVertices( -1 ),
  msCurrentSurfaceSource( "" ),
  mnCurrentSurfaceScalars( -1 ),
  maAnnotationIndicies( NULL ),
  mAnnotationTable( NULL ),
  mbShowCurvature( true ),
  mSurfaceScalarsColorMin( 2.0 ),
#if USE_MID
  mSurfaceScalarsColorMid( 2.0001 ),
#endif
  mSurfaceScalarsColorMax( 5.0 ),
  mSurfaceScalarsColorOffset( 0.0 ),
  mbSurfaceScalarsColorReverse( false ),
  mbSurfaceScalarsColorShowPositive( true ),
  mbSurfaceScalarsColorShowNegative( true ),
  mnNegativeMaxRedValue( 0.0 ),
  mnNegativeMaxGreenValue( 1.0 ),
  mnNegativeMaxBlueValue( 1.0 ),
#if USE_MID
  mnNegativeMidRedValue( 0.0 ),
  mnNegativeMidGreenValue( 0.0 ),
  mnNegativeMidBlueValue( 1.0 ),
#endif
  mnNegativeMinRedValue( 0.0 ),
  mnNegativeMinGreenValue( 0.0 ),
  mnNegativeMinBlueValue( 1.0 ),
  mnPositiveMinRedValue( 1.0 ),
  mnPositiveMinGreenValue( 0.0 ),
  mnPositiveMinBlueValue( 0.0 ),
#if USE_MID
  mnPositiveMidRedValue( 1.0 ),
  mnPositiveMidGreenValue( 0.0 ),
  mnPositiveMidBlueValue( 0.0 ),
#endif
  mnPositiveMaxRedValue( 1.0 ),
  mnPositiveMaxGreenValue( 1.0 ),
  mnPositiveMaxBlueValue( 0.0 ),
  mClusterStats( NULL ),
  mnClusters( 0 ),
  mCurrentCluster( 0 ),
  mbDrawCurvatureGreenRedIfNoScalars( false ),
  mSurfaceScalarsColorsFDRRate( 0.05 ),
  msOverlayDescription( "" ) {

  maDiscreteFactorSelection[0] = -1;
  maDiscreteFactorSelection[1] = -1;
  maContinuousFactorSelection[0] = -1;
  maContinuousFactorSelection[1] = -1;
}

vtkKWQdecWindow::~vtkKWQdecWindow () {
  
  delete mMenuLoadDataTable;
  delete mMenuLoadProjectFile;
  delete mMenuLoadLabel;
  delete mMenuLoadAnnotation;
  delete mMenuSaveProjectFile;
  delete mMenuSaveScatterPlotPostscript;
  delete mMenuSaveTIFF;
  delete mMenuQuickSnapsTIFF;
  delete mMenuSaveGDFPostscript;
  delete mMenuSaveLabel;
  delete mMenuMapLabel;
  delete mMenuClearCurvature;
  delete mMenuClearSurfaceScalars;
  delete mMenuRestoreView;
  delete mMenuZoomOut;
  delete mMenuZoomIn;
  delete mMenuShowCursor;
  delete mMenuShowCurvature;
  delete mMenuAddSelectionToROI;
  delete mMenuRemoveSelectionFromROI;
  delete mMenuClearROI;
  delete mMenuSmoothCurvatureScalars;
  delete mMenuSmoothSurfaceScalars;
  delete mMenuGraphAverageROI;
  delete mQdecProject;
  delete mStatsImportDataTable;
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

  // Show curvature
  mBtnShowCurvature.TakeReference(
    this->MakeToolbarCheckButton( toolbar, "Show curvature",
                                  this, "ShowCurvature", "ShowCurvature" )
    );

  this->AddSpacerToToolbar( toolbar );
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

  this->AddSpacerToToolbar( toolbar );
  this->AddSpacerToToolbar( toolbar );

  // Goto vertex entry.
  vtkSmartPointer<vtkKWEntryWithLabel> labeledEntryVertex =
    vtkSmartPointer<vtkKWEntryWithLabel>::New();
  labeledEntryVertex->SetParent( toolbar->GetFrame() );
  labeledEntryVertex->Create();
  labeledEntryVertex->SetLabelText( "Goto vertex number: " );
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

  // Load Annotation
  mMenuLoadAnnotation = new MenuItem();
  mMenuLoadAnnotation->
    MakeCommand( this->GetFileMenu(), nItem++,
                 "Load &Annotation...", this, "LoadAnnotationFromDlog",
                 "Ctrl+A", "LoadAnnotation" );

  // These menu items are for loading pieces of data individually, and
  // can be enabled for debugging.
#if 1
  this->GetFileMenu()->InsertCommand( nItem++, "Load Surface...",
                                      this, "LoadSurfaceFromDlog" );
//  this->GetFileMenu()->InsertCommand( nItem++, "Load GDF...",
//                                      this, "LoadGDFFromDlog" );
  this->GetFileMenu()->InsertCommand( nItem++, "Load Scalars...",
                                      this, "LoadSurfaceScalarsFromDlog" );
  this->GetFileMenu()->InsertCommand( nItem++, "Load Curvature...",
                                      this, "LoadCurvatureFromDlog" );
#endif
  this->GetFileMenu()->InsertSeparator( nItem++ );

  // Save Data Table.
  mMenuSaveDataTable = new MenuItem();
  mMenuSaveDataTable->
    MakeCommand( this->GetFileMenu(), nItem++,
                 "Save Data Table...", this, "SaveDataTableFromDlog",
                 NULL, NULL );

  // Save Project File.
  mMenuSaveProjectFile = new MenuItem();
  mMenuSaveProjectFile->
    MakeCommand( this->GetFileMenu(), nItem++,
                 "Save Project File...", this, "SaveProjectFileFromDlog",
                 NULL, NULL );

  // Save Postcript of continuous factor plot
  mMenuSaveScatterPlotPostscript = new MenuItem();
  mMenuSaveScatterPlotPostscript->
    MakeCommand( this->GetFileMenu(), nItem++,
                 "Save Scatter Plot as Postcript...", this,
                 "SaveScatterPlotPostscriptFromDlog", NULL, NULL );

  // Save TIFF
  mMenuSaveTIFF = new MenuItem();
  mMenuSaveTIFF->
    MakeCommand( this->GetFileMenu(), nItem++,
                 "Save TIFF...", this, "SaveTIFFImageFromDlog",
                 NULL, "SaveTIFF" );

  // QuickSnaps TIFF
  mMenuQuickSnapsTIFF = new MenuItem();
  mMenuQuickSnapsTIFF->
    MakeCommand( this->GetFileMenu(), nItem++,
                 "QuickSnaps TIFF...", this, "QuickSnapsTIFF",
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

  // Insert a separator,
  this->GetFileMenu()->InsertSeparator( nItem++ );

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

  // Show Curvature
  mMenuShowCurvature = new MenuItem();
  mMenuShowCurvature->
    MakeCheckButton( this->GetViewMenu(), nItem++,
                     "Show Curvature",
                     this, "SetShowCurvatureFromMenu", NULL, "ShowCurvature" );

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

  this->GetMainSplitFrame()->SetFrame1Size( 375 );

  mPanel = vtkSmartPointer<vtkKWUserInterfacePanel>::New();
  mPanel->SetUserInterfaceManager( this->GetMainUserInterfaceManager() );
  mPanel->Create();
  mPanel->AddPage(
    ksSubjectsPanelName, "Explore data from the subjects data table", NULL );
  mPanel->AddPage(
    ksDesignPanelName, "Create the GLM Design matrix", NULL );
#if 0 //HACK disable Contrast tab for now
  mPanel->AddPage(
    ksContrastPanelName, "Create the Contrast matrix, and begin the analysis",
    NULL );
#endif
  mPanel->AddPage(
    ksDisplayPanelName, "Configure the view of the results", NULL );

  //
  // Use the Subjects pane and put some sample form stuff in it
  vtkKWWidget* panelFrame = mPanel->GetPageWidget( ksSubjectsPanelName );

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
  this->Script( "pack %s -fill x", labeledEntry->GetWidgetName() );

  labeledEntry = vtkSmartPointer<vtkKWEntryWithLabel>::New();
  labeledEntry->SetParent( panelFrame );
  labeledEntry->SetLabelText( "Data Table:" );
  labeledEntry->Create();
  labeledEntry->SetLabelPositionToTop();
  mEntryDataTable = labeledEntry->GetWidget();
  mEntryDataTable->SetValue( "<not loaded>" );
  mEntryDataTable->SetReadOnly( 1 );
  this->Script( "pack %s -fill x", labeledEntry->GetWidgetName() );

  labeledEntry = vtkSmartPointer<vtkKWEntryWithLabel>::New();
  labeledEntry->SetParent( panelFrame );
  labeledEntry->SetLabelText( "Number of Subjects:" );
  labeledEntry->Create();
  mEntryNumberOfSubjects = labeledEntry->GetWidget();
  mEntryNumberOfSubjects->SetValue( "0" );
  mEntryNumberOfSubjects->SetReadOnly( 1 );
  this->Script( "pack %s -fill x", labeledEntry->GetWidgetName() );

  // Entry box for alternate fsaverage subject
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
  this->Script( "pack %s -fill x", labeledEntry->GetWidgetName() );

  // check .Qdecrc file for an alternate average subject
  const char* sAvgSubj = QdecUtilities::GetQdecrcResourceString("AVERAGE");
  if( sAvgSubj && this->mQdecProject )
  {
    mEntryAverageSubject->SetValue( sAvgSubj );
  }

  // Create the subject data scatter-plot exploration frame.
  vtkSmartPointer<vtkKWFrameWithLabel> exploreFrame = 
    vtkSmartPointer<vtkKWFrameWithLabel>::New();
  exploreFrame->SetParent( panelFrame );
  exploreFrame->Create();
  exploreFrame->SetLabelText( "Data Table View" );
  this->Script( "pack %s -fill x", exploreFrame->GetWidgetName() );

  vtkSmartPointer<vtkKWListBoxWithScrollbarsWithLabel> listBox =
    vtkSmartPointer<vtkKWListBoxWithScrollbarsWithLabel>::New();
  listBox->SetParent( exploreFrame->GetFrame() );
  listBox->SetLabelText( "Discrete and Continuous Factors:" );
  listBox->Create();
  listBox->SetLabelPositionToTop();
  mListScatterPlot = listBox->GetWidget()->GetWidget();
  mListScatterPlot->SetHeight( 5 );
  mListScatterPlot->ExportSelectionOff();
  this->Script( "pack %s -side top -fill x -expand yes",
                listBox->GetWidgetName() );
  mListScatterPlot->
    SetSelectionCommand( this, "ScatterPlotListBoxCallback" );
  mListScatterPlot->SetSelectionModeToSingle();

  vtkSmartPointer<vtkKWLabel> label = vtkSmartPointer<vtkKWLabel>::New();
  label->SetParent( exploreFrame->GetFrame() );
  label->Create();
  mLabelScatterPlotLegend = label;
  mLabelScatterPlotLegend->SetText( mScatterPlotLegend.c_str() );
  mLabelScatterPlotLegend->SetJustificationToRight();
  this->Script( "pack %s -side top -fill x -expand yes",
                label->GetWidgetName() );

  // button to remove factor from data table
  mBtnFactorRemove =  vtkSmartPointer<vtkKWPushButton>::New();
  mBtnFactorRemove->SetParent( exploreFrame->GetFrame() );
  mBtnFactorRemove->SetText( "Remove Selection from Data Table" );
  mBtnFactorRemove->Create();
  mBtnFactorRemove->SetCommand( this, "RemoveFactorFromDataTable" );
  this->Script( "pack %s", mBtnFactorRemove->GetWidgetName() );


  //
  // Create the stats data import frame.
  //
  mFrameStatsImport = vtkSmartPointer<vtkKWFrameWithLabel>::New();
  mFrameStatsImport->SetParent( panelFrame );
  mFrameStatsImport->Create();
  mFrameStatsImport->SetLabelText( "Stats Data Import" );
  this->Script( "pack %s -fill x", mFrameStatsImport->GetWidgetName() );

  // button to generate stats data
  // when the button is pressed, and after it generates the data tables,
  // then the generate button disappears and is replaced by a menu list
  // see GenerateStatsDataTables
  mBtnStatsGenerate =  vtkSmartPointer<vtkKWPushButton>::New();
  mBtnStatsGenerate->SetParent( mFrameStatsImport->GetFrame() );
  mBtnStatsGenerate->SetText( "Generate Stats Data Tables" );
  mBtnStatsGenerate->Create();
  mBtnStatsGenerate->SetCommand( this, "GenerateStatsDataTables" );

  // continue setup (this routine exists because its called if a new 
  // table or project is loaded)
  this->ResetStatsImportFrame();

  // also prepare the button used to add stats to the Data Table
  mBtnStatsAddToDataTable =  vtkSmartPointer<vtkKWPushButton>::New();
  mBtnStatsAddToDataTable->SetParent( mFrameStatsImport->GetFrame() );
  mBtnStatsAddToDataTable->SetText( "Add Selection(s) to Data Table" );
  mBtnStatsAddToDataTable->Create();
  mBtnStatsAddToDataTable->SetCommand( this, "AddStatsToDataTable" );
  // packed in SetStatsImportItem

  
  //
  // Create the subject exclusion frame.
  //
  vtkSmartPointer<vtkKWFrameWithLabel> excludeFrame = 
    vtkSmartPointer<vtkKWFrameWithLabel>::New();
  excludeFrame->SetParent( panelFrame );
  excludeFrame->Create();
  excludeFrame->SetLabelText( "Subject Exclusions" );
  this->Script( "pack %s -fill x", excludeFrame->GetWidgetName() );

  // Entry boxes for excluding subjects
  labeledEntry = vtkSmartPointer<vtkKWEntryWithLabel>::New();
  labeledEntry->SetParent( excludeFrame->GetFrame() );
  labeledEntry->SetLabelText( "Selected Factor:" );
  labeledEntry->Create();
  labeledEntry->SetLabelPositionToLeft();
  mEntryExcludeFactor = labeledEntry->GetWidget();
  mEntryExcludeFactor->SetValue( "<none selected>" );
  mEntryExcludeFactor->SetWidth( 25 );
  mEntryExcludeFactor->SetReadOnly( 1 );
  this->Script( "pack %s -fill x", labeledEntry->GetWidgetName() );

  labeledEntry = vtkSmartPointer<vtkKWEntryWithLabel>::New();
  labeledEntry->SetParent( excludeFrame->GetFrame() );
  labeledEntry->Create();
  labeledEntry->SetLabelText( "Greater Than:" );
  mEntryExcludeSubjectGT = labeledEntry->GetWidget();
  mEntryExcludeSubjectGT->SetRestrictValueToDouble();
  mEntryExcludeSubjectGT->SetValue( "" );
  mEntryExcludeSubjectGT->SetCommandTrigger ( vtkKWEntry::TriggerOnReturnKey );
  mEntryExcludeSubjectGT->RemoveBinding( "<Unmap>", 
                                         mEntryExcludeSubjectGT, 
                                         "ValueCallback" );
  mEntryExcludeSubjectGT->SetCommand ( this, "SetExcludeSubjectGT" );
  this->Script( "pack %s -fill x", labeledEntry->GetWidgetName() );

  labeledEntry = vtkSmartPointer<vtkKWEntryWithLabel>::New();
  labeledEntry->SetParent( excludeFrame->GetFrame() );
  labeledEntry->Create();
  labeledEntry->SetLabelText( "Less Than:" );
  mEntryExcludeSubjectLT = labeledEntry->GetWidget();
  mEntryExcludeSubjectLT->SetRestrictValueToDouble();
  mEntryExcludeSubjectLT->SetValue( "" );
  mEntryExcludeSubjectLT->SetCommandTrigger ( vtkKWEntry::TriggerOnReturnKey );
  mEntryExcludeSubjectLT->RemoveBinding( "<Unmap>", 
                                         mEntryExcludeSubjectLT, 
                                         "ValueCallback" );
  mEntryExcludeSubjectLT->SetCommand ( this, "SetExcludeSubjectLT" );
  this->Script( "pack %s -fill x", labeledEntry->GetWidgetName() );

  labeledEntry = vtkSmartPointer<vtkKWEntryWithLabel>::New();
  labeledEntry->SetParent( excludeFrame->GetFrame() );
  labeledEntry->Create();
  labeledEntry->SetLabelText( "Equal To:" );
  mEntryExcludeSubjectET = labeledEntry->GetWidget();
  mEntryExcludeSubjectET->SetRestrictValueToDouble();
  mEntryExcludeSubjectET->SetValue( "" );
  mEntryExcludeSubjectET->SetCommandTrigger ( vtkKWEntry::TriggerOnReturnKey );
  mEntryExcludeSubjectET->RemoveBinding( "<Unmap>", 
                                         mEntryExcludeSubjectET, 
                                         "ValueCallback" );
  mEntryExcludeSubjectET->SetCommand ( this, "SetExcludeSubjectET" );
  this->Script( "pack %s -fill x", labeledEntry->GetWidgetName() );

  // button to clear all excluded subjects
  vtkSmartPointer<vtkKWPushButton> clrButton =
    vtkSmartPointer<vtkKWPushButton>::New();
  clrButton->SetParent( excludeFrame->GetFrame() );
  clrButton->SetText( "Clear Exclusions" );
  clrButton->Create();
  clrButton->SetCommand( this, "ClearAllExcludedSubjects" );
  this->Script( "pack %s", clrButton->GetWidgetName() );

  // the Subjects panel is now complete


  // ---------------------------------------------------------------------
  //
  // Design pane gets the control for configuring the design matrix input to
  // mri_glmfit.
  panelFrame = mPanel->GetPageWidget( ksDesignPanelName );

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
  // Create the Measures frame, above the Factors frame.
  mFrameMeasures = vtkSmartPointer<vtkKWFrameWithLabel>::New();
  mFrameMeasures->SetParent( panelFrame );
  mFrameMeasures->Create();
  mFrameMeasures->SetLabelText( "Measure (Dependent variable)" );
  this->Script( "pack %s -fill x", mFrameMeasures->GetWidgetName() );

  // radio-buttons selecting measure category
  mRadBtnSetMeasure = vtkSmartPointer<vtkKWRadioButtonSet>::New();
  mRadBtnSetMeasure->SetParent(  mFrameMeasures->GetFrame() );
  mRadBtnSetMeasure->Create();
  mRadBtnSetMeasure->PackHorizontallyOn();
  this->Script( "pack %s -side top -fill both",
                mRadBtnSetMeasure->GetWidgetName() );

  int nButton = 0;
  vtkSmartPointer<vtkKWRadioButton> radBtn;
  // first category is Surface-based measures
  radBtn.TakeReference( mRadBtnSetMeasure->AddWidget( nButton++ ) );
  radBtn->SetText( "Surface-based" );
  radBtn->SelectedStateOn(); // this one is the default, so its on
  string sCmd = string( "SetCurrentMeasure " ) + radBtn->GetText();
  radBtn->SetCommand( this, sCmd.c_str() );
  // second category is Volume-based measures
  radBtn.TakeReference( mRadBtnSetMeasure->AddWidget( nButton++ ) );
  radBtn->SetText( "Volume-based" );
  radBtn->SelectedStateOff();
  sCmd = string( "SetCurrentMeasure " ) + radBtn->GetText();
  radBtn->SetCommand( this, sCmd.c_str() );
  
  //
  // Create the Factors frame.
  vtkSmartPointer<vtkKWFrameWithLabel> factorsFrame =
    vtkSmartPointer<vtkKWFrameWithLabel>::New();
  factorsFrame->SetParent( panelFrame );
  factorsFrame->Create();
  factorsFrame->SetLabelText( "Model Factors (Independent Variables)" );
  this->Script( "pack %s -fill x", factorsFrame->GetWidgetName() );

  // Discrete and continuous list boxes are inside the Factors frame.
  listBox = vtkSmartPointer<vtkKWListBoxWithScrollbarsWithLabel>::New();
  listBox->SetParent( factorsFrame->GetFrame() );
  listBox->SetLabelText( "Discrete (Fixed Factors):" );
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
  listBox->SetLabelText( "Continuous (Covariate):" );
  listBox->Create();
  listBox->SetLabelPositionToTop();
  mListContinuousFactors = listBox->GetWidget()->GetWidget();
  mListContinuousFactors->ExportSelectionOff();
  mListContinuousFactors->SetHeight ( 5 );
  this->Script( "pack %s -fill x -expand y", listBox->GetWidgetName() );
  mListContinuousFactors->SetSelectionModeToMultiple();
  mListContinuousFactors->
    SetSelectionCommand( this, "ContinuousFactorsListBoxCallback" );

  listBox = vtkSmartPointer<vtkKWListBoxWithScrollbarsWithLabel>::New();
  listBox->SetParent( factorsFrame->GetFrame() );
  listBox->SetLabelText( "Nuisance Factors:" );
  listBox->Create();
  listBox->SetLabelPositionToTop();
  mListNuisanceFactors = listBox->GetWidget()->GetWidget();
  mListNuisanceFactors->ExportSelectionOff();
  mListNuisanceFactors->SetHeight ( 5 );
  this->Script( "pack %s -fill x -expand y", listBox->GetWidgetName() );
  mListNuisanceFactors->SetSelectionModeToMultiple();
  mListNuisanceFactors->
    SetSelectionCommand( this, "NuisanceFactorsListBoxCallback" );

  // call the radio-button handler once to set it up
  this->SetCurrentMeasure( "Surface-based" );

#if 0 // HACK disable for now, until DOSS is handled properly downstream
  //
  // Create the Design Matrix Type frame.
  vtkSmartPointer<vtkKWFrameWithLabel> designMatrixTypeFrame =
    vtkSmartPointer<vtkKWFrameWithLabel>::New();
  designMatrixTypeFrame->SetParent( panelFrame );
  designMatrixTypeFrame->Create();
  designMatrixTypeFrame->SetLabelText( "Design Matrix Type" );
  this->Script( "pack %s -fill x", designMatrixTypeFrame->GetWidgetName() );
  
  // radio-buttons selecting type
  vtkSmartPointer<vtkKWRadioButtonSet> radBtnSet = 
    vtkSmartPointer<vtkKWRadioButtonSet>::New();
  radBtnSet->SetParent(  designMatrixTypeFrame->GetFrame() );
  radBtnSet->Create();
  radBtnSet->PackHorizontallyOn();
  this->Script( "pack %s -side top -fill both", radBtnSet->GetWidgetName() );

  // select between DODS and DOSS
  vtkSmartPointer<vtkKWRadioButton> radBtn2;
  radBtn2.TakeReference( radBtnSet->AddWidget( nButton++ ) );
  radBtn2->SetText( "DODS" );
  radBtn2->SelectedStateOn(); // this one is the default, so its on
  this->SetDesignMatrixType( "DODS" ); // needs to know the default
  sCmd = string( "SetDesignMatrixType " ) + radBtn2->GetText();
  radBtn2->SetCommand( this, sCmd.c_str() );
  radBtn2.TakeReference( radBtnSet->AddWidget( nButton++ ) );
  radBtn2->SetText( "DOSS" );
  radBtn2->SelectedStateOff();
  sCmd = string( "SetDesignMatrixType " ) + radBtn2->GetText();
  radBtn2->SetCommand( this, sCmd.c_str() );
#endif

  // show degrees of freedom, inside its own frame
  vtkSmartPointer<vtkKWFrameWithLabel> dofFrame =
    vtkSmartPointer<vtkKWFrameWithLabel>::New();
  dofFrame->SetParent( panelFrame );
  dofFrame->SetLabelText( "DOF" );
  dofFrame->Create();
  this->Script( "pack %s -fill x", dofFrame->GetWidgetName() );
  labeledEntry = vtkSmartPointer<vtkKWEntryWithLabel>::New();
  labeledEntry->SetParent( dofFrame->GetFrame() );
  labeledEntry->SetLabelText( "Degrees Of Freedom:" );
  labeledEntry->Create();
  mEntryDegreesOfFreedom = labeledEntry->GetWidget();
  mEntryDegreesOfFreedom->SetValue( "0" );
  mEntryDegreesOfFreedom->SetReadOnly( 1 );
  this->Script( "pack %s -fill x", labeledEntry->GetWidgetName() );

#if 0 //HACK disable Contrast tab for now
  // ---------------------------------------------------------------------
  //
  // Contrast pane gets the control for configuring the contrast matrix input
  // to mri_glmfit.
  panelFrame = mPanel->GetPageWidget( ksContrastPanelName );

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

  //
  // Create the Nuisance Factors frame.
  factorsFrame = vtkSmartPointer<vtkKWFrameWithLabel>::New();
  factorsFrame->SetParent( panelFrame );
  factorsFrame->Create();
  factorsFrame->SetLabelText( "Nuisance Factors" );
  this->Script( "pack %s -fill x", factorsFrame->GetWidgetName() );

  // List box is inside the Nuisance Factors frame.
  mListNuisanceFactors = 
    vtkSmartPointer<vtkKWListBoxWithScrollbarsWithLabel>::New();
  mListNuisanceFactors->SetParent( factorsFrame->GetFrame() );
  mListNuisanceFactors->Create();
  mListNuisanceFactors->SetLabelPositionToTop();
  mListNuisanceFactors->SetLabelText( "<none selected>" );
  mListNuisanceFactors->SetHeight( 5 );
  mListNuisanceFactors->GetWidget()->GetWidget()->ExportSelectionOff();
  mListNuisanceFactors->GetWidget()->GetWidget()->SetSelectionModeToMultiple();
  mListNuisanceFactors->GetWidget()->GetWidget()->SetSelectionCommand
    ( this, "NuisanceFactorsListBoxCallback" );
  this->Script( "pack %s -fill x -expand y", 
                mListNuisanceFactors->GetWidgetName() );
#endif

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
  panelFrame = mPanel->GetPageWidget( ksDisplayPanelName );


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

  this->Script( "pack %s %s %s -side top -fill x -pady 5",
                mFrameSurface->GetWidgetName(),
                mFrameOverlay->GetWidgetName(),
                mFrameSurfaceScalars->GetWidgetName() );

  // Create the vertex plot window (hidden until a vertex is seleted)
  // Create the vertex plot window object
//  mVertexPlot = new FsgdfPlot( this->GetApplication()->GetMainInterp() );
  mVertexPlot = new FsgdfPlot();
  mView->SetVertexPlot( mVertexPlot );


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
  callback->SetCallback( ScatterPlotGraphMouseoverEnterElementCallback );
  mGraph->AddObserver( vtkKWBltGraph::MouseoverEnterElementEvent, callback );

  callback = vtkSmartPointer<vtkCallbackCommand>::New();
  callback->SetClientData( this );
  callback->SetCallback( ScatterPlotGraphMouseoverExitElementCallback );
  mGraph->AddObserver( vtkKWBltGraph::MouseoverExitElementEvent, callback );

  callback = vtkSmartPointer<vtkCallbackCommand>::New();
  callback->SetClientData( this );
  callback->SetCallback( ScatterPlotGraphContextualMenuOpeningCallback );
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
vtkKWQdecWindow::SetCurrentMeasure( const char* isMeasure ) {

  // remove the existing frame contents
  if( mFrameSurfaceMeasures ) mFrameSurfaceMeasures->Unpack();
  if( mFrameVolumeMeasures ) mFrameVolumeMeasures->Unpack();

  // if Surface-based radio button was clicked...
  if( 0 == strcmp(isMeasure,"Surface-based") )
  {
    // Create the Surface-based frame, inside the Measures frame.
    mFrameSurfaceMeasures = vtkSmartPointer<vtkKWFrame>::New();
    mFrameSurfaceMeasures->SetParent( mFrameMeasures->GetFrame() );
    mFrameSurfaceMeasures->Create();
    this->Script( "pack %s -fill x", mFrameSurfaceMeasures->GetWidgetName() );

    // radio-buttons selecting surface measure category
    mRadBtnSetSurfaceMeasure = vtkSmartPointer<vtkKWRadioButtonSet>::New();
    mRadBtnSetSurfaceMeasure->SetParent(  mFrameSurfaceMeasures );
    mRadBtnSetSurfaceMeasure->Create();
    mRadBtnSetSurfaceMeasure->PackHorizontallyOn();
    this->Script( "pack %s -side top -fill both",
                  mRadBtnSetSurfaceMeasure->GetWidgetName() );

    int nButton = 0;
    vtkSmartPointer<vtkKWRadioButton> radBtn;
    // first category is Morphometric measures
    radBtn.TakeReference( mRadBtnSetSurfaceMeasure->AddWidget( nButton++ ) );
    radBtn->SetText( "Morphometric" );
    radBtn->SelectedStateOn(); // this one is the default, so its on
    string sCmd = string( "SetCurrentSurfaceMeasure " ) + radBtn->GetText();
    radBtn->SetCommand( this, sCmd.c_str() );
    // second category is Functional measures
    radBtn.TakeReference( mRadBtnSetSurfaceMeasure->AddWidget( nButton++ ) );
    radBtn->SetText( "Functional" );
    radBtn->SelectedStateOff();
    sCmd = string( "SetCurrentSurfaceMeasure " ) + radBtn->GetText();
    radBtn->SetCommand( this, sCmd.c_str() );
    // run once to display
    SetCurrentSurfaceMeasure( "Morphometric" );
  
  } // end if "Surface-based", else if Volume-based radio button was clicked:
  else if( 0 == strcmp(isMeasure,"Volume-based") ) {
    //
    // Create the Volume-based frame, inside the Measures frame.
    mFrameVolumeMeasures = vtkSmartPointer<vtkKWFrame>::New();
    mFrameVolumeMeasures->SetParent( mFrameMeasures->GetFrame() );
    mFrameVolumeMeasures->Create();
    this->Script( "pack %s -fill x", mFrameVolumeMeasures->GetWidgetName() );

    // temporary message:
    int nRow=0;
    vtkSmartPointer<vtkKWLabel> label = vtkSmartPointer<vtkKWLabel>::New();
    label->SetParent( mFrameVolumeMeasures );
    label->Create();
    label->SetText( "-----<not yet implemented>-----" );
    label->SetJustificationToCenter();
    this->Script( "grid %s -column 0 -row %d -sticky ne",
                  label->GetWidgetName(), nRow );

  }

  // make sure Design tab has most current info displayed
  this->UpdateDesignPage();
}

void
vtkKWQdecWindow::SetCurrentSurfaceMeasure( const char* isMeasure ) {

  if( mFrameMorphMeasures ) mFrameMorphMeasures->Unpack();
  if( mFrameFunctionalMeasures ) mFrameFunctionalMeasures->Unpack();

  // if Morphometric radio button was clicked...
  if( 0 == strcmp(isMeasure,"Morphometric") )
  {
    // The widgets inside the Surface-based Measures frame. We'll pack
    // these in a grid so we can get the actual menus as wide as
    // possible and aligned nicely.
    mFrameMorphMeasures = vtkSmartPointer<vtkKWFrame>::New();
    mFrameMorphMeasures->SetParent( mFrameSurfaceMeasures );
    mFrameMorphMeasures->Create();
    this->Script( "pack %s -fill x", mFrameMorphMeasures->GetWidgetName() );

    int nRow = 0;
    vtkSmartPointer<vtkKWLabel> label = vtkSmartPointer<vtkKWLabel>::New();
    label->SetParent( mFrameMorphMeasures );
    label->Create();
    label->SetText( "Measure: " );
    label->SetJustificationToRight();
    this->Script( "grid %s -column 0 -row %d -sticky ne",
                  label->GetWidgetName(), nRow );

    if( mMenuMorphMeasure ) mMenuMorphMeasure->Unpack();
    mMenuMorphMeasure = vtkSmartPointer<vtkKWMenuButton>::New();
    mMenuMorphMeasure->SetParent( mFrameMorphMeasures );
    mMenuMorphMeasure->Create();
    mMenuMorphMeasure->GetMenu()->AddRadioButton( "thickness" );
    mMenuMorphMeasure->GetMenu()->AddRadioButton( "area" );
    mMenuMorphMeasure->GetMenu()->AddRadioButton( "area.pial" );
    mMenuMorphMeasure->GetMenu()->AddRadioButton( "volume" );
    mMenuMorphMeasure->GetMenu()->AddRadioButton( "sulc" );
    mMenuMorphMeasure->GetMenu()->AddRadioButton( "curv" );
    mMenuMorphMeasure->GetMenu()->AddRadioButton( "white.K" );
    mMenuMorphMeasure->GetMenu()->AddRadioButton( "white.H" );
//NJS: not meaningful    mMenuMorphMeasure->GetMenu()->AddRadioButton( "jacobian_white" );
    mMenuMorphMeasure->GetMenu()->AddRadioButton( "w-g.pct.mgh" );
    mMenuMorphMeasure->GetMenu()->AddRadioButton( "intensity.deep" );
    mMenuMorphMeasure->GetMenu()->AddRadioButton( "intensity.superficial" );
    mMenuMorphMeasure->GetMenu()->AddRadioButton( "intensity.deep.mgz" );
    mMenuMorphMeasure->GetMenu()->AddRadioButton( "intensity.superficial.mgz");
    mMenuMorphMeasure->SetValue( "thickness" );

    // read .Qdecrc resource file (which can located in either ~ or in
    // $SUBJECTS_DIR/qdec, or both locations), for key = value pairs
    // consisting of MEASURE# = measure, where # is 1 to 100, and 
    // 'measure' is some user specified measure to insert in the menu.
    // user must have created this measure with 
    // recon-all -qcache -measure measure
    for (int i=1; i<=100; i++)
    {
      char cBuf[100];
      sprintf( cBuf, "MEASURE%d", i );
      const char* sUserMeasure = QdecUtilities::GetQdecrcResourceString(cBuf);
      if (sUserMeasure)
      {
        mMenuMorphMeasure->GetMenu()->AddRadioButton( sUserMeasure );
      }
    }

    this->Script( "grid %s -column 1 -row %d -sticky nw",
                  mMenuMorphMeasure->GetWidgetName(), nRow );
    nRow++;

    label = vtkSmartPointer<vtkKWLabel>::New();
    label->SetParent( mFrameMorphMeasures );
    label->Create();
    label->SetText( "Smoothing (FWHM): " );
    label->SetJustificationToRight();
    this->Script( "grid %s -column 0 -row %d -sticky ne",
                  label->GetWidgetName(), nRow );

    if( mMenuMorphSmoothness ) mMenuMorphSmoothness->Unpack();
    mMenuMorphSmoothness = vtkSmartPointer<vtkKWMenuButton>::New();
    mMenuMorphSmoothness->SetParent( mFrameMorphMeasures );
    mMenuMorphSmoothness->Create();
    mMenuMorphSmoothness->GetMenu()->AddRadioButton( "0" );
    mMenuMorphSmoothness->GetMenu()->AddRadioButton( "5" );
    mMenuMorphSmoothness->GetMenu()->AddRadioButton( "10" );
    mMenuMorphSmoothness->GetMenu()->AddRadioButton( "15" );
    mMenuMorphSmoothness->GetMenu()->AddRadioButton( "20" );
    mMenuMorphSmoothness->GetMenu()->AddRadioButton( "25" );
    mMenuMorphSmoothness->SetValue( "10" );
    this->Script( "grid %s -column 1 -row %d -sticky nw",
                  mMenuMorphSmoothness->GetWidgetName(), nRow );
    nRow++;

    label = vtkSmartPointer<vtkKWLabel>::New();
    label->SetParent( mFrameMorphMeasures );
    label->Create();
    label->SetText( "Hemisphere: " );
    label->SetJustificationToRight();
    this->Script( "grid %s -column 0 -row %d -sticky ne",
                  label->GetWidgetName(), nRow );

    if( mMenuMorphHemisphere ) mMenuMorphHemisphere->Unpack();
    mMenuMorphHemisphere = vtkSmartPointer<vtkKWMenuButton>::New();
    mMenuMorphHemisphere->SetParent( mFrameMorphMeasures );
    mMenuMorphHemisphere->Create();
    mMenuMorphHemisphere->GetMenu()->AddRadioButton( "lh" );
    mMenuMorphHemisphere->GetMenu()->AddRadioButton( "rh" );
    mMenuMorphHemisphere->SetValue( "lh" );
    this->Script( "grid %s -column 1 -row %d -sticky nw",
                  mMenuMorphHemisphere->GetWidgetName(), nRow );
    nRow++;

    // Weight our grid properly.
    this->Script( "grid columnconfigure %s 0 -weight 0",
                  mFrameMorphMeasures->GetWidgetName() );
    this->Script( "grid columnconfigure %s 1 -weight 1",
                  mFrameMorphMeasures->GetWidgetName() );
    this->Script( "grid rowconfigure %s %d -weight 1",
                  mFrameMorphMeasures->GetWidgetName(), nRow );
    for( int nRowConfigure = 0; nRowConfigure <= nRow; nRowConfigure++ )
      this->Script( "grid rowconfigure %s %d -pad 4",
                    mFrameSurfaceMeasures->GetWidgetName(), nRowConfigure );

  } // end if Morphometric, else if Functional radio button was clicked...
  else if( 0 == strcmp(isMeasure,"Functional") ) {
    //
    // Create the Functional based frame, inside the Measures frame.
    //if( mFrameFunctionalMeasures ) delete mFrameFunctionalMeasures;
    mFrameFunctionalMeasures = vtkSmartPointer<vtkKWFrame>::New();
    mFrameFunctionalMeasures->SetParent( mFrameSurfaceMeasures );
    mFrameFunctionalMeasures->Create();
    this->Script( "pack %s -fill x", 
                  mFrameFunctionalMeasures->GetWidgetName() );
    // temporary message:
    int nRow=0;
    vtkSmartPointer<vtkKWLabel> label = vtkSmartPointer<vtkKWLabel>::New();
    label->SetParent( mFrameFunctionalMeasures );
    label->Create();
    label->SetText( "-----<not yet implemented>-----" );
    label->SetJustificationToCenter();
    this->Script( "grid %s -column 0 -row %d -sticky ne",
                  label->GetWidgetName(), nRow );

  }

}


void
vtkKWQdecWindow::SetDesignMatrixType( const char* isType ) {

  assert( mQdecProject ) ;
  assert( mQdecProject->GetGlmDesign() );
  QdecGlmDesign* design =  mQdecProject->GetGlmDesign();

  if( 0 == strcmp(isType,"DODS") ) {
    design->SetDesignMatrixType( "dods" ); 
  } else if( 0 == strcmp(isType,"DOSS") ) {
    design->SetDesignMatrixType( "doss" );
  }  

  // make sure Design tab has most current info displayed
  this->UpdateDesignPage();
}


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
    try {
      this->LoadDataTable( fnDataTable.c_str() );
    } catch ( exception& e ) {
      this->GetApplication()->ErrorMessage( e.what() );
    }
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

  // Set the default dir to the fsaverage subject surf dir,
  // because more often than not the fsaverage inflated surf is desired
  string sQdecDir = this->mQdecProject->GetSubjectsDir() + 
    + "/" + this->mQdecProject->GetAverageSubject() + "/surf";
  dialog->SetLastPath( sQdecDir.c_str() );

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

  // Set the default dir to the subjects dir / qdec.
  string sQdecDir = this->mQdecProject->GetSubjectsDir() + "/qdec";
  dialog->SetLastPath( sQdecDir.c_str() );

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

  // Set the default dir to the subjects dir.
  string sQdecDir = this->mQdecProject->GetSubjectsDir();
  dialog->SetLastPath( sQdecDir.c_str() );

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

  // Set the default dir to the subjects dir
  string sQdecDir = this->mQdecProject->GetSubjectsDir();
  dialog->SetLastPath( sQdecDir.c_str() );

  if( dialog->Invoke() ) {
    dialog->SaveLastPathToRegistry( "LoadAnnotation" );
    string fn( dialog->GetFileName() );
    this->LoadAnnotation( fn.c_str() );
  }
}


void
vtkKWQdecWindow::SaveDataTableFromDlog () {

  assert( mQdecProject );

  vtkSmartPointer<vtkKWLoadSaveDialog> dialog =
    vtkSmartPointer<vtkKWLoadSaveDialog>::New();
  dialog->SetApplication( this->GetApplication() );
  dialog->SaveDialogOn();
  dialog->Create();
  dialog->SetTitle( "Save Data Table" );
  dialog->SetFileTypes( "{{Data file} {*.dat}} {{All} {*}}" );
  dialog->SetDefaultExtension( ".dat" );

  // Make the file name the currently loaded data table filename
  string fnDefault = this->mQdecProject->GetDataTable()->GetFileName();
  dialog->SetInitialFileName( fnDefault.c_str() );

  // Set the default dir to the subjects dir / qdec.
  string sQdecDir = this->mQdecProject->GetSubjectsDir() + "/qdec";
  dialog->SetLastPath( sQdecDir.c_str() );

  if( dialog->Invoke() ) {
    string fnProject( dialog->GetFileName() );
    this->mQdecProject->GetDataTable()->Save( fnProject.c_str() );
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
vtkKWQdecWindow::SaveScatterPlotPostscriptFromDlog () {

  vtkSmartPointer<vtkKWLoadSaveDialog> dialog =
    vtkSmartPointer<vtkKWLoadSaveDialog>::New();
  dialog->SetApplication( this->GetApplication() );
  dialog->SaveDialogOn();
  dialog->Create();
  dialog->SetTitle( "Save Scatter Plot" );

  dialog->SetFileTypes( "{{All} {*}}" );
  dialog->RetrieveLastPathFromRegistry( "SaveScatterPlot" );
  if( dialog->Invoke() ) {
    dialog->SaveLastPathToRegistry( "SaveScatterPlot" );
    string fnPostscript( dialog->GetFileName() );
    mGraph->SavePostscript( fnPostscript.c_str() );
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

/*
 * Time-saver 'macro' to create commonly needed snapshots
 */
void vtkKWQdecWindow::QuickSnapsTIFF () {

  // gotta have a scalar selected
  try {
    if( maSurfaceSource.find( msCurrentSurfaceSource ) ==
        maSurfaceSource.end() ||
        maSurfaceScalars.find( mnCurrentSurfaceScalars ) ==
        maSurfaceScalars.end() ) {
      throw runtime_error( "Must select from analysis results first.");
    }
  } catch ( exception& e ) {
    this->GetApplication()->ErrorMessage( e.what() );
    return;
  }

  // pop-up box stating what it will do
  string sPopUpMsg = 
    "This operation will perform the following time-saver steps:\n"
    "\n"
    " - Find the maximum value and move cursor to that vertex\n"
    " - Disable display of green cross-hair\n"
    " - Set display surface to inflated\n"
    " - Restore the camera to lateral view\n"
    " - Save TIFF image using a filename composed of relevant parameters\n"
    " - Rotate 180d to medial view\n"
    " - Save another TIFF image\n"
    " - Set display surface to pial\n"
    " - Restore the camera to lateral view\n"
    " - Save another TIFF image\n"
    " - Rotate 180d to medial view\n"
    " - Save another TIFF image\n"
    "\n"
    "Image files will be written into the directory:\n";
  sPopUpMsg += this->mQdecProject->GetWorkingDir();

  // pop-up box describing what it is going to do
  if ( vtkKWMessageDialog::PopupOkCancel ( this->GetApplication(), 
                                           this, "", sPopUpMsg.c_str() ) ) {
    // Ok was clicked, proceed...

    // - Find the maximum value and move cursor to that vertex
    this->GenerateClusterStats();

    // - Set display surface to inflated
    this->SetCurrentSurface( "inflated" );

    // - Restore the camera to lateral view
    this->RestoreView();

    // - Disable display of green cross-hair
    this->ShowCursor( false );

    // formulate a base filename composed of relevant parameters
    string sBaseFileName = this->mQdecProject->GetWorkingDir() + "/";
    if (maSurfaceScalars[mnCurrentSurfaceScalars].msLabel2.length() > 0) {
      sBaseFileName += maSurfaceScalars[mnCurrentSurfaceScalars].msLabel2;

      // append nuisance factors, if any
      QdecGlmDesign* design =  mQdecProject->GetGlmDesign();
      vector<QdecFactor*> const& lNuisanceFactors = 
        design->GetNuisanceFactors();
      for(unsigned int j=0; j < lNuisanceFactors.size(); j++) {
        if (lNuisanceFactors[j]) {
          string factorName = "_";
          factorName += lNuisanceFactors[j]->GetFactorName();
          sBaseFileName += factorName;
        }
      }
    } else {
      // doesnt have a label name (which would be the case if a scalar file
      // was loaded manually) so give it the filename (better than nothing)
      sBaseFileName = maSurfaceScalars[mnCurrentSurfaceScalars].mfnSource;
    }

    // - Save TIFF image
    string sFileName = sBaseFileName + "-inflated-lateral.tiff";
    this->SaveTIFFImage( sFileName.c_str(), 1 );
    sPopUpMsg = "\nCreated file: " + sFileName + "\n";

    // - Rotate 180d to medial view
    mView->AnimateCameraAzimuthNegative( );
    mView->AnimateCameraAzimuthNegative( );

    // - Disable display of green cross-hair
    this->ShowCursor( false );

    // - Save another TIFF image
    sFileName = sBaseFileName + "-inflated-medial.tiff";
    this->SaveTIFFImage( sFileName.c_str(), 1 );
    sPopUpMsg += "\nCreated file: " + sFileName + "\n";

    // - Set display surface to pial
    this->SetCurrentSurface( "pial" );

    // - Restore the camera to lateral view
    this->RestoreView();

    // - Disable display of green cross-hair
    this->ShowCursor( false );

    // - Save another TIFF image
    sFileName = sBaseFileName + "-pial-lateral.tiff";
    this->SaveTIFFImage( sFileName.c_str(), 1 );
    sPopUpMsg += "\nCreated file: " + sFileName + "\n";

    // - Rotate 180d to medial view
    mView->AnimateCameraAzimuthNegative( );
    mView->AnimateCameraAzimuthNegative( );

    // - Disable display of green cross-hair
    this->ShowCursor( false );

    // - Save another TIFF image
    sFileName = sBaseFileName + "-pial-medial.tiff";
    this->SaveTIFFImage( sFileName.c_str(), 1 );
    sPopUpMsg += "\nCreated file: " + sFileName + "\n";

    sPopUpMsg += "\nQuickSnaps complete!\n";
    vtkKWMessageDialog::PopupMessage ( this->GetApplication(),
                                       this, "", sPopUpMsg.c_str() );
  }
}


void
vtkKWQdecWindow::SaveGDFPostscriptFromDlog () {

  assert( mVertexPlot );
  assert( mVertexPlot->IsLoaded() );
  assert( mVertexPlot->IsWindowShowing() );

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
    string fn = dialog->GetFileName();
    mVertexPlot->SaveToPostscript( fn.c_str() );
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
      else {
        this->mQdecProject->DumpDataTable( stdout );
        cout << "Data table loading completed successfully." << endl;
      }

      // data table may have set SUBJECTS_DIR, so update our mEntry display
      string sSubjectsDir = this->mQdecProject->GetSubjectsDir();
      this->mEntrySubjectsDir->SetValue( sSubjectsDir.c_str() );
      
      // new data table means any prior stats data is invalid
      this->ResetStatsImportFrame();
      if ( mEntryExcludeFactor ) {
        this->mEntryExcludeFactor->SetValue( "<none selected>" );
      }

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

        // new data table means any prior stats data is invalid
        this->ResetStatsImportFrame();

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

      // There are three scalars to load. The first are the "Contrasts"
      // and are fnWorkingDir/<contrast>/sig.mgh. We also need to get the
      // contrast question and make it the label.
      this->SetStatusText( "Loading contrasts..." );
      this->GetProgressGauge()->SetValue( (float)nStep++ * stepIncrement );
      vector<string> fnContrastSigFiles = iGlmResults->GetContrastSigFiles();
      vector<string> sContrastQuestions = iGlmResults->GetContrastQuestions();
      vector<string> sContrastNames = iGlmResults->GetContrastNames();
      assert( fnContrastSigFiles.size() == sContrastQuestions.size() );
      assert( sContrastQuestions.size() == sContrastNames.size() );
      for( unsigned int i=0;
           i < iGlmResults->GetContrastSigFiles().size(); i++) {
        QdecUtilities::AssertFileIsReadable( fnContrastSigFiles[i] );
        int nScalar = this->LoadSurfaceScalars( fnContrastSigFiles[i].c_str(),
                                                sContrastQuestions[i].c_str(),
                                                sContrastNames[i].c_str());
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

      this->SetStatusText( "Analyzed data loaded. Select results for "
                           "display from the Scalars list." );
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
  this->RestoreView();

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

      // determine whether this is the lh or rh, and set our var
      string str (ifnSurface);
      size_t found = str.find( "rh." );
      if (found != string::npos) {
        mMenuMorphHemisphere->SetValue( "rh" );
        mQdecProject->GetGlmDesign()->SetHemi( "rh" );
      } else {
        found = str.find( "lh." );
        if (found != string::npos) {
          mMenuMorphHemisphere->SetValue( "lh" );
          mQdecProject->GetGlmDesign()->SetHemi( "lh" );
        }
      }

      // reset the view to initialized it
      mView->ResetView();

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
      if (mVertexPlot->ReadFile( ifnGDFFile )) return;
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
                                      const char* isLabel2,
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
    if( isLabel2 ) {
      maSurfaceScalars[nEntry].msLabel2 = isLabel2;
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

  // Update our menu and buttons.
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

  // Update our menu and buttons.
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

    // make sure hemi is set
    mQdecProject->GetGlmDesign()->SetHemi( mMenuMorphHemisphere->GetValue() );

    // Tell the project to run this for all our subjects.
    mQdecProject->GenerateMappedLabelForAllSubjects( ifnLabel, this );

    this->SetStatusText( "Labels mapped and saved." );

  } catch (exception& e) {
    this->GetApplication()->ErrorMessage( e.what() );
  }

  // Update our menu and buttons.
  this->UpdateCommandStatus();
}

void
vtkKWQdecWindow::SaveTIFFImage ( const char* ifnTIFF, 
                                 int iMagnificationLevel ) {

  assert( iMagnificationLevel >= 1 );

  // if there is an image on the Display page...
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

      this->SetStatusText( "Display Panel TIFF saved." );

    } catch ( exception& e ) {
      this->GetApplication()->ErrorMessage( e.what() );
    }
  }
}


void
vtkKWQdecWindow::RestoreView () {

  if( mView.GetPointer() ) {
    mView->RestoreView( mMenuMorphHemisphere->GetValue() );
  }

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
      mTableSurfaceScalars->SelectRow( inEntry );

      // Call this to display these scalars.
      this->ComposeSurfaceScalarsAndShow();

      // Histogram 
      this->UpdateSurfaceScalarsColorsEditor();

      // Update our menu.
      this->UpdateCommandStatus();

    } catch ( exception& e ) {
      this->GetApplication()->ErrorMessage( e.what() );
    }
  }

  // reset cluster info
  if( mClusterStats ) free( mClusterStats );
  mClusterStats = NULL;
  mnClusters = 0;
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

  // reset cluster info
  if( mClusterStats ) free( mClusterStats );
  mClusterStats = NULL;
  mnClusters = 0;
}


void
vtkKWQdecWindow::UnloadSurfaceScalars () {

  // unlike ClearSurfaceScalars, this just clears the currenly loaded
  // scalar (blanking the surface)
  mView->SetSurfaceScalars( NULL );
  mView->SetSurfaceScalarsColors( NULL );
  mView->SetAnnotationMessage( "" );
  mView->SetSurfaceLookupScalars( NULL );
  mView->SetSurfaceLegendColors( NULL );

  // No selected scalar.
  mnCurrentSurfaceScalars = -1;
  mTableSurfaceScalars->ClearSelection();

  // Redraw.
  this->ComposeSurfaceScalarsAndShow();

  // Update the colors editor.
  this->UpdateSurfaceScalarsColorsEditor();

  // Update our menu and page.
  this->UpdateDisplayPage();
  this->UpdateCommandStatus();

  // reset cluster info
  if( mClusterStats ) free( mClusterStats );
  mClusterStats = NULL;
  mnClusters = 0;
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

  this->ZoomBy( 1.1 );
}

void
vtkKWQdecWindow::ZoomOut () {

  this->ZoomBy( 0.9090 );
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
vtkKWQdecWindow::ShowCurvature ( int ibShow ) {

  assert( mView.GetPointer() );

  this->mbShowCurvature = ibShow;

  this->UpdateCommandStatus();

  this->ComposeSurfaceScalarsAndShow();
}

void
vtkKWQdecWindow::SetShowCurvatureFromMenu () {

  assert( mMenuShowCurvature );

  this->ShowCurvature( mMenuShowCurvature->GetSelectedState() );
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
vtkKWQdecWindow::ScatterPlotListBoxCallback () {
  assert( mListScatterPlot.GetPointer() );

  mScatterPlotSelection =
    mListScatterPlot->GetSelectionIndex();

  if (mScatterPlotSelection != -1 &&
      mEntryExcludeFactor &&
      mEntryExcludeSubjectGT &&
      mEntryExcludeSubjectLT &&
      mEntryExcludeSubjectET ) {
    // Get the name of the selected factor.
    string sFactor( mListScatterPlot->GetItem(mScatterPlotSelection) );
    // and update the Subject Exclusions box
    mEntryExcludeFactor->SetValue( sFactor.c_str() );
    mEntryExcludeSubjectGT->SetValue( "" );
    mEntryExcludeSubjectLT->SetValue( "" );
    mEntryExcludeSubjectET->SetValue( "" );

    // if the selected factor is discrete, print a legend at the bottom
    QdecFactor* lFactor = 
      this->mQdecProject->GetDataTable()->GetFactor( sFactor.c_str() );
    if( lFactor->IsDiscrete() ) {
      vector< string > lLevelNames = lFactor->GetLevelNames();
      mScatterPlotLegend = "Legend:  ";
      for( unsigned int i=0; i < lLevelNames.size(); i++ ) {
        mScatterPlotLegend += lLevelNames[i];
        mScatterPlotLegend += "=";
        stringstream val;
        val << lFactor->GetContinuousValue(lLevelNames[i].c_str());
        mScatterPlotLegend += val.str().c_str();
        mScatterPlotLegend += "  ";
        mLabelScatterPlotLegend->SetText( mScatterPlotLegend.c_str() );
      }
    } else {
      mScatterPlotLegend = "";
      mLabelScatterPlotLegend->SetText( mScatterPlotLegend.c_str() );
    }
  }
  this->UpdateScatterPlot();
}


void
vtkKWQdecWindow::DiscreteFactorsListBoxCallback () {

  assert( mListDiscreteFactors.GetPointer() );

  this->ManageFactorListBoxSelections( mListDiscreteFactors,
                                       maDiscreteFactorSelection,
                                       kMaxDiscreteFactors);

  // keep the QdecGlmDesign object up-to-date
  assert( mQdecProject ) ;
  assert( mQdecProject->GetGlmDesign() );
  QdecGlmDesign* design =  mQdecProject->GetGlmDesign();
  design->ClearDiscreteFactors();
  for( int i=0; i < kMaxDiscreteFactors; i++ ) {
    if( -1 != maDiscreteFactorSelection[i] ) {
      design->AddDiscreteFactor
        ( strdup( mListDiscreteFactors->GetItem
                  ( maDiscreteFactorSelection[i] ) ) );
    }
  }

  // update degrees of freedom value
  if( mEntryDegreesOfFreedom ) {
    mEntryDegreesOfFreedom->SetValueAsInt( design->GetDegreesOfFreedom() );
  }
}

void
vtkKWQdecWindow::ContinuousFactorsListBoxCallback () {

  assert( mListContinuousFactors.GetPointer() );

  this->ManageFactorListBoxSelections( mListContinuousFactors,
                                       maContinuousFactorSelection,
                                       kMaxContinuousFactors);

  // keep the QdecGlmDesign object up-to-date
  assert( mQdecProject ) ;
  assert( mQdecProject->GetGlmDesign() );
  QdecGlmDesign* design =  mQdecProject->GetGlmDesign();
  design->ClearContinuousFactors();
  for( int i=0; i < kMaxContinuousFactors; i++ ) {
    if( -1 != maContinuousFactorSelection[i] ) {
      design->AddContinuousFactor
        ( strdup( mListContinuousFactors->GetItem
                  ( maContinuousFactorSelection[i] ) ) );
    }
  }

  // update degrees of freedom value
  if( mEntryDegreesOfFreedom ) {
    mEntryDegreesOfFreedom->SetValueAsInt( design->GetDegreesOfFreedom() );
  }
}

void
vtkKWQdecWindow::NuisanceFactorsListBoxCallback () {

  assert( mListNuisanceFactors.GetPointer() );

  // keep the QdecGlmDesign object up-to-date
  assert( mQdecProject ) ;
  assert( mQdecProject->GetGlmDesign() );
  QdecGlmDesign* design =  mQdecProject->GetGlmDesign();
  design->ClearNuisanceFactors();
  for( int nItem = 0; nItem < mListNuisanceFactors->GetNumberOfItems(); 
       nItem++ ) {
    if(mListNuisanceFactors->GetSelectState(nItem)) {
      design->AddNuisanceFactor
        ( strdup( mListNuisanceFactors->GetItem( nItem ) ) );
    }
  }

  // update degrees of freedom value
  if( mEntryDegreesOfFreedom ) {
    mEntryDegreesOfFreedom->SetValueAsInt( design->GetDegreesOfFreedom() );
  }
}


void
vtkKWQdecWindow::ManageFactorListBoxSelections ( vtkKWListBox* iListBox,
                                                 int iaSelections[2],
                                                 int iMaxFactors) {

  // Make sure that we only have two items selected, max. If more than
  // two, unselect the first one.
  int cSelected = 0;
  for( int nItem = 0; nItem < iListBox->GetNumberOfItems(); nItem++ )
    if( iListBox->GetSelectState( nItem ) )
      cSelected++;

  if( cSelected > iMaxFactors ) {

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
  assert( mQdecProject );
  assert( mListScatterPlot.GetPointer() );

  this->UpdateNumberOfSubjects();
  this->mQdecProject->SetAverageSubject( mEntryAverageSubject->GetValue() );
  this->mQdecProject->SetSubjectsDir( mEntrySubjectsDir->GetValue() );

  mListScatterPlot->DeleteAll();

  vector< string > factors =
    this->mQdecProject->GetDataTable()->GetDiscreteFactorNames();
  for(unsigned int i=0; i < factors.size(); i++) {
    mListScatterPlot->Append( factors[i].c_str() );
  }
  factors = this->mQdecProject->GetDataTable()->GetContinuousFactorNames();
  for(unsigned int i=0; i < factors.size(); i++) {
    mListScatterPlot->Append( factors[i].c_str() );
  }
}

void
vtkKWQdecWindow::UpdateNumberOfSubjects () {

  assert( mEntryNumberOfSubjects.GetPointer() );
  assert( mQdecProject );
  assert( mQdecProject->GetDataTable() );
  assert( mQdecProject->GetGlmDesign() );

  QdecDataTable* dTable = this->mQdecProject->GetDataTable();
  QdecGlmDesign* design = this->mQdecProject->GetGlmDesign();

  mEntryDataTable->SetValue( dTable->GetFileName().c_str() );
  mEntryNumberOfSubjects->SetValueAsInt
    ( dTable->GetSubjectIDs().size() - design->GetNumberOfExcludedSubjects() );

}

void
vtkKWQdecWindow::UpdateDesignPage () {

  assert( mListDiscreteFactors.GetPointer() );
  assert( mListContinuousFactors.GetPointer() );
  assert( mListNuisanceFactors.GetPointer() );
  assert( mQdecProject ) ;
  if ( ! mQdecProject->HaveDataTable() ) return;
  assert( mQdecProject->GetDataTable() );
  assert( mQdecProject->GetGlmDesign() );

  // Fill our our entries with the values from the design.
  QdecGlmDesign* design =  mQdecProject->GetGlmDesign();
  mEntryDesignName->SetValue( design->GetName().c_str() );
  if( mFrameMorphMeasures ) {
    mMenuMorphMeasure->SetValue( design->GetMeasure().c_str() );
    mMenuMorphHemisphere->SetValue( design->GetHemi().c_str() );
    stringstream ssSmoothness;
    ssSmoothness << design->GetSmoothness();
    mMenuMorphSmoothness->SetValue( ssSmoothness.str().c_str() );
  }

  // Clear the factor lists.
  mListDiscreteFactors->DeleteAll();
  mListContinuousFactors->DeleteAll();
  mListNuisanceFactors->DeleteAll();

  // Get the current factors selected for this design.
  vector<QdecFactor*> const& lDiscreteFactors =
    design->GetDiscreteFactors();
  vector<QdecFactor*> const& lContinuousFactors = 
    design->GetContinuousFactors();
  vector<QdecFactor*> const& lNuisanceFactors = 
    design->GetNuisanceFactors();
  int nDiscreteSelection = 0;
  int nContinuousSelection = 0;

  vector< string > factors =
    this->mQdecProject->GetDataTable()->GetDiscreteFactorNames();

  for(unsigned int i=0; i < factors.size(); i++) {
    mListDiscreteFactors->Append( factors[i].c_str() );

    // If this factor is one of our chosen discrete factors from the
    // design, enter its index in our selections, and select this
    // item.
    for(unsigned int j=0; j < lDiscreteFactors.size(); j++) {
      if (lDiscreteFactors[j]) {
        string factorName = lDiscreteFactors[j]->GetFactorName();
        if (factorName == factors[i]) {
          maDiscreteFactorSelection[nDiscreteSelection++] = i;
          mListDiscreteFactors->SetSelectState( i, 1 );
        }
      }
    }
  }

  factors = this->mQdecProject->GetDataTable()->GetContinuousFactorNames();
  for(unsigned int i=0; i < factors.size(); i++) {
    mListContinuousFactors->Append( factors[i].c_str() );
    
    // If this factor is one of our chosen continuous factors from the
    // design, enter its index in our selections, and select this
    // item.
    for(unsigned int j=0; j < lContinuousFactors.size(); j++) {
      if (lContinuousFactors[j]) {
        string factorName = lContinuousFactors[j]->GetFactorName();
        if (factorName == factors[i]) {
          maContinuousFactorSelection[nContinuousSelection++] = i;
          mListContinuousFactors->SetSelectState( i, 1 );
        }
      }
    }
  }

  factors = this->mQdecProject->GetDataTable()->GetContinuousFactorNames();
  for(unsigned int i=0; i < factors.size(); i++) {
    mListNuisanceFactors->Append( factors[i].c_str() );

    // If this factor is one of our chosen nuisance factors from the
    // design, select this item.
    for(unsigned int j=0; j < lNuisanceFactors.size(); j++) {
      if (lNuisanceFactors[j]) {
        string factorName = lNuisanceFactors[j]->GetFactorName();
        if (factorName == factors[i]) {
          mListNuisanceFactors->SetSelectState( i, 1 );
        }
      }
    }
  }

  // update degrees of freedom value
  if( mEntryDegreesOfFreedom ) {
    mEntryDegreesOfFreedom->SetValueAsInt( design->GetDegreesOfFreedom() );
  }
}


#if 0 // HACK
void
vtkKWQdecWindow::UpdateContrastPage () {

  assert( mListNuisanceFactors.GetPointer() );
  assert( mQdecProject ) ;
  assert( mQdecProject->GetDataTable() );
  assert( mQdecProject->GetGlmDesign() );

  QdecGlmDesign* design =  mQdecProject->GetGlmDesign();

   // Get the current factors selected for this design.
  mListNuisanceFactors->GetWidget()->GetWidget()->DeleteAll();
  vector<QdecFactor*> const& lDiscreteFactors =
    design->GetDiscreteFactors();
  for(unsigned int j=0; j < lDiscreteFactors.size(); j++) {
    mListNuisanceFactors->GetWidget()->GetWidget()->Append
      ( lDiscreteFactors[j]->GetFactorName().c_str() );
  }
  vector<QdecFactor*> const& lContinuousFactors = 
    design->GetContinuousFactors();
  for(unsigned int j=0; j < lContinuousFactors.size(); j++) {
    mListNuisanceFactors->GetWidget()->GetWidget()->Append
      ( lContinuousFactors[j]->GetFactorName().c_str() );
  }
}
#endif

void
vtkKWQdecWindow::UpdateDisplayPage () {

  assert( mFrameSurface.GetPointer() );
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
    scalarsFrame->SetLabelText( "Analysis Results" );
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

    // Now a frame for the clear button and the reverse value checkbox
    vtkSmartPointer<vtkKWFrame> frameClearBtn =
      vtkSmartPointer<vtkKWFrame>::New();
    frameClearBtn->SetParent( scalarsFrame->GetFrame() );
    frameClearBtn->Create();

    // Button to clear scalars (ie. unload currently loaded scalar)
    vtkSmartPointer<vtkKWPushButton> buttonClearScalars = 
      vtkSmartPointer<vtkKWPushButton>::New();
    buttonClearScalars->SetParent( frameClearBtn );
    buttonClearScalars->Create();
    buttonClearScalars->SetText( "Clear" );
    buttonClearScalars->SetCommand( this, "UnloadSurfaceScalars" );

    // Checkbox for reverse values.
    vtkSmartPointer<vtkKWCheckButtonWithLabel> labeledCheckReverse =
      vtkSmartPointer<vtkKWCheckButtonWithLabel>::New();
    labeledCheckReverse->SetParent( frameClearBtn );
    labeledCheckReverse->Create();
    labeledCheckReverse->SetLabelText( "Reverse Values" );
    labeledCheckReverse->SetLabelPositionToRight();

    this->Script( "pack %s -side left -expand y -fill x -padx 5",
                  buttonClearScalars->GetWidgetName() );
    this->Script( "pack %s -side left -fill x -padx 5",
                  labeledCheckReverse->GetWidgetName() );
    this->Script( "pack %s -side top -fill x -pady 2",
                  frameClearBtn->GetWidgetName() );

    mCheckSurfaceScalarsColorReverse = labeledCheckReverse->GetWidget();
    mCheckSurfaceScalarsColorReverse->
      SetSelectedState( mbSurfaceScalarsColorReverse );
    mCheckSurfaceScalarsColorReverse->
      SetCommand( this, "SetSurfaceScalarsColorReverse" );

    // Checkboxes for positive and negative values (inside an inner
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
      mEditorSurfaceScalarColors->SetCanvasHeight( 100 );
      mEditorSurfaceScalarColors->SetLabelText("Threshold");
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

    // First do three entries with our min/mid/max/offset values. These can
    // be used to set the values directly, and will also show the
    // values in the editor.
    vtkSmartPointer<vtkKWEntryWithLabel> labeledEntryMin =
      vtkSmartPointer<vtkKWEntryWithLabel>::New();
    labeledEntryMin->SetParent( frameEntries );
    labeledEntryMin->SetLabelText( "Min: " );
    labeledEntryMin->Create();
    mEntrySurfaceScalarsColorMin = labeledEntryMin->GetWidget();
    mEntrySurfaceScalarsColorMin->SetWidth( 5 );
    mEntrySurfaceScalarsColorMin->SetRestrictValueToDouble();
    mEntrySurfaceScalarsColorMin->
      SetCommand( this, "SetSurfaceScalarsColorMin");
    mEntrySurfaceScalarsColorMin->SetValueAsDouble( mSurfaceScalarsColorMin );

#if USE_MID
    vtkSmartPointer<vtkKWEntryWithLabel> labeledEntryMid =
      vtkSmartPointer<vtkKWEntryWithLabel>::New();
    labeledEntryMid->SetParent( frameEntries );
    labeledEntryMid->SetLabelText( "Mid: " );
    labeledEntryMid->Create();
    mEntrySurfaceScalarsColorMid = labeledEntryMid->GetWidget();
    mEntrySurfaceScalarsColorMid->SetWidth( 5 );
    mEntrySurfaceScalarsColorMid->SetRestrictValueToDouble();
    mEntrySurfaceScalarsColorMid->
      SetCommand( this, "SetSurfaceScalarsColorMid");
    mEntrySurfaceScalarsColorMid->SetValueAsDouble( mSurfaceScalarsColorMid );
#endif

    vtkSmartPointer<vtkKWEntryWithLabel> labeledEntryMax = 
      vtkSmartPointer<vtkKWEntryWithLabel>::New();
    labeledEntryMax->SetParent( frameEntries );
    labeledEntryMax->SetLabelText( "Max: " );
    labeledEntryMax->Create();
    mEntrySurfaceScalarsColorMax = labeledEntryMax->GetWidget();
    mEntrySurfaceScalarsColorMax->SetWidth( 5 );
    mEntrySurfaceScalarsColorMax->SetRestrictValueToDouble();
    mEntrySurfaceScalarsColorMax->
      SetCommand( this, "SetSurfaceScalarsColorMax");
    mEntrySurfaceScalarsColorMax->SetValueAsDouble( mSurfaceScalarsColorMax );

    vtkSmartPointer<vtkKWEntryWithLabel> labeledEntryOffset =
      vtkSmartPointer<vtkKWEntryWithLabel>::New();
    labeledEntryOffset->SetParent( frameEntries );
    labeledEntryOffset->SetLabelText( "Offset: " );
    labeledEntryOffset->Create();
    mEntrySurfaceScalarsColorOffset = labeledEntryOffset->GetWidget();
    mEntrySurfaceScalarsColorOffset->SetWidth( 4 );
    mEntrySurfaceScalarsColorOffset->SetRestrictValueToDouble();
    mEntrySurfaceScalarsColorOffset->
      SetCommand( this, "SetSurfaceScalarsColorOffset");
    mEntrySurfaceScalarsColorOffset->
      SetValueAsDouble( mSurfaceScalarsColorOffset );

#if USE_MID
    this->Script( "pack %s %s %s %s -side left -padx 5 -fill x",
                  labeledEntryMin->GetWidgetName(),
                  labeledEntryMid->GetWidgetName(),
                  labeledEntryMax->GetWidgetName(),
                  labeledEntryOffset->GetWidgetName());
#else
    this->Script( "pack %s %s %s -side left -padx 5 -fill x",
                  labeledEntryMin->GetWidgetName(),
                  labeledEntryMax->GetWidgetName(),
                  labeledEntryOffset->GetWidgetName());
#endif
    this->Script( "pack %s -side top -fill x",
                  frameEntries->GetWidgetName() );


    // Buttons to generate stats on clusters and goto next cluster
    vtkSmartPointer<vtkKWFrameWithLabel> frameClusterBtns =
      vtkSmartPointer<vtkKWFrameWithLabel>::New();
    frameClusterBtns->SetParent( editorUserFrame );
    frameClusterBtns->SetLabelText("Clusters");
    frameClusterBtns->Create();
    this->Script( "pack %s -side top -fill x",
                  frameClusterBtns->GetWidgetName() );
    vtkSmartPointer<vtkKWFrame> subframeClusterBtns = 
      vtkSmartPointer<vtkKWFrame>::New();
    subframeClusterBtns->SetParent( frameClusterBtns->GetFrame() );
    subframeClusterBtns->Create();
    this->Script( "pack %s -side top -anchor center", 
                  subframeClusterBtns->GetWidgetName() );

    vtkSmartPointer<vtkKWPushButton> buttonGenerateClusterStats = 
      vtkSmartPointer<vtkKWPushButton>::New();
    buttonGenerateClusterStats->SetParent( subframeClusterBtns );
    buttonGenerateClusterStats->Create();
    buttonGenerateClusterStats->SetText("Find Clusters and Goto Max");
    buttonGenerateClusterStats->SetCommand( this, "GenerateClusterStats" );

    vtkSmartPointer<vtkKWPushButton> buttonNextCluster = 
      vtkSmartPointer<vtkKWPushButton>::New();
    buttonNextCluster->SetParent( subframeClusterBtns );
    buttonNextCluster->Create();
    buttonNextCluster->SetText("Next");
    buttonNextCluster->SetCommand( this, "GotoNextCluster" );

    vtkSmartPointer<vtkKWPushButton> buttonPrevCluster = 
      vtkSmartPointer<vtkKWPushButton>::New();
    buttonPrevCluster->SetParent( subframeClusterBtns );
    buttonPrevCluster->Create();
    buttonPrevCluster->SetText("Prev");
    buttonPrevCluster->SetCommand( this, "GotoPrevCluster" );

    this->Script( "pack %s %s %s -side left",
                  buttonGenerateClusterStats->GetWidgetName(),
                  buttonNextCluster->GetWidgetName(),
                  buttonPrevCluster->GetWidgetName() );

    // Correction for Multiple-Comparisons stuff, inside its own frame
    vtkSmartPointer<vtkKWFrameWithLabel> cmcFrame =
      vtkSmartPointer<vtkKWFrameWithLabel>::New();
    cmcFrame->SetParent(  editorUserFrame );
    cmcFrame->SetLabelText( "Correction for Multiple Comparisons" );
    cmcFrame->Create();

    // Now a frame for the button to set the color scale using FDR,
    // and its rate entry.
    vtkSmartPointer<vtkKWFrameWithLabel> frameFDR =
      vtkSmartPointer<vtkKWFrameWithLabel>::New();
    frameFDR->SetParent( cmcFrame->GetFrame() );
    frameFDR->SetLabelText( "False Discovery Rate" );
    frameFDR->Create();
    this->Script( "pack %s -side top -fill x",
                  frameFDR->GetWidgetName() );

    // sub-frame to allow clean arrangement of the stuff
    vtkSmartPointer<vtkKWFrame> subframeFDR = 
      vtkSmartPointer<vtkKWFrame>::New();
    subframeFDR->SetParent( frameFDR->GetFrame() );
    subframeFDR->Create();
    this->Script( "pack %s -side top -anchor center", 
                  subframeFDR->GetWidgetName() );

    vtkSmartPointer<vtkKWPushButton> buttonFDR = 
      vtkSmartPointer<vtkKWPushButton>::New();
    buttonFDR->SetParent( subframeFDR );
    buttonFDR->Create();
    buttonFDR->SetText( "Set Using FDR" );
    buttonFDR->SetCommand( this, "SetSurfaceScalarsColorsUsingFDR" );
    this->Script( "pack %s -side left", buttonFDR->GetWidgetName() );

    vtkSmartPointer<vtkKWEntryWithLabel> labeledEntryFDR =
      vtkSmartPointer<vtkKWEntryWithLabel>::New();
    labeledEntryFDR->SetParent( subframeFDR );
    labeledEntryFDR->SetLabelText( "Rate: " );
    labeledEntryFDR->Create();
    labeledEntryFDR->GetWidget()->SetWidth( 5 );
    labeledEntryFDR->GetWidget()->SetRestrictValueToDouble();
    labeledEntryFDR->GetWidget()->
      SetCommand( this, "SetSurfaceScalarsColorsFDRRate" );
    labeledEntryFDR->GetWidget()->
      SetValueAsDouble( mSurfaceScalarsColorsFDRRate );
    this->Script( "pack %s -side left", labeledEntryFDR->GetWidgetName() );

    //
    //
    // Now a frame for Monte Carlo Null-Z simulation 
    //
    vtkSmartPointer<vtkKWFrameWithLabel> frameSimulation =
      vtkSmartPointer<vtkKWFrameWithLabel>::New();
    frameSimulation->SetParent( cmcFrame->GetFrame() );
    frameSimulation->SetLabelText( "Monte Carlo Null-Z Simulation" );
    frameSimulation->Create();
    this->Script( "pack %s -side top -fill x",
                  frameSimulation->GetWidgetName() );

    // sub-frames to allow clean arrangement of the stuff
    vtkSmartPointer<vtkKWFrame> frameSim1 = 
      vtkSmartPointer<vtkKWFrame>::New();
    frameSim1->SetParent( frameSimulation->GetFrame() );
    frameSim1->Create();
    vtkSmartPointer<vtkKWFrame> frameSim2 = 
      vtkSmartPointer<vtkKWFrame>::New();
    frameSim2->SetParent( frameSimulation->GetFrame() );
    frameSim2->Create();
    this->Script( "pack %s %s -side top -fill x",
                  frameSim1->GetWidgetName(),
                  frameSim2->GetWidgetName() );

    // simulation threshold menu selection
    vtkSmartPointer<vtkKWLabel> label1 = vtkSmartPointer<vtkKWLabel>::New();
    label1->SetParent( frameSim1 );
    label1->Create();
    label1->SetText( "Threshold: " );
    label1->SetJustificationToRight();
    mMenuButtonSimulationThresh = vtkSmartPointer<vtkKWMenuButton>::New();
    mMenuButtonSimulationThresh->SetParent( frameSim1 );
    mMenuButtonSimulationThresh->Create();
    mMenuButtonSimulationThresh->GetMenu()->AddRadioButton( "1.3 (0.05)" );
    mMenuButtonSimulationThresh->GetMenu()->AddRadioButton( "2.0 (0.01)" );
    mMenuButtonSimulationThresh->GetMenu()->AddRadioButton( "2.3 (0.005)" );
    mMenuButtonSimulationThresh->GetMenu()->AddRadioButton( "3.0 (0.001)" );
    mMenuButtonSimulationThresh->GetMenu()->AddRadioButton( "3.3 (0.0005)" );
    mMenuButtonSimulationThresh->GetMenu()->AddRadioButton( "4.0 (0.0001)" );
    mMenuButtonSimulationThresh->SetValue( "2.0 (0.01)" );

     // simulation sign menu selection
    vtkSmartPointer<vtkKWLabel> label2 = vtkSmartPointer<vtkKWLabel>::New();
    label2->SetParent( frameSim2 );
    label2->Create();
    label2->SetText( "Sign: " );
    label2->SetJustificationToRight();
    mMenuButtonSimulationSign = vtkSmartPointer<vtkKWMenuButton>::New();
    mMenuButtonSimulationSign->SetParent( frameSim2 );
    mMenuButtonSimulationSign->Create();
    mMenuButtonSimulationSign->GetMenu()->AddRadioButton( "abs" );
    mMenuButtonSimulationSign->GetMenu()->AddRadioButton( "pos" );
    mMenuButtonSimulationSign->GetMenu()->AddRadioButton( "neg" );
    mMenuButtonSimulationSign->SetValue( "abs" );
    this->Script( "pack %s %s %s %s -side left -expand y -fill x -padx 5",
                  label1->GetWidgetName(),
                  mMenuButtonSimulationThresh->GetWidgetName(),
                  label2->GetWidgetName(),
                  mMenuButtonSimulationSign->GetWidgetName() );

    // button to run the simulation
    vtkSmartPointer<vtkKWPushButton> buttonSimulation = 
      vtkSmartPointer<vtkKWPushButton>::New();
    buttonSimulation->SetParent( frameSimulation->GetFrame() );
    buttonSimulation->Create();
    buttonSimulation->SetText( "Run" );
    buttonSimulation->SetWidth( 20 );
    buttonSimulation->SetCommand( this, "RunSimulation" );
    this->Script( "pack %s",
                  buttonSimulation->GetWidgetName() );

   // Pack the multiple comparisons items in the user frame.
    this->Script( "pack %s %s -side top -fill x -pady 2",
                  editorUserFrame->GetWidgetName(),
                  cmcFrame->GetWidgetName());

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
#if USE_MID
    mEntrySurfaceScalarsColorMid = NULL;
#endif
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
    colors->AddRGBAPoint( -mSurfaceScalarsColorMax,
                          mnNegativeMaxRedValue,
                          mnNegativeMaxGreenValue,
                          mnNegativeMaxBlueValue,
                          1 );
#if USE_MID
  mnNegativeMidValue =
    colors->AddRGBAPoint( -mSurfaceScalarsColorMid,
                          mnNegativeMidRedValue,
                          mnNegativeMidGreenValue,
                          mnNegativeMidBlueValue,
                          1 );
#endif
  mnNegativeMinValue =
    colors->AddRGBAPoint( -mSurfaceScalarsColorMin,
                          mnNegativeMinRedValue,
                          mnNegativeMinGreenValue,
                          mnNegativeMinBlueValue,
                          1 );
  mnPositiveMinValue =
    colors->AddRGBAPoint(  mSurfaceScalarsColorMin,
                           mnPositiveMinRedValue,
                           mnPositiveMinGreenValue,
                           mnPositiveMinBlueValue,
                           1 );
#if USE_MID
  mnPositiveMidValue =
    colors->AddRGBAPoint(  mSurfaceScalarsColorMid,
                           mnPositiveMidRedValue,
                           mnPositiveMidGreenValue,
                           mnPositiveMidBlueValue,
                           1 );
#endif
  mnPositiveMaxValue =
    colors->AddRGBAPoint(  mSurfaceScalarsColorMax,
                           mnPositiveMaxRedValue,
                           mnPositiveMaxGreenValue,
                           mnPositiveMaxBlueValue,
                           1 );
  colors->Build();

  // Set up the point symmetry in the colors.
  mEditorSurfaceScalarColors->SetPointCountMinimum( 6 );
  mEditorSurfaceScalarColors->SetPointCountMaximum( 6 );
  mEditorSurfaceScalarColors->
    SetPointSymmetry( mnNegativeMaxValue, mnPositiveMaxValue );
#if USE_MID
  mEditorSurfaceScalarColors->
    SetPointSymmetry( mnNegativeMidValue, mnPositiveMidValue );
#endif
  mEditorSurfaceScalarColors->
    SetPointSymmetry( mnNegativeMinValue, mnPositiveMinValue );

  // If we have scalars right now, make a histogram of them, set it in
  // the editor, and make sure the range for the editor is correct.
  if( maSurfaceScalars.find( mnCurrentSurfaceScalars ) !=
      maSurfaceScalars.end() ) {

    vtkFloatArray* surfaceScalars =
      maSurfaceScalars[mnCurrentSurfaceScalars].mValues;

    // saturate at +-100, since the vtkHistogram code has a problem with
    // extremes (like +-10000000000)
    for( int nVertex = 0; 
         nVertex < surfaceScalars->GetNumberOfTuples(); 
         nVertex++ ) {
      float val = surfaceScalars->GetValue(nVertex);
      if (val > 100) surfaceScalars->SetValue(nVertex, 100);
      else if (val < -100) surfaceScalars->SetValue(nVertex, -100);
    }
    
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
  
#if USE_MID
  mEditorSurfaceScalarColors->
    GetFunctionPointParameter( mnPositiveMidValue, &value );
  if( mSurfaceScalarsColorMid != value )
    this->SetSurfaceScalarsColorMid( value );
#endif
  
  mEditorSurfaceScalarColors->
    GetFunctionPointParameter( mnPositiveMaxValue, &value );
  if( mSurfaceScalarsColorMax != value )
    this->SetSurfaceScalarsColorMax( value );

  vtkSmartPointer<vtkRGBATransferFunction> colors = 
    mEditorSurfaceScalarColors->GetRGBATransferFunction();

  mnNegativeMaxRedValue = colors->GetRedValue( -mSurfaceScalarsColorMax );
  mnNegativeMaxGreenValue = colors->GetGreenValue( -mSurfaceScalarsColorMax );
  mnNegativeMaxBlueValue = colors->GetBlueValue( -mSurfaceScalarsColorMax );

#if USE_MID
  mnNegativeMidRedValue = colors->GetRedValue( -mSurfaceScalarsColorMid );
  mnNegativeMidGreenValue = colors->GetGreenValue( -mSurfaceScalarsColorMid );
  mnNegativeMidBlueValue = colors->GetBlueValue( -mSurfaceScalarsColorMid );
#endif

  mnNegativeMinRedValue = colors->GetRedValue( -mSurfaceScalarsColorMin );
  mnNegativeMinGreenValue = colors->GetGreenValue( -mSurfaceScalarsColorMin );
  mnNegativeMinBlueValue = colors->GetBlueValue( -mSurfaceScalarsColorMin );

  mnPositiveMaxRedValue = colors->GetRedValue( mSurfaceScalarsColorMax );
  mnPositiveMaxGreenValue = colors->GetGreenValue( mSurfaceScalarsColorMax );
  mnPositiveMaxBlueValue = colors->GetBlueValue( mSurfaceScalarsColorMax );

#if USE_MID
  mnPositiveMidRedValue = colors->GetRedValue( mSurfaceScalarsColorMid );
  mnPositiveMidGreenValue = colors->GetGreenValue( mSurfaceScalarsColorMid );
  mnPositiveMidBlueValue = colors->GetBlueValue( mSurfaceScalarsColorMid );
#endif

  mnPositiveMinRedValue = colors->GetRedValue( mSurfaceScalarsColorMin );
  mnPositiveMinGreenValue = colors->GetGreenValue( mSurfaceScalarsColorMin );
  mnPositiveMinBlueValue = colors->GetBlueValue( mSurfaceScalarsColorMin );

  // Draw with the new values.
  this->ComposeSurfaceScalarsAndShow();
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

#if USE_MID
    // Update the Mid value as well, by 0.0001, since this typically
    // what you would always do anyway.
    this->SetSurfaceScalarsColorMid( iMin + 0.0001 );
    // note that the prior call will update the editor, and
    // draw with the new values.
#else
    // Update the editor.
    this->UpdateSurfaceScalarsColorsEditor();

    // Draw with the new values.
    this->ComposeSurfaceScalarsAndShow();
#endif
  }
}

void
vtkKWQdecWindow::SetSurfaceScalarsColorMid ( double iMid ) {
#if USE_MID
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
#endif
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
#if USE_MID
  assert( mEntrySurfaceScalarsColorMid.GetPointer() );
#endif
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

#if USE_MID
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
#endif

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
vtkKWQdecWindow::GenerateClusterStats () {

  this->SetStatusText( "Generating stats on clusters..." );
  try {

    if( maSurfaceSource.find( msCurrentSurfaceSource ) ==
        maSurfaceSource.end() ||
        maSurfaceScalars.find( mnCurrentSurfaceScalars ) ==
        maSurfaceScalars.end() ) {
      throw runtime_error( "Must have a surface loaded and results selected.");
    }

    // get white surface struct, and insert current scalars into its val slots
    // but first make sure it exists.
    if( maSurfaceSource.find( "white" ) == maSurfaceSource.end() ) {
      throw runtime_error( 
        "GenerateClusterStats: white surface was not found");
    }
    MRIS* mris = maSurfaceSource["white"]->GetMRIS();
    vtkFloatArray* scalars = maSurfaceScalars[mnCurrentSurfaceScalars].mValues;
    for( int nVertex = 0; nVertex < mris->nvertices; nVertex++ ) {
      mris->vertices[nVertex].val = scalars->GetTuple1( nVertex );
      mris->vertices[nVertex].ripflag = 0;
    }

    // tal transform (not really necessary for fsaverage, since its already
    // in tal space, but not always true for other average subjects)
    MATRIX* XFM = 
      DevolveXFM(mEntryAverageSubject->GetValue(), NULL, "talairach.xfm");

    // call main cluster stats generator...
    mnClusters = 0;
    fprintf(stdout,
            "============================================================\n");
    cout << "Generating cluster stats using min threshold of " << 
      mSurfaceScalarsColorMin << "...\n";
    mClusterStats = 
      sclustMapSurfClusters(mris,
                            mSurfaceScalarsColorMin, // threshold min
                            -1,  // threshold max
                            0,   // threshold sign (0=abs)
                            0,   // min area (ignore)
                            &mnClusters, // return
                            XFM); // tal transform
    if( NULL == mClusterStats )
      throw runtime_error( "Unable to generate cluster stats\n" );
    else
      this->SetStatusText( "Completed generation of cluster stats" );
    cout << "Found " << mnClusters << " clusters" << endl;
    cout << "Contrast: '" 
         << maSurfaceScalars[mnCurrentSurfaceScalars].msLabel2
         << "', " << mQdecProject->GetGlmDesign()->GetSmoothness() 
         << "fwhm, DOF: " 
         << mQdecProject->GetGlmDesign()->GetDegreesOfFreedom() << endl;

    // screen dump cluster info
    fprintf
      (stdout,
       "ClusterNo  Max   VtxMax  Size(mm2)   "
       "TalX   TalY   TalZ NVtxs Annotation\n"
       "---------  ---   ------  ---------   "
       "----   ----   ---- ----- ----------\n");
    for (int n=0; n < mnClusters; n++) {
      fprintf(stdout," %4d  %8.4f  %6d    %6.2f  %6.1f %6.1f %6.1f %4d  %s\n",
              n+1,
              mClusterStats[n].maxval, 
              mClusterStats[n].vtxmaxval,
              mClusterStats[n].area,
              mClusterStats[n].xxfm, 
              mClusterStats[n].yxfm, 
              mClusterStats[n].zxfm,
              mClusterStats[n].nmembers,
              this->GetAnnotationForVertex(mClusterStats[n].vtxmaxval));
    }
    fprintf(stdout,
            "============================================================\n");

    // jump to max vertex
    if (mnClusters > 0) {
      this->SelectSurfaceVertex( mClusterStats[0].vtxmaxval );
    }

  } catch (exception& e) {
    stringstream ssError;
    ssError << "Error generating cluster stats\n" << e.what();
    this->GetApplication()->ErrorMessage( ssError.str().c_str() );
    this->SetStatusText( "Error generating cluster stats" );
  }
}

void
vtkKWQdecWindow::GotoNextCluster () {
  // advance index to next cluster (with wrap)
  if (++mCurrentCluster >= mnClusters) {
    mCurrentCluster = 0;
  }

  this->GotoCluster( mCurrentCluster );
}

void
vtkKWQdecWindow::GotoPrevCluster () {
    // rewind index to previous cluster (with wrap)
    if (--mCurrentCluster < 0) {
      mCurrentCluster = mnClusters-1;
    }

  this->GotoCluster( mCurrentCluster );
}

void
vtkKWQdecWindow::GotoCluster ( int iCurrentCluster ) {

  try {

    if( maSurfaceSource.find( msCurrentSurfaceSource ) ==
        maSurfaceSource.end() ||
        maSurfaceScalars.find( mnCurrentSurfaceScalars ) ==
        maSurfaceScalars.end() ) {
      throw runtime_error( "Must have a surface loaded and scalar selected.");
    }

    if( NULL == mClusterStats ) {
      throw runtime_error( "Must generate cluster stats first.");
    }

    if( 0 == mnClusters ) {
      throw runtime_error( "There are zero clusters.");
    }

    // screen dump cluster info
    fprintf
      (stdout,
       "ClusterNo  Max   VtxMax  Size(mm2)   "
       "TalX   TalY   TalZ NVtxs Annotation\n");
    int n = iCurrentCluster;
    fprintf(stdout," %4d  %8.4f  %6d    %6.2f  %6.1f %6.1f %6.1f %4d  %s\n",
            n+1,
            mClusterStats[n].maxval, 
            mClusterStats[n].vtxmaxval,
            mClusterStats[n].area,
            mClusterStats[n].xxfm, 
            mClusterStats[n].yxfm, 
            mClusterStats[n].zxfm,
            mClusterStats[n].nmembers,
            this->GetAnnotationForVertex(mClusterStats[n].vtxmaxval));

    // jump to max vertex of the current cluster
    if (mnClusters > 0) {
      this->SelectSurfaceVertex( mClusterStats[iCurrentCluster].vtxmaxval );
    }

  } catch (exception& e) {
    stringstream ssError;
    ssError << "Error advancing to next cluster\n" << e.what();
    this->GetApplication()->ErrorMessage( ssError.str().c_str() );
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
      throw runtime_error( "\nMust have a surface loaded and overlay selected.");
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

    // report the number of vertices
    int nVerticesBelowThreshold = 0;
    for( int nVertex = 0; nVertex < mris->nvertices; nVertex++ ) {
      if ( fabs( mris->vertices[nVertex].val ) > threshold ) {
        nVerticesBelowThreshold++; // found a vertex below FDR threshold
      }
    }
    cout << "Found " << nVerticesBelowThreshold << " of " << 
      mris->nvertices <<  " vertices " <<
      "above FDR threshold (of " << threshold << ")\n";

    // Set our min to the threshold, and calculate good values for mid
    // and max.
    this->SetSurfaceScalarsColors( threshold,
                                   threshold + 0.001,
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


// Creates a script to run mri_glmfit/mris_surfcluster for multiple-
// comparisons correction
void
vtkKWQdecWindow::RunSimulation () {

  this->SetStatusText( "Running Monte Carlo Simulation..." );
  try {

    if( maSurfaceSource.find( msCurrentSurfaceSource ) ==
        maSurfaceSource.end() ||
        maSurfaceScalars.find( mnCurrentSurfaceScalars ) ==
        maSurfaceScalars.end() ) {
      throw runtime_error( "\nMust have a surface loaded and "
                           "overlay selected.");
    }

    if( 0 != strcmp( this->mQdecProject->GetAverageSubject().c_str(),
                     "fsaverage") ) {
      throw runtime_error( "\nThis function only works when using fsaverage "
                           "as the common-space subject");
    }

    const char* contrast = 
      maSurfaceScalars[mnCurrentSurfaceScalars].msLabel2.c_str();

    if( 0 == strlen( contrast ) ) {
      throw runtime_error( "Select an Analysis Result corresponding to a "
                           "contrast" );
    }

    string threshold;
    if( 0 == strcmp( mMenuButtonSimulationThresh->GetValue(),
                     "1.3 (0.05)" ) ) {
      threshold = "th13";
    }
    else if( 0 == strcmp( mMenuButtonSimulationThresh->GetValue(),
                          "2.0 (0.01)" ) ) {
      threshold = "th20";
    }
    else if( 0 == strcmp( mMenuButtonSimulationThresh->GetValue(),
                          "2.3 (0.005)" ) ) {
      threshold = "th23";
    }
    else if( 0 == strcmp( mMenuButtonSimulationThresh->GetValue(),
                          "3.0 (0.001)" ) ) {
      threshold = "th30";
    }
    else if( 0 == strcmp( mMenuButtonSimulationThresh->GetValue(),
                          "3.3 (0.0005)" ) ) {
      threshold = "th33";
    }
    else if( 0 == strcmp( mMenuButtonSimulationThresh->GetValue(),
                          "4.0 (0.0001)" ) ) {
      threshold = "th40";
    }

    const char* fnClusterSigFile = NULL;
    if( this->mQdecProject->RunMonteCarloSimulation
        ( threshold.c_str(),
          mMenuButtonSimulationSign->GetValue(),
          contrast,
          &fnClusterSigFile ) ) {
      throw runtime_error( "Error running mri_surfcluster!" );
    }
    int nScalar = this->LoadSurfaceScalars( fnClusterSigFile,
                                            NULL, // will show filename
                                            contrast);
    if( nScalar < 0 ) {
      throw runtime_error( string("Could not load results file ") +
                           fnClusterSigFile );
    }

    this->SetCurrentSurfaceScalars( nScalar );

    this->SetStatusText( "Completed Monte Carlo simulation" );

  } catch (exception& e) {
    stringstream ssError;
    ssError << "Error in Monte Carlo simulation: " << e.what();
    this->GetApplication()->ErrorMessage( ssError.str().c_str() );
    this->SetStatusText( "Error during Monte Carlo simulation" );
  }
}


void
vtkKWQdecWindow::UpdateScatterPlot () {

  assert( mGraph.GetPointer() );
  assert( mQdecProject );

  this->UpdateNumberOfSubjects(); // depends on number of excluded subjects

  if( mScatterPlotSelection != -1 ) {
      
    // Start a list of points for our data to graph.
    vector<double> lPoints;
    
    // Get the name of the factor.
    string sFactor( mListScatterPlot->GetItem(mScatterPlotSelection) );
    
    // For each subject...
    vector<QdecSubject*>::iterator tSubject;
    vector<QdecSubject*> lSubjects = 
      mQdecProject->GetDataTable()->GetSubjects();
    int nIndex = 0;
    for( tSubject = lSubjects.begin();
         tSubject != lSubjects.end();
         ++tSubject, nIndex++ ) {
      
      // Get the value.
      float value = (*tSubject)->GetContinuousFactorValue( sFactor.c_str() );
      
      // dont plot it if its NaN (ie. no data for that subject)
      if (isnan(value)) {
        cout << endl << "WARNING: Excluding subject " 
             << (*tSubject)->GetId().c_str()
             << " from plot due to NaN data point!" << endl;
        mGraph->DeleteElement( (*tSubject)->GetId().c_str() );
        continue;
      }

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

      mGraph->DeleteElement( (*tSubject)->GetId().c_str() );

      // Add the data to the graph.
      mGraph->AddElement( (*tSubject)->GetId().c_str(),
                          lPoints,
                          sSymbol.c_str(),
                          0,
                          color[0], color[1], color[2] );

      // and the axis labels
      mGraph->SetYAxisTitle( sFactor.c_str() ); 
      string sXaxis = "Subject";
      mGraph->SetXAxisTitle( sXaxis.c_str() ); 

      // and un-gray the Save Factor Plot to Postscript menu item
      mMenuSaveScatterPlotPostscript->SetStateToNormal();
    }
  }
  else
  {
    mGraph->DeleteAllElements();
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
    string sLabelName = mMenuMorphHemisphere->GetValue();
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
    string sLabelName = mMenuMorphHemisphere->GetValue();
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
  assert( mVertexPlot->IsLoaded() );

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
  mVertexPlot->BeginPointList();

  // For each vertex, add it to the point list.
  vector<int>::iterator tVertex;
  for( tVertex = lVertices.begin(); tVertex != lVertices.end(); ++tVertex ) {
    int nVertex = *tVertex;
    mVertexPlot->AddPoint( nVertex );
  }

  // End the point list and draw the graph.
  mVertexPlot->EndPointList();

  // Set the info string.
  string info = "\"Average of ROI\"";
  mVertexPlot->SetInfo( info.c_str() );
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
    this->mCurrentNotebookPanelName = ksSubjectsPanelName;
    this->GetViewFrame()->UnpackChildren();
    this->Script( "pack %s -expand yes -fill both -anchor c",
                  mGraph->GetWidgetName() );
    this->SetStatusText( "Subjects selection" );
  } else if ( 0 == strcmp( isTitle, ksDesignPanelName ) ) {
    this->mCurrentNotebookPanelName = ksDesignPanelName;
    this->GetViewFrame()->UnpackChildren();
    this->SetStatusText( "Design matrix creation" );
#if 0 //HACK
  } else if ( 0 == strcmp( isTitle, ksContrastPanelName ) ) {
    this->mCurrentNotebookPanelName = ksContrastPanelName;
    this->GetViewFrame()->UnpackChildren();
    this->SetStatusText( "Contrast matrix creation and analysis launch" );
    this->UpdateContrastPage();
#endif
  } else if ( 0 == strcmp( isTitle, ksDisplayPanelName ) ) {
    this->mCurrentNotebookPanelName = ksDisplayPanelName;
    this->GetViewFrame()->UnpackChildren();
    this->Script( "pack %s -expand yes -fill both -anchor c",
                  mView->GetWidgetName() );
    this->SetStatusText( "Display of results" );
  }

  this->UpdateCommandStatus();
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
  // ScatterPlotGraphMouseoverEnterElement function.
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
vtkKWQdecWindow::ScatterPlotGraphMouseoverEnterElement 
( const char* isElement ) {

  // Just show it in our status bar.
  this->SetStatusText( isElement );
}

void 
vtkKWQdecWindow::ScatterPlotGraphMouseoverExitElement () {

  // Clear our status bar.
  this->SetStatusText( "" );
}

void
vtkKWQdecWindow::ScatterPlotGraphSetUpContextualMenu (const char* isElement,
                                                         vtkKWMenu* iMenu ) {

  assert( isElement );
  assert( iMenu );
  assert( mQdecProject );
  assert( mMenuMorphHemisphere.GetPointer() );

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
    stringstream ssTkmeditCmd;
    stringstream ssTkmeditCommand;
    // catch { exec tkmedit $s norm.mgz $hemi.white -segmentation aseg.mgz }
    ssTkmeditCmd << envVars->FREESURFER_HOME << "/tktools/tkmedit "
                 << isElement << " norm.mgz -aux T1.mgz -surfs " 
                 << "-segmentation aparc+aseg.mgz &";
    ssTkmeditCommand << "set err [catch { exec "
                     << ssTkmeditCmd.str().c_str()
                     << " }] ; if { $err != 0 } { "
                     << this->GetApplication()->GetTclName() << " "
                     << "ErrorMessage \"Exec error: " 
                     << ssTkmeditCmd.str().c_str()
                     << "\" }";
    //cout << ssTkmeditCommand.str() << endl;
    iMenu->AddCommand( "Open in tkmedit ", NULL, 
                       ssTkmeditCommand.str().c_str() );

    // tksurfer, left hemi
    stringstream ssTksurferCmd;
    stringstream ssTksurferCommand;
    // catch { exec tksurfer $s $hemi inflated -annotation aparc.annot }
    ssTksurferCmd << envVars->FREESURFER_HOME << "/tktools/tksurfer "
                  << isElement << " " 
                  << "lh inflated "
                  << "-annotation aparc.annot &";
    ssTksurferCommand << "set err [catch { exec "
                      << ssTksurferCmd.str().c_str()
                      << " }] ; ";
    //cout << ssTksurferCommand.str() << endl;
    iMenu->AddCommand( "Open lh in tksurfer ", NULL, 
                       ssTksurferCommand.str().c_str() );

    // tksurfer, right hemi
    stringstream ssTksurferCmd2;
    stringstream ssTksurferCommand2;
    ssTksurferCmd2 << envVars->FREESURFER_HOME << "/tktools/tksurfer "
                  << isElement << " " 
                  << "rh inflated "
                  << "-annotation aparc.annot &";
    ssTksurferCommand2 << "set err [catch { exec "
                       << ssTksurferCmd2.str().c_str()
                       << " }] ; ";
    //cout << ssTksurferCommand2.str() << endl;
    iMenu->AddCommand( "Open rh in tksurfer ", NULL, 
                       ssTksurferCommand2.str().c_str() );

  } else {

    // Informative disabled menu commands.
    iMenu->AddCommand( "tkmedit not found", NULL, "" );
    iMenu->SetItemStateToDisabled( 5 );
    iMenu->AddCommand( "tksurfer not found", NULL, "" );
    iMenu->SetItemStateToDisabled( 6 );
  }

}

void
vtkKWQdecWindow::ScatterPlotGraphMouseoverEnterElementCallback
( vtkObject* iCaller, unsigned long iEventId, void* iClientData,
  void* iCallData ) {

  assert( vtkKWBltGraph::MouseoverEnterElementEvent == iEventId );
  assert( iCallData );
  assert( iClientData );

  // Extract the client and call data and call the window's
  // ScatterPlotGraphMouseoverEnterElement function.
  try {

    vtkKWBltGraph::SelectedElementAndPoint* foundElement =
      static_cast<vtkKWBltGraph::SelectedElementAndPoint*>( iCallData );

    vtkKWQdecWindow* window = 
      static_cast<vtkKWQdecWindow*>( iClientData );
    
    if( window )
      window->ScatterPlotGraphMouseoverEnterElement( foundElement->msLabel);
  }
  catch(...) {
    cerr << "Invalid call or client data in "
         << "ScatterPlotGraphMouseoverEnterElementCallback" << endl;
  }
}

void
vtkKWQdecWindow::ScatterPlotGraphMouseoverExitElementCallback
( vtkObject* iCaller, unsigned long iEventId, void* iClientData,
  void* iCallData ) {

  assert( vtkKWBltGraph::MouseoverExitElementEvent == iEventId );
  assert( iClientData );

  // Extract the client data and call the window's
  // ScatterPlotGraphMouseoverExitElement function.
  try {

    vtkKWQdecWindow* window = 
      static_cast<vtkKWQdecWindow*>( iClientData );
    
    if( window )
      window->ScatterPlotGraphMouseoverExitElement();
  }
  catch(...) {
    cerr << "Invalid call or client data in "
         << "ScatterPlotGraphMouseoverExitElementCallback" << endl;
  }
}

void
vtkKWQdecWindow::ScatterPlotGraphContextualMenuOpeningCallback
( vtkObject* iCaller, unsigned long iEventId, void* iClientData,
  void* iCallData ) {

  assert( vtkKWBltGraph::ContextualMenuOpening == iEventId );
  assert( iClientData );
  assert( iCallData );

  // Extract the client and call data and call the window's
  // ScatterPlotGraphSetUpContextualMenu function.
  try {

    vtkKWBltGraph::ContextualMenuElement* clickedElement =
      static_cast<vtkKWBltGraph::ContextualMenuElement*>( iCallData );

    vtkKWQdecWindow* window = 
      static_cast<vtkKWQdecWindow*>( iClientData );
    
    if( window )
      window->ScatterPlotGraphSetUpContextualMenu(clickedElement->msElement,
                                                     clickedElement->mMenu);
  }
  catch(...) {
    cerr << "Invalid call or client data in "
         << "ScatterPlotGraphContextualMenuOpeningCallback" << endl;
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
  this->UpdateScatterPlot();
}

void
vtkKWQdecWindow::SetExcludeSubjectGT ( double inExcludeGT ) {

  if( this->mQdecProject )
  {
    string value = this->mEntryExcludeSubjectGT->GetValue();
    if ( value == "") return;

    this->mEntryExcludeSubjectGT->SetValueAsDouble( inExcludeGT );

    // Get the design from the project.
    QdecGlmDesign* design = mQdecProject->GetGlmDesign();
    assert( design );
  
    // Set the exclude flag in the design for all subjects having the
    // continuous factor value greater than iExcludeGT
    design->SetExcludeSubjectsFactorGT( this->mEntryExcludeFactor->GetValue(),
                                        inExcludeGT, 1 );

    // Redraw our graph.
    this->UpdateScatterPlot();
  }
}
void
vtkKWQdecWindow::SetExcludeSubjectGT ( const char* isExcludeGT ) {
  if ( strcmp( isExcludeGT, "" ) == 0 ) return;
}

void
vtkKWQdecWindow::SetExcludeSubjectLT ( double inExcludeLT ) {

  if( this->mQdecProject )
  {
    this->mEntryExcludeSubjectLT->SetValueAsDouble( inExcludeLT );

    // Get the design from the project.
    QdecGlmDesign* design = mQdecProject->GetGlmDesign();
    assert( design );
  
    // Set the exclude flag in the design for all subjects having the
    // continuous factor value less than iExcludeLT
    design->SetExcludeSubjectsFactorLT( this->mEntryExcludeFactor->GetValue(),
                                       inExcludeLT, 1 );

    // Redraw our graph.
    this->UpdateScatterPlot();
  }
}
void
vtkKWQdecWindow::SetExcludeSubjectLT ( const char* isExcludeLT ) {
  if ( strcmp( isExcludeLT, "" ) == 0 ) return;
}


void
vtkKWQdecWindow::SetExcludeSubjectET ( double inExcludeET ) {

  if( this->mQdecProject )
  {
    string value = this->mEntryExcludeSubjectET->GetValue();
    if ( value == "") return;

    this->mEntryExcludeSubjectET->SetValueAsDouble( inExcludeET );

    // Get the design from the project.
    QdecGlmDesign* design = mQdecProject->GetGlmDesign();
    assert( design );
  
    // Set the exclude flag in the design for all subjects having the
    // continuous factor value greater than iExcludeET
    design->SetExcludeSubjectsFactorET( this->mEntryExcludeFactor->GetValue(),
                                        inExcludeET, 1 );

    // Redraw our graph.
    this->UpdateScatterPlot();
  }
}
void
vtkKWQdecWindow::SetExcludeSubjectET ( const char* isExcludeET ) {
  if ( strcmp( isExcludeET, "" ) == 0 ) return;
}


void
vtkKWQdecWindow::ClearAllExcludedSubjects ( ) {

  if( this->mQdecProject && this->mQdecProject->HaveDataTable() )
  {
    this->mEntryExcludeSubjectGT->SetValue( "" );
    this->mEntryExcludeSubjectLT->SetValue( "" );
    this->mEntryExcludeSubjectET->SetValue( "" );

    // Get the design from the project.
    QdecGlmDesign* design = mQdecProject->GetGlmDesign();
    assert( design );
  
    // Clear the exclude flag in the design for all subjects
    design->ClearAllExcludedSubjects( );

    // Redraw our graph.
    this->UpdateScatterPlot();
  }
}


void vtkKWQdecWindow::ResetStatsImportFrame ( ) {

  // start by removing any frames that might already be displayed
  if ( mMenuStatsData ) mMenuStatsData->Unpack();
  if ( mListBoxStatsImport ) mListBoxStatsImport->Unpack();
  if ( mBtnStatsAddToDataTable ) mBtnStatsAddToDataTable->Unpack();

  // menu to hold stat table names, once those are generated
  mMenuStatsData = vtkSmartPointer<vtkKWMenuButton>::New();
  mMenuStatsData->SetParent( mFrameStatsImport->GetFrame() );
  mMenuStatsData->Create();
  mMenuStatsData->GetMenu()->AddRadioButton( "<none selected>",
                                             this,
                                             "SetStatsImportItem none");
  mMenuStatsData->SetValue( "<none selected>" );

  // prepare a new menu list to be filled with the factors from the selected
  // stats data table
  mListBoxStatsImport = 
    vtkSmartPointer<vtkKWListBoxWithScrollbarsWithLabel>::New();
  mListBoxStatsImport->SetParent( mFrameStatsImport->GetFrame() );
  mListBoxStatsImport->SetLabelText( "Continuous Factors:" );
  mListBoxStatsImport->Create();
  mListBoxStatsImport->SetLabelPositionToTop();
  mListStatsImportFactors = mListBoxStatsImport->GetWidget()->GetWidget();
  mListStatsImportFactors->ExportSelectionOff();
  mListStatsImportFactors->SetHeight ( 5 );
  mListStatsImportFactors->SetSelectionModeToMultiple();
  // notice it doesnt get packed (displayed), that is done in 
  // SetStatsImportItem

  // show just the generate button
  this->Script( "pack %s", mBtnStatsGenerate->GetWidgetName() );
  // the routine below, GenerateStatsDataTables, is what is executed when
  // that button is pressed.
}
void
vtkKWQdecWindow::GenerateStatsDataTables ( ) {

  if ( this->mQdecProject ) {

    // make sure the main data table has been loaded
    if ( ! this->mQdecProject->HaveDataTable() ) {
      vtkSmartPointer<vtkKWMessageDialog> dialog =
        vtkSmartPointer<vtkKWMessageDialog>::New();

      dialog->SetStyleToMessage();
      dialog->SetOptions( vtkKWMessageDialog::WarningIcon );
      dialog->SetApplication( this->GetApplication() );
      dialog->Create();
      dialog->SetText( "The main data table must be loaded first." );
      dialog->Invoke();
      return;
    }

    // run asegstats2table and aparcstats2table scripts
    vector< string > statsNames = this->mQdecProject->CreateStatsDataTables();

    if ( statsNames.size() ) { // it sppears to have given us some data
      // remove the generate button...
      if( mBtnStatsGenerate ) mBtnStatsGenerate->Unpack();

      // ...and replace with a menu button containing names of the generated
      // stats tables
      for( unsigned int i=0; i < statsNames.size(); i++ ) {
        string sCmd = string( "SetStatsImportItem " ) + statsNames[i].c_str();
        mMenuStatsData->GetMenu()->AddRadioButton( statsNames[i].c_str(),
                                                   this,
                                                   sCmd.c_str());
      }
      //this->Script( "pack %s -fill x", mFrameStatsImport->GetWidgetName() );
      this->Script( "pack %s", mMenuStatsData->GetWidgetName() );
    }
  }
}
void
vtkKWQdecWindow::SetStatsImportItem ( const char* isStatsImportItem ) {
  // remove prior selections
  mListStatsImportFactors->DeleteAll();

  // handle case of <none selected>
  if ( ! strcmp(isStatsImportItem,"none")) return;

  // load from file
  stringstream ssFname;
  ssFname << this->mQdecProject->GetStatsDataTablesDir() 
          << "/" << isStatsImportItem << ".stats.dat";

  // kluge: aseg.vol.stats.dat file uses word 'volume' as is fsid column name
  stringstream ssFsIdColName;
  if ( !strcmp(isStatsImportItem,"aseg.vol") ) ssFsIdColName << "volume";
  else ssFsIdColName << isStatsImportItem;

  // load the table
  mStatsImportDataTable = new QdecDataTable();
  int ret = mStatsImportDataTable->Load( ssFname.str().c_str(), 
                                         NULL, 
                                         ssFsIdColName.str().c_str() );
  if ( ret ) {
    // some error loading the table, so fall-back
    mMenuStatsData->SetValue( "<none selected>" );
    return;
  }

  // fill our menu list with the factors from this data table
  vector< string > factors = mStatsImportDataTable->GetContinuousFactorNames();
  for(unsigned int i=0; i < factors.size(); i++) {
    mListStatsImportFactors->Append( factors[i].c_str() );
  }
  this->Script( "pack %s -fill x -expand y", 
                mListBoxStatsImport->GetWidgetName() );

  // dont forget to show the button
  this->Script( "pack %s", mBtnStatsAddToDataTable->GetWidgetName() );
}
void vtkKWQdecWindow::AddStatsToDataTable ( ) {

  // this routine is called when the Add Stats To Data Table button is pressed
  int nItemsMerged=0;
  for( int nItem = 0; 
       nItem < this->mListStatsImportFactors->GetNumberOfItems(); 
       nItem++ ) {
    if( this->mListStatsImportFactors->GetSelectState( nItem ) ) {
      // Get the name of the factor.
      string sFactor( this->mListStatsImportFactors->GetItem( nItem ) );

      // and merge it into our Data Table
      this->mQdecProject->MergeFactorIntoDataTable( sFactor.c_str(), 
                                                    mStatsImportDataTable );
      nItemsMerged++;
    }
  }

  // We need to update our tabs.
  if ( nItemsMerged ) this->UpdateDesignPage();
  if ( nItemsMerged ) this->UpdateSubjectsPage();

  // so the user can see what was just added...
  if ( nItemsMerged ) this->mQdecProject->DumpDataTable( stdout );
}
void vtkKWQdecWindow::RemoveFactorFromDataTable ( ) {

  // this routine is called when the Remove Factor from Data Table 
  // button is pressed

  mScatterPlotSelection = mListScatterPlot->GetSelectionIndex();
  if ( -1 == mScatterPlotSelection ) return; // nothing selected

  // Get the name of the selected factor.
  string sFactor( mListScatterPlot->GetItem(mScatterPlotSelection) );

  string sPopUpMsg = "Remove factor '" + sFactor + "'?";
    
  // pop-up an 'are you really really really really really sure' box
  if ( vtkKWMessageDialog::PopupOkCancel ( this->GetApplication(), 
                                           this, "", sPopUpMsg.c_str() ) ) {

    // Ok was clicked, so delete factor from our Data Table
    this->mQdecProject->RemoveFactorFromDataTable( sFactor.c_str() );

    mScatterPlotSelection = -1; // default to empty factor selection

    // We need to update our tabs.
    if ( mEntryExcludeFactor ) {
      this->mEntryExcludeFactor->SetValue( "<none selected>" );
    }
    this->UpdateScatterPlot();
    this->UpdateDesignPage();
    this->UpdateSubjectsPage();

    // so the user can see what was just deleted...
    this->mQdecProject->DumpDataTable( stdout );
  }
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
    composedColors->AddRGBPoint( -mSurfaceScalarsColorMax,
                                 mnNegativeMaxRedValue,
                                 mnNegativeMaxGreenValue,
                                 mnNegativeMaxBlueValue );
#if USE_MID
    composedColors->AddRGBPoint( -mSurfaceScalarsColorMid,
                                 mnNegativeMidRedValue,
                                 mnNegativeMidGreenValue,
                                 mnNegativeMidBlueValue );
#endif
    composedColors->AddRGBPoint( -mSurfaceScalarsColorMin,
                                 mnNegativeMinRedValue,
                                 mnNegativeMinGreenValue,
                                 mnNegativeMinBlueValue );

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
    composedColors->AddRGBPoint( mSurfaceScalarsColorMax,
                                 mnPositiveMaxRedValue,
                                 mnPositiveMaxGreenValue,
                                 mnPositiveMaxBlueValue );
#if USE_MID
    composedColors->AddRGBPoint( mSurfaceScalarsColorMid,
                                 mnPositiveMidRedValue,
                                 mnPositiveMidGreenValue,
                                 mnPositiveMidBlueValue );
#endif
    composedColors->AddRGBPoint( mSurfaceScalarsColorMin,
                                 mnPositiveMinRedValue,
                                 mnPositiveMinGreenValue,
                                 mnPositiveMinBlueValue );
    composedColors->Build();

    // Set the composed scalars and colors in the view to draw on the
    // surface.
    mView->SetSurfaceScalars( composedScalars );
    mView->SetSurfaceScalarsColors( composedColors );

    // If we have surface scalars...
    if( surfaceScalars ) {

      // Show the label in the annoation (by annotation, here we mean the
      // text printed at the top of the black view screen).
      stringstream ssLabel;
      ssLabel << maSurfaceScalars[mnCurrentSurfaceScalars].msLabel;
      if ( "" != maSurfaceScalars[mnCurrentSurfaceScalars].msLabel2 ) {
        this->UpdateNumberOfSubjects();
        ssLabel << endl << "Contrast: "
                << maSurfaceScalars[mnCurrentSurfaceScalars].msLabel2
                << endl << "n=" 
                << mEntryNumberOfSubjects->GetValue() 
                << ", DOF="  
                << mQdecProject->GetGlmDesign()->GetDegreesOfFreedom();
      }
      mView->SetAnnotationMessage( ssLabel.str().c_str() );

      // Pass it to the view as the lookup scalars. This way we'll
      // only print surface scalar values and not curvature values when
      // clicking on a point.
      mView->SetSurfaceLookupScalars( surfaceScalars );

      // Make a quick color table with just the scalar colors and set
      // that in the view's legend colors.
      vtkSmartPointer<vtkColorTransferFunction> surfaceScalarsColors =
        vtkSmartPointer<vtkColorTransferFunction>::New();
      if( mbSurfaceScalarsColorReverse ) { // 'Reverse Values' check-boxed

        surfaceScalarsColors->AddRGBPoint( -mSurfaceScalarsColorMax,
                                           mnPositiveMaxRedValue,
                                           mnPositiveMaxGreenValue,
                                           mnPositiveMaxBlueValue );
#if USE_MID
        surfaceScalarsColors->AddRGBPoint( -mSurfaceScalarsColorMid,
                                           mnPositiveMidRedValue,
                                           mnPositiveMidGreenValue,
                                           mnPositiveMidBlueValue );
#endif
        surfaceScalarsColors->AddRGBPoint( -mSurfaceScalarsColorMin,
                                           mnPositiveMinRedValue,
                                           mnPositiveMinGreenValue,
                                           mnPositiveMinBlueValue );
        surfaceScalarsColors->AddRGBPoint( -mSurfaceScalarsColorMin + EPS,
                                           0.5, 0.5, 0.5 );
        surfaceScalarsColors->AddRGBPoint( mSurfaceScalarsColorMin - EPS,
                                           0.5, 0.5, 0.5 );
        surfaceScalarsColors->AddRGBPoint( mSurfaceScalarsColorMin,
                                           mnNegativeMinRedValue,
                                           mnNegativeMinGreenValue,
                                           mnNegativeMinBlueValue );
#if USE_MID
        surfaceScalarsColors->AddRGBPoint( mSurfaceScalarsColorMid,
                                           mnNegativeMidRedValue,
                                           mnNegativeMidGreenValue,
                                           mnNegativeMidBlueValue );
#endif
        surfaceScalarsColors->AddRGBPoint( mSurfaceScalarsColorMax,
                                           mnNegativeMaxRedValue,
                                           mnNegativeMaxGreenValue,
                                           mnNegativeMaxBlueValue );
      } else { // normal color scheme (not reversed)

        surfaceScalarsColors->AddRGBPoint( -mSurfaceScalarsColorMax,
                                           mnNegativeMaxRedValue,
                                           mnNegativeMaxGreenValue,
                                           mnNegativeMaxBlueValue );
#if USE_MID
        surfaceScalarsColors->AddRGBPoint( -mSurfaceScalarsColorMid,
                                           mnNegativeMidRedValue,
                                           mnNegativeMidGreenValue,
                                           mnNegativeMidBlueValue );
#endif
        surfaceScalarsColors->AddRGBPoint( -mSurfaceScalarsColorMin,
                                           mnNegativeMinRedValue,
                                           mnNegativeMinGreenValue,
                                           mnNegativeMinBlueValue );
        surfaceScalarsColors->AddRGBPoint( -mSurfaceScalarsColorMin + EPS,
                                           0.5, 0.5, 0.5 );
        surfaceScalarsColors->AddRGBPoint( mSurfaceScalarsColorMin - EPS,
                                           0.5, 0.5, 0.5 );
        surfaceScalarsColors->AddRGBPoint( mSurfaceScalarsColorMin,
                                           mnPositiveMinRedValue,
                                           mnPositiveMinGreenValue,
                                           mnPositiveMinBlueValue );
#if USE_MID
        surfaceScalarsColors->AddRGBPoint( mSurfaceScalarsColorMid,
                                           mnPositiveMidRedValue,
                                           mnPositiveMidGreenValue,
                                           mnPositiveMidBlueValue );
#endif
        surfaceScalarsColors->AddRGBPoint( mSurfaceScalarsColorMax,
                                           mnPositiveMaxRedValue,
                                           mnPositiveMaxGreenValue,
                                           mnPositiveMaxBlueValue );
      }


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
      mVertexPlot->IsLoaded() && 
      mVertexPlot->IsWindowShowing() ) {
    mMenuSaveGDFPostscript->SetStateToNormal();
  } else {
    mMenuSaveGDFPostscript->SetStateToDisabled();
  }
}

void
vtkKWQdecWindow::UpdateCommandStatus () {

  if( mPanel.GetPointer() ) {
    vtkKWWidget* panelSubjects = mPanel->GetPageWidget( "Subjects" );
    if( panelSubjects ) {
      if( this->mQdecProject->HaveDataTable() ) {
        panelSubjects->EnabledOn();
        panelSubjects->UpdateEnableState();
      } else {
        panelSubjects->EnabledOff();
        panelSubjects->UpdateEnableState();
      }
    }

    vtkKWWidget* panelDesign = mPanel->GetPageWidget( "Design" );
    if( panelDesign ) {
      if( this->mQdecProject->HaveDataTable() ) {
        panelDesign->EnabledOn();
        panelDesign->UpdateEnableState();
      } else {
        panelDesign->EnabledOff();
        panelDesign->UpdateEnableState();
      }
    }

    vtkKWWidget* panelDisplay = mPanel->GetPageWidget( "Display" );
    if( panelDisplay ) {
      if( maSurfaceSource.size() > 0 ) {
        panelDisplay->EnabledOn();
        panelDisplay->UpdateEnableState();
      } else {
        panelDisplay->EnabledOff();
        panelDisplay->UpdateEnableState();
      }
    }
  }

  assert( mMenuSaveDataTable );
  if( mQdecProject && mQdecProject->HaveDataTable() )
    mMenuSaveDataTable->SetStateToNormal();
  else
    mMenuSaveDataTable->SetStateToDisabled();

  assert( mMenuSaveProjectFile );
  if( mQdecProject && mQdecProject->GetGlmFitResults() )
    mMenuSaveProjectFile->SetStateToNormal();
  else
    mMenuSaveProjectFile->SetStateToDisabled();

  assert( mMenuSaveScatterPlotPostscript );
  if( ( mScatterPlotSelection != -1 ) && 
      mGraph &&
      ( mCurrentNotebookPanelName == ksSubjectsPanelName) )
  {
    mMenuSaveScatterPlotPostscript->SetStateToNormal();
  } else {
    mMenuSaveScatterPlotPostscript->SetStateToDisabled();
  }

  assert( mMenuSaveTIFF );
  assert( mMenuQuickSnapsTIFF );
  assert( mBtnSaveTIFF.GetPointer() );
  if( (maSurfaceSource.size() > 0) &&
      mView &&
      ( mCurrentNotebookPanelName == ksDisplayPanelName) ) {
    mBtnSaveTIFF->SetStateToNormal();
    mMenuSaveTIFF->SetStateToNormal();
    mMenuQuickSnapsTIFF->SetStateToNormal();
  } else {
    mBtnSaveTIFF->SetStateToDisabled();
    mMenuSaveTIFF->SetStateToDisabled();
    mMenuQuickSnapsTIFF->SetStateToDisabled();
  }

  assert( mMenuSaveGDFPostscript );
  if( mVertexPlot->IsLoaded() && 
      mVertexPlot->IsWindowShowing() ) {
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
  assert( mBtnShowCurvature.GetPointer() );
  assert( mMenuShowCurvature );
  assert( mEntrySelectVertex.GetPointer() );

  if( maSurfaceSource.size() > 0 ) {
    mBtnShowCursor->SetStateToNormal();
    mBtnShowCursor->SetSelectedState( mView->GetShowCursor() );
    mMenuShowCursor->SetStateToNormal();
    mMenuShowCursor->SetSelectedState( mView->GetShowCursor() );
    mBtnShowCurvature->SetStateToNormal();
    mBtnShowCurvature->SetSelectedState( mbShowCurvature );
    mMenuShowCurvature->SetStateToNormal();
    mMenuShowCurvature->SetSelectedState( mbShowCurvature );
    mEntrySelectVertex->SetStateToNormal();
    mEntrySelectVertex->SetCommand( this, "SelectSurfaceVertex" );
  } else {
    mBtnShowCursor->SetStateToDisabled();
    mMenuShowCursor->SetStateToDisabled();
    mBtnShowCurvature->SetStateToDisabled();
    mMenuShowCurvature->SetStateToDisabled();
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

  if( mVertexPlot->IsLoaded() && mROISource.GetPointer() ) {
    mMenuGraphAverageROI->SetStateToNormal();
  } else {
    mMenuGraphAverageROI->SetStateToDisabled();
  }
}

void
vtkKWQdecWindow::AnalyzeDesign () {

  assert( mListDiscreteFactors.GetPointer() );
  assert( mListContinuousFactors.GetPointer() );
  assert( mMenuMorphSmoothness.GetPointer() );
  assert( mMenuMorphMeasure.GetPointer() );
  assert( mMenuMorphHemisphere.GetPointer() );
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
    const char* measure = strdup( mMenuMorphMeasure->GetValue() );
    int smoothness = atoi( mMenuMorphSmoothness->GetValue() );
    const char* hemi = strdup( mMenuMorphHemisphere->GetValue() );
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
    int numNuisanceFactors = 0;
    const char* nuisanceFactors[1000];
    for( int nItem = 0; nItem < mListNuisanceFactors->GetNumberOfItems(); 
         nItem++ ) {
      if(mListNuisanceFactors->GetSelectState(nItem)) {
        nuisanceFactors[numNuisanceFactors] = 
          strdup( mListNuisanceFactors->GetItem( nItem ) );
        numNuisanceFactors++;
      }
      if(numNuisanceFactors>=999) break;
    }

    // Now create the design input files (checking for validity before running)
    if( this->mQdecProject->CreateGlmDesign
        (name,
         df1, // first discrete factor
         df2, // second discrete factor
         cf1, // first continuous factor
         cf2, // second continuous factor
         nuisanceFactors,
         numNuisanceFactors,
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

  // keep the QdecGlmDesign object up-to-date
  assert( mQdecProject ) ;
  assert( mQdecProject->GetGlmDesign() );
  QdecGlmDesign* design =  mQdecProject->GetGlmDesign();
  design->SetName( isDesignName );
}

bool
vtkKWQdecWindow::GetShowCurvature ( ) {

  return mbShowCurvature;
}

void
vtkKWQdecWindow::SetShowCurvature ( bool ibShow ) {

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
#if USE_MID
  double value;
  stringstream ssMid( isMid );
  ssMid >> value;

  this->SetSurfaceScalarsColorMid( value );
#endif
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

