/**
 * @file  MainWindow.cpp
 * @brief Main window.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:02 $
 *    $Revision: 1.154 $
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

#include <wx/xrc/xmlres.h>
#include <wx/config.h>
#include <wx/image.h>
#include <wx/sashwin.h>
#include <wx/laywin.h>
#include <wx/filedlg.h>
#include <wx/progdlg.h>
#include <wx/msgdlg.h>
#include <wx/splitter.h>
#include <wx/notebook.h>
#include <wx/app.h>
#include <wx/docview.h>
#include <wx/aui/auibook.h>
#include <wx/msgdlg.h>
#include <wx/textdlg.h>
#include <wx/filename.h>
#include <wx/spinctrl.h>
#include <wx/ffile.h>
#include <wx/dir.h>
#include "MainWindow.h"
#include "ControlPanel.h"
#include "PixelInfoPanel.h"
#include "RenderView2D.h"
#include "RenderView3D.h"
#include "LayerCollection.h"
#include "LayerMRI.h"
#include "LayerPropertiesMRI.h"
#include "LayerPropertiesWayPoints.h"
#include "LayerROI.h"
#include "LayerDTI.h"
#include "LayerSurface.h"
#include "LayerWayPoints.h"
#include "LayerOptimal.h"
#include "LayerPLabel.h"
#include "FSSurface.h"
#include "Interactor2DROIEdit.h"
#include "Interactor2DVoxelEdit.h"
#include "DialogPreferences.h"
#include "DialogLoadDTI.h"
#include "WindowQuickReference.h"
#include "StatusBar.h"
#include "DialogNewVolume.h"
#include "DialogNewROI.h"
#include "DialogNewWayPoints.h"
#include "DialogLoadVolume.h"
#include "WorkerThread.h"
#include "MyUtils.h"
#include "LUTDataHolder.h"
#include "BrushProperty.h"
#include "Cursor2D.h"
#include "Cursor3D.h"
#include "ToolWindowEdit.h"
#include "ToolWindowMeasure.h"
#include "DialogTransformVolume.h"
#include "DialogOptimalVolume.h"
#include "WindowHistogram.h"
#include "WindowOverlayConfiguration.h"
#include "WindowConnectivityConfiguration.h"
#include "SurfaceOverlay.h"
#include "SurfaceOverlayProperties.h"
#include "DialogSaveVolumeAs.h"
#include "VolumeFilterGradient.h"
#include "DialogGradientVolume.h"
#include "LayerPropertiesSurface.h"
#include "DialogLoadPVolumes.h"
#include "ConnectivityData.h"
#include "DialogSaveScreenshot.h"
#include "DialogSavePointSetAs.h"
#include "DialogLoadPointSet.h"
#include "DialogLoadROI.h"
#include "DialogSavePoint.h"
#include "DialogWriteMovieFrames.h"
#include "VolumeFilterMean.h"
#include "VolumeFilterMedian.h"
#include "VolumeFilterConvolve.h"
#include "chronometer.h"
#include "DialogVolumeFilter.h"
#include "DialogRepositionSurface.h"
#include "DialogCropVolume.h"
#include "VolumeCropper.h"

#define CTRL_PANEL_WIDTH 240

#define ID_TOOL_GOTO_POINT_1  (wxID_HIGHEST+100)

#define ID_TIMER_WRITE_MOVIE_FRAMES     1

// ----------------------------------------------------------------------------
// event tables and other macros for wxWindows
// ----------------------------------------------------------------------------
IMPLEMENT_CLASS(MainWindow, wxFrame)
BEGIN_EVENT_TABLE(MainWindow, wxFrame)
  EVT_MENU        ( XRCID( "ID_FILE_OPEN" ),              MainWindow::OnFileOpen )
  EVT_MENU        ( XRCID( "ID_FILE_NEW" ),               MainWindow::OnFileNew )
  EVT_UPDATE_UI   ( XRCID( "ID_FILE_NEW" ),               MainWindow::OnFileNewUpdateUI )
  EVT_MENU        ( XRCID( "ID_FILE_SAVE" ),              MainWindow::OnFileSave )
  EVT_UPDATE_UI   ( XRCID( "ID_FILE_SAVE" ),              MainWindow::OnFileSaveUpdateUI )
  EVT_MENU        ( XRCID( "ID_FILE_SAVE_AS" ),           MainWindow::OnFileSaveAs )
  EVT_UPDATE_UI   ( XRCID( "ID_FILE_SAVE_AS" ),           MainWindow::OnFileSaveAsUpdateUI )
  EVT_MENU        ( XRCID( "ID_FILE_EXIT") ,              MainWindow::OnFileExit )
  EVT_MENU_RANGE  ( wxID_FILE1, wxID_FILE1+16,            MainWindow::OnFileRecent )
  EVT_MENU        ( XRCID( "ID_FILE_NEW_ROI" ),           MainWindow::OnFileNewROI )
  EVT_UPDATE_UI   ( XRCID( "ID_FILE_NEW_ROI" ),           MainWindow::OnFileNewROIUpdateUI )
  EVT_MENU        ( XRCID( "ID_FILE_LOAD_ROI" ),          MainWindow::OnFileLoadROI )
  EVT_UPDATE_UI   ( XRCID( "ID_FILE_LOAD_ROI" ),          MainWindow::OnFileLoadROIUpdateUI )
  EVT_MENU        ( XRCID( "ID_FILE_SAVE_ROI" ),          MainWindow::OnFileSaveROI )
  EVT_UPDATE_UI   ( XRCID( "ID_FILE_SAVE_ROI" ),          MainWindow::OnFileSaveROIUpdateUI )
  EVT_MENU        ( XRCID( "ID_FILE_SAVE_ROI_AS" ),       MainWindow::OnFileSaveROIAs )
  EVT_UPDATE_UI   ( XRCID( "ID_FILE_SAVE_ROI_AS" ),       MainWindow::OnFileSaveROIAsUpdateUI )
  EVT_MENU        ( XRCID( "ID_FILE_LOAD_DTI" ),          MainWindow::OnFileLoadDTI )
  EVT_MENU        ( XRCID( "ID_FILE_LOAD_PVOLUMES" ),     MainWindow::OnFileLoadPVolumes )
  EVT_MENU        ( XRCID( "ID_FILE_NEW_WAYPOINTS" ),     MainWindow::OnFileNewWayPoints )
  EVT_UPDATE_UI   ( XRCID( "ID_FILE_NEW_WAYPOINTS" ),     MainWindow::OnFileNewWayPointsUpdateUI )
  EVT_MENU        ( XRCID( "ID_FILE_LOAD_WAYPOINTS" ),    MainWindow::OnFileLoadWayPoints )
  EVT_UPDATE_UI   ( XRCID( "ID_FILE_LOAD_WAYPOINTS" ),    MainWindow::OnFileLoadWayPointsUpdateUI )
  EVT_MENU        ( XRCID( "ID_FILE_SAVE_WAYPOINTS" ),    MainWindow::OnFileSaveWayPoints )
  EVT_UPDATE_UI   ( XRCID( "ID_FILE_SAVE_WAYPOINTS" ),    MainWindow::OnFileSaveWayPointsUpdateUI )
  EVT_MENU        ( XRCID( "ID_FILE_SAVE_WAYPOINTS_AS" ), MainWindow::OnFileSaveWayPointsAs )
  EVT_UPDATE_UI   ( XRCID( "ID_FILE_SAVE_WAYPOINTS_AS" ), MainWindow::OnFileSaveWayPointsAsUpdateUI )
  EVT_MENU        ( XRCID( "ID_FILE_SCREENSHOT" ),        MainWindow::OnFileSaveScreenshot )
  EVT_UPDATE_UI   ( XRCID( "ID_FILE_SCREENSHOT" ),        MainWindow::OnFileSaveScreenshotUpdateUI )
  EVT_MENU        ( XRCID( "ID_FILE_SAVE_MOVIE_FRAMES" ), MainWindow::OnFileSaveMovieFrames )
  EVT_UPDATE_UI   ( XRCID( "ID_FILE_SAVE_MOVIE_FRAMES" ), MainWindow::OnFileSaveMovieFramesUpdateUI )
  EVT_MENU        ( XRCID( "ID_FILE_LOAD_SURFACE" ),      MainWindow::OnFileLoadSurface )
  EVT_MENU        ( XRCID( "ID_FILE_SAVE_SURFACE" ),      MainWindow::OnFileSaveSurface )
  EVT_UPDATE_UI   ( XRCID( "ID_FILE_SAVE_SURFACE" ),      MainWindow::OnFileSaveSurfaceUpdateUI )
  EVT_MENU        ( XRCID( "ID_FILE_SAVE_SURFACE_AS" ),   MainWindow::OnFileSaveSurfaceAs )
  EVT_UPDATE_UI   ( XRCID( "ID_FILE_SAVE_SURFACE_AS" ),   MainWindow::OnFileSaveSurfaceAsUpdateUI )
  EVT_MENU        ( XRCID( "ID_MODE_NAVIGATE" ),          MainWindow::OnModeNavigate )
  EVT_UPDATE_UI   ( XRCID( "ID_MODE_NAVIGATE" ),          MainWindow::OnModeNavigateUpdateUI )
  EVT_MENU        ( XRCID( "ID_MODE_MEASURE" ),           MainWindow::OnModeMeasure )
  EVT_UPDATE_UI   ( XRCID( "ID_MODE_MEASURE" ),           MainWindow::OnModeMeasureUpdateUI )
  EVT_MENU        ( XRCID( "ID_MODE_VOXEL_EDIT" ),        MainWindow::OnModeVoxelEdit )
  EVT_UPDATE_UI   ( XRCID( "ID_MODE_VOXEL_EDIT" ),        MainWindow::OnModeVoxelEditUpdateUI )
  EVT_MENU        ( XRCID( "ID_MODE_ROI_EDIT" ),          MainWindow::OnModeROIEdit )
  EVT_UPDATE_UI   ( XRCID( "ID_MODE_ROI_EDIT" ),          MainWindow::OnModeROIEditUpdateUI )
  EVT_MENU        ( XRCID( "ID_MODE_WAYPOINTS_EDIT" ),    MainWindow::OnModeWayPointsEdit )
  EVT_UPDATE_UI   ( XRCID( "ID_MODE_WAYPOINTS_EDIT" ),    MainWindow::OnModeWayPointsEditUpdateUI )
  EVT_MENU        ( XRCID( "ID_EDIT_COPY" ),              MainWindow::OnEditCopy )
  EVT_UPDATE_UI   ( XRCID( "ID_EDIT_COPY" ),              MainWindow::OnEditCopyUpdateUI )
  EVT_MENU        ( XRCID( "ID_EDIT_COPY_STRUCTURE" ),    MainWindow::OnEditCopyStructure )
  EVT_UPDATE_UI   ( XRCID( "ID_EDIT_COPY_STRUCTURE" ),    MainWindow::OnEditCopyStructureUpdateUI )
  EVT_MENU        ( XRCID( "ID_EDIT_PASTE" ),             MainWindow::OnEditPaste )
  EVT_UPDATE_UI   ( XRCID( "ID_EDIT_PASTE" ),             MainWindow::OnEditPasteUpdateUI )
  EVT_MENU        ( XRCID( "ID_EDIT_UNDO" ),              MainWindow::OnEditUndo )
  EVT_UPDATE_UI   ( XRCID( "ID_EDIT_UNDO" ),              MainWindow::OnEditUndoUpdateUI )
  EVT_MENU        ( XRCID( "ID_EDIT_REDO" ),              MainWindow::OnEditRedo )
  EVT_UPDATE_UI   ( XRCID( "ID_EDIT_REDO" ),              MainWindow::OnEditRedoUpdateUI )
  EVT_MENU        ( XRCID( "ID_EDIT_RENAME" ),            MainWindow::OnEditRename )
  EVT_UPDATE_UI   ( XRCID( "ID_EDIT_RENAME" ),            MainWindow::OnEditRenameUpdateUI )
  EVT_MENU        ( XRCID( "ID_EDIT_PREFERENCES" ),       MainWindow::OnEditPreferences )
  EVT_MENU        ( XRCID( "ID_VIEW_LAYOUT_1X1" ),        MainWindow::OnViewLayout1X1 )
  EVT_UPDATE_UI   ( XRCID( "ID_VIEW_LAYOUT_1X1" ),        MainWindow::OnViewLayout1X1UpdateUI )
  EVT_MENU        ( XRCID( "ID_VIEW_LAYOUT_2X2" ),        MainWindow::OnViewLayout2X2 )
  EVT_UPDATE_UI   ( XRCID( "ID_VIEW_LAYOUT_2X2" ),        MainWindow::OnViewLayout2X2UpdateUI )
  EVT_MENU        ( XRCID( "ID_VIEW_LAYOUT_1N3" ),        MainWindow::OnViewLayout1N3 )
  EVT_UPDATE_UI   ( XRCID( "ID_VIEW_LAYOUT_1N3" ),        MainWindow::OnViewLayout1N3UpdateUI )
  EVT_MENU        ( XRCID( "ID_VIEW_LAYOUT_1N3_H" ),      MainWindow::OnViewLayout1N3_H )
  EVT_UPDATE_UI   ( XRCID( "ID_VIEW_LAYOUT_1N3_H" ),      MainWindow::OnViewLayout1N3_HUpdateUI )
  EVT_MENU        ( XRCID( "ID_VIEW_SAGITTAL" ),          MainWindow::OnViewSagittal )
  EVT_UPDATE_UI   ( XRCID( "ID_VIEW_SAGITTAL" ),          MainWindow::OnViewSagittalUpdateUI )
  EVT_MENU        ( XRCID( "ID_VIEW_CORONAL" ),           MainWindow::OnViewCoronal )
  EVT_UPDATE_UI   ( XRCID( "ID_VIEW_CORONAL" ),           MainWindow::OnViewCoronalUpdateUI )
  EVT_MENU        ( XRCID( "ID_VIEW_AXIAL" ),             MainWindow::OnViewAxial )
  EVT_UPDATE_UI   ( XRCID( "ID_VIEW_AXIAL" ),             MainWindow::OnViewAxialUpdateUI )
  EVT_MENU        ( XRCID( "ID_VIEW_3D" ),                MainWindow::OnView3D )
  EVT_UPDATE_UI   ( XRCID( "ID_VIEW_3D" ),                MainWindow::OnView3DUpdateUI )
  EVT_MENU        ( XRCID( "ID_VIEW_RESET" ),             MainWindow::OnViewReset )
  EVT_UPDATE_UI   ( XRCID( "ID_VIEW_RESET" ),             MainWindow::OnViewResetUpdateUI )
  EVT_MENU        ( XRCID( "ID_VIEW_SNAP_TO_AXIS" ),      MainWindow::OnViewSnapToAxis )
  EVT_UPDATE_UI   ( XRCID( "ID_VIEW_SNAP_TO_AXIS" ),      MainWindow::OnViewSnapToAxisUpdateUI )
  EVT_MENU        ( XRCID( "ID_VIEW_SLICE_FRAMES" ),      MainWindow::OnViewSliceFrames )
  EVT_UPDATE_UI   ( XRCID( "ID_VIEW_SLICE_FRAMES" ),      MainWindow::OnViewSliceFramesUpdateUI )
  EVT_MENU        ( XRCID( "ID_VIEW_SCALAR_BAR" ),        MainWindow::OnViewScalarBar )
  EVT_UPDATE_UI   ( XRCID( "ID_VIEW_SCALAR_BAR" ),        MainWindow::OnViewScalarBarUpdateUI )
  EVT_MENU        ( XRCID( "ID_VIEW_COORDINATE" ),        MainWindow::OnViewCoordinate )
  EVT_UPDATE_UI   ( XRCID( "ID_VIEW_COORDINATE" ),        MainWindow::OnViewCoordinateUpdateUI )
  EVT_MENU        ( XRCID( "ID_VIEW_CYCLE_LAYER" ),       MainWindow::OnViewCycleLayer )
  EVT_UPDATE_UI   ( XRCID( "ID_VIEW_CYCLE_LAYER" ),       MainWindow::OnViewCycleLayerUpdateUI )
  EVT_MENU        ( XRCID( "ID_VIEW_REVERSE_CYCLE_LAYER" ),   MainWindow::OnViewReverseCycleLayer )
  EVT_UPDATE_UI   ( XRCID( "ID_VIEW_REVERSE_CYCLE_LAYER" ),   MainWindow::OnViewCycleLayerUpdateUI )
  EVT_MENU        ( XRCID( "ID_VIEW_TOGGLE_VOLUME" ),     MainWindow::OnViewToggleVolumeVisibility )
  EVT_UPDATE_UI   ( XRCID( "ID_VIEW_TOGGLE_VOLUME" ),     MainWindow::OnViewToggleVolumeVisibilityUpdateUI )
  EVT_MENU        ( XRCID( "ID_VIEW_TOGGLE_ROI" ),        MainWindow::OnViewToggleROIVisibility )
  EVT_UPDATE_UI   ( XRCID( "ID_VIEW_TOGGLE_ROI" ),        MainWindow::OnViewToggleROIVisibilityUpdateUI )
  EVT_MENU        ( XRCID( "ID_VIEW_TOGGLE_SURFACE" ),    MainWindow::OnViewToggleSurfaceVisibility )
  EVT_UPDATE_UI   ( XRCID( "ID_VIEW_TOGGLE_SURFACE" ),    MainWindow::OnViewToggleSurfaceVisibilityUpdateUI )
  EVT_MENU        ( XRCID( "ID_VIEW_TOGGLE_WAYPOINTS" ),  MainWindow::OnViewToggleWayPointsVisibility )
  EVT_UPDATE_UI   ( XRCID( "ID_VIEW_TOGGLE_WAYPOINTS" ),  MainWindow::OnViewToggleWayPointsVisibilityUpdateUI )
  EVT_MENU        ( XRCID( "ID_VIEW_TOGGLE_CURSOR" ),     MainWindow::OnViewToggleCursorVisibility )
  EVT_UPDATE_UI   ( XRCID( "ID_VIEW_TOGGLE_CURSOR" ),     MainWindow::OnViewToggleCursorVisibilityUpdateUI )
  EVT_MENU        ( XRCID( "ID_VIEW_TOGGLE_VOXEL_COORDS" ), MainWindow::OnViewToggleVoxelCoordinates )
  
  EVT_MENU        ( XRCID( "ID_VIEW_SURFACE_MAIN" ),      MainWindow::OnViewSurfaceMain )
  EVT_UPDATE_UI   ( XRCID( "ID_VIEW_SURFACE_MAIN" ),      MainWindow::OnViewSurfaceMainUpdateUI )
  EVT_MENU        ( XRCID( "ID_VIEW_SURFACE_INFLATED" ),  MainWindow::OnViewSurfaceInflated )
  EVT_UPDATE_UI   ( XRCID( "ID_VIEW_SURFACE_INFLATED" ),  MainWindow::OnViewSurfaceInflatedUpdateUI )
  EVT_MENU        ( XRCID( "ID_VIEW_SURFACE_WHITE" ),     MainWindow::OnViewSurfaceWhite )
  EVT_UPDATE_UI   ( XRCID( "ID_VIEW_SURFACE_WHITE" ),     MainWindow::OnViewSurfaceWhiteUpdateUI )
  EVT_MENU        ( XRCID( "ID_VIEW_SURFACE_PIAL" ),      MainWindow::OnViewSurfacePial )
  EVT_UPDATE_UI   ( XRCID( "ID_VIEW_SURFACE_PIAL" ),      MainWindow::OnViewSurfacePialUpdateUI )
  EVT_MENU        ( XRCID( "ID_VIEW_SURFACE_ORIGINAL" ),  MainWindow::OnViewSurfaceOriginal )
  EVT_UPDATE_UI   ( XRCID( "ID_VIEW_SURFACE_ORIGINAL" ),  MainWindow::OnViewSurfaceOriginalUpdateUI )
  
  EVT_MENU        ( XRCID( "ID_VIEW_CONTROL_PANEL" ),     MainWindow::OnViewControlPanel )
  EVT_UPDATE_UI   ( XRCID( "ID_VIEW_CONTROL_PANEL" ),     MainWindow::OnViewControlPanelUpdateUI )
  
  EVT_MENU        ( XRCID( "ID_VIEW_HISTOGRAM" ),         MainWindow::OnViewHistogram )
  EVT_UPDATE_UI   ( XRCID( "ID_VIEW_HISTOGRAM" ),         MainWindow::OnViewHistogramUpdateUI )
  
  EVT_MENU        ( XRCID( "ID_LAYER_SHOW_ALL" ),         MainWindow::OnLayerShowAll )
  EVT_UPDATE_UI   ( XRCID( "ID_LAYER_SHOW_ALL" ),         MainWindow::OnLayerShowHideAllUpdateUI )
  EVT_MENU        ( XRCID( "ID_LAYER_HIDE_ALL" ),         MainWindow::OnLayerHideAll )
  EVT_UPDATE_UI   ( XRCID( "ID_LAYER_HIDE_ALL" ),         MainWindow::OnLayerShowHideAllUpdateUI )
  
  EVT_MENU        ( XRCID( "ID_TOOL_ROTATE_VOLUME" ),     MainWindow::OnToolTransformVolume )
  EVT_UPDATE_UI   ( XRCID( "ID_TOOL_ROTATE_VOLUME" ),     MainWindow::OnToolTransformVolumeUpdateUI )
  EVT_MENU        ( XRCID( "ID_TOOL_CROP_VOLUME" ),       MainWindow::OnToolCropVolume )
  EVT_UPDATE_UI   ( XRCID( "ID_TOOL_CROP_VOLUME" ),       MainWindow::OnToolCropVolumeUpdateUI )
  EVT_MENU        ( XRCID( "ID_TOOL_OPTIMAL_VOLUME" ),    MainWindow::OnToolOptimalVolume )
  EVT_UPDATE_UI   ( XRCID( "ID_TOOL_OPTIMAL_VOLUME" ),    MainWindow::OnToolOptimalVolumeUpdateUI )
  EVT_MENU        ( XRCID( "ID_TOOL_GRADIENT_VOLUME" ),   MainWindow::OnToolGradientVolume )
  EVT_UPDATE_UI   ( XRCID( "ID_TOOL_GRADIENT_VOLUME" ),   MainWindow::OnToolGradientVolumeUpdateUI )
  EVT_MENU        ( XRCID( "ID_TOOL_SAVE_POINT" ),        MainWindow::OnToolSaveGotoPoint )
  EVT_UPDATE_UI   ( XRCID( "ID_TOOL_SAVE_POINT" ),        MainWindow::OnToolSaveGotoPointUpdateUI )
  EVT_MENU        ( XRCID( "ID_TOOL_GOTO_POINT" ),        MainWindow::OnToolGotoPoint )
  EVT_UPDATE_UI   ( XRCID( "ID_TOOL_GOTO_POINT" ),        MainWindow::OnToolGotoPointUpdateUI )
  EVT_UPDATE_UI   ( XRCID( "ID_TOOL_GOTO_POINT_SUBMENU" ),        MainWindow::OnToolGotoPointUpdateUI )
  EVT_MENU_RANGE  ( ID_TOOL_GOTO_POINT_1, ID_TOOL_GOTO_POINT_1+1000,            MainWindow::OnToolMenuGotoPoint )
  EVT_MENU        ( XRCID( "ID_TOOL_REPOSITION_SURFACE" ),   MainWindow::OnToolRepositionSurface )
  EVT_UPDATE_UI   ( XRCID( "ID_TOOL_REPOSITION_SURFACE" ),   MainWindow::OnToolRepositionSurfaceUpdateUI )
  
  EVT_MENU        ( XRCID( "ID_FILTER_MEAN" ),            MainWindow::OnFilterMean )
  EVT_UPDATE_UI   ( XRCID( "ID_FILTER_MEAN" ),            MainWindow::OnFilterUpdateUI )
  EVT_MENU        ( XRCID( "ID_FILTER_MEDIAN" ),          MainWindow::OnFilterMedian )
  EVT_UPDATE_UI   ( XRCID( "ID_FILTER_MEDIAN" ),          MainWindow::OnFilterUpdateUI )
  EVT_MENU        ( XRCID( "ID_FILTER_CONVOLVE" ),        MainWindow::OnFilterConvolve )
  EVT_UPDATE_UI   ( XRCID( "ID_FILTER_CONVOLVE" ),        MainWindow::OnFilterUpdateUI )
  EVT_MENU        ( XRCID( "ID_FILTER_GRADIENT" ),        MainWindow::OnFilterGradient )
  EVT_UPDATE_UI   ( XRCID( "ID_FILTER_GRADIENT" ),        MainWindow::OnFilterUpdateUI )
  
  EVT_MENU        ( XRCID( "ID_TOOL_LABEL_STATS" ),       MainWindow::OnToolLabelStats )
  EVT_UPDATE_UI   ( XRCID( "ID_TOOL_LABEL_STATS" ),       MainWindow::OnToolLabelStatsUpdateUI )
  
  EVT_MENU  ( XRCID( "ID_HELP_QUICK_REF" ),               MainWindow::OnHelpQuickReference )
  EVT_MENU  ( XRCID( "ID_HELP_ABOUT" ),                   MainWindow::OnHelpAbout )
  /*
  EVT_SASH_DRAGGED_RANGE(ID_LOG_WINDOW, ID_LOG_WINDOW, MainWindow::OnSashDrag)
  EVT_IDLE  (MainWindow::OnIdle)
  EVT_MENU  (XRCID("ID_EVENT_LOAD_DATA"), MainWindow::OnWorkerEventLoadData)
  EVT_MENU  (XRCID("ID_EVENT_CALCULATE_MATRIX"), MainWindow::OnWorkerEventCalculateMatrix)
  */
  
  EVT_MENU      ( ID_WORKER_THREAD,                   MainWindow::OnWorkerThreadResponse )
  EVT_MENU      ( ID_THREAD_BUILD_CONTOUR,            MainWindow::OnBuildContourThreadResponse )
  EVT_SPINCTRL  ( XRCID( "ID_SPIN_BRUSH_SIZE" ),      MainWindow::OnSpinBrushSize )
  EVT_SPINCTRL  ( XRCID( "ID_SPIN_BRUSH_TOLERANCE" ), MainWindow::OnSpinBrushTolerance )
  EVT_CHECKBOX  ( XRCID( "ID_CHECK_TEMPLATE" ),       MainWindow::OnCheckBrushTemplate )
  EVT_CHOICE    ( XRCID( "ID_CHOICE_TEMPLATE" ),      MainWindow::OnChoiceBrushTemplate )
  
  EVT_ENTER_WINDOW  ( MainWindow::OnMouseEnterWindow )
  EVT_ACTIVATE      ( MainWindow::OnActivate )
  EVT_CLOSE         ( MainWindow::OnClose )
  EVT_KEY_DOWN      ( MainWindow::OnKeyDown )
  EVT_ICONIZE       ( MainWindow::OnIconize )

  EVT_TIMER         ( ID_TIMER_WRITE_MOVIE_FRAMES,    MainWindow::OnTimerWriteMovieFrames )

END_EVENT_TABLE()

// ----------------------------------------------------------------------------
// main frame
// ----------------------------------------------------------------------------

// frame constructor

void MainWindow::InitWidgetsFromXRC( wxWindow* parent)
{  
// m_bLoading = false;
// m_bSaving = false;
  m_bProcessing = false;
  m_bResampleToRAS = false;
  m_bToUpdateToolbars = false;
  m_bDoScreenshot = false;
  m_layerVolumeRef = NULL;
  m_nPrevActiveViewId = -1;
  m_nDefaultSampleMethod = SAMPLE_NEAREST;
  m_bDefaultConform = false;
  m_luts = new LUTDataHolder();
  m_propertyBrush = new BrushProperty();
  m_connectivity = new ConnectivityData();
  m_volumeCropper = new VolumeCropper();

  wxXmlResource::Get()->LoadObject( this, parent, wxT("ID_MAIN_WINDOW"), wxT("wxFrame") );

  // must be created before the following controls
  m_layerCollectionManager = new LayerCollectionManager();

//  m_panelToolbarHolder = XRCCTRL( *this, "ID_PANEL_HOLDER", wxPanel );
// wxBoxSizer* sizer = (wxBoxSizer*)panelHolder->GetSizer(); //new wxBoxSizer( wxVERTICAL );

  // create the main splitter window
// m_splitterMain = new wxSplitterWindow( this, wxID_ANY, wxDefaultPosition, wxDefaultSize, wxSP_LIVE_UPDATE );
  m_splitterMain = XRCCTRL( *this, "ID_SPLITTER_MAIN", wxSplitterWindow );
  m_splitterMain->SetMinimumPaneSize( 80 );
// sizer->Add( m_splitterMain, 1, wxEXPAND );

   m_toolbarMain = GetToolBar();

// this->SetSizer( sizer );
// sizer->Add( ( wxToolBar* )XRCCTRL( *this, "m_toolBar2", wxToolBar ), 0, wxEXPAND );

  // create the main control panel
  m_controlPanel = new ControlPanel( m_splitterMain );
  m_splitterSub = new wxSplitterWindow( m_splitterMain, wxID_ANY, wxDefaultPosition, wxDefaultSize, wxSP_LIVE_UPDATE );

  m_splitterMain->SplitVertically( m_controlPanel, m_splitterSub );

  // create the panel holder for the 4 render views
  m_renderViewHolder = new wxPanel( m_splitterSub, wxID_ANY, wxDefaultPosition, wxDefaultSize, wxTAB_TRAVERSAL );
  m_pixelInfoPanel = new PixelInfoPanel( m_splitterSub );

  m_splitterSub->SplitHorizontally( m_renderViewHolder, m_pixelInfoPanel );
  m_splitterSub->SetMinimumPaneSize( 80 );
  m_splitterSub->SetSashGravity( 1 );

  // create 4 render views
  m_viewSagittal = new RenderView2D( 0, m_renderViewHolder );
  m_viewCoronal = new RenderView2D( 1, m_renderViewHolder );
  m_viewAxial = new RenderView2D( 2, m_renderViewHolder );
  m_view3D = new RenderView3D( m_renderViewHolder );
  m_viewSagittal->SetFocusFrameColor( 1, 0, 0 );
  m_viewCoronal->SetFocusFrameColor( 0, 1, 0 );
  m_viewAxial->SetFocusFrameColor( 0, 0, 1 );
  m_viewRender[0] = m_viewSagittal;
  m_viewRender[1] = m_viewCoronal;
  m_viewRender[2] = m_viewAxial;
  m_viewRender[3] = m_view3D;
  for ( int i = 0; i < 3; i++ )
  {
    // m_viewRender[i]->Render();
    for ( int j = 0; j < 3; j++ )
    {
      if ( i != j )
        m_viewRender[i]->AddListener( m_viewRender[j] );
    }
  }

  m_layerCollectionManager->AddListener( m_viewAxial );
  m_layerCollectionManager->AddListener( m_viewSagittal );
  m_layerCollectionManager->AddListener( m_viewCoronal );
  m_layerCollectionManager->AddListener( m_view3D );
  m_layerCollectionManager->AddListener( this );
  GetLayerCollection( "MRI" )->AddListener( m_pixelInfoPanel );
  GetLayerCollection( "Surface" )->AddListener( m_pixelInfoPanel );
  
  m_connectivity->AddListener( m_view3D );
  
  m_wndQuickReference = new WindowQuickReference( this );
  m_wndQuickReference->Hide();

  m_statusBar = new StatusBar( this );
  SetStatusBar( m_statusBar );
  PositionStatusBar();

  m_toolWindowEdit = NULL;
  m_toolWindowMeasure = NULL;
  m_dlgTransformVolume = NULL;
  m_dlgGradientVolume = NULL;
  m_dlgSaveScreenshot = NULL;
  m_dlgSavePoint = NULL;
  m_dlgWriteMovieFrames = NULL;
  
  m_menuGotoPoints = NULL;

  m_wndHistogram = new WindowHistogram( this );
  m_wndHistogram->Hide();
  
  m_wndOverlayConfiguration = new WindowOverlayConfiguration( this );
  m_wndOverlayConfiguration->Hide();
  
  m_wndConnectivityConfiguration = new WindowConnectivityConfiguration( this );
  m_wndConnectivityConfiguration->Hide();
  
  m_dlgWriteMovieFrames = new DialogWriteMovieFrames( this );
  m_dlgWriteMovieFrames->Hide();
  
  m_dlgRepositionSurface = new DialogRepositionSurface( this );
  m_dlgRepositionSurface->Hide();  
  m_view3D->AddListener( m_dlgRepositionSurface );
  GetLayerCollection( "Surface" )->AddListener( m_dlgRepositionSurface );
  
  m_dlgCropVolume = new DialogCropVolume( this );
  m_dlgCropVolume->Hide();
  GetLayerCollection( "MRI" )->AddListener( m_dlgCropVolume );
  
  m_volumeCropper->AddListener( this );
  m_volumeCropper->AddListener( m_dlgCropVolume );
  m_layerCollectionManager->AddListener( m_volumeCropper );
  
  UpdateToolbars();
  
  m_nViewLayout = VL_2X2;
  m_nMainView = MV_Sagittal;
  m_nMaxRecentFiles = 10;
  wxConfigBase* config = wxConfigBase::Get();
  m_fileHistory = new wxFileHistory( m_nMaxRecentFiles, wxID_FILE1 );
  wxMenu* fileMenu = GetMenuBar()->GetMenu( 0 )->FindItem( XRCID("ID_FILE_SUBMENU_RECENT") )->GetSubMenu();
  m_settingsGeneral.SaveCopy = true;
  m_settingsMovieFrames.AngleStep = 1.0;
  m_settingsMovieFrames.OutputExtension = _("png");
  if ( config )
  {
    int x = config->Read( _T("/MainWindow/PosX"), 50L );
    int y = config->Read( _T("/MainWindow/PosY"), 50L );
    int w = config->Read( _T("/MainWindow/Width"), 950L );
    int h = config->Read( _T("/MainWindow/Height"), 750L );
    SetSize( x, y, w, h );
    m_splitterMain->SetSashPosition( config->Read( _T("/MainWindow/SplitterPosition"), 260L ) );
    // m_splitterSub->SetSashPosition( config->Read( _T("/MainWindow/SplitterPositionSub"), 50L ) );
    // m_splitterSub->UpdateSize();
    m_strLastDir = config->Read( _("/MainWindow/LastDir"), _("") );
    m_fileHistory->Load( *config );
    m_fileHistory->UseMenu( fileMenu );
    m_fileHistory->AddFilesToMenu( fileMenu );
    m_nViewLayout = config->Read( _("/MainWindow/ViewLayout"), m_nViewLayout );
    m_nMainView = config->Read( _("/MainWindow/MainView"), m_nMainView );

    m_settingsGeneral.BackgroundColor = wxColour( config->Read( _("/RenderWindow/BackgroundColor"), _("rgb(0,0,0)") ) );
    m_settingsGeneral.CursorColor = wxColour( config->Read( _("/RenderWindow/CursorColor"), _("rgb(255,0,0)") ) );
    m_settingsGeneral.CursorStyle = config->Read( _("/RenderWindow/CursorStyle"), Cursor2D::CS_Short );
    config->Read( _("/Settings/SaveCopy"), &m_settingsGeneral.SaveCopy, true );
    
    m_viewAxial   ->SetBackgroundColor( m_settingsGeneral.BackgroundColor );
    m_viewSagittal->SetBackgroundColor( m_settingsGeneral.BackgroundColor );
    m_viewCoronal ->SetBackgroundColor( m_settingsGeneral.BackgroundColor );
    m_view3D      ->SetBackgroundColor( m_settingsGeneral.BackgroundColor );
    m_viewAxial->GetCursor2D()    ->SetColor( m_settingsGeneral.CursorColor );
    m_viewSagittal->GetCursor2D() ->SetColor( m_settingsGeneral.CursorColor );
    m_viewCoronal->GetCursor2D()  ->SetColor( m_settingsGeneral.CursorColor );
    m_view3D->GetCursor3D()       ->SetColor( m_settingsGeneral.CursorColor );
    m_viewAxial->GetCursor2D()    ->SetStyle( m_settingsGeneral.CursorStyle );
    m_viewSagittal->GetCursor2D() ->SetStyle( m_settingsGeneral.CursorStyle );
    m_viewCoronal->GetCursor2D()  ->SetStyle( m_settingsGeneral.CursorStyle );

    m_nScreenshotFilterIndex = config->Read( _("/MainWindow/ScreenshotFilterIndex"), 0L );

    bool bScalarBar;
    config->Read( _("/RenderWindow/ShowScalarBar"), &bScalarBar, false );
    if ( bScalarBar )
      ShowScalarBar( bScalarBar );
    
    config->Read( _("/RenderWindow/SyncZoomFactor"), &m_settings2D.SyncZoomFactor, true );
    config->Read( _("/Screenshot/Magnification" ), &m_settingsScreenshot.Magnification, 1 );
    config->Read( _("/Screenshot/HideCursor" ), &m_settingsScreenshot.HideCursor, false );
    config->Read( _("/Screenshot/HideCoords" ), &m_settingsScreenshot.HideCoords, false );
    config->Read( _("/Screenshot/AntiAliasing" ), &m_settingsScreenshot.AntiAliasing, false );

    bool bShow = true;
    config->Read( _("/MainWindow/ShowControlPanel" ), &bShow, true );
    
    wxString tempStrg = config->Read( _("/MainWindow/GotoPoints"), _("") );
    if ( !tempStrg.IsEmpty() )
      m_strGotoPoints = MyUtils::SplitString( tempStrg, ";" );

    ShowControlPanel( bShow );
  }
  SetViewLayout( m_nViewLayout );
  
//  UpdateGotoPoints();
  
  m_timerWriteMovieFrames.SetOwner( this, ID_TIMER_WRITE_MOVIE_FRAMES );

}

// frame destructor
MainWindow::~MainWindow()
{
  for ( int i = 0; i < 4; i++ )
  {
    m_viewRender[i]->Delete();
    m_viewRender[i] = NULL;
  }
  
  delete m_luts;
  delete m_propertyBrush;
  if ( m_connectivity )
    delete m_connectivity;
  
  if ( m_volumeCropper )
    delete m_volumeCropper;
  m_volumeCropper = NULL;
  
  if ( m_menuGotoPoints )
    delete m_menuGotoPoints;
  
  delete m_layerCollectionManager;
}

MainWindow* MainWindow::GetMainWindowPointer()
{
  return wxDynamicCast( wxTheApp->GetTopWindow(), MainWindow );
}

void MainWindow::OnClose( wxCloseEvent &event )
{
  if ( IsProcessing() || IsWritingMovieFrames() )
  {
    wxMessageDialog dlg( this, _("There is on-going data processing. If you force quit, any data that is being saved can be lost. Do you really want to force quit?"), 
                         _("Force Quit"), wxYES_NO | wxNO_DEFAULT  );
    if (dlg.ShowModal() == wxID_NO )
      return;
  }

  LayerCollection* lc_mri = GetLayerCollection( "MRI" );
  LayerCollection* lc_roi = GetLayerCollection( "ROI" );
  LayerCollection* lc_wp = GetLayerCollection( "WayPoints" );
  LayerCollection* lc_surf = GetLayerCollection( "Surface" );
  wxString text = _( "" );
  for ( int i = 0; i < lc_mri->GetNumberOfLayers(); i++ )
  {
    LayerEditable* layer = ( LayerEditable* )lc_mri->GetLayer( i );
    if ( layer->IsModified() )
      text += wxString::FromAscii( layer->GetName() ) + _( "\t(" ) + 
          wxFileName( wxString::FromAscii( layer->GetFileName() ) ).GetShortPath() + _(")\n");
  }
  for ( int i = 0; i < lc_roi->GetNumberOfLayers(); i++ )
  {
    LayerEditable* layer = ( LayerEditable* )lc_roi->GetLayer( i );
    if ( layer->IsModified() )
      text += wxString::FromAscii( layer->GetName() ) + _( "\t(" ) + 
          wxFileName( wxString::FromAscii( layer->GetFileName() ) ).GetShortPath() + _(")\n");
  }
  for ( int i = 0; i < lc_surf->GetNumberOfLayers(); i++ )
  {
    LayerEditable* layer = ( LayerEditable* )lc_surf->GetLayer( i );
    if ( layer->IsModified() )
      text += wxString::FromAscii( layer->GetName() ) + _( "\t(" ) + 
          wxFileName( wxString::FromAscii( layer->GetFileName() ) ).GetShortPath() + _(")\n");
  }
  for ( int i = 0; i < lc_wp->GetNumberOfLayers(); i++ )
  {
    LayerEditable* layer = ( LayerEditable* )lc_wp->GetLayer( i );
    if ( layer->IsModified() )
      text += wxString::FromAscii( layer->GetName() ) + _( "\t(" ) + 
          wxFileName( wxString::FromAscii( layer->GetFileName() ) ).GetShortPath() + _(")\n");
  }

  if ( !text.IsEmpty() )
  {
    wxString msg = _("The following layer(s) have been modified but not saved. \n\n");
    msg += text + _("\nDo you still want to quit the program?");
    wxMessageDialog dlg( this, msg, _("Quit"), wxYES_NO | wxICON_QUESTION | wxNO_DEFAULT );
    if ( dlg.ShowModal() != wxID_YES )
      return;
  }

  wxConfigBase* config = wxConfigBase::Get();
  if ( config )
  {
    // save the frame position
    int x, y, w, h;
    if ( !IsFullScreen() && !IsIconized() && !IsMaximized() )
    {
      GetPosition( &x, &y );
      GetSize( &w, &h );
      config->Write( _("/MainWindow/PosX"), (long) x );
      config->Write( _("/MainWindow/PosY"), (long) y );
      config->Write( _("/MainWindow/Width"), (long) w );
      config->Write( _("/MainWindow/Height"), (long) h );
    }
    config->Write( _("/MainWindow/FullScreen"), IsFullScreen() );
    if ( m_controlPanel->IsShown() )
      config->Write( _("/MainWindow/SplitterPosition"), m_splitterMain->GetSashPosition() );
    config->Write( _("/MainWindow/SplitterPositionSub"), m_splitterSub->GetSashPosition() );
    config->Write( _("/MainWindow/LastDir"), m_strLastDir );
    config->Write( _("/MainWindow/ViewLayout"), m_nViewLayout );
    config->Write( _("/MainWindow/MainView"), m_nMainView );
    config->Write( _("/MainWindow/ScreenshotFilterIndex"), m_nScreenshotFilterIndex );

    config->Write( _("/RenderWindow/BackgroundColor"),  m_settingsGeneral.BackgroundColor.GetAsString( wxC2S_CSS_SYNTAX ) );
    config->Write( _("/RenderWindow/CursorColor"),      m_settingsGeneral.CursorColor.GetAsString( wxC2S_CSS_SYNTAX ) );
    config->Write( _("/RenderWindow/CursorStyle"),      m_settingsGeneral.CursorStyle );
    config->Write( _("/RenderWindow/SyncZoomFactor"),   m_settings2D.SyncZoomFactor );
    config->Write( _("/Settings/SaveCopy"),             m_settingsGeneral.SaveCopy );
    
    config->Write( _("/RenderWindow/ShowScalarBar"), 
                   ( m_nMainView >= 0 ? m_viewRender[m_nMainView]->GetShowScalarBar() : false ) );

    config->Write( _("/Screenshot/Magnification" ),  m_settingsScreenshot.Magnification );
    config->Write( _("/Screenshot/HideCursor" ),  m_settingsScreenshot.HideCursor );
    config->Write( _("/Screenshot/HideCoords" ),  m_settingsScreenshot.HideCoords );
    config->Write( _("/Screenshot/AntiAliasing" ),  m_settingsScreenshot.AntiAliasing );

    config->Write( _("/MainWindow/ShowControlPanel" ), m_controlPanel->IsShown() );
    if ( m_dlgSavePoint )
    {
      m_strGotoPoints = m_dlgSavePoint->GetGotoPoints();
      config->Write( _("/MainWindow/GotoPoints"), MyUtils::JoinStrings( m_strGotoPoints, ";" ) );
    }
    
    m_fileHistory->Save( *config );
  }

  event.Skip();
}


void MainWindow::OnFileOpen( wxCommandEvent& event )
{
  LoadVolume();
}

void MainWindow::NewVolume()
{
  // first check if there is any volume/MRI layer and if the current one is visible
  LayerCollection* col_mri = m_layerCollectionManager->GetLayerCollection( "MRI" );
  if ( !col_mri->GetActiveLayer() )
  {
    wxMessageDialog dlg( this, _("Can not create new volume without visible template volume."), 
                         _("Error"), wxOK );
    dlg.ShowModal();
    return;
  }

  // enter the name of the new ROI
  DialogNewVolume dlg( this, col_mri );
  dlg.SetVolumeName( _("New volume") );
  if ( dlg.ShowModal() != wxID_OK )
    return;

  // finally we are about to create new volume.
  LayerMRI* layer_new = new LayerMRI( dlg.GetTemplate() );
  if ( !layer_new->Create( dlg.GetTemplate(), dlg.GetCopyVoxel(), dlg.GetDataType() ) )
  {
    wxMessageDialog dlg( this, _("Can not create new volume."), 
                         _("Error"), wxOK );
    dlg.ShowModal();
    delete layer_new;
    return;
  }
  layer_new->SetName( dlg.GetVolumeName().char_str() );
  col_mri->AddLayer( layer_new );

  m_controlPanel->RaisePage( _("Volumes") );
}

wxString MainWindow::AutoSelectLastDir( wxString subdirectory )
{
  return AutoSelectLastDir( m_strLastDir, subdirectory );
}

wxString MainWindow::AutoSelectLastDir( wxString lastDir, wxString subdirectory )
{
  wxFileName fn = wxFileName::FileName( lastDir );
  fn.Normalize();
  wxArrayString dirs = fn.GetDirs();
  if ( dirs.GetCount() < 2 )
    return lastDir;
  
  wxArrayString stock;
  stock.Add( _("mri") );
  stock.Add( _("label") );
  stock.Add( _("scripts") );
  stock.Add( _("surf") );
  stock.Add( _("stats") );
  int nStop = 0;
  for ( size_t i = dirs.GetCount() - 2; i < dirs.GetCount(); i++ )
  {
    for ( size_t j = 0; j < stock.GetCount(); j++ )
    {
      if ( stock[j] == dirs[i] )
      {
        nStop = i;
        break;
      }
    }
  }
  if ( nStop == 0 )
    nStop = dirs.GetCount();
  
  wxString outDir;
  for ( int i = 0; i < nStop; i++ )
  {
    outDir += wxFileName::GetPathSeparator() + dirs[i];
  } 
  outDir += wxFileName::GetPathSeparator() + subdirectory;
  if ( wxFileName::DirExists( outDir ) )
    return outDir;
  else
    return lastDir;
}

void MainWindow::SaveVolume()
{
  // first check if there is any volume/MRI layer and if the current one is visible
  LayerCollection* col_mri = m_layerCollectionManager->GetLayerCollection( "MRI" );
  LayerMRI* layer_mri = ( LayerMRI* )col_mri->GetActiveLayer();
  if ( !layer_mri)
  {
    return;
  }
  else if ( !layer_mri->IsVisible() )
  {
    wxMessageDialog dlg( this, _("Current volume layer is not visible. Please turn it on before saving."), 
                         _("Error"), wxOK );
    dlg.ShowModal();
    return;
  }
  wxString fn = wxString::FromAscii( layer_mri->GetFileName() );
  if ( fn.IsEmpty() )
  {
    wxString name = layer_mri->GetName(); 
    name.Trim( true ).Trim( false ).Replace( _(" "), _("_") );
    wxFileDialog dlg( this, _("Save volume file"), 
                      AutoSelectLastDir( m_strLastDir, _("mri") ),
                      name + _(".mgz"),
                      _("Volume files (*.mgz;*.mgh;*.nii;*.nii.gz;*.img)|*.mgz;*.mgh;*.nii;*.nii.gz;*.img|All files (*.*)|*.*"),
                      wxFD_SAVE | wxFD_OVERWRITE_PROMPT );
    if ( dlg.ShowModal() == wxID_OK )
    {
      fn = dlg.GetPath();
    }
  }

  if ( !fn.IsEmpty() )
  {
    if ( !MyUtils::HasExtension( fn, _("nii") ) &&
          !MyUtils::HasExtension( fn, _("nii.gz") ) &&
          !MyUtils::HasExtension( fn, _("img") ) &&
          !MyUtils::HasExtension( fn, _("mgz") ) &&
          !MyUtils::HasExtension( fn, _("mgh") )
       )
    {
      fn += _(".mgz");
    }
    layer_mri->SetFileName( fn.char_str() );
    WorkerThread* thread = new WorkerThread( this );
    thread->SaveVolume( layer_mri );
  }
}

void MainWindow::SaveVolumeAs()
{
  LayerCollection* col_mri = m_layerCollectionManager->GetLayerCollection( "MRI" );
  LayerMRI* layer_mri = ( LayerMRI* )col_mri->GetActiveLayer();
  if ( !layer_mri)
  {
    return;
  }
  else if ( !layer_mri->IsVisible() )
  {
    wxMessageDialog dlg( this, _( "Current volume layer is not visible. Please turn it on before saving." ), 
                         _( "Error" ), wxOK );
    dlg.ShowModal();
    return;
  }

//  DialogSaveVolumeAs dlg( this );
//  dlg.SetFileName( layer_mri->GetFileName() );
  wxFileDialog dlg( this, _("Select file to save"), 
                    wxFileName( layer_mri->GetFileName() ).GetPath(), 
                    _(""),
                    _("Volume files (*.mgz;*.mgh;*.nii;*.nii.gz;*.img)|*.mgz;*.mgh;*.nii;*.nii.gz;*.img|All files (*.*)|*.*"),
                    wxFD_SAVE | wxFD_OVERWRITE_PROMPT );
  if ( dlg.ShowModal() == wxID_OK )
  {
    layer_mri->SetFileName( dlg.GetPath().char_str() );
    SaveVolume();
    m_controlPanel->UpdateUI();
  }
}

void MainWindow::LoadVolume()
{
  bool bHasVolume = !GetLayerCollection( "MRI" )->IsEmpty();
  bool bHasSurface = !GetLayerCollection( "Surface" )->IsEmpty();
  DialogLoadVolume dlg( this, !(bHasVolume || bHasSurface) );
  dlg.SetLastDir( AutoSelectLastDir( m_strLastDir, _("mri") ) );
  wxArrayString list;
  for ( size_t i = 0; i < m_fileHistory->GetCount(); i++ )
    list.Add( m_fileHistory->GetHistoryFile( i ) );
  dlg.SetRecentFiles( list );
  if ( dlg.ShowModal() == wxID_OK )
  {
//    this->LoadVolumeFile( dlg.GetVolumeFileName(), dlg.GetRegFileName(),
//                          ( GetLayerCollection( "MRI" )->IsEmpty() ? dlg.IsToResample() : m_bResampleToRAS ),
//                          dlg.GetSampleMethod() );
    
    wxArrayString fns = dlg.GetVolumeFileNames();
    wxString reg_fn = dlg.GetRegFileName();
    for ( size_t i = 0; i < fns.GetCount(); i++ )
    {
      wxArrayString script;
      script.Add( _("loadvolume") );
      wxString fn = fns[i];
      if ( !reg_fn.IsEmpty() )
        fn += ":reg=" + reg_fn;
      
      if ( dlg.GetSampleMethod() != SAMPLE_NEAREST )
        fn += ":sample=trilinear";
      
      fn += ":colormap=" + dlg.GetColorMap();
      if ( dlg.GetColorMap() == "lut" )
        fn += ":lut=" + dlg.GetLUT();
      
      script.Add( fn );
  
      if ( (!bHasVolume && bHasSurface) || (!bHasVolume && dlg.IsToResample()) || m_bResampleToRAS )
        script.Add( _("r") );
  
      AddScript( script );
    }
    RunScript();
  }
}

void MainWindow::LoadVolumeFile( const wxString& filename, 
                                 const wxString& reg_filename, 
                                 bool bResample, int nSampleMethod,
                                 bool bConform )
{
// cout << bResample << endl;
  m_strLastDir = MyUtils::GetNormalizedPath( filename );

  m_bResampleToRAS = bResample;
  LayerMRI* layer = new LayerMRI( m_layerVolumeRef );
  layer->SetResampleToRAS( bResample );
  layer->SetSampleMethod( nSampleMethod );
  layer->SetConform( bConform );
  layer->GetProperties()->SetLUTCTAB( m_luts->GetColorTable( 0 ) );
  wxFileName fn( filename );
  fn.Normalize();
  wxString layerName = fn.GetName();
  if ( fn.GetExt().Lower() == _("gz") )
    layerName = wxFileName( layerName ).GetName();
  layer->SetName( layerName.char_str() );
  layer->SetFileName( fn.GetFullPath().char_str() );
  if ( !reg_filename.IsEmpty() )
  {
    wxFileName reg_fn( reg_filename );
    reg_fn.Normalize( wxPATH_NORM_ALL, fn.GetPath() );
    layer->SetRegFileName( reg_fn.GetFullPath().char_str() );
  }

// if ( !bResample )
  /* 
  {
    LayerMRI* mri = (LayerMRI* )GetLayerCollection( "MRI" )->GetLayer( 0 );
    if ( mri )
    {
     layer->SetRefVolume( mri->GetSourceVolume() );
    }
   }
  */
  WorkerThread* thread = new WorkerThread( this );
  thread->LoadVolume( layer );
}


void MainWindow::RotateVolume( std::vector<RotationElement>& rotations, bool bAllVolumes )
{
//  WorkerThread* thread = new WorkerThread( this );
//  thread->RotateVolume( rotations, bAllVolumes );
  
  // rotation is much faster now so no need to do it in a separate thread 
  wxCommandEvent event( 0, 0 ); // fake event, just a place holder to call the functions
  std::vector<Layer*> layers = GetLayerCollectionManager()->GetAllLayers();
      
  // first update ROI and waypoints before their reference volume is rotated
  bool bSuccess = true;
  for ( size_t i = 0; i < layers.size(); i++ )
  {
    if ( layers[i]->IsTypeOf( "ROI" ) )
    {
      ( (LayerROI*)layers[i] )->UpdateLabelData( this, event );
    }
    else if ( layers[i]->IsTypeOf( "WayPoints" ) )
    {
      ( (LayerWayPoints*)layers[i] )->UpdateLabelData();
    }  
  }
      
  if ( bAllVolumes )
  {
    // then rotate MRI volumes
    for ( size_t i = 0; i < layers.size(); i++ )
    {
      if ( layers[i]->IsTypeOf( "MRI" ) && !layers[i]->Rotate( rotations, this, event ) )
      {
        bSuccess = false;
        break;
      }
    }
        // at last rotate others
    for ( size_t i = 0; i < layers.size() && bSuccess; i++ )
    {
      if ( !layers[i]->IsTypeOf( "MRI" ) && !layers[i]->Rotate( rotations, this, event ) )
      {
        bSuccess = false;
        break;
      }
    }
  }
  else
  {
    LayerMRI* layer = (LayerMRI*) GetActiveLayer( "MRI" );
    if ( !layer->Rotate( rotations, this, event ) )
    {
      bSuccess = false;
    }
  }
  if ( !bSuccess )
  {
    wxMessageDialog dlg( this, _("Error occured while rotating volumes."), 
                         _("Error"), wxOK );
    dlg.ShowModal(); 
  }
}


void MainWindow::OnFileExit( wxCommandEvent& event )
{
  Close();
}

void MainWindow::OnActivate( wxActivateEvent& event )
{
#if defined(__WXGTK__) || defined(__WXMAC__)
  NeedRedraw( 2 );
#endif
  event.Skip();
}

void MainWindow::OnIconize( wxIconizeEvent& event )
{
#if defined(__WXGTK__) || defined(__WXMAC__)
//#if wxCHECK_VERSION(2,9,0)
#if wxVERSION_NUMBER > 2900  
  if ( !event.IsIconized() )
#else
  if ( !event.Iconized() )
#endif
    NeedRedraw( 2 );
#endif
  event.Skip();
}


void MainWindow::OnFileRecent( wxCommandEvent& event )
{
  wxString fn( m_fileHistory->GetHistoryFile( event.GetId() - wxID_FILE1 ) );
  if ( !fn.IsEmpty() )
  {
    bool bResample = m_bResampleToRAS;
    if ( !GetLayerCollection( "Surface" )->IsEmpty() && GetLayerCollection( "MRI" )->IsEmpty() )
      bResample = true;
    this->LoadVolumeFile( fn, _(""), bResample, m_nDefaultSampleMethod, m_bDefaultConform );
  }
}

LayerCollection* MainWindow::GetLayerCollection( std::string strType )
{
  return m_layerCollectionManager->GetLayerCollection( strType );
}

Layer* MainWindow::GetActiveLayer( std::string strType )
{
  return GetLayerCollection( strType )->GetActiveLayer();
}

LayerCollectionManager* MainWindow::GetLayerCollectionManager()
{
  return m_layerCollectionManager;
}

void MainWindow::LoadROI()
{
  if ( GetLayerCollection( "MRI" )->IsEmpty() )
  {
    return;
  }
  DialogLoadROI dlg( this );
  dlg.SetLastDir( AutoSelectLastDir( m_strLastDir, _("label") ) );
  if ( dlg.ShowModal() == wxID_OK )
  {
    wxArrayString fns = dlg.GetFileNames();
    for ( size_t i = 0; i < fns.GetCount(); i++ )
    {
      wxArrayString script;
      script.Add( _("loadroi") );
      script.Add( fns[i] + _(":ref=") + dlg.GetTemplate()  );
      this->AddScript( script );
    }
    ContinueScripts();
  }
}

void MainWindow::LoadROIFile( const wxString& fn, const wxString& ref_vol )
{
  m_strLastDir = MyUtils::GetNormalizedPath( fn );

  LayerMRI* ref = NULL;
  LayerCollection* col_mri = GetLayerCollection( "MRI" );
  if ( ref_vol.IsEmpty() )
  {
    cout << "No template volume given, using current volume as template for ROI " << fn.c_str() << endl;
    ref = (LayerMRI*)col_mri->GetActiveLayer();
  }
  else
  {
    for ( int i = 0; i < col_mri->GetNumberOfLayers(); i++ )
    {
      LayerMRI* mri = ( LayerMRI* )col_mri->GetLayer( i );
      if ( ref_vol == mri->GetName() )
      {
        ref = mri;
        break;
      }
      else if ( wxFileName( mri->GetFileName() ).GetFullName() == ref_vol )
      {
        ref = mri;
        break;
      }
    }
    if ( ref == NULL )
    {
      cerr << "Can not find given template volume: " << ref_vol.c_str() 
          << ". Using current volume as template for ROI " << fn.c_str() << endl;
      ref = (LayerMRI*)col_mri->GetActiveLayer();
    }
  }
  LayerROI* roi = new LayerROI( ref );
  roi->SetName( wxFileName( fn ).GetName().char_str() );
  if ( roi->LoadROIFromFile( (const char*)fn.char_str() ) )
  {
    LayerCollection* col_roi = GetLayerCollection( "ROI" );
    if ( col_roi->IsEmpty() )
    {
      col_roi->SetWorldOrigin( col_mri->GetWorldOrigin() );
      col_roi->SetWorldSize( col_mri->GetWorldSize() );
      col_roi->SetWorldVoxelSize( col_mri->GetWorldVoxelSize() );
      col_roi->SetSlicePosition( col_mri->GetSlicePosition() );
    }
    col_roi->AddLayer( roi );

    m_controlPanel->RaisePage( _("ROIs") );
  }
  else
  {
    delete roi;
    wxMessageDialog dlg( this, _( "Can not load ROI from " ) + fn, _( "Error" ), wxOK );
    dlg.ShowModal();
  }
}

void MainWindow::NewROI()
{
  // first check if there is any volume/MRI layer and if the current one is visible
  LayerCollection* col_mri = m_layerCollectionManager->GetLayerCollection( "MRI" );
  LayerMRI* layer_mri = ( LayerMRI* )col_mri->GetActiveLayer();
  if ( !layer_mri)
  {
    wxMessageDialog dlg( this, _("Can not create new ROI without volume image."), _("Error"), wxOK );
    dlg.ShowModal();
    return;
  }

  // enter the name of the new ROI
  DialogNewROI dlg( this, col_mri );
  dlg.SetROIName( _("New ROI") );
  if ( dlg.ShowModal() != wxID_OK )
    return;

  // finally we are about to create new ROI.
  LayerCollection* col_roi = m_layerCollectionManager->GetLayerCollection( "ROI" );
  if ( col_roi->IsEmpty() )
  {
    col_roi->SetWorldOrigin( col_mri->GetWorldOrigin() );
    col_roi->SetWorldSize( col_mri->GetWorldSize() );
    col_roi->SetWorldVoxelSize( col_mri->GetWorldVoxelSize() );
    col_roi->SetSlicePosition( col_mri->GetSlicePosition() );
  }
  LayerROI* layer_roi = new LayerROI( dlg.GetTemplate() );
  layer_roi->SetName( dlg.GetROIName().char_str() );
  col_roi->AddLayer( layer_roi );

  m_controlPanel->RaisePage( _("ROIs") );

  SetMode( RenderView2D::IM_ROIEdit );
}

void MainWindow::SaveROI()
{
  LayerCollection* col_roi = m_layerCollectionManager->GetLayerCollection( "ROI" );
  LayerROI* layer_roi = ( LayerROI* )col_roi->GetActiveLayer();
  if ( !layer_roi )
  {
    return;
  }
  else if ( !layer_roi->IsVisible() )
  {
    wxMessageDialog dlg( this, _("Current ROI layer is not visible. Please turn it on before saving."), _("Error"), wxOK );
    dlg.ShowModal();
    return;
  }
  wxString fn = wxString::FromAscii( layer_roi->GetFileName() );
  if ( fn.IsEmpty() )
  {
    wxFileDialog dlg( this, _("Save ROI file"), 
                      AutoSelectLastDir( m_strLastDir, _("label") ),
                      _(""),
                      _("Label files (*.label)|*.label|All files (*.*)|*.*"),
                      wxFD_SAVE | wxFD_OVERWRITE_PROMPT );
    if ( dlg.ShowModal() == wxID_OK )
    {
      fn = dlg.GetPath();
    }
  }

  if ( !fn.IsEmpty() )
  {
    if ( !MyUtils::HasExtension( fn, _("label") ) )
    {
      fn += _(".label");
    }
    layer_roi->SetFileName( fn.char_str() );
    layer_roi->ResetModified();
    WorkerThread* thread = new WorkerThread( this );
    thread->SaveROI( layer_roi );
  }
}


void MainWindow::SaveROIAs()
{
  LayerCollection* col_roi = m_layerCollectionManager->GetLayerCollection( "ROI" );
  LayerROI* layer_roi = ( LayerROI* )col_roi->GetActiveLayer();
  if ( !layer_roi )
  {
    return;
  }
  else if ( !layer_roi->IsVisible() )
  {
    wxMessageDialog dlg( this, _("Current ROI layer is not visible. Please turn it on before saving."), 
                         _("Error"), wxOK );
    dlg.ShowModal();
    return;
  }

  wxString fn = wxString::FromAscii( layer_roi->GetFileName() );
  wxFileDialog dlg( this, _("Save ROI file as"), 
                    AutoSelectLastDir( m_strLastDir, _("label") ),
                    fn,
                    _("Label files (*.label)|*.label|All files (*.*)|*.*"),
                    wxFD_SAVE | wxFD_OVERWRITE_PROMPT );
  if ( dlg.ShowModal() == wxID_OK )
  {
    layer_roi->SetFileName( dlg.GetPath().char_str() );
    SaveROI();
    m_controlPanel->UpdateUI();
  }
}

void MainWindow::OnFileNewROI( wxCommandEvent& event )
{
  NewROI();
}

void MainWindow::OnFileNewROIUpdateUI( wxUpdateUIEvent& event )
{
  event.Enable( !GetLayerCollection( "MRI" )->IsEmpty() );
}


void MainWindow::OnFileNewWayPoints( wxCommandEvent& event )
{
  NewWayPoints();
}

void MainWindow::OnFileNewWayPointsUpdateUI( wxUpdateUIEvent& event )
{
  event.Enable( !GetLayerCollection( "MRI" )->IsEmpty() );
}

/*
void MainWindow::LoadWayPoints()
{
  if ( GetLayerCollection( "MRI" )->IsEmpty() )
  {
    return;
  }
  wxFileDialog dlg( this, _("Open Way Points file"), m_strLastDir, _(""),
                    _("Way Points files (*.label)|*.label|All files (*.*)|*.*"),
                    wxFD_OPEN );
  if ( dlg.ShowModal() == wxID_OK )
  {
    this->LoadWayPointsFile( dlg.GetPath() );
  }
}
*/

void MainWindow::LoadWayPoints()
{
  if ( GetLayerCollection( "MRI" )->IsEmpty() )
  {
    return;
  }
  DialogLoadPointSet dlg( this );
  if ( dlg.ShowModal() == wxID_OK )
  {
    int nType = dlg.GetPointSetType();
    if ( nType == -1 )  // auto
    {
      if ( FSWayPoints::IsLabelFormat( dlg.GetFileName().c_str() ) )
        nType = LayerPropertiesWayPoints::WayPoints;
      else 
        nType = LayerPropertiesWayPoints::ControlPoints;
    }
    if ( nType == LayerPropertiesWayPoints::WayPoints )
      this->LoadWayPointsFile( dlg.GetFileName() );
    else if ( nType == LayerPropertiesWayPoints::ControlPoints )
      this->LoadControlPointsFile( dlg.GetFileName() );
  }
}

void MainWindow::LoadWayPointsFile( const wxString& fn )
{
  this->LoadPointSetFile( fn, LayerPropertiesWayPoints::WayPoints );
}

void MainWindow::LoadControlPointsFile( const wxString& fn )
{
  this->LoadPointSetFile( fn, LayerPropertiesWayPoints::ControlPoints );
}

void MainWindow::LoadPointSetFile( const wxString& fn, int type )
{
  m_strLastDir = MyUtils::GetNormalizedPath( fn );

  LayerCollection* col_mri = GetLayerCollection( "MRI" );
  LayerMRI* mri = ( LayerMRI* )col_mri->GetActiveLayer();
  LayerWayPoints* wp = new LayerWayPoints( mri, type );
  wp->SetName( wxFileName( fn ).GetName().char_str() );
  if ( wp->LoadFromFile( (const char*)fn.char_str() ) )
  {
    LayerCollection* col_wp = GetLayerCollection( "WayPoints" );
    if ( col_wp->IsEmpty() )
    {
      col_wp->SetWorldOrigin( col_mri->GetWorldOrigin() );
      col_wp->SetWorldSize( col_mri->GetWorldSize() );
      col_wp->SetWorldVoxelSize( col_mri->GetWorldVoxelSize() );
      col_wp->SetSlicePosition( col_mri->GetSlicePosition() );
    }
    col_wp->AddLayer( wp );

    m_controlPanel->RaisePage( _("Point Sets") );
  }
  else
  {
    delete wp;
    wxMessageDialog dlg( this, _( "Can not load Way Points from " ) + fn, _("Error"), wxOK );
    dlg.ShowModal();
  }
}

void MainWindow::NewWayPoints()
{
  // first check if there is any volume/MRI layer and if the current one is visible
  LayerCollection* col_mri = m_layerCollectionManager->GetLayerCollection( "MRI" );
  LayerMRI* layer_mri = ( LayerMRI* )col_mri->GetActiveLayer();
  if ( !layer_mri)
  {
    wxMessageDialog dlg( this, _("Can not create new way points without volume image."), _("Error"), wxOK );
    dlg.ShowModal();
    return;
  }

  // enter the name of the new point set
  DialogNewWayPoints dlg( this, col_mri );
  dlg.SetWayPointsName( _("New Point Set") );
  if ( dlg.ShowModal() != wxID_OK )
    return;

  // finally we are about to create new point set.
  LayerCollection* col_wp = m_layerCollectionManager->GetLayerCollection( "WayPoints" );
  if ( col_wp->IsEmpty() )
  {
    col_wp->SetWorldOrigin( col_mri->GetWorldOrigin() );
    col_wp->SetWorldSize( col_mri->GetWorldSize() );
    col_wp->SetWorldVoxelSize( col_mri->GetWorldVoxelSize() );
    col_wp->SetSlicePosition( col_mri->GetSlicePosition() );
  }
  LayerWayPoints* layer_wp = new LayerWayPoints( dlg.GetTemplate(), dlg.GetType() );
  layer_wp->SetName( dlg.GetWayPointsName().char_str() );
  col_wp->AddLayer( layer_wp );

  m_controlPanel->RaisePage( _("Point Sets") );

  SetMode( RenderView2D::IM_WayPointsEdit );
}

void MainWindow::SaveWayPoints()
{
  LayerCollection* col_wp = m_layerCollectionManager->GetLayerCollection( "WayPoints" );
  LayerWayPoints* layer_wp = ( LayerWayPoints* )col_wp->GetActiveLayer();
  if ( !layer_wp )
  {
    return;
  }
  else if ( !layer_wp->IsVisible() )
  {
    wxMessageDialog dlg( this, _("Current Way Points layer is not visible. Please turn it on before saving."), 
                         _("Error"), wxOK );
    dlg.ShowModal();
    return;
  }
  wxString fn = wxString::FromAscii( layer_wp->GetFileName() );
  if ( fn.IsEmpty() )
  {
    /*
    wxFileDialog dlg( this, _("Save Way Points file"), m_strLastDir, _(""),
                      _("Way Points files (*.label)|*.label|All files (*.*)|*.*"),
                      wxFD_SAVE | wxFD_OVERWRITE_PROMPT );
    if ( dlg.ShowModal() == wxID_OK )
    {
      fn = dlg.GetPath();
    }
    */
    SaveWayPointsAs();
  }
  else
  {
    if ( layer_wp->GetProperties()->GetType() == LayerPropertiesWayPoints::WayPoints &&
         !MyUtils::HasExtension( fn, _("label") ) )
    {
      fn += _(".label");
    }
    if ( layer_wp->GetProperties()->GetType() == LayerPropertiesWayPoints::ControlPoints &&
         !MyUtils::HasExtension( fn, _("dat") ) )
    {
      fn += _(".dat");
    }
    layer_wp->SetFileName( fn.char_str() );
    layer_wp->ResetModified();
    WorkerThread* thread = new WorkerThread( this );
    thread->SaveWayPoints( layer_wp );
  }
}

void MainWindow::SaveWayPointsAs()
{
  LayerCollection* col_wp = m_layerCollectionManager->GetLayerCollection( "WayPoints" );
  LayerWayPoints* layer_wp = ( LayerWayPoints* )col_wp->GetActiveLayer();
  if ( !layer_wp )
  {
    return;
  }
  else if ( !layer_wp->IsVisible() )
  {
    wxMessageDialog dlg( this, _("Current Point Set layer is not visible. Please turn it on before saving."), 
                         _("Error"), wxOK );
    dlg.ShowModal();
    return;
  }

  wxString fn = wxString::FromAscii( layer_wp->GetFileName() );
  DialogSavePointSetAs dlg( this );
  dlg.SetPointSetType( layer_wp->GetProperties()->GetType() );
  if ( dlg.ShowModal() == wxID_OK )
  {
    layer_wp->SetFileName( dlg.GetFileName().char_str() );
    layer_wp->GetProperties()->SetType( dlg.GetPointSetType() );
    SaveWayPoints();
    m_controlPanel->UpdateUI();
  }
}

void MainWindow::OnKeyDown( wxKeyEvent& event )
{
// cout << "test" << endl;
// m_viewAxial->GetEventHandler()->ProcessEvent( event );

  event.Skip();
}

void MainWindow::UpdateToolbars()
{
  m_bToUpdateToolbars = true;
}

void MainWindow::DoUpdateToolbars()
{
  if ( !IsShown() )
    return;
  
  if ( !m_toolWindowEdit )
    m_toolWindowEdit = new ToolWindowEdit( this );

  if ( !m_toolWindowMeasure )
  {
    m_toolWindowMeasure = new ToolWindowMeasure( this );
    for ( int i = 0; i < 4; i++ )
      m_viewRender[i]->AddListener( m_toolWindowMeasure );
  }
    
  m_toolWindowEdit->Show( m_viewAxial->GetInteractionMode() == RenderView::IM_VoxelEdit ||
      m_viewAxial->GetInteractionMode() == RenderView::IM_ROIEdit );
  
  m_toolWindowMeasure->Show( m_viewAxial->GetInteractionMode() == RenderView::IM_Measure );

  m_toolWindowEdit->UpdateTools();
  
  m_dlgCropVolume->Show( m_viewAxial->GetInteractionMode() == RenderView::IM_VolumeCrop );

  m_bToUpdateToolbars = false;
}

void MainWindow::SetMode( int nMode )
{
  m_viewAxial->SetInteractionMode( nMode );
  m_viewCoronal->SetInteractionMode( nMode );
  m_viewSagittal->SetInteractionMode( nMode );
  m_view3D->SetInteractionMode( nMode );

  UpdateToolbars();
}

void MainWindow::OnModeNavigate( wxCommandEvent& event )
{
  SetMode( RenderView2D::IM_Navigate );
}

void MainWindow::OnModeNavigateUpdateUI( wxUpdateUIEvent& event)
{
  event.Check( m_viewAxial->GetInteractionMode() == RenderView2D::IM_Navigate );
}

void MainWindow::OnModeMeasure( wxCommandEvent& event )
{
  SetMode( RenderView2D::IM_Measure );
}

void MainWindow::OnModeMeasureUpdateUI( wxUpdateUIEvent& event)
{
  event.Check( m_viewAxial->GetInteractionMode() == RenderView2D::IM_Measure );
  event.Enable( m_layerCollectionManager->HasLayer( "MRI" ) );
}

void MainWindow::OnModeVoxelEdit( wxCommandEvent& event )
{
  SetMode( RenderView2D::IM_VoxelEdit );
}

void MainWindow::OnModeVoxelEditUpdateUI( wxUpdateUIEvent& event)
{
  event.Check( m_viewAxial->GetInteractionMode() == RenderView2D::IM_VoxelEdit );
  event.Enable( m_layerCollectionManager->HasLayer( "MRI" ) );
}

void MainWindow::OnModeROIEdit( wxCommandEvent& event )
{
  SetMode( RenderView2D::IM_ROIEdit );
}

void MainWindow::OnModeROIEditUpdateUI( wxUpdateUIEvent& event)
{
  event.Check( m_viewAxial->GetInteractionMode() == RenderView2D::IM_ROIEdit );
  event.Enable( m_layerCollectionManager->HasLayer( "ROI" ) );
}

void MainWindow::OnModeWayPointsEdit( wxCommandEvent& event )
{
  SetMode( RenderView2D::IM_WayPointsEdit );
}

void MainWindow::OnModeWayPointsEditUpdateUI( wxUpdateUIEvent& event)
{
  event.Check( m_viewAxial->GetInteractionMode() == RenderView2D::IM_WayPointsEdit );
  event.Enable( m_layerCollectionManager->HasLayer( "WayPoints" ) );
}

void MainWindow::SetAction( int nAction )
{
  m_viewAxial->SetAction( nAction );
  m_viewCoronal->SetAction( nAction );
  m_viewSagittal->SetAction( nAction );
  m_view3D->SetAction( nAction );
  UpdateToolbars();
}

void MainWindow::OnEditUndo( wxCommandEvent& event )
{
  if ( m_viewAxial->GetInteractionMode() == RenderView2D::IM_ROIEdit )
  {
    LayerROI* roi = ( LayerROI* )GetLayerCollection( "ROI" )->GetActiveLayer();
    if ( roi )
      roi->Undo();
  }
  else if ( m_viewAxial->GetInteractionMode() == RenderView2D::IM_VoxelEdit )
  {
    LayerMRI* mri = ( LayerMRI* )GetLayerCollection( "MRI" )->GetActiveLayer();
    if ( mri )
      mri->Undo();
  }
  else if ( m_viewAxial->GetInteractionMode() == RenderView2D::IM_WayPointsEdit )
  {
    LayerWayPoints* wp = ( LayerWayPoints* )GetLayerCollection( "WayPoints" )->GetActiveLayer();
    if ( wp )
      wp->Undo();
  }
}

void MainWindow::OnEditUndoUpdateUI( wxUpdateUIEvent& event )
{
  if ( m_viewAxial->GetInteractionMode() == RenderView2D::IM_ROIEdit )
  {
    LayerROI* roi = ( LayerROI* )GetLayerCollection( "ROI" )->GetActiveLayer();
    event.Enable( roi && roi->IsVisible() && roi->HasUndo() );
  }
  else if ( m_viewAxial->GetInteractionMode() == RenderView2D::IM_VoxelEdit )
  {
    LayerMRI* mri = ( LayerMRI* )GetLayerCollection( "MRI" )->GetActiveLayer();
    event.Enable( mri && mri->IsVisible() && mri->HasUndo() );
  }
  else if ( m_viewAxial->GetInteractionMode() == RenderView2D::IM_WayPointsEdit )
  {
    LayerWayPoints* wp = ( LayerWayPoints* )GetLayerCollection( "WayPoints" )->GetActiveLayer();
    event.Enable( wp && wp->IsVisible() && wp->HasUndo() );
  }
  else
    event.Enable( false );
}

void MainWindow::OnEditRedo( wxCommandEvent& event )
{
  if ( m_viewAxial->GetInteractionMode() == RenderView2D::IM_ROIEdit )
  {
    LayerROI* roi = ( LayerROI* )GetLayerCollection( "ROI" )->GetActiveLayer();
    if ( roi )
      roi->Redo();
  }
  else if ( m_viewAxial->GetInteractionMode() == RenderView2D::IM_VoxelEdit )
  {
    LayerMRI* mri = ( LayerMRI* )GetLayerCollection( "MRI" )->GetActiveLayer();
    if ( mri )
      mri->Redo();
  }
  else if ( m_viewAxial->GetInteractionMode() == RenderView2D::IM_WayPointsEdit )
  {
    LayerWayPoints* wp = ( LayerWayPoints* )GetLayerCollection( "WayPoints" )->GetActiveLayer();
    if ( wp )
      wp->Redo();
  }
}

void MainWindow::OnEditRedoUpdateUI( wxUpdateUIEvent& event )
{
  if ( m_viewAxial->GetInteractionMode() == RenderView2D::IM_ROIEdit )
  {
    LayerROI* roi = ( LayerROI* )GetLayerCollection( "ROI" )->GetActiveLayer();
    event.Enable( roi && roi->IsVisible() && roi->HasRedo() );
  }
  else if ( m_viewAxial->GetInteractionMode() == RenderView2D::IM_VoxelEdit )
  {
    LayerMRI* mri = ( LayerMRI* )GetLayerCollection( "MRI" )->GetActiveLayer();
    event.Enable( mri && mri->IsVisible() && mri->HasRedo() );
  }
  else if ( m_viewAxial->GetInteractionMode() == RenderView2D::IM_WayPointsEdit )
  {
    LayerWayPoints* wp = ( LayerWayPoints* )GetLayerCollection( "WayPoints" )->GetActiveLayer();
    event.Enable( wp && wp->IsVisible() && wp->HasRedo() );
  }
  else
    event.Enable( false );
}


void MainWindow::OnEditCopy( wxCommandEvent& event )
{
  int nWnd = GetActiveViewId();
  if ( m_viewAxial->GetInteractionMode() == RenderView2D::IM_ROIEdit )
  {
    LayerROI* roi = ( LayerROI* )GetLayerCollection( "ROI" )->GetActiveLayer();
    if ( roi && nWnd >= 0 && nWnd < 3 )
      roi->Copy( nWnd );
  }
  else if ( m_viewAxial->GetInteractionMode() == RenderView2D::IM_VoxelEdit )
  {
    LayerMRI* mri = ( LayerMRI* )GetLayerCollection( "MRI" )->GetActiveLayer();
    if ( mri && nWnd >= 0 && nWnd < 3 )
      mri->Copy( nWnd );
  }
}

void MainWindow::OnEditCopyUpdateUI( wxUpdateUIEvent& event )
{
  int nWnd = GetActiveViewId();
  if ( m_viewAxial->GetInteractionMode() == RenderView2D::IM_ROIEdit )
  {
    LayerROI* roi = ( LayerROI* )GetLayerCollection( "ROI" )->GetActiveLayer();
    event.Enable( roi && roi->IsVisible() && nWnd >= 0 && nWnd < 3 );
  }
  else if ( m_viewAxial->GetInteractionMode() == RenderView2D::IM_VoxelEdit )
  {
    LayerMRI* mri = ( LayerMRI* )GetLayerCollection( "MRI" )->GetActiveLayer();
    event.Enable( mri && mri->IsVisible() && nWnd >= 0 && nWnd < 3 );
  }
  else
    event.Enable( false );
}


void MainWindow::OnEditCopyStructure( wxCommandEvent& event )
{
  int nWnd = GetActiveViewId();
  if ( m_viewAxial->GetInteractionMode() == RenderView2D::IM_VoxelEdit )
  {
    LayerMRI* mri = ( LayerMRI* )GetLayerCollection( "MRI" )->GetActiveLayer();
    if ( mri && nWnd >= 0 && nWnd < 3 )
    {
      double* pos = mri->GetSlicePosition();
      if ( !mri->CopyStructure( nWnd, pos ) )
      {
        wxMessageDialog dlg( this, 
                             _( "Please move the cursor to the structure you want to copy and try again." ), 
                             _( "Copy Structure Failed" ),
                             wxOK );
        dlg.ShowModal();                            
      }
    }
  }
}

void MainWindow::OnEditCopyStructureUpdateUI( wxUpdateUIEvent& event )
{
  int nWnd = GetActiveViewId();
  if ( m_viewAxial->GetInteractionMode() == RenderView2D::IM_VoxelEdit )
  {
    LayerMRI* mri = ( LayerMRI* )GetLayerCollection( "MRI" )->GetActiveLayer();
    event.Enable( mri && mri->IsVisible() && nWnd >= 0 && nWnd < 3 );
  }
  else
    event.Enable( false );
}

void MainWindow::OnEditPaste( wxCommandEvent& event )
{
  int nWnd = GetActiveViewId();
  if ( m_viewAxial->GetInteractionMode() == RenderView2D::IM_ROIEdit )
  {
    LayerROI* roi = ( LayerROI* )GetLayerCollection( "ROI" )->GetActiveLayer();
    if ( roi && nWnd >= 0 && nWnd < 3 )
      roi->Paste( nWnd );
  }
  else if ( m_viewAxial->GetInteractionMode() == RenderView2D::IM_VoxelEdit )
  {
    LayerMRI* mri = ( LayerMRI* )GetLayerCollection( "MRI" )->GetActiveLayer();
    if ( mri && nWnd >= 0 && nWnd < 3 )
      mri->Paste( nWnd );
  }
}

void MainWindow::OnEditPasteUpdateUI( wxUpdateUIEvent& event )
{
  int nWnd = GetActiveViewId();
  if ( m_viewAxial->GetInteractionMode() == RenderView2D::IM_ROIEdit )
  {
    LayerROI* roi = ( LayerROI* )GetLayerCollection( "ROI" )->GetActiveLayer();
    event.Enable( roi && roi->IsVisible() && roi->IsEditable() && nWnd >= 0 && nWnd < 3 && roi->IsValidToPaste( nWnd ) );
  }
  else if ( m_viewAxial->GetInteractionMode() == RenderView2D::IM_VoxelEdit )
  {
    LayerMRI* mri = ( LayerMRI* )GetLayerCollection( "MRI" )->GetActiveLayer();
    event.Enable( mri && mri->IsVisible() && mri->IsEditable() && nWnd >= 0 && nWnd < 3 && mri->IsValidToPaste( nWnd ) );
  }
  else
    event.Enable( false );
}

void MainWindow::OnViewToggleVoxelCoordinates( wxCommandEvent& event )
{
  m_pixelInfoPanel->ToggleShowVoxelCoordinates();
}


void MainWindow::SetViewLayout( int nLayout )
{
  m_nViewLayout = nLayout;

  RenderView* view[4] = { 0 };
  switch ( m_nMainView )
  {
  case MV_Coronal:
    view[0] = m_viewCoronal;
    view[1] = m_viewSagittal;
    view[2] = m_viewAxial;
    view[3] = m_view3D;
    break;
  case MV_Axial:
    view[0] = m_viewAxial;
    view[1] = m_viewSagittal;
    view[2] = m_viewCoronal;
    view[3] = m_view3D;
    break;
  case MV_3D:
    view[0] = m_view3D;
    view[1] = m_viewSagittal;
    view[2] = m_viewCoronal;
    view[3] = m_viewAxial;
    break;
  default:
    view[0] = m_viewSagittal;
    view[1] = m_viewCoronal;
    view[2] = m_viewAxial;
    view[3] = m_view3D;
    break;
  }

  if ( m_nViewLayout == VL_1X1 )
  {
    view[0]->Show();
    view[1]->Hide();
    view[2]->Hide();
    view[3]->Hide();
    wxBoxSizer* sizer = new wxBoxSizer( wxVERTICAL );
    m_renderViewHolder->SetSizer( sizer );
    sizer->Add( view[0], 1, wxEXPAND );
  }
  else if ( m_nViewLayout == VL_1N3 )
  {
    for ( int i = 0; i < 4; i++ )
      view[i]->Show();

    wxBoxSizer* sizer = new wxBoxSizer( wxVERTICAL );
    m_renderViewHolder->SetSizer( sizer );

    sizer->Add( view[0], 2, wxEXPAND );
    sizer->AddSpacer( 1 );
    wxBoxSizer* sizer2 = new wxBoxSizer( wxHORIZONTAL );
    sizer->Add( sizer2, 1, wxEXPAND );
    sizer2->Add( view[1], 1, wxEXPAND );
    sizer2->AddSpacer( 1 );
    sizer2->Add( view[2], 1, wxEXPAND );
    sizer2->AddSpacer( 1 );
    sizer2->Add( view[3], 1, wxEXPAND );
  }
  else if ( m_nViewLayout == VL_1N3_H )
  {
    for ( int i = 0; i < 4; i++ )
      view[i]->Show();

    wxBoxSizer* sizer = new wxBoxSizer( wxHORIZONTAL );
    m_renderViewHolder->SetSizer( sizer );

    sizer->Add( view[0], 2, wxEXPAND );
    sizer->AddSpacer( 1 );
    wxBoxSizer* sizer2 = new wxBoxSizer( wxVERTICAL );
    sizer->Add( sizer2, 1, wxEXPAND );
    sizer2->Add( view[1], 1, wxEXPAND );
    sizer2->AddSpacer( 1 );
    sizer2->Add( view[2], 1, wxEXPAND );
    sizer2->AddSpacer( 1 );
    sizer2->Add( view[3], 1, wxEXPAND );
  }
  else
  {
    wxGridSizer* grid = new wxGridSizer( 2, 2, 1, 1 );
    m_renderViewHolder->SetSizer( grid );
    for ( int i = 0; i < 4; i++ )
    {
      view[i]->Show();
      grid->Add( view[i], 1, wxEXPAND );
    }
  }

  m_renderViewHolder->Layout();
  view[0]->SetFocus();

  NeedRedraw( 2 );
}

void MainWindow::SetMainView( int nView )
{
  if ( m_nMainView >= 0 )
  {
    // udpate scalar bar
    bool bScalarBar = m_viewRender[m_nMainView]->GetShowScalarBar();
    if ( bScalarBar )
    {
      m_viewRender[m_nMainView]->ShowScalarBar( false );
      m_viewRender[nView]->ShowScalarBar( true );
    }
  }
  m_nMainView = nView;

  SetViewLayout( m_nViewLayout );
  this->SendBroadcast( "MainViewChanged", this );
}

void MainWindow::OnViewLayout1X1( wxCommandEvent& event )
{
  SetViewLayout( VL_1X1 );
}

void MainWindow::OnViewLayout1X1UpdateUI( wxUpdateUIEvent& event )
{
  event.Check( m_nViewLayout == VL_1X1 );
}

void MainWindow::OnViewLayout2X2( wxCommandEvent& event )
{
  SetViewLayout( VL_2X2 );
}

void MainWindow::OnViewLayout2X2UpdateUI( wxUpdateUIEvent& event )
{
  event.Check( m_nViewLayout == VL_2X2 );
}

void MainWindow::OnViewLayout1N3( wxCommandEvent& event )
{
  SetViewLayout( VL_1N3 );
}

void MainWindow::OnViewLayout1N3UpdateUI( wxUpdateUIEvent& event )
{
  event.Check( m_nViewLayout == VL_1N3 );
}


void MainWindow::OnViewLayout1N3_H( wxCommandEvent& event )
{
  SetViewLayout( VL_1N3_H );
}

void MainWindow::OnViewLayout1N3_HUpdateUI( wxUpdateUIEvent& event )
{
  event.Check( m_nViewLayout == VL_1N3_H );
}

void MainWindow::OnViewSagittal( wxCommandEvent& event )
{
  SetMainView( MV_Sagittal );
}

void MainWindow::OnViewSagittalUpdateUI( wxUpdateUIEvent& event )
{
  event.Check( m_nMainView == MV_Sagittal );
}

void MainWindow::OnViewCoronal( wxCommandEvent& event )
{
  SetMainView( MV_Coronal );
}

void MainWindow::OnViewCoronalUpdateUI( wxUpdateUIEvent& event )
{
  event.Check( m_nMainView == MV_Coronal );
}

void MainWindow::OnViewAxial( wxCommandEvent& event )
{
  SetMainView( MV_Axial );
}

void MainWindow::OnViewAxialUpdateUI( wxUpdateUIEvent& event )
{
  event.Check( m_nMainView == MV_Axial );
}

void MainWindow::OnView3D( wxCommandEvent& event )
{
  SetMainView( MV_3D );
}

void MainWindow::OnView3DUpdateUI( wxUpdateUIEvent& event )
{
  event.Check( m_nMainView == MV_3D );
}

void MainWindow::OnViewReset( wxCommandEvent& event )
{
  wxWindow* wnd = FindFocus();
  if ( wnd == m_viewAxial )
    m_viewAxial->ResetView();
  else if ( wnd == m_viewSagittal )
    m_viewSagittal->ResetView();
  else if ( wnd == m_viewCoronal )
    m_viewCoronal->ResetView();
  else if ( wnd == m_view3D )
    m_view3D->ResetView();
}

void MainWindow::OnViewResetUpdateUI( wxUpdateUIEvent& event )
{
  event.Enable( !GetLayerCollection( "MRI" )->IsEmpty() || !GetLayerCollection( "Surface" )->IsEmpty() );
}

void MainWindow::OnViewSnapToAxis( wxCommandEvent& event )
{
  m_view3D->SnapToNearestAxis();
}

void MainWindow::OnViewSnapToAxisUpdateUI( wxUpdateUIEvent& event )
{
  event.Enable( (m_view3D->IsShown() && !GetLayerCollection( "MRI" )->IsEmpty()) || !GetLayerCollection( "Surface" )->IsEmpty() );
}

void MainWindow::OnViewSliceFrames( wxCommandEvent& event )
{
  m_view3D->SetShowSliceFrames( !m_view3D->GetShowSliceFrames() );
}

void MainWindow::OnViewSliceFramesUpdateUI( wxUpdateUIEvent& event )
{
  event.Enable( m_view3D->IsShown() );
  event.Check( m_view3D->GetShowSliceFrames() );
}

void MainWindow::ShowControlPanel( bool bShow )
{
  wxConfigBase* config = wxConfigBase::Get();
  if ( bShow )
  {
    m_splitterMain->SplitVertically( m_controlPanel, m_splitterSub,
                                     config ? config->Read( _T("/MainWindow/SplitterPosition"), 260L ) : 260 );
  }
  else
  {
    if ( config )
      config->Write( _("/MainWindow/SplitterPosition"), m_splitterMain->GetSashPosition() );
    m_splitterMain->Unsplit( m_controlPanel );
  }
}

void MainWindow::OnViewControlPanel( wxCommandEvent& event )
{
  ShowControlPanel( event.IsChecked() );
}

void MainWindow::OnViewControlPanelUpdateUI( wxUpdateUIEvent& event )
{
  event.Check( m_controlPanel->IsShown() );
}

void MainWindow::OnViewScalarBar( wxCommandEvent& event )
{
  ShowScalarBar( event.IsChecked() );
}

void MainWindow::ShowScalarBar( bool bShow )
{
  for ( int i = 0; i < 4; i++ )
  {
    if ( i != m_nMainView )
      m_viewRender[i]->ShowScalarBar( false );
  }
  if ( m_nMainView >= 0 && m_nMainView < 4 )
  {
    m_viewRender[m_nMainView]->ShowScalarBar( bShow );
    m_viewRender[m_nMainView]->UpdateScalarBar();
  }

  NeedRedraw( 1 );
}

void MainWindow::OnViewScalarBarUpdateUI( wxUpdateUIEvent& event )
{
//  LayerMRI* mri = (LayerMRI*)MainWindow::GetMainWindowPointer()->GetLayerCollection( "MRI" )->GetActiveLayer();
//  LayerSurface* surf = (LayerSurface*)MainWindow::GetMainWindowPointer()->GetLayerCollection( "Surface" )->GetActiveLayer();
//  event.Enable( mri || ( m_nMainView == 3 && surf && surf->GetActiveOverlay() ) );
  if ( m_nMainView >= 0 && m_nMainView < 4 )
    event.Check( m_viewRender[m_nMainView]->GetShowScalarBar() );
}


void MainWindow::OnViewCoordinate( wxCommandEvent& event )
{
  for ( int i = 0; i < 3; i++ )
  {
    ( (RenderView2D*)m_viewRender[i] )->ShowCoordinateAnnotation( event.IsChecked() );
  }

  NeedRedraw( 1 );
}

void MainWindow::OnViewCoordinateUpdateUI( wxUpdateUIEvent& event )
{
  event.Enable( MainWindow::GetMainWindowPointer()->GetLayerCollection( "MRI" )->GetActiveLayer() );
  event.Check( m_viewAxial->GetShowCoordinateAnnotation() );
}

LayerCollection* MainWindow::GetCurrentLayerCollection()
{
  LayerCollection* lc = NULL;
  wxString name = m_controlPanel->GetCurrentLayerCollectionName();
  if ( name == _("Volumes") )
    lc = GetLayerCollection( "MRI" );
  else if ( name == _("ROIs") )
    lc = GetLayerCollection( "ROI" );
  else if ( name == _("Surfaces") )
    lc = GetLayerCollection( "Surface" );
  else if ( name == _("Point Sets") )
    lc = GetLayerCollection( "WayPoints" );
  
  return lc;
}

void MainWindow::OnViewCycleLayer( wxCommandEvent& event )
{
  LayerCollection* lc = GetCurrentLayerCollection();
  if ( lc )
  {
    lc->CycleLayer();
  }
}

void MainWindow::OnViewReverseCycleLayer( wxCommandEvent& event )
{
  LayerCollection* lc = GetCurrentLayerCollection();

  if ( lc )
  {
    lc->CycleLayer( false );
  }
}

void MainWindow::OnViewCycleLayerUpdateUI( wxUpdateUIEvent& event )
{
  LayerCollection* lc = GetCurrentLayerCollection();
  event.Enable( lc && lc->GetNumberOfLayers() > 1 );
}

void MainWindow::OnViewToggleVolumeVisibility( wxCommandEvent& event )
{
  LayerCollection* lc = GetLayerCollection( "MRI" );
  if ( !lc->IsEmpty() )
  {
    Layer* layer = lc->GetActiveLayer();
    if ( layer )
      layer->SetVisible( !layer->IsVisible() );
  }
}

void MainWindow::OnViewToggleVolumeVisibilityUpdateUI( wxUpdateUIEvent& event )
{
  event.Enable( !GetLayerCollection( "MRI" )->IsEmpty() );
}


void MainWindow::OnViewToggleROIVisibility( wxCommandEvent& event )
{
  LayerCollection* lc = GetLayerCollection( "ROI" );
  if ( !lc->IsEmpty() )
  {
    Layer* layer = lc->GetActiveLayer();
    if ( layer )
      layer->SetVisible( !layer->IsVisible() );
  }
}

void MainWindow::OnViewToggleROIVisibilityUpdateUI( wxUpdateUIEvent& event )
{
  event.Enable( !GetLayerCollection( "ROI" )->IsEmpty() );
}


void MainWindow::OnViewToggleSurfaceVisibility( wxCommandEvent& event )
{
  LayerCollection* lc = GetLayerCollection( "Surface" );
  if ( !lc->IsEmpty() )
  {
    Layer* layer = lc->GetActiveLayer();  
    if ( layer )
      layer->SetVisible( !layer->IsVisible() );
  }
}

void MainWindow::OnViewToggleSurfaceVisibilityUpdateUI( wxUpdateUIEvent& event )
{
  event.Enable( !GetLayerCollection( "Surface" )->IsEmpty() );
}


void MainWindow::OnViewToggleWayPointsVisibility( wxCommandEvent& event )
{
  LayerCollection* lc = GetLayerCollection( "WayPoints" );
  if ( !lc->IsEmpty() )
  {
    Layer* layer = lc->GetActiveLayer();
    if ( layer )
      layer->SetVisible( !layer->IsVisible() );
  }
}

void MainWindow::OnViewToggleWayPointsVisibilityUpdateUI( wxUpdateUIEvent& event )
{
  event.Enable( !GetLayerCollection( "WayPoints" )->IsEmpty() );
}

void MainWindow::OnViewToggleCursorVisibility( wxCommandEvent& event )
{
  bool bCur = m_viewAxial->GetCursor2D()->IsShown();
  m_viewAxial->GetCursor2D()->Show( !bCur );
  m_viewSagittal->GetCursor2D()->Show( !bCur );
  m_viewCoronal->GetCursor2D()->Show( !bCur );
  m_view3D->GetCursor3D()->Show( !bCur );
  NeedRedraw( 1 );
}

void MainWindow::OnViewToggleCursorVisibilityUpdateUI( wxUpdateUIEvent& event )
{}

void MainWindow::OnViewSurfaceMain( wxCommandEvent& event )
{
  LayerSurface* layer = (LayerSurface*)GetLayerCollection( "Surface" )->GetActiveLayer();
  if ( layer )
    layer->SetActiveSurface( FSSurface::SurfaceMain );
}

void MainWindow::OnViewSurfaceMainUpdateUI( wxUpdateUIEvent& event )
{
  LayerSurface* layer = (LayerSurface*)GetLayerCollection( "Surface" )->GetActiveLayer();
  event.Enable( layer );
  event.Check( layer && layer->GetActiveSurface() == FSSurface::SurfaceMain );
}

void MainWindow::OnViewSurfaceInflated( wxCommandEvent& event )
{
  LayerSurface* layer = (LayerSurface*)GetLayerCollection( "Surface" )->GetActiveLayer();
  if ( layer )
    layer->SetActiveSurface( FSSurface::SurfaceInflated );
}

void MainWindow::OnViewSurfaceInflatedUpdateUI( wxUpdateUIEvent& event )
{
  LayerSurface* layer = (LayerSurface*)GetLayerCollection( "Surface" )->GetActiveLayer();
  FSSurface* surf = ( layer ? layer->GetSourceSurface() : NULL );
  event.Enable( layer && surf && surf->IsSurfaceLoaded( FSSurface::SurfaceInflated ) );
  event.Check( layer && layer->GetActiveSurface() == FSSurface::SurfaceInflated );
}

void MainWindow::OnViewSurfaceWhite( wxCommandEvent& event )
{
  LayerSurface* layer = (LayerSurface*)GetLayerCollection( "Surface" )->GetActiveLayer();
  if ( layer )
    layer->SetActiveSurface( FSSurface::SurfaceWhite );
}

void MainWindow::OnViewSurfaceWhiteUpdateUI( wxUpdateUIEvent& event )
{
  LayerSurface* layer = (LayerSurface*)GetLayerCollection( "Surface" )->GetActiveLayer();
  FSSurface* surf = ( layer ? layer->GetSourceSurface() : NULL );
  event.Enable( layer && surf && surf->IsSurfaceLoaded( FSSurface::SurfaceWhite ) );
  event.Check( layer && layer->GetActiveSurface() == FSSurface::SurfaceWhite );
}

void MainWindow::OnViewSurfacePial( wxCommandEvent& event )
{
  LayerSurface* layer = (LayerSurface*)GetLayerCollection( "Surface" )->GetActiveLayer();
  if ( layer )
    layer->SetActiveSurface( FSSurface::SurfacePial );
}

void MainWindow::OnViewSurfacePialUpdateUI( wxUpdateUIEvent& event )
{
  LayerSurface* layer = (LayerSurface*)GetLayerCollection( "Surface" )->GetActiveLayer();
  FSSurface* surf = ( layer ? layer->GetSourceSurface() : NULL );
  event.Enable( layer && surf && surf->IsSurfaceLoaded( FSSurface::SurfacePial ) );
  event.Check( layer && layer->GetActiveSurface() == FSSurface::SurfacePial );
}

void MainWindow::OnViewSurfaceOriginal( wxCommandEvent& event )
{
  LayerSurface* layer = (LayerSurface*)GetLayerCollection( "Surface" )->GetActiveLayer();
  if ( layer )
    layer->SetActiveSurface( FSSurface::SurfaceOriginal );
}

void MainWindow::OnViewSurfaceOriginalUpdateUI( wxUpdateUIEvent& event )
{
  LayerSurface* layer = (LayerSurface*)GetLayerCollection( "Surface" )->GetActiveLayer();
  FSSurface* surf = ( layer ? layer->GetSourceSurface() : NULL );
  event.Enable( layer && surf && surf->IsSurfaceLoaded( FSSurface::SurfaceOriginal ) );
  event.Check( layer && layer->GetActiveSurface() == FSSurface::SurfaceOriginal );
}

void MainWindow::NeedRedraw( int nCount )
{
  m_nRedrawCount = nCount;
}

void MainWindow::ForceRedraw()
{
  for ( int i = 0; i < 4; i++ )
  {
    if ( m_viewRender[i]->IsShown() )
      m_viewRender[i]->Render();
  }
}

void MainWindow::OnInternalIdle()
{
  wxFrame::OnInternalIdle();

#if defined(__WXGTK__) || defined(__WXMAC__)
  if ( IsShown() && m_nRedrawCount > 0 )
  {
    for ( int i = 0; i < 4; i++ )
      if ( m_viewRender[i]->IsShown() )
        m_viewRender[i]->Render();
    m_nRedrawCount--;
  }
#endif

  if ( m_bToUpdateToolbars )
    DoUpdateToolbars();

  // hack to restore previously saved sub-splitter postion.
  static bool run_once = false;
  if ( !run_once && IsShown() )
  {
    wxConfigBase* config = wxConfigBase::Get();
    if ( config )
    {
      if ( config->Exists( _("/MainWindow/SplitterPositionSub") ) )
      {
        m_splitterSub->SetSashPosition( config->Read( _("/MainWindow/SplitterPositionSub"), 80L ) );
        m_splitterSub->UpdateSize();
      }
    }
    
    run_once = true;
  }
  
  if ( m_bDoScreenshot )
    DoSaveScreenshot();
}

void MainWindow::OnEditPreferences( wxCommandEvent& event )
{
  DialogPreferences dlg( this );
  dlg.SetGeneralSettings( m_settingsGeneral );
  dlg.Set2DSettings( m_settings2D );
  dlg.SetScreenshotSettings( m_settingsScreenshot );

  if ( dlg.ShowModal() == wxID_OK )
  {
    m_settingsGeneral = dlg.GetGeneralSettings();
    m_settings2D = dlg.Get2DSettings();
    m_settingsScreenshot = dlg.GetScreenshotSettings();

    m_viewAxial   ->SetBackgroundColor( m_settingsGeneral.BackgroundColor );
    m_viewSagittal->SetBackgroundColor( m_settingsGeneral.BackgroundColor );
    m_viewCoronal ->SetBackgroundColor( m_settingsGeneral.BackgroundColor );
    m_view3D      ->SetBackgroundColor( m_settingsGeneral.BackgroundColor );
    m_viewAxial->GetCursor2D()    ->SetColor( m_settingsGeneral.CursorColor );
    m_viewSagittal->GetCursor2D() ->SetColor( m_settingsGeneral.CursorColor );
    m_viewCoronal->GetCursor2D()  ->SetColor( m_settingsGeneral.CursorColor );
    m_view3D->GetCursor3D()       ->SetColor( m_settingsGeneral.CursorColor );
    m_viewAxial->GetCursor2D()    ->SetStyle( m_settingsGeneral.CursorStyle );
    m_viewSagittal->GetCursor2D() ->SetStyle( m_settingsGeneral.CursorStyle );
    m_viewCoronal->GetCursor2D()  ->SetStyle( m_settingsGeneral.CursorStyle );
  }
}

void MainWindow::OnHelpQuickReference( wxCommandEvent& event )
{
  m_wndQuickReference->Show();
  wxConfigBase* config = wxConfigBase::Get();
  if ( config )
  {
    int x = config->Read( _("/QuickRefWindow/PosX"), 280L );
    int y = config->Read( _("/QuickRefWindow/PosY"), 30L );
    int x1, y1;
    GetPosition( &x1, &y1 );
    m_wndQuickReference->Move( x1 + x, y1 + y );
  }
}

void MainWindow::OnHelpAbout( wxCommandEvent& event )
{
  wxString msg = _( "freeview 1.0 (internal) \r\nbuild " ) + MyUtils::GetDateAndTime();

  wxMessageDialog dlg( this, msg, _("About freeview"), wxOK | wxICON_INFORMATION );
  dlg.ShowModal();
}

void MainWindow::OnFileNew( wxCommandEvent& event )
{
  NewVolume();
}

void MainWindow::OnFileNewUpdateUI( wxUpdateUIEvent& event )
{
  event.Enable( !GetLayerCollection( "MRI" )->IsEmpty() );
}

void MainWindow::OnFileSave( wxCommandEvent& event )
{
  SaveVolume();
}

void MainWindow::OnFileSaveUpdateUI( wxUpdateUIEvent& event )
{
  LayerMRI* layer = ( LayerMRI* )( GetLayerCollection( "MRI" )->GetActiveLayer() );
  event.Enable( layer && layer->IsModified() && !IsProcessing() && !IsWritingMovieFrames() );
}


void MainWindow::OnFileSaveAs( wxCommandEvent& event )
{
  SaveVolumeAs();
}

void MainWindow::OnFileSaveAsUpdateUI( wxUpdateUIEvent& event )
{
  LayerMRI* layer = ( LayerMRI* )( GetLayerCollection( "MRI" )->GetActiveLayer() );
  event.Enable( layer && layer->IsEditable() && !IsProcessing() && !IsWritingMovieFrames() );
}

void MainWindow::OnFileLoadROI( wxCommandEvent& event )
{
  LoadROI();
}

void MainWindow::OnFileLoadROIUpdateUI( wxUpdateUIEvent& event )
{
  event.Enable( !GetLayerCollection( "MRI" )->IsEmpty() && !IsProcessing() && !IsWritingMovieFrames() );
}


void MainWindow::OnFileSaveROI( wxCommandEvent& event )
{
  SaveROI();
}

void MainWindow::OnFileSaveROIUpdateUI( wxUpdateUIEvent& event )
{
  LayerROI* layer = ( LayerROI* )( GetLayerCollection( "ROI" )->GetActiveLayer() );
  event.Enable( layer && layer->IsModified() && !IsProcessing() && !IsWritingMovieFrames() );
}


void MainWindow::OnFileSaveROIAs( wxCommandEvent& event )
{
  SaveROIAs();
}

void MainWindow::OnFileSaveROIAsUpdateUI( wxUpdateUIEvent& event )
{
  LayerROI* layer = ( LayerROI* )( GetLayerCollection( "ROI" )->GetActiveLayer() );
  event.Enable( layer && !IsProcessing() && !IsWritingMovieFrames() );
}


void MainWindow::OnFileLoadWayPoints( wxCommandEvent& event )
{
  LoadWayPoints();
}

void MainWindow::OnFileLoadWayPointsUpdateUI( wxUpdateUIEvent& event )
{
  event.Enable( !GetLayerCollection( "MRI" )->IsEmpty() && !IsProcessing() && !IsWritingMovieFrames() );
}


void MainWindow::OnFileSaveWayPoints( wxCommandEvent& event )
{
  SaveWayPoints();
}

void MainWindow::OnFileSaveWayPointsUpdateUI( wxUpdateUIEvent& event )
{
  LayerWayPoints* layer = ( LayerWayPoints* )( GetLayerCollection( "WayPoints" )->GetActiveLayer() );
  event.Enable( layer && layer->IsModified() && !IsProcessing() && !IsWritingMovieFrames() );
}


void MainWindow::OnFileSaveWayPointsAs( wxCommandEvent& event )
{
  SaveWayPointsAs();
}

void MainWindow::OnFileSaveWayPointsAsUpdateUI( wxUpdateUIEvent& event )
{
  LayerWayPoints* layer = ( LayerWayPoints* )( GetLayerCollection( "WayPoints" )->GetActiveLayer() );
  event.Enable( layer && !IsProcessing() && !IsWritingMovieFrames() );
}


void MainWindow::OnWorkerThreadResponse( wxCommandEvent& event )
{
  wxString strg = event.GetString();

  if ( strg.Left( 6 ) == _("Failed") )
  {
    m_statusBar->m_gaugeBar->Hide();
    m_bProcessing = false;
    EnableControls( true );
    m_controlPanel->UpdateUI();
    strg = strg.Mid( 6 );
    if ( strg.IsEmpty() )
      strg = _("Operation failed. See console for more information.");
    wxMessageDialog dlg( this, strg, _("Error"), wxOK | wxICON_ERROR );
    dlg.ShowModal();
    return;
  }

  if ( event.GetInt() == -1 )  // successfully finished
  {
    m_bProcessing = false;
    EnableControls( true );
    m_statusBar->m_gaugeBar->Hide();

    Layer* layer = ( Layer* )(void*)event.GetClientData();
    wxASSERT( layer != NULL || strg != _("Load") || strg != _("Save") );
    LayerCollection* lc_mri = GetLayerCollection( "MRI" );
    LayerCollection* lc_surface = GetLayerCollection( "Surface" );

    // loading operation finished
    if ( strg == _("Load" ) )
    {
      // volume loaded
      if ( layer->IsTypeOf( "MRI" ) )
      {
        LayerMRI* mri = (LayerMRI*)layer;
        if ( lc_mri->IsEmpty() )
        {
          double worigin[3], wsize[3];
          mri->GetWorldOrigin( worigin );
          mri->GetWorldSize( wsize );
          if ( lc_surface->IsEmpty() )
          {
            mri->SetSlicePositionToWorldCenter();
            m_viewAxial->SetWorldCoordinateInfo( worigin, wsize );
            m_viewSagittal->SetWorldCoordinateInfo( worigin, wsize );
            m_viewCoronal->SetWorldCoordinateInfo( worigin, wsize );
            m_view3D->SetWorldCoordinateInfo( worigin, wsize );
          }
          else
          {
            mri->SetSlicePosition( lc_surface->GetSlicePosition() );
            lc_surface->SetWorldVoxelSize( mri->GetWorldVoxelSize() );
            lc_surface->SetWorldOrigin( mri->GetWorldOrigin() );
            lc_surface->SetWorldSize( mri->GetWorldSize() );
          }

          lc_mri->AddLayer( mri, true );
          lc_mri->SetCursorRASPosition( lc_mri->GetSlicePosition() );
          m_layerVolumeRef = mri;

          mri->AddListener( this );
        }
        else
          lc_mri->AddLayer( layer );

        m_fileHistory->AddFileToHistory( MyUtils::GetNormalizedFullPath( wxString::FromAscii( mri->GetFileName() ) ) );

        m_controlPanel->RaisePage( _("Volumes") );
      }
      // surface loaded
      else if ( layer->IsTypeOf( "Surface" ) )
      {
        LayerSurface* sf = (LayerSurface*)layer;
        if ( lc_surface->IsEmpty() )
        {
          double worigin[3], wsize[3];
          sf->GetWorldOrigin( worigin );
          sf->GetWorldSize( wsize );
          sf->SetSlicePositionToWorldCenter();
          if ( lc_mri->IsEmpty() )
          {
            m_viewAxial->SetWorldCoordinateInfo( worigin, wsize );
            m_viewSagittal->SetWorldCoordinateInfo( worigin, wsize );
            m_viewCoronal->SetWorldCoordinateInfo( worigin, wsize );
            m_view3D->SetWorldCoordinateInfo( worigin, wsize );
            lc_surface->AddLayer( sf, true );
          }
          else
          {
            lc_surface->SetWorldOrigin( lc_mri->GetWorldOrigin() );
            lc_surface->SetWorldSize( lc_mri->GetWorldSize() );
            lc_surface->SetWorldVoxelSize( lc_mri->GetWorldVoxelSize() );
            lc_surface->SetSlicePosition( lc_mri->GetSlicePosition() );
            lc_surface->AddLayer( sf );
          }
          // lc_surface->SetCursorRASPosition( lc_surface->GetSlicePosition() );
        }
        else
          lc_surface->AddLayer( layer );
        
        if ( !sf->HasVolumeGeometry() )
        {
          wxMessageDialog dlg( this, _("Surface does not contain valid volume geometry information. It may not align with volumes and other surfaces."), 
                               _("Warning"), wxOK );
          dlg.ShowModal();
        }
        // m_fileHistory->AddFileToHistory( MyUtils::GetNormalizedFullPath( layer->GetFileName() ) );

        m_controlPanel->RaisePage( _("Surfaces") );
      }

      SetMode( RenderView2D::IM_Navigate );
    }

    // Saving operation finished
    else if ( strg == _("Save") )
    {
      cout << ( (LayerEditable*)layer )->GetFileName() << " saved successfully." << endl;
      if ( layer->IsTypeOf( "Surface" ) )
          m_dlgRepositionSurface->UpdateUI();
    }

    else if ( strg == _("Rotate") )
    {
      m_bResampleToRAS = false;
      m_layerCollectionManager->RefreshSlices();
    }
   
    m_controlPanel->UpdateUI( true ); 

#ifdef __WXMAC__
    // On Mac OS X (Carbon) with wxWidgets 2.8.10 we are seeing a problem where
    // after loading a volume, the PanelVolume widgets do not get laid out     
    // properly.  Resizing the window fixes the problem.  This *HACK* fix
    // is to set the size to a different size, and then set it back to
    // its old size, which seems to force a resize event and works around the
    // problem - DRG
    int width, height;
    GetSize(&width, &height);
    SetSize(width, height-1);
    SetSize(width, height);
#endif

    RunScript();
  }
  else if ( event.GetInt() == 0 )  // just started
  {
    m_bProcessing = true;
    if ( strg == _("Rotate") )
    {
      EnableControls( false );
    }

    m_statusBar->ActivateProgressBar();
    NeedRedraw();
    m_controlPanel->UpdateUI( true );
  }
  else
  {
    // if ( event.GetInt() > m_statusBar->m_gaugeBar->GetValue() )
    int val = event.GetInt();
    if ( val > m_statusBar->m_gaugeBar->GetRange() )
      val = m_statusBar->m_gaugeBar->GetRange();
    m_statusBar->m_gaugeBar->SetValue( val );
  }
}

void MainWindow::OnBuildContourThreadResponse( wxCommandEvent& event )
{
  // if the thread id returned (by event.GetInt() ) is not the same as the current thread id,
  // it means there is a new thread currently processing the contour, so this event should 
  // be ignored. 
  LayerMRI* mri = (LayerMRI*)( (void*)event.GetClientData() );
  if ( GetLayerCollection( "MRI" )->GetLayerIndex( mri ) < 0 ||
       mri->GetBuildContourThreadID() != event.GetInt() )
    return;
  
  mri->RealizeContourActor();
  if ( m_volumeCropper->GetEnabled() )
  {
    m_volumeCropper->UpdateProps();
    m_volumeCropper->Apply();
  }
  
  if ( mri->GetProperties()->GetShowAsContour() )
    NeedRedraw();
}

void MainWindow::DoListenToMessage ( std::string const iMsg, void* iData, void* sender )
{
  if ( iMsg == "LayerAdded" || iMsg == "LayerRemoved" || iMsg == "LayerMoved" )
  {
    UpdateToolbars();
  }
  else if ( iMsg == "ActiveLayerChanged" )
  {
    if ( !GetLayerCollection( "MRI" )->IsEmpty() )
    {
      LayerMRI* mri = (LayerMRI*)GetLayerCollection( "MRI" )->GetActiveLayer();
      wxString fn;
      if ( mri->IsTypeOf( "DTI" ) )
        fn = wxString::FromAscii( ((LayerDTI*)mri)->GetVectorFileName() );
      else
        fn = wxString::FromAscii( mri->GetFileName() );
      SetTitle( wxString::FromAscii( mri->GetName() ) + _(" (") + fn + _(") - freeview") );
    }
    else if ( !GetLayerCollection( "Surface" )->IsEmpty() )
    {
      LayerSurface* surf = (LayerSurface*)GetLayerCollection( "Surface" )->GetActiveLayer();
      SetTitle( wxString::FromAscii( surf->GetName() ) + _(" (") + surf->GetFileName() + _(") - freeview") );
    }
    else
    {
      SetTitle( _("freeview") );
    }
  }
  else if ( iMsg == "LayerObjectDeleted" )
  {
    if ( m_layerVolumeRef == iData )
    {
      m_layerVolumeRef = (LayerMRI*)GetActiveLayer( "MRI" );
    }
  }
  else if ( iMsg == "MRINotVisible" )
  {
    wxMessageDialog dlg( this, _("Active volume is not visible. Please turn it on before editing."), 
                         _("Error"), wxOK | wxICON_ERROR );
    dlg.ShowModal();
  }
  else if ( iMsg == "MRINotEditable" )
  {
    wxMessageDialog dlg( this, _("Active volume is not editable."), _("Error"), wxOK | wxICON_ERROR );
    dlg.ShowModal();
  }
  else if ( iMsg == "MRINotEditableForRotation" )
  {
    wxMessageDialog dlg( this, _("Active volume has been transformed. It is not a good idea to directly edit on transformed volume. Because partial volume effect may cause \"what you see is NOT what you get\". Please save, close and reload the volume to edit. This is a temporary and safe solution."), _("Error"), wxOK | wxICON_ERROR );
    dlg.ShowModal();
  }
  else if ( iMsg == "MRIReferenceNotSet" )
  {
    wxMessageDialog dlg( this, _("Reference volume is not set."), _("Error"), wxOK | wxICON_ERROR );
    dlg.ShowModal();
  }
  else if ( iMsg == "ROINotVisible" )
  {
    wxMessageDialog dlg( this, _("Active ROI is not visible. Please turn it on before editing."), 
                         _("Error"), wxOK | wxICON_ERROR );
    dlg.ShowModal();
  }
  else if ( iMsg == "SurfacePositionChanged" )
  {
    m_view3D->UpdateConnectivityDisplay();
  }
  else if ( iMsg == "CropBoundChanged" )
  {
    NeedRedraw();
  }
  
  // Update world geometry
  LayerCollection* mri_col = GetLayerCollection( "MRI" );
  LayerCollection* surf_col = GetLayerCollection( "Surface" );
  if ( !mri_col || !surf_col )
    return;
  Layer* mri = mri_col->GetActiveLayer();
  double new_origin[3], new_ext[3];
  if ( mri )
  {
    mri->GetWorldOrigin( new_origin);
    double* ws = mri->GetWorldSize();
    new_ext[0] = new_origin[0] + ws[0];
    new_ext[1] = new_origin[1] + ws[1];
    new_ext[2] = new_origin[2] + ws[2];
  }
  for ( int i = 0; i < mri_col->GetNumberOfLayers(); i++ )
  {
    Layer* layer = mri_col->GetLayer( i );
    double* orig = layer->GetWorldOrigin();
    double ext[3];
    double* ws = layer->GetWorldSize(); 
    for ( int i = 0; i < 3; i++ )
    {   
      ext[i] = orig[i] + ws[i];
      if ( orig[i] < new_origin[i] )
        new_origin[i] = orig[i];
      if ( ext[i] > new_ext[i] )
        new_ext[i] = ext[i];
    }
  }
  if ( mri )
  {
    double new_ws[3];
    double* voxel_size = mri->GetWorldVoxelSize();
    double* origin = mri->GetWorldOrigin();
    for ( int i = 0; i < 3; i++ )
    {
      new_origin[i] = origin[i] - ( (int) ( (origin[i] - new_origin[i])/voxel_size[i] + 0.5 ) ) * voxel_size[i];
      new_ws[i] = ((int)((new_ext[i]-new_origin[i])/voxel_size[i]+0.5)) * voxel_size[i];
    }
    mri_col->SetWorldOrigin( new_origin );
    mri_col->SetWorldSize( new_ws );
    mri_col->SetWorldVoxelSize( voxel_size );
    surf_col->SetWorldOrigin( new_origin );
    surf_col->SetWorldSize( new_ws );
    surf_col->SetWorldVoxelSize( voxel_size );
  }
}

void MainWindow::OnFileLoadDTI( wxCommandEvent& event )
{
  DialogLoadDTI dlg( this );
  dlg.Initialize( m_bResampleToRAS, GetLayerCollection( "MRI" )->IsEmpty() );
  dlg.SetLastDir( m_strLastDir );

  wxArrayString list;
  for ( int i = 0; i < m_fileHistory->GetMaxFiles(); i++ )
    list.Add( m_fileHistory->GetHistoryFile( i ) );
  dlg.SetRecentFiles( list );

  if ( dlg.ShowModal() != wxID_OK )
    return;

  this->LoadDTIFile( dlg.GetVectorFileName(), dlg.GetFAFileName(), dlg.GetRegFileName(), dlg.IsToResample() );
}

void MainWindow::LoadDTIFile( const wxString& fn_vector,
                              const wxString& fn_fa,
                              const wxString& reg_filename,
                              bool bResample )
{
  m_strLastDir = MyUtils::GetNormalizedPath( fn_fa );
  m_bResampleToRAS = bResample;

  LayerDTI* layer = new LayerDTI( m_layerVolumeRef );
  layer->SetResampleToRAS( bResample );
  wxString layerName = wxFileName( fn_vector ).GetName();
  if ( wxFileName( fn_fa ).GetExt().Lower() == _("gz") )
    layerName = wxFileName( layerName ).GetName();
  layer->SetName( layerName.char_str() );
  layer->SetFileName( fn_fa.char_str() );
  layer->SetVectorFileName( fn_vector.char_str() );
  if ( !reg_filename.IsEmpty() )
  {
    wxFileName reg_fn( reg_filename );
    reg_fn.Normalize( wxPATH_NORM_ALL, m_strLastDir );
    layer->SetRegFileName( reg_fn.GetFullPath().char_str() );
  }

  /* LayerMRI* mri = (LayerMRI* )GetLayerCollection( "MRI" )->GetLayer( 0 );
   if ( mri )
   {
    layer->SetRefVolume( mri->GetSourceVolume() );
   }
  */
  WorkerThread* thread = new WorkerThread( this );
  thread->LoadVolume( layer );
}

void MainWindow::OnFileLoadPVolumes( wxCommandEvent& event )
{
  DialogLoadPVolumes dlg( this );

  if ( dlg.ShowModal() != wxID_OK )
    return;
  
  wxArrayString filenames = dlg.GetVolumeFileNames();
  wxString prefix = dlg.GetFileNamePrefix();
  wxString lut = dlg.GetLUT();
  LoadPVolumeFiles( filenames, prefix, lut );
}

void MainWindow::LoadPVolumeFiles( const wxArrayString& filenames, const wxString& prefix, const wxString& lut )
{
  LayerPLabel* layer = new LayerPLabel( m_layerVolumeRef );
  layer->SetVolumeFileNames( filenames );
  layer->SetFileNamePrefix( prefix );
  layer->SetLUT( lut );
  
  WorkerThread* thread = new WorkerThread( this );
  thread->LoadVolume( layer );
}

void MainWindow::LoadConnectivityDataFile( const wxString& data_file, const wxString& lut )
{
  LayerCollection* lc = GetLayerCollection( "Surface" );
  if ( lc->GetNumberOfLayers() < 2 )
  {
    cerr << "Can not load connectivity data. Please load both hemisphere surfaces first!" << endl;
    return;
  }
  
  m_connectivity->Load( data_file, lut, (LayerSurface*)lc->GetLayer( 0 ), (LayerSurface*)lc->GetLayer( 1 ) );
  if ( m_connectivity->IsValid() )
    m_wndConnectivityConfiguration->ShowWindow();
}

void MainWindow::AddScript( const wxArrayString& script )
{
  m_scripts.push_back( script );
}

void MainWindow::RunScript()
{
  if ( m_scripts.size() == 0 )
    return;

  wxArrayString sa = m_scripts[0];
  m_scripts.erase( m_scripts.begin() );
  if ( sa[0] == _("loadvolume") )
  {
    CommandLoadVolume( sa );
  }
  else if ( sa[0] == _("loaddti") )
  {
    CommandLoadDTI( sa );
  }
  else if ( sa[0] == _("loadsurface") )
  {
    CommandLoadSurface( sa );
  }
  else if ( sa[0] == _("loadsurfacevector") )
  {
    CommandLoadSurfaceVector( sa );
  }
  else if ( sa[0] == _("loadsurfacecurvature") )
  {
    CommandLoadSurfaceCurvature( sa );
  }
  else if ( sa[0] == _("loadsurfaceoverlay") )
  {
    CommandLoadSurfaceOverlay( sa );
  }
  else if ( sa[0] == _("loadsurfaceannotation") )
  {
    CommandLoadSurfaceAnnotation( sa );
  }
  else if ( sa[0] == _("loadsurfacelabel") )
  {
    CommandLoadSurfaceLabel( sa );
  }
  else if ( sa[0] == _("loadconnectivity") )
  {
    CommandLoadConnectivityData( sa );
  }
  else if ( sa[0] == _("loadroi") || sa[0] == _("loadlabel") )
  {
    CommandLoadROI( sa );
  }
  else if ( sa[0] == _("loadwaypoints") )
  {
    CommandLoadWayPoints( sa );
  }
  else if ( sa[0] == _("loadcontrolpoints") )
  {
    CommandLoadControlPoints( sa );
  }
  else if ( sa[0] == _("loadpvolumes") )
  {
    CommandLoadPVolumes( sa );
  }
  else if ( sa[0] == _("screencapture") )
  {
    CommandScreenCapture( sa );
  }
  else if ( sa[0] == _("quit") || sa[0] == _("exit") )
  {
    Close();
  }
  else if ( sa[0] == _("setviewport") )
  {
    CommandSetViewport( sa );
  }
  else if ( sa[0] == _("setviewsize") )
  {
    CommandSetViewSize( sa );
  }
  else if ( sa[0] == _("zoom") )
  {
    CommandZoom( sa );
  }
  else if ( sa[0] == _("ras") )
  {
    CommandSetRAS( sa );
  }
  else if ( sa[0] == _("slice") )
  {
    CommandSetSlice( sa );
  }
  else if ( sa[0] == _("setcolormap") )
  {
    CommandSetColorMap( sa );
  }
  else if ( sa[0] == _("setheatscaleoptions") )
  {
    CommandSetHeadScaleOptions( sa );
  }
  else if ( sa[0] == _("setlut") )
  {
    CommandSetLUT( sa );
  }
  else if ( sa[0] == _("setopacity") )
  {
    CommandSetOpacity( sa );
  }
  else if ( sa[0] == _("setdisplayisosurface") )
  {
    CommandSetDisplayIsoSurface( sa );
  }
  else if ( sa[0] == _("setisosurfacecolor") )
  {
    CommandSetIsoSurfaceColor( sa );
  }
  else if ( sa[0] == _("loadisosurfaceregion") )
  {
    CommandLoadIsoSurfaceRegion( sa );
  }
  else if ( sa[0] == _("setsurfaceoverlaymethod") )
  {
    CommandSetSurfaceOverlayMethod( sa );
  }
  else if ( sa[0] == _("setsurfaceoffset") )
  {
    CommandSetSurfaceOffset( sa );
  }
  else if ( sa[0] == _("setwaypointscolor") )
  {
    CommandSetWayPointsColor( sa );
  }
  else if ( sa[0] == _("setwaypointsradius") )
  {
    CommandSetWayPointsRadius( sa );
  }
  else if ( sa[0] == _("setdisplayvector") )
  {
    CommandSetDisplayVector( sa );
  }
  else if ( sa[0] == _("setdisplaytensor") )
  {
    CommandSetDisplayTensor( sa );
  }
  else if ( sa[0] == _("setsurfacecolor") )
  {
    CommandSetSurfaceColor( sa );
  }
  else if ( sa[0] == _("setsurfaceedgecolor") )
  {
    CommandSetSurfaceEdgeColor( sa );
  }
  else if ( sa[0] == _("setsurfaceedgethickness") )
  {
    CommandSetSurfaceEdgeThickness( sa );
  }
  else if ( sa[0] == _("setlayername") )
  {
    CommandSetLayerName( sa );
  }
  else if ( sa[0] == _("locklayer") )
  {
    CommandLockLayer( sa );
  }
  else if ( sa[0] == _("showlayer") )
  {
    CommandShowLayer( sa );
  }
}

void MainWindow::CommandLoadVolume( const wxArrayString& sa )
{
  wxArrayString sa_vol = MyUtils::SplitString( sa[1], _(":") );
  wxString fn = sa_vol[0];
  wxString reg_fn;
  wxArrayString scales;
  wxString colormap = _("grayscale");
  wxString colormap_scale = _("grayscale");
  wxString lut_name;
  wxString vector_display = _("no"), 
           vector_inversion = _("none"), 
           vector_render = _("line"),
           tensor_display = _("no"),
           tensor_render = _("boxoid");
  int nSampleMethod = m_nDefaultSampleMethod;
  bool bConform = m_bDefaultConform;
  for ( size_t i = 1; i < sa_vol.GetCount(); i++ )
  {
    wxString strg = sa_vol[i];
    int n = strg.Find( _("=") );
    if ( n != wxNOT_FOUND )
    {
      wxString subOption = strg.Left( n ).Lower();
      wxString subArgu = strg.Mid( n + 1 );
      if ( subOption == _("colormap") )
      {
        colormap = subArgu.Lower();
      }
      else if ( subOption == _("grayscale") || 
                subOption == _("heatscale") || 
                subOption == _("colorscale") ) 
      {
        colormap_scale = subOption;    // colormap scale might be different from colormap!
        scales = MyUtils::SplitString( subArgu, _(",") );
      }
      else if ( subOption == _("heatscaleoption") || 
                subOption == _("heatscaleoptions") )
      {
        wxArrayString script;
        script.Add( _("setheatscaleoptions") );
        wxArrayString opts = MyUtils::SplitString( subArgu, _(",") );
        for ( size_t i = 0; i < opts.size(); i++ )
          script.Add( opts[i] );
    
        m_scripts.insert( m_scripts.begin(), script );
      }
      else if ( subOption == _("lut") )
      {
        lut_name = subArgu;
        if ( lut_name.IsEmpty() )
        {
          cerr << "Missing lut name." << endl;
        }
      }
      else if ( subOption == _("vector") )
      {
        vector_display = subArgu.Lower();
        if ( vector_display.IsEmpty() )
        {
          cerr << "Missing vector display argument." << endl;
        }
      }
      else if ( subOption == _("tensor") )
      {
        tensor_display = subArgu.Lower();
        if ( tensor_display.IsEmpty() )
        {
          cerr << "Missing tensor display argument." << endl;
        }
      }
      else if ( subOption == _("inversion") || 
                subOption == _("invert") )
      {
        vector_inversion = subArgu.Lower();
        if ( vector_inversion.IsEmpty() )
        {
          cerr << "Missing inversion argument." << endl;
          vector_inversion = _("none");
        }
      }
      else if ( subOption == _("render") )
      {
        vector_render = subArgu.Lower();
        tensor_render = vector_render;
        if ( vector_render.IsEmpty() )
        {
          cerr << "Missing render argument." << endl;
          vector_render = _("line");
          tensor_render = _("boxoid");
        }
      }
      else if ( subOption == _("reg") )
      {
        reg_fn = subArgu;
      }
      else if ( subOption == _("sample") )
      {
        if ( subArgu.Lower() == _("nearest") )
          nSampleMethod = SAMPLE_NEAREST;
        else if ( subArgu.Lower() == _("trilinear") )
          nSampleMethod = SAMPLE_TRILINEAR;
      }
      else if ( subOption == _("opacity") )
      {
        wxArrayString script;
        script.Add( _("setopacity") );
        script.Add( subArgu );
        
        m_scripts.insert( m_scripts.begin(), script );
      }
      else if ( subOption == _("isosurface") )
      {
        wxArrayString script;
        script.Add( _("setdisplayisosurface") );
        wxArrayString argus = MyUtils::SplitString( subArgu, _(",") );
        if ( argus.size() > 0 && argus[0].size() > 0 )
        {
          script.Add( argus[0] );
        }
        if ( argus.size() > 1 && argus[1].size() > 0 )
        {
          script.Add( argus[1] );
        }
        m_scripts.insert( m_scripts.begin(), script );
      }
      else if (subOption == _("color"))
      {
        wxArrayString script;
        script.Add( _("setisosurfacecolor") );
        script.Add(subArgu);
        m_scripts.insert( m_scripts.begin(), script );
      }
      else if ( subOption == _("surface_region") || subOption == _("surface_regions") )
      {
        wxArrayString script;
        script.Add( _("loadisosurfaceregion") );
        wxFileName fn( subArgu );
        fn.Normalize();
        script.Add( fn.GetFullPath() );       
        m_scripts.insert( m_scripts.begin(), script );
      }
      else if ( subOption == _("name") )
      {
        wxArrayString script;
        script.Add( _("setlayername") );
        script.Add( _("MRI") );
        script.Add( subArgu );
        m_scripts.insert( m_scripts.begin(), script );
      }
      else if ( subOption == _("lock") )
      {
        wxArrayString script;
        script.Add( _("locklayer") );
        script.Add( _("MRI") );
        script.Add( subArgu );
        m_scripts.insert( m_scripts.begin(), script );
      }
      else if ( subOption == _("visible") )
      {
        wxArrayString script;
        script.Add( _("showlayer") );
        script.Add( _("MRI") );
        script.Add( subArgu );
        m_scripts.insert( m_scripts.begin(), script );
      }
      else
      {
        cerr << "Unrecognized sub-option flag '" << strg << "'." << endl;
        return;
      }
    }
    else
    {
      cerr << "Unrecognized sub-option flag '" << strg << "'." << endl;
      return;
    }
  }
  bool bResample = false;
  if ( sa[ sa.GetCount()-1 ] == _("r") )
    bResample = true;

  if ( scales.size() > 0 || colormap != _("grayscale") )
  {
    wxArrayString script;
    script.Add( _("setcolormap") );
    script.Add( colormap );
    script.Add( colormap_scale ); 
    for ( size_t i = 0; i < scales.size(); i++ )
      script.Add( scales[i] );
    
    m_scripts.insert( m_scripts.begin(), script );
  }
  
  if ( !lut_name.IsEmpty() )
  {
    wxArrayString script;
    script.Add( _("setlut") );
    script.Add( lut_name );
        
    m_scripts.insert( m_scripts.begin(), script );
  }
  
  if ( !tensor_display.IsEmpty() && tensor_display != _("no") )
  {
    wxArrayString script;
    script.Add( _("setdisplaytensor") );
    script.Add( tensor_display );
    script.Add( tensor_render );
    script.Add( vector_inversion );
  
    m_scripts.insert( m_scripts.begin(), script );  
  }
  else if ( !vector_display.IsEmpty() && vector_display != _("no") )
  {
    wxArrayString script;
    script.Add( _("setdisplayvector") );
    script.Add( vector_display );
    script.Add( vector_render );
    script.Add( vector_inversion );
  
    m_scripts.insert( m_scripts.begin(), script );  
  }
  
  LoadVolumeFile( fn, reg_fn, bResample, nSampleMethod, bConform );
}

void MainWindow::CommandSetColorMap( const wxArrayString& sa )
{
  int nColorMap = LayerPropertiesMRI::Grayscale;;
  wxString strg = sa[1];
  if ( strg == _("heat") || strg == _("heatscale") )
    nColorMap = LayerPropertiesMRI::Heat;
  else if ( strg == _("jet") || strg == _("jetscale") )
    nColorMap = LayerPropertiesMRI::Jet;
  else if ( strg == _("lut") )
    nColorMap = LayerPropertiesMRI::LUT;
  else if ( strg == _("gecolor" ) || strg == _("ge_color" ) )
    nColorMap = LayerPropertiesMRI::GEColor;
  else if ( strg == _("nih") )
    nColorMap = LayerPropertiesMRI::NIH;
  else if ( strg != _("grayscale") )
    cerr << "Unrecognized colormap name '" << strg << "'." << endl;
  
  int nColorMapScale = LayerPropertiesMRI::Grayscale;
  strg = sa[2];
  if ( strg == _("heatscale") )
    nColorMapScale = LayerPropertiesMRI::Heat;
  else if ( strg == _("colorscale") )
    nColorMapScale = LayerPropertiesMRI::Jet;
  else if ( strg == _("lut") )
    nColorMapScale = LayerPropertiesMRI::LUT;
  
  std::vector<double> pars;
  for ( size_t i = 3; i < sa.size(); i++ )
  {
    double dValue;
    if ( !sa[i].ToDouble( &dValue ) )
    {
      cerr << "Invalid color scale value(s). " << endl;
      break;
    }
    else
    {
      pars.push_back( dValue );
    }
  }
  
  SetVolumeColorMap( nColorMap, nColorMapScale, pars );
  
  ContinueScripts();
}

void MainWindow::CommandSetHeadScaleOptions( const wxArrayString& sa )
{
  if ( GetLayerCollection( "MRI" )->GetActiveLayer() )
  {
    LayerPropertiesMRI* p = ( (LayerMRI*)GetLayerCollection( "MRI" )->GetActiveLayer() )->GetProperties();
    for ( size_t i = 1; i < sa.size(); i++ )
    {
      if ( sa[i] == _("invert") )
        p->SetHeatScaleInvert( true );
      else if ( sa[i] == _("truncate") )
        p->SetHeatScaleTruncate( true );
    }
  }
  
  ContinueScripts();
}

void MainWindow::CommandSetLayerName( const wxArrayString& cmd )
{
  if ( cmd.size() > 2 )
  {
    LayerCollection* lc = GetLayerCollection( (const char*)cmd[1].c_str() );
    if ( lc && !lc->IsEmpty() )
    {
      lc->GetActiveLayer()->SetName( (const char*)cmd[2].c_str() );
    }
  }
  ContinueScripts();
}

void MainWindow::CommandLockLayer( const wxArrayString& cmd )
{
  if ( cmd.size() > 2 && ( cmd[2] == _("1") || cmd[2] == _("true") ) )
  {
    LayerCollection* lc = GetLayerCollection( (const char*)cmd[1].c_str() );
    if ( lc && !lc->IsEmpty() )
    {
      lc->GetActiveLayer()->Lock( true );
    }
  }
  ContinueScripts();
}

void MainWindow::CommandShowLayer( const wxArrayString& cmd )
{
  if ( cmd.size() > 2 && ( cmd[2] == _("0") || cmd[2] == _("false") ) )
  {
    LayerCollection* lc = GetLayerCollection( (const char*)cmd[1].c_str() );
    if ( lc && !lc->IsEmpty() )
    {
      lc->GetActiveLayer()->SetVisible( false );
    }
  }
  ContinueScripts();
}


void MainWindow::CommandSetDisplayVector( const wxArrayString& cmd )
{
  if ( cmd[1].Lower() == _("yes") || cmd[1].Lower() == _("true") || cmd[1].Lower() == _("1") )
  {
    LayerMRI* mri = (LayerMRI*)GetLayerCollection( "MRI" )->GetActiveLayer();
    if ( mri )
    {
      if ( !mri->IsTypeOf( "DTI" ) && mri->GetNumberOfFrames() < 3 )
      {
        cerr << "Volume has less than 3 frames. Can not display as vectors." << endl;
      }
      else
      {
        mri->GetProperties()->SetDisplayVector( true );  
      
        if ( cmd[2].Lower() == _("line") )
        {
          mri->GetProperties()->SetVectorRepresentation( LayerPropertiesMRI::VR_Line );
        }
        else if ( cmd[2].Lower() == _("bar") )
        {
          mri->GetProperties()->SetVectorRepresentation( LayerPropertiesMRI::VR_Bar );
        }
        else
        {
          cerr << "Unrecognized argument '" << cmd[2] << "' for vector rendering." << endl;
        }
        
        if ( cmd[3].Lower() != _("none") )
        {
          if ( cmd[3].Lower() == _("x") )
            mri->GetProperties()->SetVectorInversion( LayerPropertiesMRI::VI_X );
          else if ( cmd[3].Lower() == _("y") )
            mri->GetProperties()->SetVectorInversion( LayerPropertiesMRI::VI_Y );
          else if ( cmd[3].Lower() == _("z") )
            mri->GetProperties()->SetVectorInversion( LayerPropertiesMRI::VI_Z );
          else
            cerr << "Unknown inversion flag '" << cmd[2] << "'." << endl;
        }
      }
    }
  }
  
  ContinueScripts();
}


void MainWindow::CommandSetDisplayTensor( const wxArrayString& cmd )
{
  if ( cmd[1].Lower() == _("yes") || cmd[1].Lower() == _("true") || cmd[1].Lower() == _("1") )
  {
    LayerMRI* mri = (LayerMRI*)GetLayerCollection( "MRI" )->GetActiveLayer();
    if ( mri )
    {
      if ( mri->GetNumberOfFrames() < 9 )
      {
        cerr << "Volume has less than 9 frames. Can not display as tensor." << endl;
      }
      else
      {
        mri->GetProperties()->SetDisplayTensor( true );  
      
        if ( cmd[2].Lower().Find( _("box") ) != wxNOT_FOUND )
        {
          mri->GetProperties()->SetTensorRepresentation( LayerPropertiesMRI::TR_Boxoid );
        }
        else if ( cmd[2].Lower().Find( _("ellips") ) != wxNOT_FOUND )
        {
          mri->GetProperties()->SetTensorRepresentation( LayerPropertiesMRI::TR_Ellipsoid );
        }
        else
        {
          cerr << "Unrecognized argument '" << cmd[2] << "' for tensor rendering." << endl;
        }
        
        if ( cmd[3].Lower() != _("none") )
        {
          if ( cmd[3].Lower() == _("x") )
            mri->GetProperties()->SetTensorInversion( LayerPropertiesMRI::VI_X );
          else if ( cmd[3].Lower() == _("y") )
            mri->GetProperties()->SetTensorInversion( LayerPropertiesMRI::VI_Y );
          else if ( cmd[3].Lower() == _("z") )
            mri->GetProperties()->SetTensorInversion( LayerPropertiesMRI::VI_Z );
          else
            cerr << "Unknown inversion flag '" << cmd[2] << "'." << endl;
        }
      }
    }
  }
  
  ContinueScripts();
}


void MainWindow::CommandSetLUT( const wxArrayString& sa )
{
  LayerMRI* mri = (LayerMRI*)GetLayerCollection( "MRI" )->GetActiveLayer();
  if ( mri )
  {  
    COLOR_TABLE* ct = m_luts->LoadColorTable( sa[1].char_str() );

    if ( ct )
    {
      mri->GetProperties()->SetLUTCTAB( ct );
    }
  }
  
  ContinueScripts();
}

void MainWindow::CommandSetOpacity( const wxArrayString& sa )
{
  LayerMRI* mri = (LayerMRI*)GetLayerCollection( "MRI" )->GetActiveLayer();
  if ( mri )
  {  
    double dValue;
    if ( sa[1].ToDouble( &dValue ) )
    {
      mri->GetProperties()->SetOpacity( dValue );
    }
    else
    {
      cerr << "Opacity value is not valid." << endl;
    }
  }
  
  ContinueScripts();
}

void MainWindow::CommandSetDisplayIsoSurface( const wxArrayString& sa )
{
  LayerMRI* mri = (LayerMRI*)GetLayerCollection( "MRI" )->GetActiveLayer();
  if ( mri )
  {  
    double dValue;
    if ( sa.size() > 1 )
    {
      if ( sa[1].ToDouble( &dValue ) )
      {
        mri->GetProperties()->SetContourMinThreshold( dValue );
      }
      else if ( sa[1].Lower() != "on" )
      {
        cerr << "Isosurface threshold value is not valid." << endl;
      }
    }
    if ( sa.size() > 2 )
    {
      if ( sa[2].ToDouble( &dValue ) )
      {
        mri->GetProperties()->SetContourMaxThreshold( dValue );
      }
      else
      {
        cerr << "Isosurface threshold value is not valid." << endl;
      }
    }
    mri->GetProperties()->SetShowAsContour( true );
  }
  
  ContinueScripts();
}

void MainWindow::CommandSetIsoSurfaceColor( const wxArrayString& cmd )
{
  LayerMRI* mri = (LayerMRI*)GetLayerCollection( "MRI" )->GetActiveLayer();
  if ( mri )
  { 
    wxColour color( cmd[1] );
    if ( !color.IsOk() )      
    {
      long rgb[3];
      wxArrayString rgb_strs = MyUtils::SplitString( cmd[1], _(",") );
      if ( rgb_strs.GetCount() < 3 || 
           !rgb_strs[0].ToLong( &(rgb[0]) ) ||
           !rgb_strs[1].ToLong( &(rgb[1]) ) ||
           !rgb_strs[2].ToLong( &(rgb[2]) ) )
      {
        cerr << "Invalid isosurface color name or value " << cmd[1] << endl;
      }
      else
      {
        color.Set( rgb[0], rgb[1], rgb[2] );
      }
    }
      
    if ( color.IsOk() )
      mri->GetProperties()->SetContourColor( color.Red()/255.0, color.Green()/255.0, color.Blue()/255.0 );
    else
      cerr << "Invalid isosurface color name or value " << cmd[1] << endl;  
  }
  ContinueScripts();
}

void MainWindow::CommandLoadIsoSurfaceRegion( const wxArrayString& sa )
{
  LayerMRI* mri = (LayerMRI*)GetLayerCollection( "MRI" )->GetActiveLayer();
  if ( mri )
  {  
    if ( sa.size() > 1 )
    {
      if ( !mri->LoadSurfaceRegions( sa[1] ) )
      {
        cerr << "Can not load surfacer region(s) from " << sa[1].c_str() << endl;
      }
    }
  }
  
  ContinueScripts();
}

void MainWindow::ContinueScripts()
{
  if ( m_bProcessing )
    return;
  
  // create a fake worker event to notify end of operation 
  // so scripts in queue will continue on
  wxCommandEvent event( wxEVT_COMMAND_MENU_SELECTED, ID_WORKER_THREAD );
  event.SetInt( -1 );
  event.SetString( _("ByPass") );
  wxPostEvent( this, event );
}

void MainWindow::CommandLoadDTI( const wxArrayString& sa )
{
  bool bResample = false;
  if ( sa.GetCount() > 3 && sa[3] == _("r") )
    bResample = true;
  
  if ( sa.GetCount() > 2 )
  {
    wxArrayString sa_vol = MyUtils::SplitString( sa[1], _(":") );
    wxString fn = sa_vol[0];
    wxString strg, reg_fn;
    wxString vector_display = _("no"), 
             vector_inversion = _("none"),
             vector_render = _("line");
    
    for ( size_t i = 1; i < sa_vol.GetCount(); i++ )
    {
      wxString strg = sa_vol[i];
      int n = strg.Find( _("=") );
      if ( n != wxNOT_FOUND )
      {  
        if ( strg.Left( n ).Lower() == _("vector") )
        {
          vector_display = strg.Mid( n + 1 ).Lower();
          if ( vector_display.IsEmpty() )
          {
            cerr << "Missing vector display argument." << endl;
          }
        }
        else if ( strg.Left( n ).Lower() == _("inversion") || 
                  strg.Left( n ).Lower() == _("invert") )
        {
          vector_inversion = strg.Mid( n + 1 ).Lower();
          if ( vector_inversion.IsEmpty() )
          {
            cerr << "Missing vector inversion argument." << endl;
            vector_inversion = _("none");
          }
        }
        else if ( strg.Left( n ).Lower() == _("render") )
        {
          vector_render = strg.Mid( n + 1 ).Lower();
          {
            if ( vector_render.IsEmpty() )
            {
              cerr << "Missing vector render argument." << endl;
              vector_render = _("line");
            }
          }
        }
        else if ( strg.Left( n ).Lower() == _("reg") )
        {
          reg_fn = strg.Mid( n + 1 );
        }
      }
    }
    
    if ( !vector_display.IsEmpty() && vector_display != _("no") )
    {
      wxArrayString script;
      script.Add( _("setdisplayvector") );
      script.Add( vector_display );
      script.Add( vector_render );
      script.Add( vector_inversion );
  
      m_scripts.insert( m_scripts.begin(), script );  
    }
    
    // cout << reg_fn.c_str() << endl;
    this->LoadDTIFile( fn, sa[2], reg_fn, bResample );
  }
}

void MainWindow::CommandLoadPVolumes( const wxArrayString& cmd )
{
  wxArrayString files = MyUtils::SplitString( cmd[1], _(";") );  
  wxString lut = _("");
  if ( cmd.size() > 3 )
  {
    lut = cmd[3];
    COLOR_TABLE* ct = m_luts->LoadColorTable( lut.char_str() );
    if ( !ct )
    {
      cerr << "Can not load look up table " << lut.c_str() << endl;
      return;
    }
  }
  this->LoadPVolumeFiles( files, cmd[2], lut );
}

void MainWindow::CommandLoadConnectivityData( const wxArrayString& cmd )
{
  if ( cmd.size() > 2 )
    this->LoadConnectivityDataFile( cmd[1], cmd[2] );
  
  ContinueScripts();
}

void MainWindow::CommandLoadROI( const wxArrayString& cmd )
{
  wxArrayString options = MyUtils::SplitString( cmd[1], _(":") );
  wxString fn = options[0];
  wxString ref;
  for ( size_t i = 1; i < options.GetCount(); i++ )
  {
    wxString strg = options[i];
    int n = strg.Find( _("=") );
    if ( n != wxNOT_FOUND )
    {
      wxString option = strg.Left( n ).Lower();
      wxString argu = strg.Mid( n+1 );
      if ( option == _("ref") || option == _("template") )
      {
        ref = argu;
      }
      else
      {
        cerr << "Unrecognized sub-option flag '" << strg << "'." << endl;
      }
    }
  }
  
  LoadROIFile( fn, ref );

  ContinueScripts();
}

void MainWindow::CommandLoadSurface( const wxArrayString& cmd )
{
  wxString fullfn = cmd[1];
  int nIgnoreStart = fullfn.find( _("#") );
  int nIgnoreEnd = fullfn.find( _("#"), nIgnoreStart+1 );
  wxArrayString sa_fn = MyUtils::SplitString( fullfn, _(":"), nIgnoreStart, nIgnoreEnd - nIgnoreStart + 1 );
  wxString fn = sa_fn[0];
  wxString fn_patch = _("");
  wxString fn_target = _("");
  for ( size_t k = 1; k < sa_fn.size(); k++ )
  {
    int n = sa_fn[k].Find( _("=") );
    if ( n != wxNOT_FOUND  )
    {
      wxString subOption = sa_fn[k].Left( n ).Lower(); 
      wxString subArgu = sa_fn[k].Mid( n+1 );
      if ( subOption == _( "color" ) )
      {
        wxArrayString script;
        script.Add( _("setsurfacecolor") );
        script.Add( subArgu );
            
        m_scripts.insert( m_scripts.begin(), script );
      }
      else if ( subOption == _( "edgecolor" ) )
      {
        wxArrayString script;
        script.Add( _("setsurfaceedgecolor") );
        script.Add( subArgu );
            
        m_scripts.insert( m_scripts.begin(), script );
      }
      else if ( subOption == _( "edgethickness" ) )
      {
        wxArrayString script;
        script.Add( _("setsurfaceedgethickness") );
        script.Add( subArgu );
            
        m_scripts.insert( m_scripts.begin(), script );
      }
      else if ( subOption == _("curv") || subOption == _("curvature") )
      {
        // add script to load surface curvature file
        wxArrayString script;
        script.Add( _("loadsurfacecurvature") );
        script.Add( subArgu );
        m_scripts.insert( m_scripts.begin(), script );
      }
      else if ( subOption == _("overlay") )
      {
        // add script to load surface overlay files
        nIgnoreStart = subArgu.find( _("#") );
        nIgnoreEnd = subArgu.find( _("#"), nIgnoreStart+1 );
        wxArrayString overlay_fns = MyUtils::SplitString( subArgu, _(","), nIgnoreStart, nIgnoreEnd - nIgnoreStart + 1 );
        for ( int i = overlay_fns.size() - 1 ; i >= 0 ; i-- )
        {
          wxArrayString script;
          script.Add( _("loadsurfaceoverlay") );          
          int nSubStart = overlay_fns[i].find( _("#") );      
          int nSubEnd = overlay_fns[i].find( _("#"), nSubStart+1 );
          if ( nSubEnd == wxNOT_FOUND )
            nSubEnd = overlay_fns[i].Length() - 1;
          if ( nSubStart != wxNOT_FOUND )
            script.Add( overlay_fns[i].Left( nSubStart ) );
          else
            script.Add( overlay_fns[i] );
          m_scripts.insert( m_scripts.begin(), script );
          
          // if there are sub-options attached with overlay file, parse them
          wxString opt_strg;
          if ( nSubStart != wxNOT_FOUND )
            opt_strg = overlay_fns[i].Mid( nSubStart+1, nSubEnd - nSubStart - 1 ); 

          wxArrayString overlay_opts = MyUtils::SplitString( opt_strg, _(":") );
          
          wxString method = _("linearopaque");
          wxArrayString thresholds;
          for ( size_t j = 0; j < overlay_opts.GetCount(); j++ )
          {
            wxString strg = overlay_opts[j];
            if ( ( n = strg.Find( _("=") ) ) != wxNOT_FOUND && strg.Left( n ).Lower() == _("method") )
            {
              method = strg.Mid( n+1 ).Lower();
            }
            else if ( ( n = strg.Find( _("=") ) ) != wxNOT_FOUND && strg.Left( n ).Lower() == _("threshold") )
            {
              thresholds = MyUtils::SplitString( strg.Mid( n+1 ), _(",") );
            }
          }
          
          if ( method != _("linearopaque") || !thresholds.IsEmpty() )
          {
            script.Clear();
            script.Add( _("setsurfaceoverlaymethod") );
            script.Add( method );
            for ( size_t j = 0; j < thresholds.GetCount(); j++ )
              script.Add( thresholds[j] );
            
            // insert right AFTER loadsurfaceoverlay command
            m_scripts.insert( m_scripts.begin()+1, script );
          }
        }
      }
      else if ( subOption == _("annot") || subOption == _("annotation") )
      {
        // add script to load surface annotation files
        wxArrayString annot_fns = MyUtils::SplitString( subArgu, _(",") );
        for ( int i = annot_fns.size()-1; i >= 0; i-- )
        {
          wxArrayString script;
          script.Add( _("loadsurfaceannotation") );
          script.Add( annot_fns[i] );
          m_scripts.insert( m_scripts.begin(), script );
        }
      }
      else if ( subOption == _("label") )
      {
        // add script to load surface label files
        wxArrayString fns = MyUtils::SplitString( subArgu, _(",") );
        for ( int i = fns.size()-1; i >= 0; i-- )
        {
          wxArrayString script;
          script.Add( _("loadsurfacelabel") );
          script.Add( fns[i] );
          m_scripts.insert( m_scripts.begin(), script );
        }
      }
      else if ( subOption == _("vector") )
      {
        // add script to load surface vector files
        wxArrayString vector_fns = MyUtils::SplitString( subArgu, _(",") );
        for ( int i = vector_fns.size() - 1 ; i >= 0 ; i-- )
        {
          wxArrayString script;
          script.Add( _("loadsurfacevector") );
          script.Add( vector_fns[i] );
          m_scripts.insert( m_scripts.begin(), script );
        }
      }    
      else if ( subOption == _("patch") )
      {
        if ( subArgu.Contains( _("/") ) )
          subArgu = MyUtils::GetNormalizedFullPath( subArgu );
        fn_patch = subArgu;
      }  
      else if ( subOption == _("target") || subOption == _("target_surf"))
      {
        if ( subArgu.Contains( _("/") ) )
          subArgu = MyUtils::GetNormalizedFullPath( subArgu );
        fn_target = subArgu;
      }
      else if ( subOption == _("name") )
      {
        wxArrayString script;
        script.Add( _("setlayername") );
        script.Add( _("Surface") );
        script.Add( subArgu );
        m_scripts.insert( m_scripts.begin(), script );
      }
      else if ( subOption == _("lock") )
      {
        wxArrayString script;
        script.Add( _("locklayer") );
        script.Add( _("Surface") );
        script.Add( subArgu );
        m_scripts.insert( m_scripts.begin(), script );
      }
      else if ( subOption == _("visible") )
      {
        wxArrayString script;
        script.Add( _("showlayer") );
        script.Add( _("Surface") );
        script.Add( subArgu );
        m_scripts.insert( m_scripts.begin(), script );
      }
      else if ( subOption == _( "offset" ) )
      {
        wxArrayString script;
        script.Add( _("setsurfaceoffset") );
        wxArrayString values = MyUtils::SplitString( subArgu, _(",") );
        for ( size_t i = 0; i < values.size(); i++ )
          script.Add( values[i] );
        m_scripts.insert( m_scripts.begin(), script );
      }
      else
      {
        cerr << "Unrecognized sub-option flag '" << subOption << "'." << endl;
        return;
      }
    }
  }
  LoadSurfaceFile( fn, fn_patch, fn_target );
}

void MainWindow::CommandSetSurfaceOverlayMethod( const wxArrayString& cmd )
{
  LayerSurface* surf = (LayerSurface*)GetLayerCollection( "Surface" )->GetActiveLayer();
  if ( surf )
  {
    SurfaceOverlay* overlay = surf->GetActiveOverlay();
    if ( overlay )
    {
      int nMethod = SurfaceOverlayProperties::CM_LinearOpaque;
      if ( cmd[1] == _("linear") )
        nMethod = SurfaceOverlayProperties::CM_Linear;
      else if ( cmd[1] == _("piecewise") )
        nMethod = SurfaceOverlayProperties::CM_Piecewise;
      else if ( cmd[1] != _("linearopaque") )
      {
        cerr << "Unrecognized overlay method name '" << cmd[1] << "'." << endl;
        ContinueScripts();
        return;
      }
      
      overlay->GetProperties()->SetColorMethod( nMethod );
      
      double values[3];
      if ( cmd.GetCount() - 2 >= 3 )   // 3 values
      {
        if ( cmd[2].ToDouble( &(values[0]) ) &&
             cmd[3].ToDouble( &(values[1]) ) &&
             cmd[4].ToDouble( &(values[2]) ) )
        {
          overlay->GetProperties()->SetMinPoint( values[0] );
          overlay->GetProperties()->SetMidPoint( values[1] );
          overlay->GetProperties()->SetMaxPoint( values[2] );
        }
        else
        {
          cerr << "Invalid input for overlay threshold." << endl;
        }
      }   
      else if ( cmd.GetCount() - 2 == 2 )   // 2 values
      {
        if ( cmd[2].ToDouble( &(values[0]) ) &&
             cmd[3].ToDouble( &(values[1]) ) )
        {
          overlay->GetProperties()->SetMinPoint( values[0] );
          overlay->GetProperties()->SetMaxPoint( values[1] );
          overlay->GetProperties()->SetMidPoint( ( values[0] + values[1] ) / 2 );
        }
        else
        {
          cerr << "Invalid input for overlay threshold." << endl;
        }
      }  
      surf->UpdateOverlay();    
    }
  }
  
  ContinueScripts();
}


void MainWindow::CommandSetSurfaceColor( const wxArrayString& cmd )
{
  LayerSurface* surf = (LayerSurface*)GetLayerCollection( "Surface" )->GetActiveLayer();
  if ( surf && cmd[1] != _("null") )
  {
    wxColour color( cmd[1] );
    if ( !color.IsOk() )      
    {
      long rgb[3];
      wxArrayString rgb_strs = MyUtils::SplitString( cmd[1], _(",") );
      if ( rgb_strs.GetCount() < 3 || 
           !rgb_strs[0].ToLong( &(rgb[0]) ) ||
           !rgb_strs[1].ToLong( &(rgb[1]) ) ||
           !rgb_strs[2].ToLong( &(rgb[2]) ) )
      {
        cerr << "Invalid color name or value " << cmd[1] << endl;
      }
      else
      {
        color.Set( rgb[0], rgb[1], rgb[2] );
      }
    }
      
    if ( color.IsOk() )
      surf->GetProperties()->SetBinaryColor( color.Red()/255.0, color.Green()/255.0, color.Blue()/255.0 );
    else
      cerr << "Invalid color name or value " << cmd[1] << endl;
  }
  
  ContinueScripts();
}

void MainWindow::CommandSetSurfaceEdgeColor( const wxArrayString& cmd )
{
  LayerSurface* surf = (LayerSurface*)GetLayerCollection( "Surface" )->GetActiveLayer();
  if ( surf && cmd[1] != _("null") )
  {
    wxColour color( cmd[1] );
    if ( !color.IsOk() )      
    {
      long rgb[3];
      wxArrayString rgb_strs = MyUtils::SplitString( cmd[1], _(",") );
      if ( rgb_strs.GetCount() < 3 || 
           !rgb_strs[0].ToLong( &(rgb[0]) ) ||
           !rgb_strs[1].ToLong( &(rgb[1]) ) ||
           !rgb_strs[2].ToLong( &(rgb[2]) ) )
      {
        cerr << "Invalid color name or value " << cmd[1] << endl;
      }
      else
      {
        color.Set( rgb[0], rgb[1], rgb[2] );
      }
    }
      
    if ( color.IsOk() )
      surf->GetProperties()->SetEdgeColor( color.Red()/255.0, color.Green()/255.0, color.Blue()/255.0 );
    else
      cerr << "Invalid color name or value " << cmd[1] << endl;
  }
  
  ContinueScripts();
}

void MainWindow::CommandSetSurfaceEdgeThickness( const wxArrayString& cmd )
{
  LayerSurface* surf = (LayerSurface*)GetLayerCollection( "Surface" )->GetActiveLayer();
  if ( surf )
  {
    long thickness;
    if ( !cmd[1].ToLong( &thickness ) )
      cerr << "Invalid edge thickness value. Must be a integer." << endl;
    else
      surf->GetProperties()->SetEdgeThickness( thickness );
  }
  
  ContinueScripts();
}


void MainWindow::CommandSetSurfaceOffset( const wxArrayString& cmd )
{
  LayerSurface* surf = (LayerSurface*)GetLayerCollection( "Surface" )->GetActiveLayer();
  if ( surf )
  {
    double pos[3];
    if ( cmd.size() < 4 || 
         !cmd[1].ToDouble( &(pos[0]) ) ||
         !cmd[2].ToDouble( &(pos[1]) ) ||
         !cmd[3].ToDouble( &(pos[2]) ) )
    {
      cerr << "Invalid surface offset inputs. Need 3 numbers." << endl;
    }
    else
    {
      surf->GetProperties()->SetPosition( pos ); 
    } 
  }
  ContinueScripts();
}


void MainWindow::CommandLoadSurfaceVector( const wxArrayString& cmd )
{
  LoadSurfaceVectorFile( cmd[1] );
  
  // do not run ContinueScripts() here because LoadSurfaceVectorFile calls worker thread.
}

void MainWindow::CommandLoadSurfaceCurvature( const wxArrayString& cmd )
{
  LoadSurfaceCurvatureFile( cmd[1] );
 
  ContinueScripts();
}

void MainWindow::CommandLoadSurfaceOverlay( const wxArrayString& cmd )
{
  LoadSurfaceOverlayFile( cmd[1] );

  ContinueScripts();
}

void MainWindow::CommandLoadSurfaceAnnotation( const wxArrayString& cmd )
{
  LoadSurfaceAnnotationFile( cmd[1] );
   
  ContinueScripts();
}

void MainWindow::CommandLoadSurfaceLabel( const wxArrayString& cmd )
{
  LoadSurfaceLabelFile( cmd[1] );
   
  ContinueScripts();
}

void MainWindow::CommandLoadWayPoints( const wxArrayString& cmd )
{
  wxArrayString options = MyUtils::SplitString( cmd[1], _(":") );
  wxString fn = options[0];
  wxString color = _("null");
  wxString spline_color = _("null");
  wxString radius = _("0");
  wxString spline_radius = _("0");
  for ( size_t i = 1; i < options.GetCount(); i++ )
  {
    wxString strg = options[i];
    int n = strg.Find( _("=") );
    if ( n != wxNOT_FOUND )
    {
      wxString option = strg.Left( n ).Lower();
      wxString argu = strg.Mid( n+1 );
      if ( option == _("color") )
      {
        color = argu;
      }
      else if ( option == _("splinecolor") )
      {
        spline_color = argu;
      }
      else if ( option == _("radius") )
      {
        radius = argu;
      }
      else if ( option == _("splineradius") )
      {
        spline_radius = argu;
      }
      else if ( option == _("name") )
      {
        wxArrayString script;
        script.Add( _("setlayername") );
        script.Add( _("WayPoints") );
        script.Add( argu );
        m_scripts.insert( m_scripts.begin(), script );
      }
      else if ( option == _("visible") )
      {
        wxArrayString script;
        script.Add( _("showlayer") );
        script.Add( _("WayPoints") );
        script.Add( argu );
        m_scripts.insert( m_scripts.begin(), script );
      }
      else
      {
        cerr << "Unrecognized sub-option flag '" << strg << "'." << endl;
      }
    }
  }
  
  if ( color != _("null") || spline_color != _("null") )
  {
    wxArrayString script;
    script.Add( _("setwaypointscolor") );
    script.Add( color );
    script.Add( spline_color ); 
    
    m_scripts.insert( m_scripts.begin(), script );
  }
  
  if ( radius != _("0") || spline_radius != _("0") )
  {
    wxArrayString script;
    script.Add( _("setwaypointsradius") );
    script.Add( radius );
    script.Add( spline_radius ); 
    
    m_scripts.insert( m_scripts.begin(), script );
  }
  
  LoadWayPointsFile( fn );

  ContinueScripts();
}


void MainWindow::CommandLoadControlPoints( const wxArrayString& cmd )
{
  wxArrayString options = MyUtils::SplitString( cmd[1], _(":") );
  wxString fn = options[0];
  wxString color = _("null");
  wxString radius = _("0");
  for ( size_t i = 1; i < options.GetCount(); i++ )
  {
    wxString strg = options[i];
    int n = strg.Find( _("=") );
    if ( n != wxNOT_FOUND )
    {
      wxString option = strg.Left( n ).Lower();
      wxString argu = strg.Mid( n+1 );
      if ( option == _("color") )
      {
        color = argu;
      }
      else if ( option == _("radius") )
      {
        radius = argu;
      }
      else if ( option == _("name") )
      {
        wxArrayString script;
        script.Add( _("setlayername") );
        script.Add( _("WayPoints") );
        script.Add( argu );
        m_scripts.insert( m_scripts.begin(), script );
      }
      else if ( option == _("visible") )
      {
        wxArrayString script;
        script.Add( _("showlayer") );
        script.Add( _("WayPoints") );
        script.Add( argu );
        m_scripts.insert( m_scripts.begin(), script );
      }
      else
      {
        cerr << "Unrecognized sub-option flag '" << strg << "'." << endl;
      }
    }
  }
  
  if ( color != _("null") )
  {
    wxArrayString script;
    script.Add( _("setwaypointscolor") );
    script.Add( color );
    
    m_scripts.insert( m_scripts.begin(), script );
  }
  
  if ( radius != _("0") )
  {
    wxArrayString script;
    script.Add( _("setwaypointsradius") );
    script.Add( radius );
    
    m_scripts.insert( m_scripts.begin(), script );
  }
  
  LoadControlPointsFile( fn );

  ContinueScripts();
}

void MainWindow::CommandSetWayPointsColor( const wxArrayString& cmd )
{
  LayerWayPoints* wp = (LayerWayPoints*)GetLayerCollection( "WayPoints" )->GetActiveLayer();
  if ( wp )
  {
    if ( cmd[1] != _("null") )
    {
      wxColour color( cmd[1] );
      if ( !color.IsOk() )      
      {
        long rgb[3];
        wxArrayString rgb_strs = MyUtils::SplitString( cmd[1], _(",") );
        if ( rgb_strs.GetCount() < 3 || 
             !rgb_strs[0].ToLong( &(rgb[0]) ) ||
             !rgb_strs[1].ToLong( &(rgb[1]) ) ||
             !rgb_strs[2].ToLong( &(rgb[2]) ) )
        {
          cerr << "Invalid color name or value " << cmd[1] << endl;
        }
        else
        {
          color.Set( rgb[0], rgb[1], rgb[2] );
        }
      }
      
      if ( color.IsOk() )
        wp->GetProperties()->SetColor( color.Red()/255.0, color.Green()/255.0, color.Blue()/255.0 );
      else
        cerr << "Invalid color name or value " << cmd[1] << endl;
    }
    
    if ( cmd.size() > 2 && cmd[2] != _("null") )
    {
      wxColour color( cmd[2] );
      if ( !color.IsOk() )      
      {
        long rgb[3];
        wxArrayString rgb_strs = MyUtils::SplitString( cmd[2], _(",") );
        if ( rgb_strs.GetCount() < 3 || 
             !rgb_strs[0].ToLong( &(rgb[0]) ) ||
             !rgb_strs[1].ToLong( &(rgb[1]) ) ||
             !rgb_strs[2].ToLong( &(rgb[2]) ) )
        {
          cerr << "Invalid color name or value " << cmd[2] << endl;
        }
        else
        {
          color.Set( rgb[0], rgb[1], rgb[2] );
        }
      }
      
      if ( color.IsOk() )
        wp->GetProperties()->SetSplineColor( color.Red()/255.0, color.Green()/255.0, color.Blue()/255.0 );
      else
        cerr << "Invalid color name or value " << cmd[1] << endl;
    }
  }
  
  ContinueScripts();
}

void MainWindow::CommandSetWayPointsRadius( const wxArrayString& cmd )
{
  LayerWayPoints* wp = (LayerWayPoints*)GetLayerCollection( "WayPoints" )->GetActiveLayer();
  if ( wp )
  {
    if ( cmd[1] != _("0") )
    {
      double dvalue;
      if ( cmd[1].ToDouble( &dvalue ) ) 
        wp->GetProperties()->SetRadius( dvalue );
      else
        cerr << "Invalid way points radius." << endl;
    }
    
    if ( cmd.size() > 2 && cmd[2] != _("0") )
    {
      double dvalue;
      if ( cmd[2].ToDouble( &dvalue ) )
        wp->GetProperties()->SetSplineRadius( dvalue );
      else
        cerr << "Invalid spline radius." << endl;
    }
  }
  
  ContinueScripts();
}

void MainWindow::CommandScreenCapture( const wxArrayString& cmd )
{
  m_viewRender[m_nMainView]->SaveScreenshot( cmd[1].c_str(), 
                                             m_settingsScreenshot.Magnification, 
                                             m_settingsScreenshot.AntiAliasing );

  ContinueScripts();
}

void MainWindow::CommandSetViewport( const wxArrayString& cmd )
{
  if ( cmd[1] == _("x") )
    SetMainView( MV_Sagittal );
  else if ( cmd[1] == _("y") )
    SetMainView( MV_Coronal );
  else if ( cmd[1] == _("z") )
    SetMainView( MV_Axial );
  else if ( cmd[1] == _("3d") )
    SetMainView( MV_3D );

  ContinueScripts();
}

void MainWindow::CommandSetViewSize( const wxArrayString& cmd )
{
  long x, y;
  if ( !cmd[1].ToLong( &x ) || !cmd[2].ToLong( &y ) )
  {
    cerr << "Invalid view size." << endl;
  }
  
  int vx, vy;
  ((wxWindow*)m_viewRender[m_nMainView])->GetSize( &vx, &vy );
  wxConfigBase* config = wxConfigBase::Get();
  if ( config )
  {
    if ( config->Exists( _("/MainWindow/SplitterPositionSub") ) )
    {
      vy = config->Read( _("/MainWindow/SplitterPositionSub"), 80L );
      switch( m_nViewLayout )
      {
        case VL_2X2:
          vy = vy/2;
          break;
        case VL_1N3:
          vy = vy/3*2;
          break;
      }
    }
  }
  
  int offsetx = x - vx, offsety = y - vy;
  switch( m_nViewLayout )
  {
    case VL_2X2:
      offsetx *= 2;
      offsety *= 2;
      break;
    case VL_1N3:
      offsety += offsety/2;
      break;
    case VL_1N3_H:
      offsetx += offsetx/2;
      break;
    default:
      break;
  }
  int mx, my;
  this->GetSize( &mx, &my );
  this->SetSize( mx + offsetx, my + offsety ); 
  if ( config )
  {
    if ( config->Exists( _("/MainWindow/SplitterPositionSub") ) )
    {
      int npos = config->Read( _("/MainWindow/SplitterPositionSub"), 80L );
      npos += offsety;
      m_splitterSub->SetSashPosition( npos );
      m_splitterSub->UpdateSize();      
      config->Write( _("/MainWindow/SplitterPositionSub"), npos );
    }
  }
  
  ContinueScripts();
}


void MainWindow::CommandZoom( const wxArrayString& cmd )
{
  double dValue;
  cmd[1].ToDouble( &dValue );
  if ( m_nMainView >= 0 )
  {
    m_viewRender[m_nMainView]->Zoom( dValue );
  }
  
  ContinueScripts();
}

void MainWindow::CommandSetRAS( const wxArrayString& cmd )
{
  double ras[3];
  if ( cmd[1].ToDouble( &(ras[0]) ) &&
       cmd[2].ToDouble( &(ras[1]) ) &&
       cmd[3].ToDouble( &(ras[2]) ) )
  {
    LayerCollection* lc = GetLayerCollection( "MRI" );
    LayerMRI* layer = (LayerMRI*)lc->GetLayer( 0 );
    if ( layer )
    {
      if ( cmd.GetCount() > 4 && cmd[4] == "tkreg" )
        layer->TkRegToNativeRAS( ras, ras );
      layer->RASToTarget( ras, ras );
    }
    lc->SetCursorRASPosition( ras );
    m_layerCollectionManager->SetSlicePosition( ras );
  }
  else
  {
    cerr << "Invalid input values for RAS coordinates. " << endl;
  }

  ContinueScripts();
}


void MainWindow::CommandSetSlice( const wxArrayString& cmd )
{
  LayerCollection* lc_mri = GetLayerCollection( "MRI" );
  if ( !lc_mri->IsEmpty() )
  {
    LayerMRI* mri = (LayerMRI*)lc_mri->GetLayer( lc_mri->GetNumberOfLayers()-1 );;
    long x, y, z;
    if ( cmd[1].ToLong( &x ) && cmd[2].ToLong( &y ) && cmd[3].ToLong( &z ) )
    {
      int slice[3] = { x, y, z };
      double ras[3];
      mri->OriginalIndexToRAS( slice, ras );
      mri->RASToTarget( ras, ras );
      
      lc_mri->SetCursorRASPosition( ras );
      m_layerCollectionManager->SetSlicePosition( ras );
    }
    else
    {
      cerr << "Invalide slice number(s). " << endl;
    }
  }
  else
  {
    cerr << "No volume was loaded. Set slice failed." << endl;
  }

  ContinueScripts();
}

void MainWindow::OnSpinBrushSize( wxSpinEvent& event )
{
  m_propertyBrush->SetBrushSize( event.GetInt() );
}

void MainWindow::OnSpinBrushTolerance( wxSpinEvent& event )
{
  m_propertyBrush->SetBrushTolerance( event.GetInt() );
}

void MainWindow::OnChoiceBrushTemplate( wxCommandEvent& event )
{
// wxChoice* choiceTemplate = XRCCTRL( *m_toolbarBrush, "ID_CHOICE_TEMPLATE", wxChoice );
// LayerEditable* layer = (LayerEditable*)(void*)choiceTemplate->GetClientData( event.GetSelection() );
// if ( layer )
//  m_propertyBrush->SetReferenceLayer( layer );

  UpdateToolbars();
}

void MainWindow::OnCheckBrushTemplate( wxCommandEvent& event )
{
  UpdateToolbars();
}

int MainWindow::GetActiveViewId()
{
  wxWindow* wnd = FindFocus();
  if ( !wnd )
    return m_nPrevActiveViewId;
  
  int nId = -1;
  if ( wnd == m_viewSagittal || wnd->GetParent() == m_viewSagittal )
    nId = 0;
  else if ( wnd == m_viewCoronal || wnd->GetParent() == m_viewCoronal )
    nId = 1;
  else if ( wnd == m_viewAxial || wnd->GetParent() == m_viewAxial )
    nId = 2;
  else if ( wnd == m_view3D || wnd->GetParent() == m_view3D )
    nId = 3;

  if ( nId >= 0 )
    m_nPrevActiveViewId = nId;
  if ( nId < 0 )
    nId = m_nPrevActiveViewId;

  return nId;
}

RenderView* MainWindow::GetActiveView()
{
  int nId = GetActiveViewId();
  if ( nId >= 0 )
    return m_viewRender[nId];
  else
    return NULL;
}

RenderView* MainWindow::GetPreviousActiveView()
{
  if ( m_nPrevActiveViewId < 0 )
    return NULL;
  else
    return m_viewRender[ m_nPrevActiveViewId ];
}

void MainWindow::OnToolSaveGotoPoint( wxCommandEvent& event )
{
  /*
  if ( !m_dlgSavePoint )
  {
    m_dlgSavePoint = new DialogSavePoint( this );
    m_dlgSavePoint->SetGotoPoints( m_strGotoPoints );
  } 
  m_dlgSavePoint->Show();
  */
 
  wxString fn;
  LayerCollection* lc = GetLayerCollection( "MRI" );
  for ( int i = 0; i < lc->GetNumberOfLayers(); i++ )
  {
    fn = ( (LayerMRI*)lc->GetLayer( i ) )->GetFileName();
    if ( !fn.IsEmpty() )
      break;
  }
  if ( fn.IsEmpty() )
  {
    lc = GetLayerCollection( "Surface" );
    for ( int i = 0; i < lc->GetNumberOfLayers(); i++ )
    {
      fn = ( (LayerSurface*)lc->GetLayer( i ) )->GetFileName();
      if ( !fn.IsEmpty() )
        break;
    }
  }
  
  bool bError = false;
  wxString msg;
  if ( !fn.IsEmpty() )
  {
    wxFileName wfn( fn );
    wfn.Normalize();
    wxString dir = AutoSelectLastDir( wfn.GetPath(), _("tmp") );
    if ( wxFileName::DirExists( dir ) )
    { 
      double ras[3];
      GetCursorRAS( ras, true );  // in tkReg coordinate
      wxFFile file( dir + wxFileName::GetPathSeparator() + _("edit.dat"), "w" );
      wxString strg;
      strg << ras[0] << " " << ras[1] << " " << ras[2] << "\n";
      bError = !file.Write( strg );
      file.Close();
      if ( bError )
        msg = _("Can not write to file."); 
    }
    else
    {
      bError = true;
      msg = _("Directory ") + dir + _(" does not exist.");
    }
  }
  else
  {
    bError = true;
    msg = _("Layer file name is empty. Can not decide where to save.");
  }
  if ( bError )
  {
    wxMessageDialog dlg( this, msg, 
                         _("Error"), wxOK );
    dlg.ShowModal();
  }
}

void MainWindow::OnToolSaveGotoPointUpdateUI( wxUpdateUIEvent& event )
{
  event.Enable( !GetLayerCollection("MRI")->IsEmpty() || !GetLayerCollection( "Surface" )->IsEmpty() );
}

void MainWindow::UpdateGotoPoints()
{  
  if ( m_dlgSavePoint )
    m_strGotoPoints = m_dlgSavePoint->GetGotoPoints();
  
  if ( m_strGotoPoints.size() == 0 )
    return;
  
  wxMenuBar* menubar = GetMenuBar();
  if ( menubar )
  {
    int nId = menubar->FindMenuItem( _("Tools"), _("Goto Saved Points") );
    if ( nId != wxNOT_FOUND )
    {
      wxMenuItem* item = menubar->FindItem( nId );
      wxMenu* menu = item->GetMenu();
      menu->Remove( item );
      wxMenu* submenu = new wxMenu;
      BuildGotoPointMenu( submenu );
      if ( item->GetSubMenu() )
        delete item->GetSubMenu();
      item->SetSubMenu( submenu );
      menu->Insert( 1, item );
    }
  }
}

void MainWindow::BuildGotoPointMenu( wxMenu* menu )
{
  for ( size_t i = 0; i < m_strGotoPoints.size(); i++ )
  {
    wxArrayString ar = MyUtils::SplitString( m_strGotoPoints[i], "," );
    int nItemId = menu->FindItem( ar[0] );
    wxMenu* subMenu = NULL;
    if ( nItemId == wxNOT_FOUND )
    {
      subMenu = new wxMenu;
      menu->AppendSubMenu( subMenu, ar[0] );
    }
    else
    {
      subMenu = menu->FindItem( nItemId )->GetSubMenu();
    }
    
    if ( subMenu )
    {
      subMenu->Append( ID_TOOL_GOTO_POINT_1 + (int)i, ar[1] );
    }
  }
}

void MainWindow::OnToolGotoPoint( wxCommandEvent& event )
{
  /*
  if( !m_menuGotoPoints )
    m_menuGotoPoints = new wxMenu;

  BuildGotoPointMenu( m_menuGotoPoints );
  PopupMenu( m_menuGotoPoints, 800, 0 );
  delete m_menuGotoPoints;
  m_menuGotoPoints = NULL;
  */
  LayerCollection* lc = GetLayerCollection( "MRI" );
  wxString fn;
  for ( int i = 0; i < lc->GetNumberOfLayers(); i++ )
  {
    fn = ( (LayerMRI*)lc->GetLayer( i ) )->GetFileName();
    wxFileName wfn( fn );
    wfn.Normalize();
    fn = AutoSelectLastDir( wfn.GetPath(), _("tmp") ) + wxFileName::GetPathSeparator() + _("edit.dat");
    if ( wxFileName::FileExists( fn ) )
      break;
    else  
      fn = _("");
  }
  
  if ( !fn.IsEmpty() )
  {
    wxFFile file( fn );
    wxString strg;
    file.ReadAll( &strg );
    wxArrayString ar_line = MyUtils::SplitString( strg, "\n" );
    wxArrayString ar = MyUtils::SplitString( ar_line[0], " " );
    ar.insert( ar.begin(), _("setras") );
    ar.push_back( _("tkreg") );
    CommandSetRAS( ar ); 
  }
  else
  {
    wxMessageDialog dlg( this, _("Could not find saved point."), 
                         _("Error"), wxOK );
    dlg.ShowModal();
  }
}

void MainWindow::OnToolGotoPointUpdateUI( wxUpdateUIEvent& event )
{
  /*
  event.Enable( ( !GetLayerCollection("MRI")->IsEmpty() || !GetLayerCollection( "Surface" )->IsEmpty() ) &&
      !m_strGotoPoints.IsEmpty() );
  event.Check( m_menuGotoPoints );*/
  event.Enable( !GetLayerCollection("MRI")->IsEmpty() || !GetLayerCollection( "Surface" )->IsEmpty() );
}

void MainWindow::OnToolMenuGotoPoint( wxCommandEvent& event )
{
  int n = event.GetId() - ID_TOOL_GOTO_POINT_1;
  if ( n < 0 )
    return; 

  wxArrayString ar = MyUtils::SplitString( m_strGotoPoints[n], "," );
  ar.RemoveAt( 0 );
  CommandSetRAS( ar );  
}

void MainWindow::OnToolLabelStats( wxCommandEvent& event )
{
  LayerMRI* mri = (LayerMRI*)GetLayerCollection("MRI")->GetActiveLayer();
  float fLabel, fArea = 0;
  int nCount = 0;
  mri->GetCurrentLabelStats( m_nMainView, &fLabel, &nCount, &fArea );
  wxString strg = ( wxString() << "Label: " << (int)fLabel << "\nCount: " << nCount << "\nArea: " << fArea << " mm2" ); 
  wxMessageDialog dlg( this, strg, 
                       _("Label Stats"), wxOK );
  dlg.ShowModal();
}

void MainWindow::OnToolLabelStatsUpdateUI( wxUpdateUIEvent& event )
{
  LayerMRI* mri = (LayerMRI*)GetLayerCollection("MRI")->GetActiveLayer();
  event.Enable( m_nMainView != MV_3D && mri && mri->GetProperties()->GetColorMap() == LayerPropertiesMRI::LUT );
}

void MainWindow::OnFileSaveScreenshot( wxCommandEvent& event )
{
  if ( !m_dlgSaveScreenshot )
    m_dlgSaveScreenshot = new DialogSaveScreenshot( this );
  
  if ( !m_dlgSaveScreenshot->IsVisible() )
  {
    m_dlgSaveScreenshot->SetSettings( m_settingsScreenshot );
    m_dlgSaveScreenshot->SetFileName( m_strLastDir );
    m_dlgSaveScreenshot->SetFilterIndex( m_nScreenshotFilterIndex );
    m_dlgSaveScreenshot->Show();
  }
}
  
bool MainWindow::SaveScreenshot()
{
  int nId = GetActiveViewId();
  if ( nId < 0 )
    nId = m_nPrevActiveViewId;

  if ( nId < 0 )
  {
    wxMessageDialog dlg( this, _("Can not identify active view."), _("Error"), wxOK );
    dlg.ShowModal();
    return false;
  }
  
  wxString fn = m_dlgSaveScreenshot->GetFileName();
  if ( fn.IsEmpty() )
    return false;
  
  m_viewRender[nId]->Render();
  m_bDoScreenshot = true;
  return true;
}

void MainWindow::DoSaveScreenshot()
{
  wxString fn = m_dlgSaveScreenshot->GetFileName();
  m_settingsScreenshot = m_dlgSaveScreenshot->GetSettings();
  int nId = GetActiveViewId();
  if ( nId < 0 )
    nId = m_nPrevActiveViewId;
   
  if ( !fn.IsEmpty() )
  {
    m_strLastDir = MyUtils::GetNormalizedPath( fn );
    m_nScreenshotFilterIndex = m_dlgSaveScreenshot->GetFilterIndex();
    if ( !m_viewRender[nId]->SaveScreenshot( fn, 
          m_settingsScreenshot.Magnification,
          m_settingsScreenshot.AntiAliasing ) )
    {
      wxMessageDialog dlg( this, _("Error occured writing to file. Please make sure you have right permission and the disk is not full."), _("Error"), wxOK );
      dlg.ShowModal();
    }
  }
  m_bDoScreenshot = false;
}

void MainWindow::OnFileSaveScreenshotUpdateUI( wxUpdateUIEvent& event )
{
  event.Enable( ( GetActiveViewId() >= 0 || m_nPrevActiveViewId >= 0 ) && m_layerCollectionManager->HasAnyLayer() );
}


void MainWindow::OnFileLoadSurface( wxCommandEvent& event )
{
  LoadSurface();
}

void MainWindow::LoadSurface()
{
  {
    wxFileDialog dlg( this, _("Open surface file"), 
                      AutoSelectLastDir( m_strLastDir, _("surf") ),
                      _(""),
                      _("Surface files (*.*)|*.*"),
                      wxFD_OPEN | wxFD_MULTIPLE );
    if ( dlg.ShowModal() == wxID_OK )
    {
      wxArrayString fns;
      dlg.GetPaths( fns );
      for ( size_t i = 0; i < fns.GetCount(); i++ )
      {
        wxArrayString script;
        script.Add( _("loadsurface") );
        script.Add( fns[i] );
        this->AddScript( script );
      }
      ContinueScripts();
    }
  }
}

void MainWindow::OnFileSaveSurface( wxCommandEvent& event )
{
  SaveSurface();
}

void MainWindow::OnFileSaveSurfaceUpdateUI( wxUpdateUIEvent& event )
{
  LayerSurface* layer = ( LayerSurface* )GetActiveLayer( "Surface" );
  event.Enable( layer && layer->IsModified() && !IsProcessing() && !IsWritingMovieFrames() );
}


void MainWindow::OnFileSaveSurfaceAs( wxCommandEvent& event )
{
  SaveSurfaceAs();
}

void MainWindow::OnFileSaveSurfaceAsUpdateUI( wxUpdateUIEvent& event )
{
  LayerSurface* layer = ( LayerSurface* )GetActiveLayer( "Surface" );
  event.Enable( layer && layer->IsEditable() && !IsProcessing() && !IsWritingMovieFrames() );
}

void MainWindow::LoadSurfaceFile( const wxString& filename, const wxString& fn_patch, const wxString& fn_target )
{
  m_strLastDir = MyUtils::GetNormalizedPath( filename );

  LayerSurface* layer = new LayerSurface( m_layerVolumeRef );
  wxFileName fn( filename );
  wxString layerName = fn.GetFullName();
// if ( fn.GetExt().Lower() == "gz" )
//  layerName = wxFileName( layerName ).GetName();
  layer->SetName( layerName.char_str() );
  layer->SetFileName( fn.GetFullPath().char_str() );
  layer->SetPatchFileName( fn_patch.char_str() );
  layer->SetTargetFileName( fn_target.char_str() );

  WorkerThread* thread = new WorkerThread( this );
  thread->LoadSurface( layer );
}

void MainWindow::SaveSurface()
{
  // first check if there is any volume/MRI layer and if the current one is visible
  LayerSurface* layer_surf = ( LayerSurface* )GetActiveLayer( "Surface" );
  if ( !layer_surf)
  {
    return;
  }
  else if ( !layer_surf->IsVisible() )
  {
    wxMessageDialog dlg( this, _("Current surface layer is not visible. Please turn it on before saving."), 
                         _("Error"), wxOK );
    dlg.ShowModal();
    return;
  }
  wxString fn = wxString::FromAscii( layer_surf->GetFileName() );
  if ( fn.IsEmpty() )
  {
    wxString name = layer_surf->GetName(); 
    name.Trim( true ).Trim( false ).Replace( _(" "), _("_") );
    wxFileDialog dlg( this, _("Save surface file"), 
                      AutoSelectLastDir( m_strLastDir, _("surf") ),
                      name,
                      _("Surface files (*.*)|*.*"),
                      wxFD_SAVE | wxFD_OVERWRITE_PROMPT );
    if ( dlg.ShowModal() == wxID_OK )
    {
      fn = dlg.GetPath();
    }
  }

  if ( !fn.IsEmpty() )
  {
    layer_surf->SetFileName( fn.char_str() );
    WorkerThread* thread = new WorkerThread( this );
    thread->SaveSurface( layer_surf );
  }
}

void MainWindow::SaveSurfaceAs()
{
  LayerSurface* layer_surf = ( LayerSurface* )GetActiveLayer( "Surface" );
  if ( !layer_surf)
  {
    return;
  }
  else if ( !layer_surf->IsVisible() )
  {
    wxMessageDialog dlg( this, _( "Current surface layer is not visible. Please turn it on before saving." ), 
                         _( "Error" ), wxOK );
    dlg.ShowModal();
    return;
  }

  wxFileDialog dlg( this, _("Save surface file"), 
                    AutoSelectLastDir( m_strLastDir, _("surf") ),
                    layer_surf->GetFileName(),
                    _("Surface files (*.*)|*.*"),
                    wxFD_SAVE | wxFD_OVERWRITE_PROMPT );
  if ( dlg.ShowModal() == wxID_OK )
  {
    layer_surf->SetFileName( dlg.GetPath().char_str() );
    SaveSurface();
    m_controlPanel->UpdateUI();
  }
}  

void MainWindow::LoadSurfaceVector()
{
  wxFileDialog dlg( this, _("Open surface file as vector"), 
                    AutoSelectLastDir( m_strLastDir, _("surf") ),
                    _(""),
                    _("Surface files (*.*)|*.*"),
                    wxFD_OPEN );
  if ( dlg.ShowModal() == wxID_OK )
  {
    this->LoadSurfaceVectorFile( dlg.GetPath() );
//    m_strLastDir = MyUtils::GetNormalizedFullPath( dlg.GetPath() );
  }
}

void MainWindow::LoadSurfaceVectorFile( const wxString& filename )
{
  wxString fn = filename;
  if ( fn.Contains( _("/") ) )
    fn = MyUtils::GetNormalizedFullPath( filename );

  LayerSurface* layer = ( LayerSurface* )GetLayerCollection( "Surface" )->GetActiveLayer();
  if ( layer )
  {
    layer->SetVectorFileName( fn.char_str() );

    WorkerThread* thread = new WorkerThread( this );
    thread->LoadSurfaceVector( layer );   
  }
}

void MainWindow::LoadSurfaceCurvature()
{
  wxFileDialog dlg( this, _("Open curvature file"), 
                    AutoSelectLastDir( m_strLastDir, _("surf") ), 
                    _(""),
                    _("Curvature files (*.*)|*.*"),
                    wxFD_OPEN );
  if ( dlg.ShowModal() == wxID_OK )
  {
    this->LoadSurfaceCurvatureFile( dlg.GetPath() );
  }
}

void MainWindow::LoadSurfaceCurvatureFile( const wxString& filename )
{
  wxString fn = filename;
  if ( fn.Contains( _("/") ) )
    fn = MyUtils::GetNormalizedFullPath( filename );
  LayerSurface* layer = ( LayerSurface* )GetLayerCollection( "Surface" )->GetActiveLayer();
  if ( layer )
  {
    layer->LoadCurvatureFromFile( fn.char_str() );
//    m_strLastDir = MyUtils::GetNormalizedPath( filename );
  }
}


void MainWindow::LoadSurfaceOverlay()
{
  wxFileDialog dlg( this, _("Open overlay file"), 
                    AutoSelectLastDir( m_strLastDir, _("surf") ), 
                    _(""),
                    _("Overlay files (*.*)|*.*"),
                    wxFD_OPEN );
  if ( dlg.ShowModal() == wxID_OK )
  {
    this->LoadSurfaceOverlayFile( dlg.GetPath() );
  }
}

void MainWindow::LoadSurfaceOverlayFile( const wxString& filename )
{
  wxString fn = filename;
  if ( fn.Contains( _("/") ) )
    fn = MyUtils::GetNormalizedFullPath( filename );
  LayerSurface* layer = ( LayerSurface* )GetLayerCollection( "Surface" )->GetActiveLayer();
  if ( layer )
  {
    if ( MyUtils::HasExtension( fn, "mgh" ) || MyUtils::HasExtension( fn, "mgz" ) ||
         MyUtils::HasExtension( fn, "nii" ) || MyUtils::HasExtension( fn, "nii.gz" ) ||
         MyUtils::HasExtension( fn, "mnc" ) || MyUtils::HasExtension( fn, "img" ) )
      layer->LoadCorrelationFromFile( fn.char_str() );
    else
      layer->LoadOverlayFromFile( fn.char_str() );
//    m_strLastDir = MyUtils::GetNormalizedPath( filename );
  }
}

void MainWindow::LoadSurfaceAnnotation()
{
  wxFileDialog dlg( this, _("Open annotation file"), 
                    AutoSelectLastDir( m_strLastDir, _("label") ), 
                    _(""),
                    _("Annotation files (*.*)|*.*"),
                    wxFD_OPEN );
  if ( dlg.ShowModal() == wxID_OK )
  {
    this->LoadSurfaceAnnotationFile( dlg.GetPath() );
  }  
}


void MainWindow::LoadSurfaceAnnotationFile( const wxString& filename )
{
  wxString fn = filename;
  if ( fn.Contains( _("/") ) )
    fn = MyUtils::GetNormalizedFullPath( filename );
  LayerSurface* layer = ( LayerSurface* )GetLayerCollection( "Surface" )->GetActiveLayer();
  if ( layer )
  {
    layer->LoadAnnotationFromFile( fn.char_str() );
//  m_strLastDir = MyUtils::GetNormalizedPath( filename );
  }
}

void MainWindow::LoadSurfaceLabel()
{
  wxFileDialog dlg( this, _("Open label file"), 
                    AutoSelectLastDir( m_strLastDir, _("label") ), 
                    _(""),
                    _("Label files (*.*)|*.*"),
                    wxFD_OPEN );
  if ( dlg.ShowModal() == wxID_OK )
  {
    this->LoadSurfaceLabelFile( dlg.GetPath() );
  }  
}

void MainWindow::LoadSurfaceLabelFile( const wxString& filename )
{
  wxString fn = filename;
  if ( fn.Contains( _("/") ) )
    fn = MyUtils::GetNormalizedFullPath( filename );
  LayerSurface* layer = ( LayerSurface* )GetLayerCollection( "Surface" )->GetActiveLayer();
  if ( layer )
  {
    layer->LoadLabelFromFile( fn.char_str() );
    //m_strLastDir = MyUtils::GetNormalizedPath( filename );
  }
}

void MainWindow::OnToolTransformVolume( wxCommandEvent& event )
{
  if ( !m_dlgTransformVolume )
  {
    m_dlgTransformVolume = new DialogTransformVolume( this );
    GetLayerCollection( "MRI" )->AddListener( m_dlgTransformVolume );
  }

  if ( !m_dlgTransformVolume->IsVisible() )
  {
    wxMessageDialog dlg( this, _("Transformation can only apply to volumes for now. If you data includes ROI/Surface/Way Points, please do not use this feature yet."), _("Warning"), wxOK );
    dlg.ShowModal();
    m_dlgTransformVolume->Show();
    m_dlgTransformVolume->UpdateUI();
  }
}

void MainWindow::OnToolTransformVolumeUpdateUI( wxUpdateUIEvent& event )
{
// event.Check( m_dlgTransformVolume && m_dlgTransformVolume->IsShown() );
  event.Enable( !GetLayerCollection( "MRI" )->IsEmpty() && !IsProcessing() && !IsWritingMovieFrames() );
}

void MainWindow::OnToolOptimalVolume( wxCommandEvent& event )
{
  DialogOptimalVolume dlg( this, GetLayerCollection( "MRI" ) );
  if ( dlg.ShowModal() == wxID_OK )
  {
    std::vector<LayerMRI*> layers = dlg.GetSelectedLayers();
    LayerOptimal* layer_new = new LayerOptimal( layers[0] );
    layer_new->Create( dlg.GetLabelVolume(), layers );

    layer_new->SetName( dlg.GetVolumeName().char_str() );
    GetLayerCollection( "MRI" )->AddLayer( layer_new );

    m_controlPanel->RaisePage( _("Volumes") );
  }
}

void MainWindow::OnToolOptimalVolumeUpdateUI( wxUpdateUIEvent& event )
{
// event.Check( m_dlgTransformVolume && m_dlgTransformVolume->IsShown() );
  event.Enable( GetLayerCollection( "MRI" )->GetNumberOfLayers() > 1 && !IsProcessing() && !IsWritingMovieFrames() );
}

void MainWindow::OnToolGradientVolume( wxCommandEvent& event )
{
  LayerCollection* col_mri = m_layerCollectionManager->GetLayerCollection( "MRI" );  
  LayerMRI* mri = (LayerMRI*)col_mri->GetActiveLayer();
  if ( mri )
  {
    LayerMRI* layer_new = new LayerMRI( mri );
    layer_new->Create( mri, false );
    wxString name = mri->GetName();
    name += _("_gradient");
    layer_new->SetName( name.c_str() );
    VolumeFilterGradient filter( mri, layer_new );
    filter.Update();
    layer_new->ResetWindowLevel();
    col_mri->AddLayer( layer_new );
  
    m_controlPanel->RaisePage( _("Volumes") );
    
    if ( !m_dlgGradientVolume )
      m_dlgGradientVolume = new DialogGradientVolume( mri, layer_new, this );
    else
      m_dlgGradientVolume->SetVolumes( mri, layer_new );
    m_dlgGradientVolume->Show();
  }
}

void MainWindow::OnToolGradientVolumeUpdateUI( wxUpdateUIEvent& event )
{
// event.Check( m_dlgTransformVolume && m_dlgTransformVolume->IsShown() );
  event.Enable( !GetLayerCollection( "MRI" )->IsEmpty() && !IsProcessing() && !IsWritingMovieFrames() );
}


void MainWindow::EnableControls( bool bEnable )
{
  m_controlPanel->Enable( bEnable );
  if ( m_dlgTransformVolume )
    m_dlgTransformVolume->Enable( bEnable );
}

void MainWindow::OnMouseEnterWindow( wxMouseEvent& event )
{
  if ( m_viewAxial->GetInteractionMode() != RenderView2D::IM_Navigate && FindFocus() == m_toolWindowEdit )
  {
    this->Raise();
    SetFocus();
  }
}

void MainWindow::OnViewHistogram( wxCommandEvent& event )
{
  if ( m_wndHistogram->IsVisible() )
    m_wndHistogram->Hide();
  else
	  m_wndHistogram->Show( !m_wndHistogram->IsVisible() );
}

void MainWindow::OnViewHistogramUpdateUI( wxUpdateUIEvent& event )
{
  event.Check( m_wndHistogram->IsVisible() );
  event.Enable( !GetLayerCollection( "MRI" )->IsEmpty() || !GetLayerCollection( "Surface" )->IsEmpty() );
}

void MainWindow::ConfigureOverlay()
{
  m_wndOverlayConfiguration->ShowWindow( (LayerSurface*)GetLayerCollection( "Surface" )->GetActiveLayer() );
}

void MainWindow::SetVolumeColorMap( int nColorMap, int nColorMapScale, std::vector<double>& scales )
{
  if ( GetLayerCollection( "MRI" )->GetActiveLayer() )
  {
    LayerPropertiesMRI* p = ( (LayerMRI*)GetLayerCollection( "MRI" )->GetActiveLayer() )->GetProperties();
    p->SetColorMap( (LayerPropertiesMRI::ColorMapType) nColorMap );
    switch ( nColorMapScale )
    {
      case LayerPropertiesMRI::Grayscale:
        if ( scales.size() >= 2 )
        {
          p->SetMinMaxGrayscaleWindow( scales[0], scales[1] );
        }
        else if ( !scales.empty() )
          cerr << "Need 2 values for grayscale." << endl;
        break;
      case LayerPropertiesMRI::Heat:
        if ( scales.size() >= 3 )
        {
          p->SetHeatScaleMinThreshold( scales[0] );
          p->SetHeatScaleMidThreshold( scales[1] );
          p->SetHeatScaleMaxThreshold( scales[2] );
        }
        else if ( !scales.empty() )
          cerr << "Need 3 values for heatscale." << endl;
        break;
      case LayerPropertiesMRI::LUT:
        if ( scales.size() >= 1 )
        {
        }
        else if ( !scales.empty() )
          cerr << "Need a value for lut." << endl;
        break;    
      default:
        if ( scales.size() >= 2 )
        {
          p->SetMinMaxGenericThreshold( scales[0], scales[1] );
        }
        else if ( !scales.empty() )
          cerr << "Need 2 values for colorscale." << endl;
        break;
    }
  }
}

void MainWindow::LoadLUT()
{
  wxFileDialog dlg( this, _("Load lookup table file"), m_strLastDir, _(""),
                  _("LUT files (*.*)|*.*"),
                  wxFD_OPEN );
  if ( dlg.ShowModal() == wxID_OK )
  {
    m_scripts.clear();
    wxArrayString sa;
    sa.Add( _("setlut") );
    sa.Add( dlg.GetPath().c_str() );
    CommandSetLUT( sa );
  }
}

bool MainWindow::GetCursorRAS( double* ras_out, bool tkReg )
{
  LayerCollection* lc_mri = GetLayerCollection( "MRI" );
  LayerCollection* lc_surf = GetLayerCollection( "Surface" );
  
  if ( !lc_mri->IsEmpty() )
  {
    LayerMRI* mri = ( LayerMRI* )lc_mri->GetLayer( 0 );
    mri->RemapPositionToRealRAS( lc_mri->GetCursorRASPosition(), ras_out );
    if ( tkReg )
      mri->NativeRASToTkReg( ras_out, ras_out );
    return true;
  }
  else if ( !lc_surf->IsEmpty() )
  {
    lc_surf->GetCurrentRASPosition( ras_out );
    return true;
  }
  else
    return false;
}

void MainWindow::OnFileSaveMovieFrames( wxCommandEvent& event )
{  
  if ( !m_dlgWriteMovieFrames->IsVisible() )
  {
    m_dlgWriteMovieFrames->Show();
  }
}

void MainWindow::OnFileSaveMovieFramesUpdateUI( wxUpdateUIEvent& event )
{
  event.Enable( !m_bProcessing && 
      ( !GetLayerCollection( "MRI" )->IsEmpty() || !GetLayerCollection( "Surface" )->IsEmpty() ) ); 
}

void MainWindow::OnTimerWriteMovieFrames( wxTimerEvent& event )
{
  double angle_step = m_settingsMovieFrames.AngleStep;
  wxString fn;
  fn.Printf( "%sframe%03d.%s", 
             m_settingsMovieFrames.OutputLocation.c_str(),
             m_settingsMovieFrames.StepCount,
             m_settingsMovieFrames.OutputExtension.c_str() );
  if ( m_nMainView == 3 ) // 3D view
  {
    m_view3D->SaveScreenshot( fn, 
                              m_settingsScreenshot.Magnification,
                              m_settingsScreenshot.AntiAliasing );
    m_view3D->Azimuth( angle_step );
    ForceRedraw( );
    m_settingsMovieFrames.StepCount ++;
    double angle = m_settingsMovieFrames.StepCount * fabs(angle_step);
    m_statusBar->m_gaugeBar->SetValue( (int)(100*angle/360) );
    if ( angle > 360 )
      StopWriteMovieFrames();
  }
  else
  {
    m_viewRender[m_nMainView]->SaveScreenshot( fn, 
                                               m_settingsScreenshot.Magnification,
                                               m_settingsScreenshot.AntiAliasing );
    m_settingsMovieFrames.StepCount ++;
    if ( ((RenderView2D*)m_viewRender[m_nMainView])->SetSliceNumber( m_settingsMovieFrames.StepCount ) )
    {
      ForceRedraw();
      LayerMRI* mri = (LayerMRI*)GetActiveLayer( "MRI" );
      if ( mri )
      {
        double nMaxSlice = (mri->GetWorldSize())[m_nMainView] / (mri->GetWorldVoxelSize())[m_nMainView];
        m_statusBar->m_gaugeBar->SetValue( (int)(100*m_settingsMovieFrames.StepCount / nMaxSlice) );
      }
    }
    else
      StopWriteMovieFrames();
  }
}

void MainWindow::StartWriteMovieFrames()
{
  m_statusBar->m_gaugeBar->Show();
  m_settingsMovieFrames.StepCount = 0;
  if ( m_nMainView != 3 )
  {
    ((RenderView2D*)m_viewRender[m_nMainView])->SetSliceNumber( 0 );
    ForceRedraw();
  }
  if ( !m_dlgWriteMovieFrames->GetOutputLocation().IsEmpty() )
  {
    m_settingsMovieFrames.OutputLocation = m_dlgWriteMovieFrames->GetOutputLocation(); 
    if ( m_settingsMovieFrames.OutputLocation[m_settingsMovieFrames.OutputLocation.Len()-1] != '/' )
      m_settingsMovieFrames.OutputLocation += _("/");
  }
  if ( m_dlgWriteMovieFrames->GetAngleStep() != 0 )
  {
    m_settingsMovieFrames.AngleStep = m_dlgWriteMovieFrames->GetAngleStep();
  }
  if ( !wxDir::Exists( m_settingsMovieFrames.OutputLocation ) )
  {
    wxString cmd = "mkdir ";
    cmd += m_settingsMovieFrames.OutputLocation;
    system( cmd.c_str() );
  }
  m_settingsMovieFrames.OutputExtension = m_dlgWriteMovieFrames->GetOutputExtension();
  m_statusBar->m_gaugeBar->SetValue( 0 );
  m_timerWriteMovieFrames.Start( 1000 ); 
}

void MainWindow::StopWriteMovieFrames()
{
  m_timerWriteMovieFrames.Stop();
  m_settingsMovieFrames.StepCount = 0;
  m_statusBar->m_gaugeBar->Hide();
  m_dlgWriteMovieFrames->UpdateUI();
}

void MainWindow::OnFilterMean( wxCommandEvent& event )
{  
  LayerMRI* mri = (LayerMRI*)GetLayerCollection( "MRI" )->GetActiveLayer();
  if ( mri )
  {
    VolumeFilterMean* filter = new VolumeFilterMean( mri, mri );
    DialogVolumeFilter* dlg = new DialogVolumeFilter( this );
    dlg->SetFilter( filter );
    dlg->ShowSigma( false );
    if ( dlg->ShowModal() == wxID_OK )
    {
      filter->SetKernelSize( dlg->GetKernelSize() );
      filter->Update();
    }
    delete filter;
  }
}

void MainWindow::OnFilterUpdateUI( wxUpdateUIEvent& event )
{
  Layer* layer = GetLayerCollection( "MRI" )->GetActiveLayer();
  event.Enable( !m_bProcessing && layer && layer->IsVisible() ); 
}

void MainWindow::OnFilterMedian( wxCommandEvent& event )
{  
  LayerMRI* mri = (LayerMRI*)GetLayerCollection( "MRI" )->GetActiveLayer();
  if ( mri )
  {
    /*
    Chronometer tSample;
    InitChronometer( &tSample );
    StartChronometer( &tSample );
    */
    
    VolumeFilterMedian* filter = new VolumeFilterMedian( mri, mri );
    DialogVolumeFilter* dlg = new DialogVolumeFilter( this );
    dlg->ShowSigma( false );
    dlg->SetFilter( filter );
    if ( dlg->ShowModal() == wxID_OK )
    {
      filter->SetKernelSize( dlg->GetKernelSize() );
      filter->Update();
    }
    delete filter;
    
    /*
    StopChronometer( &tSample );
    printf( "  tSample : %9.3f ms\n", GetChronometerValue( &tSample ) );
    fflush( 0 );
    */
  }  
}

void MainWindow::OnFilterConvolve( wxCommandEvent& event )
{  
  LayerMRI* mri = (LayerMRI*)GetLayerCollection( "MRI" )->GetActiveLayer();
  if ( mri )
  {
    VolumeFilterConvolve* filter = new VolumeFilterConvolve( mri, mri );
    DialogVolumeFilter* dlg = new DialogVolumeFilter( this );
    dlg->SetFilter( filter );
    dlg->SetSigma( filter->GetSigma() );
    if ( dlg->ShowModal() == wxID_OK )
    {
      filter->SetKernelSize( dlg->GetKernelSize() );
      filter->SetSigma( dlg->GetSigma() );
      filter->Update();
    }
    delete filter;
  }  
}

void MainWindow::OnFilterGradient( wxCommandEvent& event )
{  
  LayerMRI* mri = (LayerMRI*)GetLayerCollection( "MRI" )->GetActiveLayer();
  if ( mri )
  {
    VolumeFilterGradient* filter = new VolumeFilterGradient( mri, mri );
    filter->Update();
    mri->ResetWindowLevel();
    delete filter;
  }  
}

void MainWindow::OnToolRepositionSurface( wxCommandEvent& event )
{  
  m_dlgRepositionSurface->Show();
}

void MainWindow::OnToolRepositionSurfaceUpdateUI( wxUpdateUIEvent& event )
{
  event.Enable( !m_bProcessing && GetLayerCollection( "Surface" )->GetActiveLayer() && GetLayerCollection( "MRI" )->GetActiveLayer() ); 
}

void MainWindow::OnEditRename( wxCommandEvent& event )
{
  LayerCollection* lc = GetCurrentLayerCollection();
  if ( lc && lc->GetActiveLayer() )
  {
    Layer* layer = lc->GetActiveLayer();
    wxString new_name = ::wxGetTextFromUser( _("Enter new layer name"), _("Rename Layer"), layer->GetName() );
    if ( !new_name.IsEmpty() )
      layer->SetName( new_name.c_str() );
  }  
}

void MainWindow::OnEditRenameUpdateUI( wxUpdateUIEvent& event )
{
  LayerCollection* lc = GetCurrentLayerCollection();
  event.Enable( lc && lc->GetActiveLayer() );
}

void MainWindow::OnToolCropVolume( wxCommandEvent& event )
{
  LayerMRI* mri = (LayerMRI*)GetLayerCollection( "MRI" )->GetActiveLayer();
  m_dlgCropVolume->SetVolume( mri );
  m_dlgCropVolume->Show();
  m_volumeCropper->SetEnabled( true );
  m_volumeCropper->SetVolume( mri );
  m_volumeCropper->Show();
  SetMode( RenderView::IM_VolumeCrop );
}

void MainWindow::OnToolCropVolumeUpdateUI( wxUpdateUIEvent& event )
{
  Layer* layer = GetLayerCollection( "MRI" )->GetActiveLayer();
  event.Enable( !IsProcessing() && layer && layer->IsVisible() );
}

void MainWindow::SaveRegistrationAs()
{
  LayerMRI* layer_mri = ( LayerMRI* )GetActiveLayer( "MRI" );
  if ( !layer_mri)
  {
    return;
  }
  
  wxFileDialog dlg( this, _("Select file to save"), 
                    wxFileName( layer_mri->GetFileName() ).GetPath(), 
                    _(""),
                    _("LTA files (*.lta)|*.lta|All files (*.*)|*.*"),
                    wxFD_SAVE | wxFD_OVERWRITE_PROMPT );
  if ( dlg.ShowModal() == wxID_OK )
  {
    layer_mri->SaveRegistration( dlg.GetPath().char_str() );
  }
}

void MainWindow::OnLayerShowAll( wxCommandEvent& event )
{
  LayerCollection* lc = GetCurrentLayerCollection();
  if ( lc )
  {
    for ( int i = 0; i < lc->GetNumberOfLayers(); i++ )
      lc->GetLayer( i )->Show();
  }
}

void MainWindow::OnLayerHideAll( wxCommandEvent& event )
{
  LayerCollection* lc = GetCurrentLayerCollection();
  if ( lc )
  {
    for ( int i = 0; i < lc->GetNumberOfLayers(); i++ )
      lc->GetLayer( i )->Hide();
  }
}

void MainWindow::OnLayerShowHideAllUpdateUI( wxUpdateUIEvent& event )
{
  LayerCollection* lc = GetCurrentLayerCollection();
  event.Enable( lc && lc->GetNumberOfLayers() > 0 );
}
