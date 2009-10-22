/**
 * @file  MainWindow.cpp
 * @brief Main window.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2009/10/22 18:29:45 $
 *    $Revision: 1.75 $
 *
 * Copyright (C) 2008-2009,
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
#include "DialogRotateVolume.h"
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

#define CTRL_PANEL_WIDTH 240

#define ID_FILE_RECENT1     10001

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
  EVT_MENU        ( XRCID( "ID_FILE_LOAD_SURFACE" ),      MainWindow::OnFileLoadSurface )
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
  EVT_MENU        ( XRCID( "ID_VIEW_SCALAR_BAR" ),        MainWindow::OnViewScalarBar )
  EVT_UPDATE_UI   ( XRCID( "ID_VIEW_SCALAR_BAR" ),        MainWindow::OnViewScalarBarUpdateUI )
  EVT_MENU        ( XRCID( "ID_VIEW_COORDINATE" ),        MainWindow::OnViewCoordinate )
  EVT_UPDATE_UI   ( XRCID( "ID_VIEW_COORDINATE" ),        MainWindow::OnViewCoordinateUpdateUI )
  EVT_MENU        ( XRCID( "ID_VIEW_CYCLE_LAYER" ),       MainWindow::OnViewCycleLayer )
  EVT_UPDATE_UI   ( XRCID( "ID_VIEW_CYCLE_LAYER" ),       MainWindow::OnViewCycleLayerUpdateUI )
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
  
  EVT_MENU        ( XRCID( "ID_TOOL_ROTATE_VOLUME" ),     MainWindow::OnToolRotateVolume )
  EVT_UPDATE_UI   ( XRCID( "ID_TOOL_ROTATE_VOLUME" ),     MainWindow::OnToolRotateVolumeUpdateUI )
  EVT_MENU        ( XRCID( "ID_TOOL_OPTIMAL_VOLUME" ),    MainWindow::OnToolOptimalVolume )
  EVT_UPDATE_UI   ( XRCID( "ID_TOOL_OPTIMAL_VOLUME" ),    MainWindow::OnToolOptimalVolumeUpdateUI )
  EVT_MENU        ( XRCID( "ID_TOOL_GRADIENT_VOLUME" ),   MainWindow::OnToolGradientVolume )
  EVT_UPDATE_UI   ( XRCID( "ID_TOOL_GRADIENT_VOLUME" ),   MainWindow::OnToolGradientVolumeUpdateUI )
  
  EVT_MENU  ( XRCID( "ID_HELP_QUICK_REF" ),               MainWindow::OnHelpQuickReference )
  EVT_MENU  ( XRCID( "ID_HELP_ABOUT" ),                   MainWindow::OnHelpAbout )
  /*
  EVT_SASH_DRAGGED_RANGE(ID_LOG_WINDOW, ID_LOG_WINDOW, MainWindow::OnSashDrag)
  EVT_IDLE  (MainWindow::OnIdle)
  EVT_MENU  (XRCID("ID_EVENT_LOAD_DATA"), MainWindow::OnWorkerEventLoadData)
  EVT_MENU  (XRCID("ID_EVENT_CALCULATE_MATRIX"), MainWindow::OnWorkerEventCalculateMatrix)
  */
  
  EVT_MENU      ( ID_WORKER_THREAD,                   MainWindow::OnWorkerThreadResponse )
  EVT_SPINCTRL  ( XRCID( "ID_SPIN_BRUSH_SIZE" ),      MainWindow::OnSpinBrushSize )
  EVT_SPINCTRL  ( XRCID( "ID_SPIN_BRUSH_TOLERANCE" ), MainWindow::OnSpinBrushTolerance )
  EVT_CHECKBOX  ( XRCID( "ID_CHECK_TEMPLATE" ),       MainWindow::OnCheckBrushTemplate )
  EVT_CHOICE    ( XRCID( "ID_CHOICE_TEMPLATE" ),      MainWindow::OnChoiceBrushTemplate )
  
  EVT_ENTER_WINDOW  ( MainWindow::OnMouseEnterWindow )
  EVT_ACTIVATE      ( MainWindow::OnActivate )
  EVT_CLOSE         ( MainWindow::OnClose )
  EVT_KEY_DOWN      ( MainWindow::OnKeyDown )
  EVT_ICONIZE       ( MainWindow::OnIconize )
END_EVENT_TABLE()

// ----------------------------------------------------------------------------
// main frame
// ----------------------------------------------------------------------------

// frame constructor
MainWindow::MainWindow() : Listener( "MainWindow" ), Broadcaster( "MainWindow" )
{
// m_bLoading = false;
// m_bSaving = false;
  m_bProcessing = false;
  m_bResampleToRAS = false;
  m_bToUpdateToolbars = false;
  m_layerVolumeRef = NULL;
  m_nPrevActiveViewId = -1;
  m_luts = new LUTDataHolder();
  m_propertyBrush = new BrushProperty();
  m_connectivity = new ConnectivityData();

  wxXmlResource::Get()->LoadFrame( this, NULL, wxT("ID_MAIN_WINDOW") );

  // must be created before the following controls
  m_layerCollectionManager = new LayerCollectionManager();

  m_panelToolbarHolder = XRCCTRL( *this, "ID_PANEL_HOLDER", wxPanel );
// wxBoxSizer* sizer = (wxBoxSizer*)panelHolder->GetSizer(); //new wxBoxSizer( wxVERTICAL );

  // create the main splitter window
// m_splitterMain = new wxSplitterWindow( this, wxID_ANY, wxDefaultPosition, wxDefaultSize, wxSP_LIVE_UPDATE );
  m_splitterMain = XRCCTRL( *this, "ID_SPLITTER_MAIN", wxSplitterWindow );
  m_splitterMain->SetMinimumPaneSize( 80 );
// sizer->Add( m_splitterMain, 1, wxEXPAND );

  m_toolbarVoxelEdit = XRCCTRL( *this, "ID_TOOLBAR_VOXEL_EDIT", wxToolBar );
  m_toolbarROIEdit = XRCCTRL( *this, "ID_TOOLBAR_ROI_EDIT", wxToolBar );
  m_toolbarBrush = XRCCTRL( *this, "ID_TOOLBAR_BRUSH", wxToolBar );
  m_toolbarBrush->Show( false );
  m_toolbarVoxelEdit->Show( false );
  m_toolbarROIEdit->Show( false );

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
  m_viewCoronal = new RenderView2D( m_renderViewHolder, wxID_ANY );
  m_viewSagittal = new RenderView2D( m_renderViewHolder, wxID_ANY );
  m_viewAxial = new RenderView2D( m_renderViewHolder, wxID_ANY );
  m_view3D = new RenderView3D( m_renderViewHolder, wxID_ANY );
  m_viewSagittal->SetViewPlane( 0 );
  m_viewCoronal->SetViewPlane( 1 );
  m_viewAxial->SetViewPlane( 2 );
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
  m_layerCollectionManager->GetLayerCollection( "MRI" )->AddListener( m_pixelInfoPanel );
  m_layerCollectionManager->GetLayerCollection( "Surface" )->AddListener( m_pixelInfoPanel );
  m_layerCollectionManager->AddListener( this );

  m_connectivity->AddListener( m_view3D );
  
  m_wndQuickReference = new WindowQuickReference( this );
  m_wndQuickReference->Hide();

  m_statusBar = new StatusBar( this );
  SetStatusBar( m_statusBar );
  PositionStatusBar();

  m_toolWindowEdit = NULL;
  m_dlgRotateVolume = NULL;
  m_dlgGradientVolume = NULL;

  m_wndHistogram = new WindowHistogram( this );
  m_wndHistogram->Hide();
  
  m_wndOverlayConfiguration = new WindowOverlayConfiguration( this );
  m_wndOverlayConfiguration->Hide();
  
  m_wndConnectivityConfiguration = new WindowConnectivityConfiguration( this );
  m_wndConnectivityConfiguration->Hide();
  
  UpdateToolbars();
  
  m_nViewLayout = VL_2X2;
  m_nMainView = MV_Sagittal;
  m_nMaxRecentFiles = 10;
  wxConfigBase* config = wxConfigBase::Get();
  m_fileHistory = new wxFileHistory( m_nMaxRecentFiles, wxID_FILE1 );
  wxMenu* fileMenu = GetMenuBar()->GetMenu( 0 )->FindItem( XRCID("ID_FILE_SUBMENU_RECENT") )->GetSubMenu();
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

    ShowControlPanel( bShow );
  }
  SetViewLayout( m_nViewLayout );
}

// frame destructor
MainWindow::~MainWindow()
{
  m_viewAxial->Delete();
  m_viewSagittal->Delete();
  m_viewCoronal->Delete();
  m_view3D->Delete();

  delete m_layerCollectionManager;
  delete m_luts;
  delete m_propertyBrush;
  if ( m_connectivity )
    delete m_connectivity;
}

MainWindow* MainWindow::GetMainWindowPointer()
{
  return wxDynamicCast( wxTheApp->GetTopWindow(), MainWindow );
}

void MainWindow::OnClose( wxCloseEvent &event )
{
  if ( IsProcessing() )
  {
    wxMessageDialog dlg( this, _("There is on-going data processing. If you force quit, any data that is being saved can be lost. Do you really want to force quit?"), 
                         _("Force Quit"), wxYES_NO | wxNO_DEFAULT  );
    if (dlg.ShowModal() == wxID_NO )
      return;
  }

  LayerCollection* lc_mri = GetLayerCollection( "MRI" );
  LayerCollection* lc_roi = GetLayerCollection( "ROI" );
  LayerCollection* lc_wp = GetLayerCollection( "WayPoints" );
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
  for ( int i = 0; i < lc_wp->GetNumberOfLayers(); i++ )
  {
    LayerEditable* layer = ( LayerEditable* )lc_wp->GetLayer( i );
    if ( layer->IsModified() )
      text += wxString::FromAscii( layer->GetName() ) + _( "\t(" ) + 
          wxFileName( wxString::FromAscii( layer->GetFileName() ) ).GetShortPath() + _(")\n");
  }

  if ( !text.IsEmpty() )
  {
    wxString msg = _("The following volume(s) and/or label(s) have been modified but not saved. \n\n");
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
    
    config->Write( _("/RenderWindow/ShowScalarBar"), 
                   ( m_nMainView >= 0 ? m_viewRender[m_nMainView]->GetShowScalarBar() : false ) );

    config->Write( _("/Screenshot/Magnification" ),  m_settingsScreenshot.Magnification );
    config->Write( _("/Screenshot/HideCursor" ),  m_settingsScreenshot.HideCursor );
    config->Write( _("/Screenshot/HideCoords" ),  m_settingsScreenshot.HideCoords );
    config->Write( _("/Screenshot/AntiAliasing" ),  m_settingsScreenshot.AntiAliasing );

    config->Write( _("/MainWindow/ShowControlPanel" ), m_controlPanel->IsShown() );

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
  layer_new->Create( dlg.GetTemplate(), dlg.GetCopyVoxel() );
  layer_new->SetName( dlg.GetVolumeName().char_str() );
  col_mri->AddLayer( layer_new );

  m_controlPanel->RaisePage( _("Volumes") );

// m_viewAxial->SetInteractionMode( RenderView2D::IM_ROIEdit );
// m_viewCoronal->SetInteractionMode( RenderView2D::IM_ROIEdit );
// m_viewSagittal->SetInteractionMode( RenderView2D::IM_ROIEdit );
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
    wxFileDialog dlg( this, _("Save volume file"), m_strLastDir, _(""),
                      _("Volume files (*.nii;*.nii.gz;*.img;*.mgz)|*.nii;*.nii.gz;*.img;*.mgz|All files (*.*)|*.*"),
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
          !MyUtils::HasExtension( fn, _("mgz") )
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

  DialogSaveVolumeAs dlg( this );
  dlg.SetFileName( layer_mri->GetFileName() );
  if ( dlg.ShowModal() == wxID_OK )
  {
    layer_mri->SetFileName( dlg.GetFileName().char_str() );
    layer_mri->SetReorient( dlg.GetReorient() );
    SaveVolume();
    m_controlPanel->UpdateUI();
  }
}

void MainWindow::LoadVolume()
{
  DialogLoadVolume dlg( this, GetLayerCollection( "MRI" )->IsEmpty() );
  dlg.SetLastDir( m_strLastDir );
  wxArrayString list;
  for ( int i = 0; i < m_fileHistory->GetMaxFiles(); i++ )
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
  
      if ( ( GetLayerCollection( "MRI" )->IsEmpty() && dlg.IsToResample() ) || m_bResampleToRAS )
        script.Add( _("r") );
  
      AddScript( script );
    }
    RunScript();
  }
}

void MainWindow::LoadVolumeFile( const wxString& filename, 
                                 const wxString& reg_filename, 
                                 bool bResample, int nSampleMethod )
{
// cout << bResample << endl;
  m_strLastDir = MyUtils::GetNormalizedPath( filename );

  m_bResampleToRAS = bResample;
  LayerMRI* layer = new LayerMRI( m_layerVolumeRef );
  layer->SetResampleToRAS( bResample );
  layer->SetSampleMethod( nSampleMethod );
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
  WorkerThread* thread = new WorkerThread( this );
  thread->RotateVolume( rotations, bAllVolumes );
}


void MainWindow::OnFileExit( wxCommandEvent& event )
{
  Close();
}

void MainWindow::OnActivate( wxActivateEvent& event )
{
#ifdef __WXGTK__
  NeedRedraw();
#endif
  event.Skip();
}

void MainWindow::OnIconize( wxIconizeEvent& event )
{
#ifdef __WXGTK__
  if ( !event.Iconized() )
    NeedRedraw();
#endif
  event.Skip();
}


void MainWindow::OnFileRecent( wxCommandEvent& event )
{
  wxString fn( m_fileHistory->GetHistoryFile( event.GetId() - wxID_FILE1 ) );
  if ( !fn.IsEmpty() )
    this->LoadVolumeFile( fn, _(""), m_bResampleToRAS );
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
  wxFileDialog dlg( this, _("Open ROI file"), m_strLastDir, _(""),
                    _("Label files (*.label)|*.label|All files (*.*)|*.*"),
                    wxFD_OPEN );
  if ( dlg.ShowModal() == wxID_OK )
  {
    this->LoadROIFile( dlg.GetPath() );
  }
}

void MainWindow::LoadROIFile( const wxString& fn )
{
  m_strLastDir = MyUtils::GetNormalizedPath( fn );

  LayerCollection* col_mri = GetLayerCollection( "MRI" );
  LayerMRI* mri = ( LayerMRI* )col_mri->GetActiveLayer();
  LayerROI* roi = new LayerROI( mri );
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

  m_viewAxial->SetInteractionMode( RenderView2D::IM_ROIEdit );
  m_viewCoronal->SetInteractionMode( RenderView2D::IM_ROIEdit );
  m_viewSagittal->SetInteractionMode( RenderView2D::IM_ROIEdit );
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
    wxFileDialog dlg( this, _("Save ROI file"), m_strLastDir, _(""),
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
  wxFileDialog dlg( this, _("Save ROI file as"), m_strLastDir, fn,
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

void MainWindow::LoadWayPointsFile( const wxString& fn )
{
  m_strLastDir = MyUtils::GetNormalizedPath( fn );

  LayerCollection* col_mri = GetLayerCollection( "MRI" );
  LayerMRI* mri = ( LayerMRI* )col_mri->GetActiveLayer();
  LayerWayPoints* wp = new LayerWayPoints( mri );
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

    m_controlPanel->RaisePage( _("Way Points") );
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

  // enter the name of the new ROI
  DialogNewWayPoints dlg( this, col_mri );
  dlg.SetWayPointsName( _("New Way Points") );
  if ( dlg.ShowModal() != wxID_OK )
    return;

  // finally we are about to create new ROI.
  LayerCollection* col_wp = m_layerCollectionManager->GetLayerCollection( "WayPoints" );
  if ( col_wp->IsEmpty() )
  {
    col_wp->SetWorldOrigin( col_mri->GetWorldOrigin() );
    col_wp->SetWorldSize( col_mri->GetWorldSize() );
    col_wp->SetWorldVoxelSize( col_mri->GetWorldVoxelSize() );
    col_wp->SetSlicePosition( col_mri->GetSlicePosition() );
  }
  LayerWayPoints* layer_wp = new LayerWayPoints( dlg.GetTemplate() );
  layer_wp->SetName( dlg.GetWayPointsName().char_str() );
  col_wp->AddLayer( layer_wp );

  m_controlPanel->RaisePage( _("Way Points") );

  m_viewAxial->SetInteractionMode( RenderView2D::IM_WayPointsEdit );
  m_viewCoronal->SetInteractionMode( RenderView2D::IM_WayPointsEdit );
  m_viewSagittal->SetInteractionMode( RenderView2D::IM_WayPointsEdit );
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
    wxFileDialog dlg( this, _("Save Way Points file"), m_strLastDir, _(""),
                      _("Way Points files (*.label)|*.label|All files (*.*)|*.*"),
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
    wxMessageDialog dlg( this, _("Current Way Points layer is not visible. Please turn it on before saving."), 
                         _("Error"), wxOK );
    dlg.ShowModal();
    return;
  }

  wxString fn = wxString::FromAscii( layer_wp->GetFileName() );
  wxFileDialog dlg( this, _("Save Way Points file as"), m_strLastDir, fn,
                    _("Way Points files (*.label)|*.label|All files (*.*)|*.*"),
                    wxFD_SAVE | wxFD_OVERWRITE_PROMPT );
  if ( dlg.ShowModal() == wxID_OK )
  {
    layer_wp->SetFileName( dlg.GetPath().char_str() );
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
  if ( !m_toolWindowEdit )
    m_toolWindowEdit = new ToolWindowEdit( this );
    
  if ( m_viewAxial->GetInteractionMode() == RenderView2D::IM_VoxelEdit ||
       m_viewAxial->GetInteractionMode() == RenderView2D::IM_ROIEdit )
  {
    m_toolWindowEdit->Show();
  }
  else
    m_toolWindowEdit->Hide();

  m_toolWindowEdit->UpdateTools();

  m_bToUpdateToolbars = false;
}

void MainWindow::SetMode( int nMode )
{
  m_viewAxial->SetInteractionMode( nMode );
  m_viewCoronal->SetInteractionMode( nMode );
  m_viewSagittal->SetInteractionMode( nMode );

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

  NeedRedraw();
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

void MainWindow::OnViewCycleLayer( wxCommandEvent& event )
{
  LayerCollection* lc = NULL;
  switch ( m_controlPanel->GetCurrentLayerCollectionIndex() )
  {
  case 0:   // Volume
    lc = GetLayerCollection( "MRI" );
    break;
  case 1:   // ROI
    lc = GetLayerCollection( "ROI" );
    break;
  case 2:
    lc = GetLayerCollection( "Surface" );
  }

  if ( lc )
  {
    lc->CycleLayer();
  }
}

void MainWindow::OnViewCycleLayerUpdateUI( wxUpdateUIEvent& event )
{
  LayerCollection* lc = NULL;
  switch ( m_controlPanel->GetCurrentLayerCollectionIndex() )
  {
  case 0:   // Volume
    lc = GetLayerCollection( "MRI" );
    break;
  case 1:   // ROI
    lc = GetLayerCollection( "ROI" );
    break;
  }

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

void MainWindow::OnInternalIdle()
{
  wxFrame::OnInternalIdle();

#ifdef __WXGTK__
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
  event.Enable( layer && layer->IsModified() && !IsProcessing() );
}


void MainWindow::OnFileSaveAs( wxCommandEvent& event )
{
  SaveVolumeAs();
}

void MainWindow::OnFileSaveAsUpdateUI( wxUpdateUIEvent& event )
{
  LayerMRI* layer = ( LayerMRI* )( GetLayerCollection( "MRI" )->GetActiveLayer() );
  event.Enable( layer && layer->IsEditable() && !IsProcessing() );
}

void MainWindow::OnFileLoadROI( wxCommandEvent& event )
{
  LoadROI();
}

void MainWindow::OnFileLoadROIUpdateUI( wxUpdateUIEvent& event )
{
  event.Enable( !GetLayerCollection( "MRI" )->IsEmpty() && !IsProcessing() );
}


void MainWindow::OnFileSaveROI( wxCommandEvent& event )
{
  SaveROI();
}

void MainWindow::OnFileSaveROIUpdateUI( wxUpdateUIEvent& event )
{
  LayerROI* layer = ( LayerROI* )( GetLayerCollection( "ROI" )->GetActiveLayer() );
  event.Enable( layer && layer->IsModified() && !IsProcessing() );
}


void MainWindow::OnFileSaveROIAs( wxCommandEvent& event )
{
  SaveROIAs();
}

void MainWindow::OnFileSaveROIAsUpdateUI( wxUpdateUIEvent& event )
{
  LayerROI* layer = ( LayerROI* )( GetLayerCollection( "ROI" )->GetActiveLayer() );
  event.Enable( layer && !IsProcessing() );
}


void MainWindow::OnFileLoadWayPoints( wxCommandEvent& event )
{
  LoadWayPoints();
}

void MainWindow::OnFileLoadWayPointsUpdateUI( wxUpdateUIEvent& event )
{
  event.Enable( !GetLayerCollection( "MRI" )->IsEmpty() && !IsProcessing() );
}


void MainWindow::OnFileSaveWayPoints( wxCommandEvent& event )
{
  SaveWayPoints();
}

void MainWindow::OnFileSaveWayPointsUpdateUI( wxUpdateUIEvent& event )
{
  LayerWayPoints* layer = ( LayerWayPoints* )( GetLayerCollection( "WayPoints" )->GetActiveLayer() );
  event.Enable( layer && layer->IsModified() && !IsProcessing() );
}


void MainWindow::OnFileSaveWayPointsAs( wxCommandEvent& event )
{
  SaveWayPointsAs();
}

void MainWindow::OnFileSaveWayPointsAsUpdateUI( wxUpdateUIEvent& event )
{
  LayerWayPoints* layer = ( LayerWayPoints* )( GetLayerCollection( "WayPoints" )->GetActiveLayer() );
  event.Enable( layer && !IsProcessing() );
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

        // m_fileHistory->AddFileToHistory( MyUtils::GetNormalizedFullPath( layer->GetFileName() ) );

        m_controlPanel->RaisePage( _("Surfaces") );
      }

      m_viewAxial->SetInteractionMode( RenderView2D::IM_Navigate );
      m_viewCoronal->SetInteractionMode( RenderView2D::IM_Navigate );
      m_viewSagittal->SetInteractionMode( RenderView2D::IM_Navigate );
    }

    // Saving operation finished
    else if ( strg == _("Save") )
    {
      cout << ( (LayerEditable*)layer )->GetFileName() << " saved successfully." << endl;
    }

    else if ( strg == _("Rotate") )
    {
      m_bResampleToRAS = false;
      m_layerCollectionManager->RefreshSlices();
    }

    m_controlPanel->UpdateUI();

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
    m_statusBar->m_gaugeBar->SetValue( event.GetInt() );
  }
}

void MainWindow::DoListenToMessage ( std::string const iMsg, void* iData, void* sender )
{
  if ( iMsg == "LayerAdded" || iMsg == "LayerRemoved" || iMsg == "LayerMoved" )
  {
    UpdateToolbars();
    if ( !GetLayerCollection( "MRI" )->IsEmpty() )
    {
      SetTitle( wxString::FromAscii( GetLayerCollection( "MRI" )->GetLayer( 0 )->GetName() ) + _(" - freeview") );
    }
    else if ( !GetLayerCollection( "Surface" )->IsEmpty() )
    {
      SetTitle( wxString::FromAscii( GetLayerCollection( "Surface" )->GetLayer( 0 )->GetName() ) + _(" - freeview") );
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
  int nSampleMethod = SAMPLE_NEAREST;
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
                subOption == _("jetscale") ) 
      {
        colormap_scale = subOption;    // colormap scale might be different from colormap!
        scales = MyUtils::SplitString( subArgu, _(",") );
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
        if ( subArgu.Lower() == _("trilinear") )
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
      else if ( subOption == _("name") )
      {
        wxArrayString script;
        script.Add( _("setlayername") );
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
  
  LoadVolumeFile( fn, reg_fn, bResample, nSampleMethod );
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
  else if ( strg != _("grayscale") )
    cerr << "Unrecognized colormap name '" << strg << "'." << endl;
  
  int nColorMapScale = LayerPropertiesMRI::Grayscale;
  strg = sa[2];
  if ( strg == _("heatscale") )
    nColorMapScale = LayerPropertiesMRI::Heat;
  else if ( strg == _("jetscale") )
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

void MainWindow::CommandSetLayerName( const wxArrayString& cmd )
{
  if ( cmd.size() > 2 )
  {
    LayerCollection* lc = GetLayerCollection( cmd[1].c_str() );
    if ( lc && !lc->IsEmpty() )
    {
      lc->GetActiveLayer()->SetName( cmd[2].c_str() );
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
      else
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

void MainWindow::ContinueScripts()
{
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
  LoadROIFile( cmd[1] );

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
        fn_patch = subArgu;
      }  
      else if ( subOption == _("name") )
      {
        wxArrayString script;
        script.Add( _("setlayername") );
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
  LoadSurfaceFile( fn, fn_patch );
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
    
    if ( cmd[2] != _("null") )
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
    
    if ( cmd[2] != _("0") )
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
      layer->RASToTarget( ras, ras );
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

void MainWindow::OnFileSaveScreenshot( wxCommandEvent& event )
{
  int nId = GetActiveViewId();
  if ( nId < 0 )
    nId = m_nPrevActiveViewId;

  if ( nId < 0 )
    return;

  wxString fn;
  wxFileDialog dlg( this, _("Save screenshot as"), m_strLastDir, _(""),
                    _("PNG files (*.png)|*.png|JPEG files (*.jpg;*.jpeg)|*.jpg;*.jpeg|TIFF files (*.tif;*.tiff)|*.tif;*.tiff|Bitmap files (*.bmp)|*.bmp|PostScript files (*.ps)|*.ps|VRML files (*.wrl)|*.wrl|All files (*.*)|*.*"),
                    wxFD_SAVE | wxFD_OVERWRITE_PROMPT );
  dlg.SetFilterIndex( m_nScreenshotFilterIndex );
  if ( dlg.ShowModal() == wxID_OK )
  {
    fn = dlg.GetPath();
  }

  if ( !fn.IsEmpty() )
  {
    m_strLastDir = MyUtils::GetNormalizedPath( fn );
    m_nScreenshotFilterIndex = dlg.GetFilterIndex();
    if ( !m_viewRender[nId]->SaveScreenshot( fn, 
                                       m_settingsScreenshot.Magnification,
                                       m_settingsScreenshot.AntiAliasing ) )
    {
      wxMessageDialog dlg( this, _("Error occured writing to file. Please make sure you have right permission and the disk is not full."), _("Error"), wxOK );
      dlg.ShowModal();
    }
  }
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
    wxFileDialog dlg( this, _("Open surface file"), m_strLastDir, _(""),
                      _("Surface files (*.*)|*.*"),
                      wxFD_OPEN );
    if ( dlg.ShowModal() == wxID_OK )
    {
      this->LoadSurfaceFile( dlg.GetPath() );
    }
  }
}

void MainWindow::LoadSurfaceFile( const wxString& filename, const wxString& fn_patch )
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

  WorkerThread* thread = new WorkerThread( this );
  thread->LoadSurface( layer );
}


void MainWindow::LoadSurfaceVector()
{
  wxFileDialog dlg( this, _("Open surface file as vector"), m_strLastDir, _(""),
                    _("Surface files (*.*)|*.*"),
                    wxFD_OPEN );
  if ( dlg.ShowModal() == wxID_OK )
  {
    this->LoadSurfaceVectorFile( dlg.GetPath() );
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
    
    m_strLastDir = MyUtils::GetNormalizedPath( filename );
  }
}

void MainWindow::LoadSurfaceCurvature()
{
  wxFileDialog dlg( this, _("Open curvature file"), m_strLastDir, _(""),
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
    if ( layer->LoadCurvatureFromFile( fn.char_str() ) )
      m_strLastDir = MyUtils::GetNormalizedPath( filename );
  }
}


void MainWindow::LoadSurfaceOverlay()
{
  wxFileDialog dlg( this, _("Open overlay file"), m_strLastDir, _(""),
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
    if ( layer->LoadOverlayFromFile( fn.char_str() ) )
       m_strLastDir = MyUtils::GetNormalizedPath( filename );
  }
}

void MainWindow::LoadSurfaceAnnotation()
{
  wxFileDialog dlg( this, _("Open annotation file"), m_strLastDir, _(""),
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
    if ( layer->LoadAnnotationFromFile( fn.char_str() ) )
      m_strLastDir = MyUtils::GetNormalizedPath( filename );
  }
}

void MainWindow::OnToolRotateVolume( wxCommandEvent& event )
{
  if ( !m_dlgRotateVolume )
    m_dlgRotateVolume = new DialogRotateVolume( this );

  if ( !m_dlgRotateVolume->IsVisible() )
  {
    wxMessageDialog dlg( this, _("Rotation can only apply to volume for now. If you data includes ROI/Surface/Way Points, please do not use this feature yet."), _("Warning"), wxOK );
    dlg.ShowModal();
    m_dlgRotateVolume->Show();
  }
}

void MainWindow::OnToolRotateVolumeUpdateUI( wxUpdateUIEvent& event )
{
// event.Check( m_dlgRotateVolume && m_dlgRotateVolume->IsShown() );
  event.Enable( !GetLayerCollection( "MRI" )->IsEmpty() && !IsProcessing() );
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
// event.Check( m_dlgRotateVolume && m_dlgRotateVolume->IsShown() );
  event.Enable( GetLayerCollection( "MRI" )->GetNumberOfLayers() > 1 && !IsProcessing() );
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
// event.Check( m_dlgRotateVolume && m_dlgRotateVolume->IsShown() );
  event.Enable( !GetLayerCollection( "MRI" )->IsEmpty() && !IsProcessing() );
}


void MainWindow::EnableControls( bool bEnable )
{
  m_controlPanel->Enable( bEnable );
  if ( m_dlgRotateVolume )
    m_dlgRotateVolume->Enable( bEnable );
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
      case LayerPropertiesMRI::Jet:
        if ( scales.size() >= 2 )
        {
          p->SetMinMaxJetScaleWindow( scales[0], scales[1] );
        }
        else if ( !scales.empty() )
          cerr << "Need 2 values for jetscale." << endl;
        break;
      case LayerPropertiesMRI::LUT:
        if ( scales.size() >= 1 )
        {
        }
        else if ( !scales.empty() )
          cerr << "Need a value for lut." << endl;
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
