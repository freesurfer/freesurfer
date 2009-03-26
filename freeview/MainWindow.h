/**
 * @file  MainWindow.h
 * @brief Main window.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2009/03/26 21:13:38 $
 *    $Revision: 1.31 $
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

#ifndef MainWindow_h
#define MainWindow_h

#include <wx/wx.h>
#include <wx/frame.h>
#include <wx/dynarray.h>
#include "Listener.h"
#include "Broadcaster.h"
#include "LayerCollectionManager.h"
#include "CommonDataStruct.h"
#include <vector>

#define ID_WORKER_THREAD wxID_HIGHEST + 1

class ControlPanel;
class wxPanel;
class wxSplitterWindow;
class wxFileHistory;
class wxSpinEvent;
class wxMouseEvent;
class RenderView;
class RenderView2D;
class RenderView3D;
class PixelInfoPanel;
class WindowQuickReference;
class StatusBar;
class LayerCollectionManager;
class LUTDataHolder;
class wxToolBar;
class BrushProperty;
class ToolWindowEdit;
class DialogRotateVolume;
class LayerMRI;
class WindowHistogram;
class WindowOverlayConfiguration;

class MainWindow : public wxFrame, public Listener, public Broadcaster
{
public:
  MainWindow();
  virtual ~MainWindow();

  enum ViewLayout { VL_1X1 = 0, VL_2X2, VL_1N3, VL_1N3_H };
  enum MainView  { MV_Sagittal = 0, MV_Coronal, MV_Axial, MV_3D };

  // event handlers (these functions should _not_ be virtual)
  void OnFileOpen    ( wxCommandEvent& event );
  void OnFileNew    ( wxCommandEvent& event );
  void OnFileNewUpdateUI  ( wxUpdateUIEvent& event );
  void OnFileSave    ( wxCommandEvent& event );
  void OnFileSaveUpdateUI  ( wxUpdateUIEvent& event );
  void OnFileSaveAs   ( wxCommandEvent& event );
  void OnFileSaveAsUpdateUI ( wxUpdateUIEvent& event );
  void OnFileExit    ( wxCommandEvent& event );
  void OnFileRecent   ( wxCommandEvent& event );
  void OnFileNewROI   ( wxCommandEvent& event );
  void OnFileNewROIUpdateUI ( wxUpdateUIEvent& event );
  void OnFileLoadROI   ( wxCommandEvent& event );
  void OnFileLoadROIUpdateUI ( wxUpdateUIEvent& event );
  void OnFileSaveROI   ( wxCommandEvent& event );
  void OnFileSaveROIUpdateUI ( wxUpdateUIEvent& event );
  void OnFileSaveROIAs  ( wxCommandEvent& event );
  void OnFileSaveROIAsUpdateUI( wxUpdateUIEvent& event );
  void OnFileLoadDTI   ( wxCommandEvent& event );
  void OnFileSaveScreenshot ( wxCommandEvent& event );
  void OnFileSaveScreenshotUpdateUI( wxUpdateUIEvent& event );
  void OnFileLoadSurface  ( wxCommandEvent& event );

  void OnFileNewWayPoints   ( wxCommandEvent& event );
  void OnFileNewWayPointsUpdateUI ( wxUpdateUIEvent& event );
  void OnFileLoadWayPoints  ( wxCommandEvent& event );
  void OnFileLoadWayPointsUpdateUI ( wxUpdateUIEvent& event );
  void OnFileSaveWayPoints   ( wxCommandEvent& event );
  void OnFileSaveWayPointsUpdateUI ( wxUpdateUIEvent& event );
  void OnFileSaveWayPointsAs   ( wxCommandEvent& event );
  void OnFileSaveWayPointsAsUpdateUI ( wxUpdateUIEvent& event );

  void OnViewLayout1X1  ( wxCommandEvent& event );
  void OnViewLayout1X1UpdateUI( wxUpdateUIEvent& event );
  void OnViewLayout2X2  ( wxCommandEvent& event );
  void OnViewLayout2X2UpdateUI( wxUpdateUIEvent& event );
  void OnViewLayout1N3  ( wxCommandEvent& event );
  void OnViewLayout1N3UpdateUI( wxUpdateUIEvent& event );
  void OnViewLayout1N3_H  ( wxCommandEvent& event );
  void OnViewLayout1N3_HUpdateUI( wxUpdateUIEvent& event );
  void OnViewSagittal   ( wxCommandEvent& event );
  void OnViewSagittalUpdateUI ( wxUpdateUIEvent& event );
  void OnViewCoronal   ( wxCommandEvent& event );
  void OnViewCoronalUpdateUI ( wxUpdateUIEvent& event );
  void OnViewAxial   ( wxCommandEvent& event );
  void OnViewAxialUpdateUI ( wxUpdateUIEvent& event );
  void OnView3D    ( wxCommandEvent& event );
  void OnView3DUpdateUI  ( wxUpdateUIEvent& event );
  void OnViewReset   ( wxCommandEvent& event );
  void OnViewResetUpdateUI ( wxUpdateUIEvent& event );

  void OnViewControlPanel   ( wxCommandEvent& event );
  void OnViewControlPanelUpdateUI ( wxUpdateUIEvent& event );

  void OnViewScalarBar   ( wxCommandEvent& event );
  void OnViewScalarBarUpdateUI ( wxUpdateUIEvent& event );
  void OnViewCoordinate   ( wxCommandEvent& event );
  void OnViewCoordinateUpdateUI ( wxUpdateUIEvent& event );

  void OnViewHistogram   ( wxCommandEvent& event );
  void OnViewHistogramUpdateUI ( wxUpdateUIEvent& event );

  void OnViewCycleLayer     ( wxCommandEvent& event );
  void OnViewCycleLayerUpdateUI   ( wxUpdateUIEvent& event );
  void OnViewToggleVoxelCoordinates  ( wxCommandEvent& event );
  void OnViewToggleVolumeVisibility  ( wxCommandEvent& event );
  void OnViewToggleVolumeVisibilityUpdateUI( wxUpdateUIEvent& event );
  void OnViewToggleROIVisibility   ( wxCommandEvent& event );
  void OnViewToggleROIVisibilityUpdateUI ( wxUpdateUIEvent& event );
  void OnViewToggleSurfaceVisibility  ( wxCommandEvent& event );
  void OnViewToggleSurfaceVisibilityUpdateUI( wxUpdateUIEvent& event );
  void OnViewToggleWayPointsVisibility ( wxCommandEvent& event );
  void OnViewToggleWayPointsVisibilityUpdateUI( wxUpdateUIEvent& event );
  void OnViewToggleCursorVisibility   ( wxCommandEvent& event );
  void OnViewToggleCursorVisibilityUpdateUI ( wxUpdateUIEvent& event );

  void OnViewSurfaceMain   ( wxCommandEvent& event );
  void OnViewSurfaceMainUpdateUI ( wxUpdateUIEvent& event );
  void OnViewSurfaceInflated  ( wxCommandEvent& event );
  void OnViewSurfaceInflatedUpdateUI( wxUpdateUIEvent& event );
  void OnViewSurfaceWhite   ( wxCommandEvent& event );
  void OnViewSurfaceWhiteUpdateUI ( wxUpdateUIEvent& event );
  void OnViewSurfacePial   ( wxCommandEvent& event );
  void OnViewSurfacePialUpdateUI ( wxUpdateUIEvent& event );
  void OnViewSurfaceOriginal  ( wxCommandEvent& event );
  void OnViewSurfaceOriginalUpdateUI( wxUpdateUIEvent& event );

  void OnModeNavigate   ( wxCommandEvent& event );
  void OnModeNavigateUpdateUI ( wxUpdateUIEvent& event);
  void OnModeVoxelEdit  ( wxCommandEvent& event );
  void OnModeVoxelEditUpdateUI( wxUpdateUIEvent& event);
  void OnModeROIEdit   ( wxCommandEvent& event );
  void OnModeROIEditUpdateUI ( wxUpdateUIEvent& event);
  void OnModeWayPointsEdit   ( wxCommandEvent& event );
  void OnModeWayPointsEditUpdateUI ( wxUpdateUIEvent& event);

  void OnEditCopy( wxCommandEvent& event );
  void OnEditCopyUpdateUI( wxUpdateUIEvent& event );
  void OnEditPaste( wxCommandEvent& event );
  void OnEditPasteUpdateUI( wxUpdateUIEvent& event );
  void OnEditUndo( wxCommandEvent& event );
  void OnEditUndoUpdateUI( wxUpdateUIEvent& event );
  void OnEditRedo( wxCommandEvent& event );
  void OnEditRedoUpdateUI( wxUpdateUIEvent& event );
  void OnEditPreferences( wxCommandEvent& event );

  void OnHelpQuickReference( wxCommandEvent& event );
  void OnHelpAbout( wxCommandEvent& event );

  void OnToolRotateVolume( wxCommandEvent& event );
  void OnToolRotateVolumeUpdateUI( wxUpdateUIEvent& event );
  void OnToolCreateOptimalVolume( wxCommandEvent& event );
  void OnToolCreateOptimalVolumeUpdateUI( wxUpdateUIEvent& event );

  void OnMouseEnterWindow( wxMouseEvent& event );

  void OnWorkerThreadResponse( wxCommandEvent& event );

  void OnSpinBrushSize( wxSpinEvent& event );
  void OnSpinBrushTolerance( wxSpinEvent& event );
  void OnCheckBrushTemplate( wxCommandEvent& event );
  void OnChoiceBrushTemplate( wxCommandEvent& event );

  void OnActivate ( wxActivateEvent& event );
  void OnIconize ( wxIconizeEvent& event );
  void OnClose ( wxCloseEvent &event );
  void OnKeyDown ( wxKeyEvent &event );

  void LoadVolume();
  void NewVolume();
  void SaveVolume();
  void SaveVolumeAs();

  void RotateVolume( std::vector<RotationElement>& rotations );

  void LoadROI();
  void NewROI();
  void SaveROI();
  void SaveROIAs();

  void LoadWayPoints();
  void NewWayPoints();
  void SaveWayPoints();
  void SaveWayPointsAs();

  void LoadSurface();
  void LoadSurfaceVector();
  void LoadSurfaceCurvature();
  void LoadSurfaceOverlay();

  void LoadVolumeFile( const wxString& fn, 
		       const wxString& reg_fn, 
		       bool bResample = true, 
		       int nColorMap = 0 );
  void LoadDTIFile( const wxString& fn_vector, 
		    const wxString& fn_fa, 
		    const wxString& reg_fn, 
		    bool Resample = true );
  void LoadROIFile( const wxString& fn );
  void LoadSurfaceFile( const wxString& fn, 
			const wxString& fn_vector = _("") );
  void LoadSurfaceVectorFile( const wxString& fn );
  void LoadWayPointsFile( const wxString& fn );
  void LoadSurfaceCurvatureFile( const wxString& fn );
  void LoadSurfaceOverlayFile( const wxString& fn );

// bool IsSaving()
//  { return m_bSaving; }

// bool IsLoading()
//  { return m_bLoading; }

  bool IsProcessing()
  {
    return m_bProcessing;
  }

  LayerCollection* GetLayerCollection( std::string strType );
  LayerCollectionManager* GetLayerCollectionManager();

  LUTDataHolder* GetLUTData()
  {
    return m_luts;
  }

  static MainWindow* GetMainWindowPointer();

  void NeedRedraw( int nCount = 2 );

  void AddScript( const wxArrayString& script );
  void RunScript();

  int GetActiveViewId();

  RenderView* GetActiveView();

  RenderView* GetPreviousActiveView();

  RenderView* GetRenderView( int nIndex )
  {
    return m_viewRender[nIndex];
  }

  BrushProperty* GetBrushProperty()
  {
    return m_propertyBrush;
  }

  Settings2D Get2DSettings()
  {
    return m_settings2D;
  }

  SettingsScreenshot GetScreenshotSettings()
  {
    return m_settingsScreenshot;
  }

  void SetAction( int nAction );

  void SetMode( int nMode );

  void EnableControls( bool bEnable );

  wxString GetLastDir()
  {
    return m_strLastDir;
  }
  
  void ConfigureOverlay();

protected:
  void CommandLoadVolume( const wxArrayString& cmd );
  void CommandLoadDTI( const wxArrayString& cmd );
  void CommandLoadROI( const wxArrayString& cmd );
  void CommandLoadSurface( const wxArrayString& cmd );
  void CommandLoadSurfaceVector( const wxArrayString& cmd );
  void CommandLoadSurfaceCurvature( const wxArrayString& cmd );
  void CommandLoadSurfaceOverlay( const wxArrayString& cmd );
  void CommandLoadWayPoints( const wxArrayString& cmd );
  void CommandScreenCapture( const wxArrayString& cmd );
  void CommandSetViewport( const wxArrayString& cmd );
  void CommandZoom( const wxArrayString& cmd );
  void CommandSetRAS( const wxArrayString& cmd );

  void OnInternalIdle();
  virtual void DoListenToMessage ( std::string const iMsg, void* iData, void* sender );

private:
  void SetViewLayout( int nLayout );
  void SetMainView( int nView );
  void ShowControlPanel( bool bShow );
  void UpdateToolbars();
  void DoUpdateToolbars();

  ControlPanel*     m_controlPanel;
  PixelInfoPanel*   m_pixelInfoPanel;
  wxSplitterWindow* m_splitterMain;
  wxSplitterWindow* m_splitterSub;
  wxPanel*          m_renderViewHolder;
  WindowQuickReference* m_wndQuickReference;
  StatusBar*    m_statusBar;
  wxToolBar*    m_toolbarVoxelEdit;
  wxToolBar*    m_toolbarROIEdit;
  wxToolBar*    m_toolbarBrush;
  wxPanel*      m_panelToolbarHolder;
  ToolWindowEdit*     m_toolWindowEdit;
  DialogRotateVolume* m_dlgRotateVolume;
  WindowHistogram*    m_wndHistogram;
  WindowOverlayConfiguration* m_wndOverlayConfiguration;

  RenderView2D*  m_viewAxial;
  RenderView2D*  m_viewSagittal;
  RenderView2D*  m_viewCoronal;
  RenderView3D*  m_view3D;
  RenderView*   m_viewRender[4];
  int     m_nViewLayout;
  int     m_nMainView;
  int     m_nPrevActiveViewId;

  LayerCollectionManager* m_layerCollectionManager;
  LayerMRI*   m_layerVolumeRef;

  Settings2D   m_settings2D;
  SettingsScreenshot m_settingsScreenshot;

  wxString   m_strLastDir;
  wxFileHistory*  m_fileHistory;
  int     m_nMaxRecentFiles;
  bool    m_bResampleToRAS;
  int     m_nScreenshotFilterIndex;

  LUTDataHolder*  m_luts;

  int    m_nRedrawCount;
  bool   m_bToUpdateToolbars;

// bool   m_bSaving;
// bool   m_bLoading;
  bool   m_bProcessing;

  std::vector<wxArrayString> m_scripts;

  BrushProperty* m_propertyBrush;

  DECLARE_CLASS( MainWindow )
  // any class wishing to process wxWindows events must use this macro
  DECLARE_EVENT_TABLE()

};

#endif


