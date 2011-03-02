/**
 * @file  MainWindow.h
 * @brief Main window.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:02 $
 *    $Revision: 1.85 $
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

#ifndef MainWindow_h
#define MainWindow_h

#include <wx/wx.h>
#include <wx/frame.h>
#include <wx/dynarray.h>
#include <wx/timer.h>
#include "Listener.h"
#include "Broadcaster.h"
#include "LayerCollectionManager.h"
#include "CommonDataStruct.h"
#include <vector>

#define ID_WORKER_THREAD          wxID_HIGHEST + 1
#define ID_THREAD_BUILD_CONTOUR   wxID_HIGHEST + 2

class ControlPanel;
class wxPanel;
class wxSplitterWindow;
class wxFileHistory;
class wxSpinEvent;
class wxMouseEvent;
class wxToolBar;
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
class ToolWindowMeasure;
class DialogTransformVolume;
class LayerMRI;
class WindowHistogram;
class WindowOverlayConfiguration;
class WindowConnectivityConfiguration;
class DialogGradientVolume;
class DialogSaveScreenshot;
class DialogSavePoint;
class DialogWriteMovieFrames;
class DialogRepositionSurface;
class DialogCropVolume;
class ConnectivityData;
class VolumeCropper;

class MainWindow : public wxFrame, public Listener, public Broadcaster
{
public:
  MainWindow(wxWindow* parent = NULL) : Listener( "MainWindow" ), Broadcaster( "MainWindow" )
  {
    InitWidgetsFromXRC(parent);
  }
  virtual ~MainWindow();
  
private:
  void InitWidgetsFromXRC(wxWindow* parent);

public:
  enum ViewLayout { VL_1X1 = 0, VL_2X2, VL_1N3, VL_1N3_H };
  enum MainView  { MV_Sagittal = 0, MV_Coronal, MV_Axial, MV_3D };

  // event handlers (these functions should _not_ be virtual)
  void OnFileOpen             ( wxCommandEvent& event );
  void OnFileNew              ( wxCommandEvent& event );
  void OnFileNewUpdateUI      ( wxUpdateUIEvent& event );
  void OnFileSave             ( wxCommandEvent& event );
  void OnFileSaveUpdateUI     ( wxUpdateUIEvent& event );
  void OnFileSaveAs           ( wxCommandEvent& event );
  void OnFileSaveAsUpdateUI   ( wxUpdateUIEvent& event );
  void OnFileExit             ( wxCommandEvent& event );
  void OnFileRecent           ( wxCommandEvent& event );
  void OnFileNewROI           ( wxCommandEvent& event );
  void OnFileNewROIUpdateUI   ( wxUpdateUIEvent& event );
  void OnFileLoadROI          ( wxCommandEvent& event );
  void OnFileLoadROIUpdateUI  ( wxUpdateUIEvent& event );
  void OnFileSaveROI          ( wxCommandEvent& event );
  void OnFileSaveROIUpdateUI  ( wxUpdateUIEvent& event );
  void OnFileSaveROIAs        ( wxCommandEvent& event );
  void OnFileSaveROIAsUpdateUI( wxUpdateUIEvent& event );
  void OnFileLoadDTI          ( wxCommandEvent& event );
  void OnFileLoadPVolumes     ( wxCommandEvent& event );
  void OnFileSaveScreenshot   ( wxCommandEvent& event );
  void OnFileSaveScreenshotUpdateUI   ( wxUpdateUIEvent& event );
  void OnFileSaveMovieFrames  ( wxCommandEvent& event );
  void OnFileSaveMovieFramesUpdateUI  ( wxUpdateUIEvent& event );
  void OnFileLoadSurface      ( wxCommandEvent& event );
  void OnFileSaveSurface      ( wxCommandEvent& event );
  void OnFileSaveSurfaceUpdateUI   ( wxUpdateUIEvent& event );
  void OnFileSaveSurfaceAs    ( wxCommandEvent& event );
  void OnFileSaveSurfaceAsUpdateUI   ( wxUpdateUIEvent& event );
  void OnFileNewWayPoints     ( wxCommandEvent& event );
  void OnFileNewWayPointsUpdateUI ( wxUpdateUIEvent& event );
  void OnFileLoadWayPoints    ( wxCommandEvent& event );
  void OnFileLoadWayPointsUpdateUI ( wxUpdateUIEvent& event );
  void OnFileSaveWayPoints    ( wxCommandEvent& event );
  void OnFileSaveWayPointsUpdateUI ( wxUpdateUIEvent& event );
  void OnFileSaveWayPointsAs  ( wxCommandEvent& event );
  void OnFileSaveWayPointsAsUpdateUI ( wxUpdateUIEvent& event );

  void OnViewLayout1X1        ( wxCommandEvent& event );
  void OnViewLayout1X1UpdateUI( wxUpdateUIEvent& event );
  void OnViewLayout2X2        ( wxCommandEvent& event );
  void OnViewLayout2X2UpdateUI( wxUpdateUIEvent& event );
  void OnViewLayout1N3        ( wxCommandEvent& event );
  void OnViewLayout1N3UpdateUI( wxUpdateUIEvent& event );
  void OnViewLayout1N3_H      ( wxCommandEvent& event );
  void OnViewLayout1N3_HUpdateUI( wxUpdateUIEvent& event );
  void OnViewSagittal         ( wxCommandEvent& event );
  void OnViewSagittalUpdateUI ( wxUpdateUIEvent& event );
  void OnViewCoronal          ( wxCommandEvent& event );
  void OnViewCoronalUpdateUI  ( wxUpdateUIEvent& event );
  void OnViewAxial            ( wxCommandEvent& event );
  void OnViewAxialUpdateUI    ( wxUpdateUIEvent& event );
  void OnView3D               ( wxCommandEvent& event );
  void OnView3DUpdateUI       ( wxUpdateUIEvent& event );
  void OnViewReset            ( wxCommandEvent& event );
  void OnViewResetUpdateUI    ( wxUpdateUIEvent& event );
  void OnViewSnapToAxis       ( wxCommandEvent& event );
  void OnViewSnapToAxisUpdateUI ( wxUpdateUIEvent& event );

  void OnViewControlPanel     ( wxCommandEvent& event );
  void OnViewControlPanelUpdateUI ( wxUpdateUIEvent& event );

  void OnViewScalarBar        ( wxCommandEvent& event );
  void OnViewScalarBarUpdateUI( wxUpdateUIEvent& event );
  void OnViewCoordinate       ( wxCommandEvent& event );
  void OnViewCoordinateUpdateUI   ( wxUpdateUIEvent& event );
  void OnViewSliceFrames      ( wxCommandEvent& event );
  void OnViewSliceFramesUpdateUI  ( wxUpdateUIEvent& event );

  void OnViewHistogram        ( wxCommandEvent& event );
  void OnViewHistogramUpdateUI( wxUpdateUIEvent& event );

  void OnViewCycleLayer                     ( wxCommandEvent& event );
  void OnViewCycleLayerUpdateUI             ( wxUpdateUIEvent& event );
  void OnViewReverseCycleLayer              ( wxCommandEvent& event );
  void OnViewToggleVoxelCoordinates         ( wxCommandEvent& event );
  void OnViewToggleVolumeVisibility         ( wxCommandEvent& event );
  void OnViewToggleVolumeVisibilityUpdateUI ( wxUpdateUIEvent& event );
  void OnViewToggleROIVisibility            ( wxCommandEvent& event );
  void OnViewToggleROIVisibilityUpdateUI    ( wxUpdateUIEvent& event );
  void OnViewToggleSurfaceVisibility        ( wxCommandEvent& event );
  void OnViewToggleSurfaceVisibilityUpdateUI( wxUpdateUIEvent& event );
  void OnViewToggleWayPointsVisibility      ( wxCommandEvent& event );
  void OnViewToggleWayPointsVisibilityUpdateUI( wxUpdateUIEvent& event );
  void OnViewToggleCursorVisibility         ( wxCommandEvent& event );
  void OnViewToggleCursorVisibilityUpdateUI ( wxUpdateUIEvent& event );

  void OnViewSurfaceMain            ( wxCommandEvent& event );
  void OnViewSurfaceMainUpdateUI    ( wxUpdateUIEvent& event );
  void OnViewSurfaceInflated        ( wxCommandEvent& event );
  void OnViewSurfaceInflatedUpdateUI( wxUpdateUIEvent& event );
  void OnViewSurfaceWhite           ( wxCommandEvent& event );
  void OnViewSurfaceWhiteUpdateUI   ( wxUpdateUIEvent& event );
  void OnViewSurfacePial            ( wxCommandEvent& event );
  void OnViewSurfacePialUpdateUI    ( wxUpdateUIEvent& event );
  void OnViewSurfaceOriginal        ( wxCommandEvent& event );
  void OnViewSurfaceOriginalUpdateUI( wxUpdateUIEvent& event );
  
  void OnLayerShowAll               ( wxCommandEvent& event );
  void OnLayerHideAll               ( wxCommandEvent& event );
  void OnLayerShowHideAllUpdateUI   ( wxUpdateUIEvent& event );

  void OnModeNavigate               ( wxCommandEvent& event );
  void OnModeNavigateUpdateUI       ( wxUpdateUIEvent& event);
  void OnModeMeasure                ( wxCommandEvent& event );
  void OnModeMeasureUpdateUI        ( wxUpdateUIEvent& event);
  void OnModeVoxelEdit              ( wxCommandEvent& event );
  void OnModeVoxelEditUpdateUI      ( wxUpdateUIEvent& event);
  void OnModeROIEdit                ( wxCommandEvent& event );
  void OnModeROIEditUpdateUI        ( wxUpdateUIEvent& event);
  void OnModeWayPointsEdit          ( wxCommandEvent& event );
  void OnModeWayPointsEditUpdateUI  ( wxUpdateUIEvent& event);

  void OnEditCopy                   ( wxCommandEvent& event );
  void OnEditCopyUpdateUI           ( wxUpdateUIEvent& event );
  void OnEditCopyStructure          ( wxCommandEvent& event );
  void OnEditCopyStructureUpdateUI  ( wxUpdateUIEvent& event );
  void OnEditPaste                  ( wxCommandEvent& event );
  void OnEditPasteUpdateUI          ( wxUpdateUIEvent& event );
  void OnEditUndo                   ( wxCommandEvent& event );
  void OnEditUndoUpdateUI           ( wxUpdateUIEvent& event );
  void OnEditRedo                   ( wxCommandEvent& event );
  void OnEditRedoUpdateUI           ( wxUpdateUIEvent& event );
  void OnEditRename                 ( wxCommandEvent& event );
  void OnEditRenameUpdateUI         ( wxUpdateUIEvent& event );
  void OnEditPreferences            ( wxCommandEvent& event );

  void OnHelpQuickReference         ( wxCommandEvent& event );
  void OnHelpAbout                  ( wxCommandEvent& event );

  void OnToolTransformVolume           ( wxCommandEvent& event );
  void OnToolTransformVolumeUpdateUI   ( wxUpdateUIEvent& event );
  void OnToolCropVolume             ( wxCommandEvent& event );
  void OnToolCropVolumeUpdateUI     ( wxUpdateUIEvent& event );
  void OnToolOptimalVolume          ( wxCommandEvent& event );
  void OnToolOptimalVolumeUpdateUI  ( wxUpdateUIEvent& event );
  void OnToolGradientVolume         ( wxCommandEvent& event );
  void OnToolGradientVolumeUpdateUI ( wxUpdateUIEvent& event );
  void OnToolSaveGotoPoint          ( wxCommandEvent& event );
  void OnToolSaveGotoPointUpdateUI  ( wxUpdateUIEvent& event );
  void OnToolGotoPoint              ( wxCommandEvent& event );
  void OnToolGotoPointUpdateUI      ( wxUpdateUIEvent& event );
  void OnToolMenuGotoPoint          ( wxCommandEvent& event );
  void OnToolRepositionSurface      ( wxCommandEvent& event );
  void OnToolRepositionSurfaceUpdateUI  ( wxUpdateUIEvent& event );
  
  void OnToolLabelStats             ( wxCommandEvent& event );
  void OnToolLabelStatsUpdateUI     ( wxUpdateUIEvent& event );
  
  void OnFilterMean                 ( wxCommandEvent& event );
  void OnFilterMedian               ( wxCommandEvent& event );
  void OnFilterConvolve             ( wxCommandEvent& event );
  void OnFilterGradient             ( wxCommandEvent& event );
  void OnFilterUpdateUI             ( wxUpdateUIEvent& event );

  void OnMouseEnterWindow           ( wxMouseEvent& event );

  void OnWorkerThreadResponse       ( wxCommandEvent& event );
  void OnBuildContourThreadResponse ( wxCommandEvent& event );

  void OnSpinBrushSize              ( wxSpinEvent& event );
  void OnSpinBrushTolerance         ( wxSpinEvent& event );
  void OnCheckBrushTemplate         ( wxCommandEvent& event );
  void OnChoiceBrushTemplate        ( wxCommandEvent& event );

  void OnActivate           ( wxActivateEvent& event );
  void OnIconize            ( wxIconizeEvent& event );
  void OnClose              ( wxCloseEvent &event );
  void OnKeyDown            ( wxKeyEvent &event );
  void OnTimerWriteMovieFrames  ( wxTimerEvent& event );
  
  void LoadVolume();
  void NewVolume();
  void SaveVolume();
  void SaveVolumeAs();
  void SaveRegistrationAs();

  void RotateVolume( std::vector<RotationElement>& rotations, bool bAllVolumes );

  void LoadLUT();
  
  void LoadROI();
  void NewROI();
  void SaveROI();
  void SaveROIAs();

  void LoadWayPoints();
  void NewWayPoints();
  void SaveWayPoints();
  void SaveWayPointsAs();

  void LoadSurface();
  void SaveSurface();
  void SaveSurfaceAs();
  void LoadSurfaceVector();
  void LoadSurfaceCurvature();
  void LoadSurfaceOverlay();
  void LoadSurfaceAnnotation();
  void LoadSurfaceLabel();

  void LoadVolumeFile ( const wxString& fn, 
		       const wxString& reg_fn, 
		       bool bResample = false, 
           int nSampleMethod = 0,
           bool bConform = false );
  void LoadDTIFile    ( const wxString& fn_vector, 
		    const wxString& fn_fa, 
		    const wxString& reg_fn, 
		    bool Resample = true );
  void LoadPVolumeFiles( const wxArrayString& filenames, const wxString& prefix, const wxString& lut );
  void LoadROIFile    ( const wxString& fn, const wxString& ref_vol = "" );
  void LoadSurfaceFile( const wxString& fn, 
			const wxString& fn_patch = _(""),
      const wxString& fn_target = _("")
      );
  void LoadSurfaceVectorFile    ( const wxString& fn );
  void LoadWayPointsFile        ( const wxString& fn );  
  void LoadControlPointsFile    ( const wxString& fn );
  void LoadSurfaceCurvatureFile ( const wxString& fn );
  void LoadSurfaceOverlayFile   ( const wxString& fn );
  void LoadSurfaceAnnotationFile( const wxString& fn );
  void LoadSurfaceLabelFile     ( const wxString& fn );
  void LoadConnectivityDataFile ( const wxString& data_file, const wxString& lut );

// bool IsSaving()
//  { return m_bSaving; }

// bool IsLoading()
//  { return m_bLoading; }
  
  void ShowScalarBar( bool bShow = true );

  bool IsProcessing()
  {
    return m_bProcessing;
  }

  LayerCollection* GetLayerCollection( std::string strType );
  Layer* GetActiveLayer( std::string strType );
  LayerCollectionManager* GetLayerCollectionManager();

  LUTDataHolder* GetLUTData()
  {
    return m_luts;
  }

  static MainWindow* GetMainWindowPointer();

  void NeedRedraw( int nCount = 1 );
  
  void ForceRedraw();

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
  
  bool SaveScreenshot();

  void SetAction( int nAction );

  void SetMode( int nMode );

  void EnableControls( bool bEnable );

  wxString GetLastDir()
  {
    return m_strLastDir;
  }
  
  void ConfigureOverlay();
  
  void SetVolumeColorMap( int nColorMap, int nColorMapScale, std::vector<double>& scales );
  
  ConnectivityData* GetConnectivityData()
  {
    return m_connectivity;
  }
  
  bool GetCursorRAS( double* ras_out, bool bTkreg = false );
  
  void UpdateGotoPoints();
  
  bool GetSaveCopy()
  {
    return m_settingsGeneral.SaveCopy;
  }
  
  void StartWriteMovieFrames();
  
  void StopWriteMovieFrames();
  
  bool IsWritingMovieFrames()
  {
    return m_timerWriteMovieFrames.IsRunning();
  }
  
  ToolWindowMeasure* GetToolWindowMeasure()
  {
    return m_toolWindowMeasure;
  }
  
  VolumeCropper* GetVolumeCropper()
  {
    return m_volumeCropper;
  }
  
  int GetMainViewId()
  {
    return m_nMainView;
  }
  
  void SetDefaultSampleMethod( int nMethod )
  {
    m_nDefaultSampleMethod = nMethod;
  }
  
  void SetDefaultConform( bool bConform )
  {
    m_bDefaultConform = bConform;
  }
  
  wxString AutoSelectLastDir( wxString subdirectory );
  
  static wxString AutoSelectLastDir( wxString lastDir, wxString subdirectory );
  
protected:
  void CommandLoadVolume        ( const wxArrayString& cmd );
  void CommandLoadDTI           ( const wxArrayString& cmd );
  void CommandLoadROI           ( const wxArrayString& cmd );
  void CommandLoadSurface       ( const wxArrayString& cmd );
  void CommandLoadSurfaceVector ( const wxArrayString& cmd );
  void CommandLoadSurfaceCurvature  ( const wxArrayString& cmd );
  void CommandLoadSurfaceOverlay( const wxArrayString& cmd );
  void CommandLoadSurfaceAnnotation ( const wxArrayString& cmd );
  void CommandLoadSurfaceLabel  ( const wxArrayString& cmd );
  void CommandLoadConnectivityData  ( const wxArrayString& cmd );
  void CommandLoadWayPoints     ( const wxArrayString& cmd );
  void CommandLoadControlPoints ( const wxArrayString& cmd );
  void CommandLoadPVolumes      ( const wxArrayString& cmd );
  void CommandScreenCapture     ( const wxArrayString& cmd );
  void CommandSetViewport       ( const wxArrayString& cmd );
  void CommandSetViewSize       ( const wxArrayString& cmd );
  void CommandZoom              ( const wxArrayString& cmd );
  void CommandSetRAS            ( const wxArrayString& cmd );
  void CommandSetSlice          ( const wxArrayString& cmd );
  void CommandSetColorMap       ( const wxArrayString& cmd );
  void CommandSetLUT            ( const wxArrayString& cmd );
  void CommandSetHeadScaleOptions( const wxArrayString& sa );
  void CommandSetOpacity        ( const wxArrayString& cmd ); 
  void CommandSetSurfaceOverlayMethod     ( const wxArrayString& cmd );
  void CommandSetSurfaceColor   ( const wxArrayString& cmd );
  void CommandSetSurfaceEdgeColor ( const wxArrayString& cmd );
  void CommandSetSurfaceEdgeThickness ( const wxArrayString& cmd );
  void CommandSetSurfaceOffset  ( const wxArrayString& cmd );
  void CommandSetWayPointsColor ( const wxArrayString& cmd );
  void CommandSetWayPointsRadius( const wxArrayString& cmd );
  void CommandSetDisplayVector  ( const wxArrayString& cmd );
  void CommandSetDisplayTensor  ( const wxArrayString& cmd );
  void CommandSetDisplayIsoSurface  ( const wxArrayString& cmd );  
  void CommandSetIsoSurfaceColor( const wxArrayString& cmd );
  void CommandLoadIsoSurfaceRegion  ( const wxArrayString& cmd );
  void CommandLockLayer         ( const wxArrayString& cmd );
  void CommandShowLayer         ( const wxArrayString& cmd );
  
  void CommandSetLayerName      ( const wxArrayString& cmd );

  void OnInternalIdle();

  virtual void DoListenToMessage ( std::string const iMsg, void* iData, void* sender );

private:
  void LoadPointSetFile( const wxString& fn, int type );
  void ContinueScripts();
  void SetViewLayout( int nLayout );
  void SetMainView( int nView );
  void ShowControlPanel( bool bShow );
  void UpdateToolbars();
  void DoUpdateToolbars();
  void DoSaveScreenshot(); 
  void BuildGotoPointMenu( wxMenu* menu );
  LayerCollection* GetCurrentLayerCollection();

  ControlPanel*       m_controlPanel;
  PixelInfoPanel*     m_pixelInfoPanel;
  wxSplitterWindow*   m_splitterMain;
  wxSplitterWindow*   m_splitterSub;
  wxPanel*            m_renderViewHolder;
  WindowQuickReference* m_wndQuickReference;
  StatusBar*          m_statusBar;
  wxToolBar*          m_toolbarMain;
  wxPanel*            m_panelToolbarHolder;
  ToolWindowEdit*     m_toolWindowEdit;
  ToolWindowMeasure*  m_toolWindowMeasure;
  DialogTransformVolume* m_dlgTransformVolume;
  WindowHistogram*    m_wndHistogram;
  WindowOverlayConfiguration*       m_wndOverlayConfiguration;
  WindowConnectivityConfiguration*  m_wndConnectivityConfiguration; 
  DialogGradientVolume*       m_dlgGradientVolume;
  DialogSaveScreenshot*       m_dlgSaveScreenshot;
  DialogWriteMovieFrames*     m_dlgWriteMovieFrames;
  DialogSavePoint*            m_dlgSavePoint;
  DialogRepositionSurface*    m_dlgRepositionSurface;
  DialogCropVolume*           m_dlgCropVolume;
  wxMenu*           m_menuGotoPoints;

  RenderView2D*   m_viewAxial;
  RenderView2D*   m_viewSagittal;
  RenderView2D*   m_viewCoronal;
  RenderView3D*   m_view3D;
  RenderView*     m_viewRender[4];
  int             m_nViewLayout;
  int             m_nMainView;
  int             m_nPrevActiveViewId;

  LayerCollectionManager* m_layerCollectionManager;
  LayerMRI*       m_layerVolumeRef;
  
  VolumeCropper*  m_volumeCropper;
  
  ConnectivityData*   m_connectivity;

  SettingsGeneral     m_settingsGeneral;
  Settings2D          m_settings2D;
  SettingsScreenshot  m_settingsScreenshot;
  SettingsMovieFrames m_settingsMovieFrames;

  wxString        m_strLastDir;
  wxFileHistory*  m_fileHistory;
  int             m_nMaxRecentFiles;
  bool            m_bResampleToRAS;
  int             m_nScreenshotFilterIndex;

  LUTDataHolder*  m_luts;
  
  wxArrayString   m_strGotoPoints;

  int             m_nRedrawCount;
  bool            m_bToUpdateToolbars;

  bool            m_bProcessing;
  bool            m_bDoScreenshot;
  
  int             m_nDefaultSampleMethod;
  bool            m_bDefaultConform;

  std::vector<wxArrayString> m_scripts;

  BrushProperty*  m_propertyBrush;

  wxTimer         m_timerWriteMovieFrames;
  
  DECLARE_CLASS( MainWindow )
  // any class wishing to process wxWindows events must use this macro
  DECLARE_EVENT_TABLE()

};

#endif


