/*
 * Original Author: Ruopeng Wang
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
#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QVariantMap>
#include <QStringList>
#include "CommonDataStruct.h"

class QAction;
class Layer;
class LayerCollection;
class LayerMRI;
class BrushProperty;
class RenderView;
class LUTDataHolder;
class ThreadIOWorker;
class VolumeCropper;
class ToolWindowMeasure;
class ToolWindowEdit;
class ToolWindowROIEdit;
class DialogTransformVolume;
class DialogCropVolume;
class DialogSaveScreenshot;
class DialogPreferences;
class WindowQuickReference;
class FloatingStatusBar;
class TermWidget;
class MyCmdLineParser;
class LayerSurface;
class DialogWriteMovieFrames;
class DialogRepositionSurface;
class DialogSmoothSurface;
class QMessageBox;
class WindowTimeCourse;
class WindowGroupPlot;
class DialogLabelStats;
class VolumeFilterWorkerThread;
class VolumeFilter;
class DialogLineProfile;
class LayerFCD;
class DialogSetCamera;
class DialogThresholdVolume;
class DialogVolumeSegmentation;
class BinaryTreeView;
class WindowLayerInfo;
class QFileSystemWatcher;
class DialogTransformSurface;
class DialogMovePoint;

#define MAX_RECENT_FILES    10

namespace Ui
{
class MainWindow;
}

class MainWindow : public QMainWindow
{
  Q_OBJECT

public:
  explicit MainWindow( QWidget *parent, MyCmdLineParser* cmdParser );
  ~MainWindow();

  enum ViewLayout { VL_1x1 = 0, VL_2x2, VL_1n3, VL_1n3h };
  enum MainView  { MV_Sagittal = 0, MV_Coronal, MV_Axial, MV_3D };

  static MainWindow* GetMainWindow();

  static void WriteLog(const QString& str, const QString& filename = "freeview_log.txt", bool bOverwrite = false);

  bool HadError()
  {
    return m_bHadError;
  }

  void SetHadError(bool b)
  {
    m_bHadError = b;
  }

  BrushProperty* GetBrushProperty()
  {
    return m_propertyBrush;
  }

  bool GetSaveCopy()
  {
    if (m_settings.contains("SaveCopy"))
    {
      return m_settings["SaveCopy"].toBool();
    }
    else
    {
      return true;
    }
  }

  bool IsBusy();
  bool IsEmpty();

  int GetActiveViewId();

  LayerCollection* GetLayerCollection( const QString& strType );
  Layer* GetActiveLayer( const QString& strType );
  Layer* GetTopVisibleLayer( const QString& strType );
  QList<Layer*> GetLayers( const QString& strType );
  QList<Layer*> GetVisibleLayers( const QString& strType );

  bool SetSlicePosition( int nPlane, double dPos, bool bRoundToGrid = true );
  bool SetSlicePosition( double* pos );
  bool OffsetSlicePosition( int nPlane, double dPosDiff, bool bRoundToGrid = true );

  LUTDataHolder* GetLUTData()
  {
    return m_luts;
  }

  VolumeCropper* GetVolumeCropper()
  {
    return m_volumeCropper;
  }

  void LoadSurfaceCurvatureFile( const QString& filename );
  void LoadSurfaceOverlayFile( const QString& filename, const QString& reg_file = "", bool bCorrelation = false, bool bSecondHalfData = false );
  void LoadSurfaceAnnotationFile( const QString& filename );
  bool LoadSurfaceLabelFile( const QString& filename );
  void LoadSurfaceVectorFile( const QString& filename );
  void LoadSurfaceSplineFile( const QString& filename );
  bool LoadSurfaceRGBMap(const QString& filename = "");

  void SetDefaultSampleMethod( int nMethod )
  {
    m_nDefaultSampleMethod = nMethod;
  }

  void SetDefaultColorMapType(const QString& colormap)
  {
    m_strDefaultColorMapType = colormap;
  }

  void SetDefaultConform( bool bConform )
  {
    m_bDefaultConform = bConform;
  }

  RenderView* GetRenderView( int n )
  {
    return m_views[n];
  }

  void SetAction( int nAction );

  int GetMainViewId()
  {
    return m_nMainView;
  }

  RenderView* GetMainView()
  {
    return m_views[m_nMainView];
  }

  void RotateVolume( std::vector<RotationElement>& rotations, bool bAllVolumes );

  void TransformVolume(double* mat, int sample_method);

  void AddScript(const QStringList& command);
  void AddScripts(const QList<QStringList>& cmds);

  QString AutoSelectLastDir( const QString& subdir );
  static QString AutoSelectLastDir( const QString& lastdir, const QString& subdir );

  SettingsScreenshot GetScreenShotSettings()
  {
    return m_settingsScreenshot;
  }

  void SetScreenShotSettings(SettingsScreenshot s)
  {
    m_settingsScreenshot = s;
  }

  bool ParseCommand(int argc, char* argv[], bool bAutoQuit = false);

  TermWidget* GetCommandConsole()
  {
    return m_term;
  }

  Layer* GetSupplementLayer(const QString& type);

  bool IsRepositioningSurface();

  bool GetSplinePicking()
  {
    return m_bSplinePicking;
  }

  QVariantMap GetGeneralSettings()
  {
    return m_settings;
  }

  QVariant GetSetting(const QString& key);

  void SetSetting(const QString& key, const QVariant& value );

  QList<Layer*> GetSelectedLayers(const QString& layerType);

  int GetMode();

  QVariantMap GetDefaultSettings()
  {
    return m_defaultSettings;
  }

  QString GetCurrentLayerType();

  void SaveLayers(const QList<Layer*>& layers);

  Layer* FindSupplementLayer(const QString& name);

Q_SIGNALS:
  void MainViewChanged( int n );
  void ViewLayoutChanged( int n );
  void SlicePositionChanged(bool bCenterView = false);
  void SurfaceRepositionVertexChanged();
  void SurfaceRepositionIntensityChanged();
  void NewVolumeCreated();
  void CycleOverlayRequested();
  void CycleAnnotationRequested();
  void SupplementLayerChanged();
  void OverlayMaskRequested(const QString& fn);
  void RefreshLookUpTableRequested();
  void LinkVolumeRequested(LayerMRI* mri);

public slots:
  void SetMode( int nMode );
  void SetViewLayout( int n );
  void SetMainView  ( int n );
  void LoadSurfaceCurvature();
  void LoadSurfaceOverlay(bool bCorrelation = false);
  void LoadSurfaceAnnotation();
  void LoadSurfaceLabel();
  void LoadSurfaceVector();
  void LoadSurfaceSpline();
  void LoadLUT();
  void RequestRedraw();
  bool SaveVolumeAs();
  void SaveVolumeAsAndReload();
  void SetSaveCopy(bool bSaveCopy)
  {
    m_settings["SaveCopy"] = bSaveCopy;
  }
  void SyncZoom(bool bSync);
  void SetUseCommandControl(bool b);
  void SetUnifiedTitleAndToolBar(bool b);
  //  void ShowAllLayers();
  //  void HideAllLayers();
  bool ParseCommand(const QString& cmd, bool bAutoQuit = false);
  bool ParseCommand(MyCmdLineParser* parser, const QString& cmd, bool bAutoQuit = false);

  void SetProgress(int n);

  void SaveSurface();
  void SaveSurfaceAs();

  void ToggleSplinePicking();
  void SetSplinePicking(bool b);

  void UpdateSettings();

  void GoToContralateralPoint();

  void UpdateSurfaceContralateralInfo();

  void OnApplyVolumeTransform();

  void CenterAtWorldPosition(double* pos, bool mainview_only = false);

  void OnLoadSurfaceParameterization();

  void OnStereoRender(bool bOn);

  void AbortScripts()
  {
    ClearScripts();
  }

  void OnExportLabelStats();

  bool ExportLineProfileThickness(const QString& filename, const QVariantMap& options);

  void OnShowControlPanel(bool bShow);

  void OnFloatPanels(bool bFloat);

  void OnPointSetToLabel();

  DialogMovePoint* GetMovePointDlg()
  {
    return m_dlgMovePoint;
  }

protected:
  void closeEvent   ( QCloseEvent * event );
  void resizeEvent  (QResizeEvent * event);
  void moveEvent    (QMoveEvent * event);
  void showEvent    (QShowEvent * event);
  void LoadVolumeFile(  const QString& filename,
                        const QString& reg_filename = "",
                        bool bResample = false,
                        int nSampleMethod = 0,
                        bool bConform = false,
                        int nGotoLabelOrientation = -1,
                        const QString& strGotoLabelName = "",
                        const QVariantMap& sup_data = QVariantMap());

  void LoadDTIFile( const QString& fn_vector,
                    const QString& fn_fa,
                    const QString& fn_scale = "",
                    const QString& reg_fn = "",
                    bool Resample = true );
  void LoadVolumeTrackFile( const QString& fn,
                            bool Resample = false );

  void LoadSurfaceFile( const QString& filename,
                        const QString& fn_patch = "",
                        const QString& fn_target = "",
                        const QStringList& sup_files = QStringList(), const QVariantMap& sup_options = QVariantMap());
  void LoadPVolumeFiles( const QStringList& filenames, const QString& prefix, const QString& lut );
  void LoadROIFile( const QString& fn, const QString& ref_vol, const QVariantMap& args = QVariantMap() );
  void LoadWayPointsFile        ( const QString& fn, const QVariantMap& args = QVariantMap() );
  void LoadControlPointsFile    ( const QString& fn, const QVariantMap& args = QVariantMap() );
  void LoadTrackFile            ( const QString& fn );
  void LoadFCD        ( const QString& subdir, const QString& subject, const QString& suffix = "");
  void LoadSurfaceParameterization(const QString& filename);
  void LoadSurfaceCoordsFromParameterization(const QString& filename);
  void SetVolumeColorMap( int nColorMap, int nColorMapScale, const QList<double>& scales );
  bool GetCursorRAS( double* ras_out, bool tkReg );

  void RunScript();
  void ClearScripts();
  void CommandLoadCommand( const QStringList& sa );
  void CommandLoadSubject( const QStringList& sa );
  void CommandHideLayer( const QStringList& sa);
  void CommandUnloadLayer( const QStringList& sa);
  void CommandLoadVolume( const QStringList& sa );
  void CommandLoadDTI           ( const QStringList& cmd );
  void CommandLoadVolumeTrack   ( const QStringList& cmd );
  void CommandLoadROI           ( const QStringList& cmd );
  void CommandLoadTrack         ( const QStringList& cmd );
  void CommandLoadSurface       ( const QStringList& cmd );
  void CommandLoadSurfaceVector ( const QStringList& cmd );
  void CommandLoadSurfaceCurvature  ( const QStringList& cmd );
  void CommandSetSurfaceCurvatureMap(const QStringList& cmd);
  void CommandLoadSurfaceOverlay( const QStringList& cmd );
  void CommandLoadSurfaceAnnotation ( const QStringList& cmd );
  void CommandLoadSurfaceLabel  ( const QStringList& cmd );
  void CommandLoadSurfaceSpline ( const QStringList& cmd );
  void CommandLoadSurfaceCoordsFromParameterization ( const QStringList& cmd );
  void CommandLoadConnectomeMatrix  ( const QStringList& cmd );
  void CommandLoadODF           ( const QStringList& cmd );
  void CommandLoadFCD           ( const QStringList& cmd );
  void CommandLoadWayPoints     ( const QStringList& cmd );
  void CommandLoadControlPoints ( const QStringList& cmd );
  void CommandLoadPVolumes      ( const QStringList& cmd );
  void CommandScreenCapture     ( const QStringList& cmd );
  void CommandFlyThrough        ( const QStringList& cmd );
  void CommandSetViewport       ( const QStringList& cmd );
  void CommandSetViewSize       ( const QStringList& cmd );
  void CommandZoom              ( const QStringList& cmd );
  void CommandSetRAS            ( const QStringList& cmd );
  void CommandSetSlice          ( const QStringList& cmd );
  void CommandWriteSurfaceIntersection( const QStringList& cmd);
  void CommandSetColorMap       ( const QStringList& cmd );
  void CommandSetLUT            ( const QStringList& cmd );
  void CommandSetHeatScaleOptions( const QStringList& sa );
  void CommandSetHeatScaleOffset ( const QStringList& sa );
  void CommandSetOpacity        ( const QStringList& cmd );
  void CommandSetLabelOutline   ( const QStringList& cmd );
  void CommandSetSelectedLabels (const QStringList& cmd);
  void CommandSetSurfaceOverlayMethod     ( const QStringList& cmd );
  void CommandSetSurfaceOverlayColormap   ( const QStringList& cmd );
  void CommandSetSurfaceOverlayOpacity    ( const QStringList& cmd );
  void CommandSetSurfaceOverlayFrame      ( const QStringList& cmd );
  void CommandSetSurfaceOverlaySmooth     ( const QStringList& cmd );
  void CommandSetSurfaceOverlayMask      ( const QStringList& cmd );
  void CommandSetSurfaceOverlayCustom     ( const QStringList& cmd );
  void CommandSetSurfaceColor   ( const QStringList& cmd );
  void CommandSetSurfaceEdgeColor ( const QStringList& cmd );
  void CommandSetSurfaceEdgeThickness ( const QStringList& cmd );
  void CommandSetSurfaceOpacity ( const QStringList& cmd );
  void CommandSetSurfaceOffset  ( const QStringList& cmd );
  void CommandSetSurfaceLabelOutline   ( const QStringList& cmd );
  void CommandSetSurfaceLabelOpacity   ( const QStringList& cmd );
  void CommandSetSurfaceAnnotationOutline   ( const QStringList& cmd );
  void CommandGoToSurfaceVertex        ( const QStringList& cmd );
  void CommandSetDisplaySurfaceVertex  ( const QStringList& cmd );
  void CommandHideSurfaceIn3D       ( const QStringList &cmd );
  void CommandSetSurfaceVertexColor ( const QStringList& cmd );
  void CommandSetSurfaceLabelColor  ( const QStringList& cmd );
  void CommandSetSurfaceLabelThreshold  ( const QStringList& cmd );
  void CommandHideSurfaceLabel (const QStringList& cmd );
  void CommandSetPointSetColor ( const QStringList& cmd );
  void CommandSetPointSetRadius( const QStringList& cmd );
  void CommandSetPointSetHeatmap( const QStringList& cmd );
  void CommandSetDisplayVector  ( const QStringList& cmd );
  void CommandSetDisplayTensor  ( const QStringList& cmd );
  void CommandSetDisplayIsoSurface  ( const QStringList& cmd );
  void CommandSetIsoSurfaceColor( const QStringList& cmd );
  void CommandSetIsoSurfaceUpsample ( const QStringList& cmd );
  void CommandSetIsoSurfaceSmooth ( const QStringList& cmd );
  void CommandSetExtractAllRegions ( const QStringList& cmd );
  void CommandLoadIsoSurfaceRegion  ( const QStringList& cmd );
  void CommandLockLayer         ( const QStringList& cmd );
  void CommandShowLayer         ( const QStringList& cmd );
  void CommandSetLayerName      ( const QStringList& cmd );
  void CommandSetVolumeMask     ( const QStringList& cmd );
  void CommandSetSmoothed       ( const QStringList& cmd );
  void CommandSetRgb            ( const QStringList& cmd );
  void CommandGoToLabel         ( const QStringList& cmd );
  void CommandSaveLayer         ( const QStringList& cmd );
  void CommandSetTrackColor     ( const QStringList& cmd );
  void CommandSetTrackRender    ( const QStringList& cmd );
  void CommandLoadTractCluster  ( const QStringList& cmd );
  void CommandReorderLayers   ( const QStringList& cmd );
  void CommandUnloadLayers    ( const QStringList& cmd );
  void CommandSetActiveFrame    ( const QStringList& cmd );
  void CommandSetAutoAdjustFrameContrast ( const QStringList& cmd );
  void CommandSetActiveLayer    ( const QStringList& cmd );
  void CommandExportLineProfileThickness  (const QStringList& cmd);
  void CommandSetVolumeTrackFrame   ( const QStringList& cmd );
  void CommandLinkVolume        ( const QStringList& cmd );

public:
  void CommandSetCamera         ( const QStringList& cmd );
  void SetViewSize  (int x, int y);

protected slots:
  void OnIdle();
  void OnAbout();
  void OnRecentVolumeFile();
  void OnRecentSurfaceFile();
  void OnSetViewLayout( QAction* );
  void OnSetMainView  ( QAction* );
  void OnNewVolume();
  void OnLoadVolume();
  bool OnCloseVolume(const QList<Layer*>& layers = QList<Layer*>());
  void OnSaveVolume();
  void OnReloadVolume();
  void OnLoadDTI();
  void OnLoadTrackVolume();
  void OnLoadSurface();
  void OnCloseSurface(const QList<Layer*>& layers = QList<Layer*>());
  void OnReloadSurface();
  void OnLoadPatch();
  void OnSavePatchAs();
  void OnNewROI();
  void OnLoadROI();
  void OnSaveROI();
  void OnSaveROIAs();
  void OnCloseROI(const QList<Layer*>& layers = QList<Layer*>());
  void OnNewPointSet(bool bSilent = false);
  void OnLoadPointSet();
  void OnSavePointSet(bool bForce = false);
  void OnSavePointSetAs();
  void OnClosePointSet(const QList<Layer*>& layers = QList<Layer*>());
  void OnLoadTrack();
  void OnCloseTrack();
  void OnIOError( Layer* layer, int ntype );
  void OnIOFinished( Layer* layer, int ntype );
  void OnFCDLoadFinished(LayerFCD* layer);
  void OnCycleLayer();
  void OnReverseCycleLayer();
  void OnSetModeNavigate();
  void OnSetModeMeasure();
  void OnSetModeVoxelEdit();
  void OnSetModeReconEdit();
  void OnSetModeROIEdit();
  void OnSetModePointSetEdit();
  void OnPreferences();
  void OnEditUndo();
  void OnEditRedo();
  void OnTransformVolume();
  void OnTransformSurface();
  void OnCropVolume();
  void OnThresholdVolume();
  void OnSegmentVolume();
  void OnSaveScreenshot();
  void OnVolumeFilterMean();
  void OnVolumeFilterMedian();
  void OnVolumeFilterConvolve();
  void OnVolumeFilterGradient();
  void OnVolumeFilterSobel();
  void OnVolumeFilterErode();
  void OnVolumeFilterDilate();
  void OnVolumeFilterOpen();
  void OnVolumeFilterClose();
  void OnVolumeFilterThreshold();
  void OnVolumeFilterBoundary();
  void OnResetView();
  void OnSavePoint();
  void OnGoToPoint();
  void OnToolSaveCamera();
  void OnToolLoadCamera(const QString& fn = "");
  void OnShowAnnotation(bool bShow);
  void OnShowColorBar(bool bShow);
  void OnCopy();
  void OnCopyStructure();
  void OnCopyView();
  void OnPaste();
  void OnToggleShowVolume();
  void OnToggleShowSurface();
  void OnToggleShowROI();
  void OnToggleShowPointSet();
  void OnToggleAllSurfaces();
  void OnToggleWm();
  void OnToggleAseg();
  void OnToggleBrainmask();
  void OnLoadCommand();
  void OnWriteMovieFrames();
  void OnIncreaseOpacity();
  void OnDecreaseOpacity();
  void OnToggleCursorVisibility(bool bShow);
  void OnRepositionSurface();
  void OnSmoothSurface();
  void OnRemoveIntersectionsFromSurface();
  void OnShowLabelStats();
  void OnSaveIsoSurface(const QString& fn = "");
  void OnPlot();
  void OnLineProfile();
  void OnCycleSurfaceLabel();
  void OnGoToROI(bool center = false);
  void OnGoToPointSet(bool center = false);
  void OnLoadFCD();
  void OnCloseFCD();
  void OnGoToSurfaceLabel(bool center = true);
  void OnReloadROI();
  void OnReloadPointSet();

  void OnViewSetCamera();

  void OnLoadConnectomeMatrix();
  void OnCloseConnectomeMatrix();

  void OnActiveLayerChanged(Layer*);

  void SetSlicePosition(double x, double y, double z)
  {
    double ras[3] = {x, y, z};
    SetSlicePosition(ras);
  }

  void SetProcessing( bool bProcessing = true )
  {
    m_bProcessing = bProcessing;
  }

  void SetProcessingFinished()
  {
    SetProcessing(false);
  }

  void ReassureGeometry();

  void OnVolumeFilterFinished(VolumeFilter* filter);

  void SlotActivateWindow()
  {
    this->activateWindow();
  }

  void UpdateInfoPanel();

  void OnSurfaceVertexClicked(LayerSurface* surf);

  void On2DCursorClicked();

  void ReorderLayers(const QList<Layer*>& layers);

  void OnLoadSurfaceLabelRequested(const QString& fn);

  void OnLoadTractCluster();

  void OnTractClusterLoaded(const QVariantMap& data);

  void ShowTractClusterMap();

  void OnLoadVolumeTransform();

  void OnUnloadVolumeTransform();

  void SetCurrentTimeCourseFrame(int nFrame);

  void OnViewLayerInfo();

  void UpdateLayerInfo(Layer* layer);

  void OnSyncInstances(bool bChecked);

  void OnSyncFileChanged(const QString& fn);

  void UpdateSyncCoord();

  void OnTileSyncedWindows();

  void OnLoadODF();

  void OnCloseODF();

  void OnSaveLabelAsVolume();

  void OnCreateOptimalVolume();

  void OnDeleteLayer();

private:
  bool DoParseCommand(MyCmdLineParser* parser, bool bAutoQuit);
  void SaveSettings();
  void LoadSettings();
  void SetCurrentFile( const QString &fileName, int type = 0 );
  void LoadPointSetFile( const QString& fn, int type, const QVariantMap& args = QVariantMap() );
  void UpdateRecentFileActions();
  void ToggleShowLayer(const QString& type );
  void ToggleSpecialVolume(const QString& name);
  bool UpdateSurfaceCorrelation(LayerSurface* layer);
  void ShowNonModalMessage(const QString& title, const QString& msg);
  void LoadConnectomeMatrixFile(const QString& fn_cmat, const QString& fn_parcel, const QString& fn_ctab);
  void GoToContralateralPoint(LayerSurface* layer);
  void ConnectMRILayer(LayerMRI* mri);

  QColor ParseColorInput(const QString& cmd);

  void LoadSphereLeftRightIfNeeded(LayerSurface* layer);

  void UpdateSyncIds(bool bAdd = true);

  void TileWindow(int n);

  int m_nViewLayout;
  int m_nMainView;

public:
  Ui::MainWindow *ui;

private:
  RenderView*   m_views[4];
  //  LayerCollectionManager*   m_layerCollectionManager;
  QMap<QString, LayerCollection*> m_layerCollections;
  BrushProperty*    m_propertyBrush;
  bool              m_bResampleToRAS;
  int               m_nDefaultSampleMethod;
  bool              m_bDefaultConform;
  QString           m_strDefaultColorMapType;
  LayerMRI*         m_layerVolumeRef;
  LUTDataHolder*    m_luts;
  ThreadIOWorker*   m_threadIOWorker;
  bool              m_bProcessing;
  VolumeCropper*    m_volumeCropper;

  QString           m_strLastDir;
  QString           m_strLastFsgdDir;
  QList<QAction*>   m_actionRecentVolumes;
  QList<QAction*>   m_actionRecentSurfaces;

  QList<QStringList>  m_scripts;
  bool              m_bScriptRunning;

  bool              m_bSplinePicking;

  MyCmdLineParser*  m_cmdParser;

  ToolWindowEdit*       m_toolWindowEdit;
  ToolWindowMeasure*    m_toolWindowMeasure;
  ToolWindowROIEdit*    m_toolWindowROIEdit;
  DialogTransformVolume*    m_dlgTransformVolume;
  DialogTransformSurface*   m_dlgTransformSurface;
  DialogCropVolume*     m_dlgCropVolume;
  DialogSaveScreenshot* m_dlgSaveScreenshot;
  DialogWriteMovieFrames*   m_dlgWriteMovieFrames;
  DialogPreferences*    m_dlgPreferences;
  DialogRepositionSurface*  m_dlgRepositionSurface;
  DialogSmoothSurface*  m_dlgSmoothSurface;
  FloatingStatusBar*    m_statusBar;
  TermWidget*           m_term;
  WindowGroupPlot*      m_wndGroupPlot;
  DialogLabelStats*     m_dlgLabelStats;
  DialogLineProfile*    m_dlgLineProfile;
  DialogSetCamera*      m_dlgSetCamera;
  DialogThresholdVolume* m_dlgThresholdVolume;
  DialogVolumeSegmentation* m_dlgVolumeSegmentation;
  BinaryTreeView*       m_wndTractCluster;
  WindowQuickReference* m_wndQuickRef;
  WindowTimeCourse*     m_wndTimeCourse;
  WindowLayerInfo*      m_wndLayerInfo;
  QWidget*              m_widgetFloatControlPanel;
  QWidget*              m_widgetFloatInfoPanel;
  DialogMovePoint*      m_dlgMovePoint;

  VolumeFilterWorkerThread* m_threadVolumeFilter;

  SettingsScreenshot    m_settingsScreenshot;
  QVariantMap           m_settings;
  QPoint                m_ptBackUpPos;      // for X11 geometry hack
  QMessageBox*          m_dlgMessage;

  QMap<int, QVariantMap>  m_layerSettings;
  QVariantMap           m_defaultSettings;
  bool                  m_bShowTransformWindow;

  bool                  m_bVerbose;
  bool                  m_bContinue;

  bool                  m_bHadError;
  QString               m_sTitle;

  QFileSystemWatcher*   m_syncFileWatcher;
  QString               m_sSyncFilePath;
};

#endif // MAINWINDOW_H
