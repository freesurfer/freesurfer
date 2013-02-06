/**
 * @file  MainWindow.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2013/02/06 18:35:43 $
 *    $Revision: 1.125 $
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

  LayerCollection* GetCurrentLayerCollection();
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
  void LoadSurfaceOverlayFile( const QString& filename, const QString& reg_file = "", bool bCorrelation = false );
  void LoadSurfaceAnnotationFile( const QString& filename );
  void LoadSurfaceLabelFile( const QString& filename );
  void LoadSurfaceVectorFile( const QString& filename );
  void LoadSurfaceSplineFile( const QString& filename );

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

  void AddScript(const QString& command);

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

Q_SIGNALS:
  void MainViewChanged( int n );
  void ViewLayoutChanged( int n );
  void SlicePositionChanged();
  void SurfaceRepositionVertexChanged();
  void SurfaceRepositionIntensityChanged();

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
  void SaveVolumeAs();
  void SetSaveCopy(bool bSaveCopy)
  {
    m_settings["SaveCopy"] = bSaveCopy;
  }
  void SyncZoom(bool bSync);
  void SetUseCommandControl(bool b);
  void SetUnifiedTitleAndToolBar(bool b);
  void ShowAllLayers();
  void HideAllLayers();
  bool ParseCommand(const QString& cmd, bool bAutoQuit = false);

  void SetProgress(int n);

  void SaveSurface();
  void SaveSurfaceAs();

  void ToggleSplinePicking();
  void SetSplinePicking(bool b);

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
                        const QString& strGotoLabelName = "" );

  void LoadDTIFile( const QString& fn_vector,
                    const QString& fn_fa,
                    const QString& reg_fn = "",
                    bool Resample = true );
  void LoadVolumeTrackFile( const QString& fn,
                    bool Resample = false );

  void LoadSurfaceFile( const QString& filename,
                        const QString& fn_patch = "",
                        const QString& fn_target = "",
                        bool bAllSurfaces = false );
  void LoadPVolumeFiles( const QStringList& filenames, const QString& prefix, const QString& lut );
  void LoadROIFile( const QString& fn, const QString& ref_vol, const QColor& color = Qt::yellow );
  void LoadWayPointsFile        ( const QString& fn );
  void LoadControlPointsFile    ( const QString& fn );
  void LoadTrackFile            ( const QString& fn );
  void SetVolumeColorMap( int nColorMap, int nColorMapScale, const QList<double>& scales );
  bool GetCursorRAS( double* ras_out, bool tkReg );

  void RunScript();
  void CommandLoadCommand( const QStringList& sa );
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
  void CommandLoadSurfaceOverlay( const QStringList& cmd );
  void CommandLoadSurfaceAnnotation ( const QStringList& cmd );
  void CommandLoadSurfaceLabel  ( const QStringList& cmd );
  void CommandLoadSurfaceSpline ( const QStringList& cmd );
  void CommandLoadConnectivityData  ( const QStringList& cmd );
  void CommandLoadWayPoints     ( const QStringList& cmd );
  void CommandLoadControlPoints ( const QStringList& cmd );
  void CommandLoadPVolumes      ( const QStringList& cmd );
  void CommandScreenCapture     ( const QStringList& cmd );
  void CommandSetViewport       ( const QStringList& cmd );
  void CommandSetViewSize       ( const QStringList& cmd );
  void CommandZoom              ( const QStringList& cmd );
  void CommandSetRAS            ( const QStringList& cmd );
  void CommandSetSlice          ( const QStringList& cmd );
  void CommandSetColorMap       ( const QStringList& cmd );
  void CommandSetLUT            ( const QStringList& cmd );
  void CommandSetHeadScaleOptions( const QStringList& sa );
  void CommandSetOpacity        ( const QStringList& cmd );
  void CommandSetLabelOutline   ( const QStringList& cmd );
  void CommandSetSurfaceOverlayMethod     ( const QStringList& cmd );
  void CommandSetSurfaceColor   ( const QStringList& cmd );
  void CommandSetSurfaceEdgeColor ( const QStringList& cmd );
  void CommandSetSurfaceEdgeThickness ( const QStringList& cmd );
  void CommandSetSurfaceOffset  ( const QStringList& cmd );
  void CommandSetSurfaceLabelOutline   ( const QStringList& cmd );
  void CommandSetDisplaySurfaceVertex  ( const QStringList& cmd );
  void CommandSetSurfaceVertexColor ( const QStringList& cmd );
  void CommandSetPointSetColor ( const QStringList& cmd );
  void CommandSetPointSetRadius( const QStringList& cmd );
  void CommandSetDisplayVector  ( const QStringList& cmd );
  void CommandSetDisplayTensor  ( const QStringList& cmd );
  void CommandSetDisplayIsoSurface  ( const QStringList& cmd );
  void CommandSetIsoSurfaceColor( const QStringList& cmd );
  void CommandSetIsoSurfaceUpsample ( const QStringList& cmd );
  void CommandLoadIsoSurfaceRegion  ( const QStringList& cmd );
  void CommandLockLayer         ( const QStringList& cmd );
  void CommandShowLayer         ( const QStringList& cmd );
  void CommandSetLayerName      ( const QStringList& cmd );
  void CommandSetCamera         ( const QStringList& cmd );
  void CommandGotoLabel         ( const QStringList& cmd );

protected slots:
  void OnIdle();
  void OnAbout();
  void OnRecentVolumeFile();
  void OnRecentSurfaceFile();
  void OnSetViewLayout( QAction* );
  void OnSetMainView  ( QAction* );
  void OnNewVolume();
  void OnLoadVolume();
  bool OnCloseVolume();
  void OnSaveVolume();
  void OnReloadVolume();
  void OnLoadDTI();
  void OnLoadTrackVolume();
  void OnLoadSurface();
  void OnCloseSurface();
  void OnReloadSurface();
  void OnNewROI();
  void OnLoadROI();
  void OnSaveROI();
  void OnSaveROIAs();
  void OnCloseROI();
  void OnNewPointSet();
  void OnLoadPointSet();
  void OnSavePointSet();
  void OnSavePointSetAs();
  void OnClosePointSet();
  void OnLoadTrack();
  void OnCloseTrack();
  void OnIOError( Layer* layer, int ntype );
  void OnIOFinished( Layer* layer, int ntype );
  void OnCycleLayer();
  void OnReverseCycleLayer();
  void OnSetModeNavigate();
  void OnSetModeMeasure();
  void OnSetModeVoxelEdit();
  void OnSetModeROIEdit();
  void OnSetModePointSetEdit();
  void OnPreferences();
  void OnEditUndo();
  void OnEditRedo();
  void OnTransformVolume();
  void OnCropVolume();
  void OnSaveScreenshot();
  void OnVolumeFilterMean();
  void OnVolumeFilterMedian();
  void OnVolumeFilterConvolve();
  void OnVolumeFilterGradient();
  void OnVolumeFilterSobel();
  void OnResetView();
  void OnSavePoint();
  void OnGoToPoint();
  void OnShowAnnotation(bool bShow);
  void OnShowColorBar(bool bShow);
  void OnCopy();
  void OnCopyStructure();
  void OnPaste();
  void OnToggleShowVolume();
  void OnToggleShowSurface();
  void OnToggleShowROI();
  void OnToggleShowPointSet();
  void OnLoadCommand();
  void OnWriteMovieFrames();
  void OnIncreaseOpacity();
  void OnDecreaseOpacity();
  void OnToggleCursorVisibility(bool bShow);
  void OnRepositionSurface();
  void OnSmoothSurface();
  void OnRemoveIntersectionsFromSurface();
  void OnShowLabelStats();
  void OnSaveIsoSurface();
  void OnPlot();
  void OnLineProfile();

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

  void ReassureGeometry();

  void OnVolumeFilterFinished(VolumeFilter* filter);

  void SlotActivateWindow()
  {
    this->activateWindow();
  }

  void UpdateInfoPanel();

private:
  bool DoParseCommand(bool bAutoQuit);
  void SaveSettings();
  void LoadSettings();
  void SetCurrentFile( const QString &fileName, int type = 0 );
  void LoadPointSetFile( const QString& fn, int type );
  void UpdateRecentFileActions();
  void ToggleShowLayer(const QString& type );
  bool UpdateSurfaceCorrelation(LayerSurface* layer);
  void ShowNonModalMessage(const QString& title, const QString& msg);

  QColor ParseColorInput(const QString& cmd);

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
  QList<QAction*>   m_actionRecentVolumes;
  QList<QAction*>   m_actionRecentSurfaces;

  QStringList       m_scripts;
  bool              m_bScriptRunning;

  bool              m_bSplinePicking;

  MyCmdLineParser*  m_cmdParser;

  ToolWindowEdit*       m_toolWindowEdit;
  ToolWindowMeasure*    m_toolWindowMeasure;
  ToolWindowROIEdit*    m_toolWindowROIEdit;
  DialogTransformVolume*    m_dlgTransformVolume;
  DialogCropVolume*     m_dlgCropVolume;
  DialogSaveScreenshot* m_dlgSaveScreenshot;
  DialogWriteMovieFrames*   m_dlgWriteMovieFrames;
  DialogPreferences*    m_dlgPreferences;
  DialogRepositionSurface*  m_dlgRepositionSurface;
  DialogSmoothSurface*  m_dlgSmoothSurface;
  WindowQuickReference* m_wndQuickRef;
  FloatingStatusBar*    m_statusBar;
  TermWidget*           m_term;
  WindowTimeCourse*     m_wndTimeCourse;
  WindowGroupPlot*      m_wndGroupPlot;
  DialogLabelStats*     m_dlgLabelStats;
  DialogLineProfile*    m_dlgLineProfile;

  VolumeFilterWorkerThread* m_threadVolumeFilter;

  SettingsScreenshot    m_settingsScreenshot;
  QVariantMap           m_settings;
  QPoint                m_ptBackUpPos;      // for X11 geometry hack
  QMessageBox*          m_dlgMessage;

  QVariantMap           m_volumeSettings;
  QVariantMap           m_surfaceSettings;
};

#endif // MAINWINDOW_H
