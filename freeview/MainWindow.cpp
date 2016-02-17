/**
 * @file  MainWindow.cpp
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: zkaufman $
 *    $Date: 2016/02/17 20:36:46 $
 *    $Revision: 1.313 $
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
#include "MainWindow.h"
#include "ui_MainWindow.h"
#include "LayerMRI.h"
#include "LayerPropertyMRI.h"
#include "LayerSurface.h"
#include "LayerROI.h"
#include "LayerTrack.h"
#include "LayerDTI.h"
#include "LayerVolumeTrack.h"
#include "LayerCollection.h"
#include "BrushProperty.h"
#include <QtGui>
#include <QtCore>
#include <QFileInfo>
#include <QMessageBox>
#include "LUTDataHolder.h"
#include "DialogLoadVolume.h"
#include "ThreadIOWorker.h"
#include "VolumeCropper.h"
#include "DialogPreferences.h"
#include "ToolWindowMeasure.h"
#include "ToolWindowEdit.h"
#include "ToolWindowROIEdit.h"
#include "Interactor2DVoxelEdit.h"
#include "DialogNewVolume.h"
#include "DialogNewROI.h"
#include "DialogNewPointSet.h"
#include "DialogLoadDTI.h"
#include "MyUtils.h"
#include "SurfaceOverlay.h"
#include "SurfaceOverlayProperty.h"
#include "LayerPLabel.h"
#include "LayerPointSet.h"
#include "LayerPropertySurface.h"
#include "LayerPropertyPointSet.h"
#include "LayerPropertyROI.h"
#include "FSPointSet.h"
#include "DialogLoadPointSet.h"
#include "DialogTransformVolume.h"
#include "DialogCropVolume.h"
#include "DialogSaveScreenshot.h"
#include "VolumeFilterConvolve.h"
#include "VolumeFilterMean.h"
#include "VolumeFilterGradient.h"
#include "VolumeFilterMedian.h"
#include "VolumeFilterSobel.h"
#include "VolumeFilterErode.h"
#include "VolumeFilterDilate.h"
#include "VolumeFilterOpen.h"
#include "VolumeFilterClose.h"
#include "VolumeFilterThreshold.h"
#include "DialogVolumeFilter.h"
#include "DialogGradientFilter.h"
#include "Cursor2D.h"
#include "Cursor3D.h"
#include "DialogAbout.h"
#include "WindowQuickReference.h"
#include "FloatingStatusBar.h"
#include "TermWidget.h"
#include "MyCmdLineParser.h"
#include "DialogSavePointSet.h"
#include "DialogSaveVolume.h"
#include "DialogWriteMovieFrames.h"
#include "LayerLandmarks.h"
#include "Interactor2DNavigate.h"
#include "MainApplication.h"
#include "DialogRepositionSurface.h"
#include "WindowTimeCourse.h"
#include "DialogLabelStats.h"
#include "VolumeFilterWorkerThread.h"
#include "FSGroupDescriptor.h"
#include "WindowGroupPlot.h"
#include "DialogLoadSurfaceOverlay.h"
#include "DialogReloadLayer.h"
#include "DialogSmoothSurface.h"
#include "DialogLineProfile.h"
#include "LayerLineProfile.h"
#include "DialogLoadConnectome.h"
#include "LayerConnectomeMatrix.h"
#include "DialogLoadSurface.h"
#include "LayerFCD.h"
#include "LayerPropertyFCD.h"
#include "DialogSetCamera.h"
#include "DialogThresholdVolume.h"
#include "DialogVolumeSegmentation.h"
#include <QProcessEnvironment>
#include "Json.h"
#include "DialogThresholdFilter.h"

MainWindow::MainWindow( QWidget *parent, MyCmdLineParser* cmdParser ) :
  QMainWindow( parent ),
  ui(new Ui::MainWindow),
  m_bResampleToRAS(false),
  m_nDefaultSampleMethod( SAMPLE_NEAREST ),
  m_strDefaultColorMapType( "grayscale" ),
  m_bDefaultConform( false ),
  m_layerVolumeRef( NULL ),
  m_bScriptRunning(false),
  m_bProcessing(false),
  m_bSplinePicking(true),
  m_cmdParser(cmdParser)
{ 
  // must create layer collections first before setupui()
  m_layerCollections["MRI"] = new LayerCollection( "MRI", this );
  m_layerCollections["ROI"] = new LayerCollection( "ROI", this );
  m_layerCollections["Surface"] = new LayerCollection( "Surface", this );
  m_layerCollections["PointSet"] = new LayerCollection( "PointSet", this );
  m_layerCollections["Track"] = new LayerCollection( "Track", this );
  m_layerCollections["CMAT"] = new LayerCollection("CMAT", this);
  m_layerCollections["FCD"] = new LayerCollection("FCD", this);

  // supplemental layers will not show on control panel
  m_layerCollections["Supplement"] = new LayerCollection( "Supplement", this);
  LayerLandmarks* landmarks = new LayerLandmarks(this);
  m_layerCollections["Supplement"]->AddLayer(landmarks);  

  m_luts = new LUTDataHolder();
  m_propertyBrush = new BrushProperty();
  m_volumeCropper = new VolumeCropper( this );
  connect(m_volumeCropper, SIGNAL(CropBoundChanged(LayerMRI*)), this, SLOT(RequestRedraw()));
  connect(m_layerCollections["MRI"], SIGNAL(LayerRemoved(Layer*)),
          m_propertyBrush, SLOT(OnLayerRemoved(Layer*)));

  ui->setupUi(this);
#ifndef DEVELOPMENT
//  ui->tabWidgetControlPanel->removeTab(ui->tabWidgetControlPanel->indexOf(ui->tabTrack));
#endif

  this->addAction(ui->actionIncreaseOpacity);
  this->addAction(ui->actionDecreaseOpacity);
  this->addAction(ui->actionCycleSurfaceLabel);

  m_statusBar = new FloatingStatusBar(this);
  m_statusBar->hide();

  ui->viewSagittal->SetViewPlane( 0 );
  ui->viewCoronal ->SetViewPlane( 1 );
  ui->viewAxial   ->SetViewPlane( 2 );
  ui->viewSagittal->SetFocusFrameColor( 1, 0, 0 );
  ui->viewCoronal ->SetFocusFrameColor( 0, 1, 0 );
  ui->viewAxial   ->SetFocusFrameColor( 0, 0, 1 );
  m_views[0] = ui->viewSagittal;
  m_views[1] = ui->viewCoronal;
  m_views[2] = ui->viewAxial;
  m_views[3] = ui->view3D;
  m_nMainView = MV_Sagittal;
  m_nViewLayout = VL_2x2;

  m_toolWindowMeasure = new ToolWindowMeasure( this );
  m_toolWindowMeasure->hide();
  m_toolWindowEdit = new ToolWindowEdit( this );
  m_toolWindowEdit->hide();

  m_toolWindowROIEdit = new ToolWindowROIEdit( this );
  m_toolWindowROIEdit->hide();
  m_dlgTransformVolume = new DialogTransformVolume( this );
  m_dlgTransformVolume->hide();
  for (int i = 0; i < 3; i++)
  {
    connect(m_dlgTransformVolume, SIGNAL(CurrentLandmarkChanged(int)),
            ((RenderView2D*)m_views[i])->GetInteractorNavigate(),
            SLOT(SetCurrentLandmark(int)));
  }

  m_dlgCropVolume = new DialogCropVolume(this);
  m_dlgCropVolume->hide();
  connect(m_volumeCropper, SIGNAL(CropBoundChanged(LayerMRI*)),
          m_dlgCropVolume, SLOT(OnCropBoundChanged(LayerMRI*)));
  connect(m_layerCollections["MRI"], SIGNAL(LayerRemoved(Layer*)),
          m_dlgCropVolume, SLOT(OnLayerRemoved(Layer*)));
  m_dlgSaveScreenshot = NULL;
  m_dlgPreferences = NULL;

  m_dlgThresholdVolume = new DialogThresholdVolume(this);
  m_dlgThresholdVolume->hide();
  connect(m_layerCollections["MRI"], SIGNAL(LayerAdded(Layer*)), m_dlgThresholdVolume, SLOT(UpdateVolumes()));
  connect(m_layerCollections["MRI"], SIGNAL(LayerRemoved(Layer*)), m_dlgThresholdVolume, SLOT(UpdateVolumes()));

  m_dlgVolumeSegmentation = new DialogVolumeSegmentation(this);
  m_dlgVolumeSegmentation->hide();
  connect(m_layerCollections["MRI"], SIGNAL(LayerAdded(Layer*)), m_dlgVolumeSegmentation, SLOT(UpdateVolumes()));
  connect(m_layerCollections["MRI"], SIGNAL(LayerRemoved(Layer*)), m_dlgVolumeSegmentation, SLOT(UpdateVolumes()));

  m_wndQuickRef = new WindowQuickReference(this);
  m_wndQuickRef->hide();

  m_dlgWriteMovieFrames = new DialogWriteMovieFrames(this);
  m_dlgWriteMovieFrames->hide();
  connect(m_dlgWriteMovieFrames, SIGNAL(Started()),
          m_statusBar, SLOT(ShowProgress()));
  connect(m_dlgWriteMovieFrames, SIGNAL(Stopped()),
          m_statusBar, SLOT(HideProgress()));
  connect(m_dlgWriteMovieFrames, SIGNAL(Progress(int)),
          m_statusBar, SLOT(SetProgress(int)));

  m_term = new TermWidget(this);
  m_term->hide();
  connect(ui->actionShowCommandConsole, SIGNAL(toggled(bool)),
          m_term, SLOT(setVisible(bool)));
  connect(ui->actionRunCommand, SIGNAL(triggered()),
          m_term, SLOT(show()));

  m_dlgMessage = new QMessageBox(QMessageBox::Warning, "Warning", "", QMessageBox::Ok,
                                 this, Qt::Tool | Qt::MSWindowsFixedSizeDialogHint);
  m_dlgMessage->setModal(false);
  m_dlgMessage->hide();

  m_dlgRepositionSurface = new DialogRepositionSurface(this);
  m_dlgRepositionSurface->hide();
  connect(m_layerCollections["Surface"], SIGNAL(LayerModified()),
          m_dlgRepositionSurface, SLOT(UpdateUI()));
//  connect(ui->view3D, SIGNAL(SurfaceVertexClicked()),
//          m_dlgRepositionSurface, SLOT(OnSurfaceVertexClicked()));
  connect(this, SIGNAL(SurfaceRepositionVertexChanged()),
          m_dlgRepositionSurface, SLOT(UpdateVertex()), Qt::QueuedConnection);
  connect(this, SIGNAL(SurfaceRepositionIntensityChanged()),
          m_dlgRepositionSurface, SLOT(UpdateIntensity()), Qt::QueuedConnection);

  m_dlgSmoothSurface = new DialogSmoothSurface(this);
  m_dlgSmoothSurface->hide();

  m_wndTimeCourse = new WindowTimeCourse(this);
  m_wndTimeCourse->hide();
  connect(this, SIGNAL(SlicePositionChanged()), m_wndTimeCourse, SLOT(UpdateData()));
  connect(ui->actionTimeCourse, SIGNAL(toggled(bool)),
          m_wndTimeCourse, SLOT(setVisible(bool)));

  m_wndGroupPlot = new WindowGroupPlot(this);
  m_wndGroupPlot->hide();

  m_dlgLabelStats = new DialogLabelStats(this);
  m_dlgLabelStats->hide();
  connect(this, SIGNAL(SlicePositionChanged()), m_dlgLabelStats, SLOT(UpdateStats()), Qt::QueuedConnection);
  connect(m_layerCollections["MRI"], SIGNAL(ActiveLayerChanged(Layer*)), m_dlgLabelStats, SLOT(UpdateStats()), Qt::QueuedConnection);
  connect(m_layerCollections["ROI"], SIGNAL(ActiveLayerChanged(Layer*)), m_dlgLabelStats, SLOT(UpdateStats()), Qt::QueuedConnection);

  m_dlgLineProfile = new DialogLineProfile(this);
  m_dlgLineProfile->hide();
  connect(m_layerCollections["PointSet"], SIGNAL(LayerAdded(Layer*)), m_dlgLineProfile, SLOT(UpdatePointSetList()));
  connect(m_layerCollections["PointSet"], SIGNAL(LayerRemoved(Layer*)), m_dlgLineProfile, SLOT(UpdatePointSetList()));
  connect(m_layerCollections["PointSet"], SIGNAL(LayerNameChanged()), m_dlgLineProfile, SLOT(UpdatePointSetList()));
  for (int i = 0; i < 3; i++)
  {
    connect(this->m_views[i], SIGNAL(LineProfileIdPicked(LayerLineProfile*, int)),
            m_dlgLineProfile, SLOT(OnLineProfileIdPicked(LayerLineProfile*,int)));
  }

  m_dlgSetCamera = new DialogSetCamera(this);
  m_dlgSetCamera->hide();

  QStringList keys = m_layerCollections.keys();
  for ( int i = 0; i < keys.size(); i++ )
  {
    for ( int j = 0; j < 4; j++ )
    {
      connect( m_layerCollections[keys[i]], SIGNAL(LayerAdded(Layer*)),
               m_views[j], SLOT(RefreshAllActors()) );
      connect( m_layerCollections[keys[i]], SIGNAL(LayerRemoved(Layer*)),
               m_views[j], SLOT(RefreshAllActors()) );
      connect( m_layerCollections[keys[i]], SIGNAL(LayerMoved(Layer*)),
               m_views[j], SLOT(RefreshAllActors()) );
      connect( m_layerCollections[keys[i]], SIGNAL(LayerActorChanged()),
               m_views[j], SLOT(RefreshAllActors()) );
      connect( m_layerCollections[keys[i]], SIGNAL(LayerActorUpdated()),
               m_views[j], SLOT(RequestRedraw()) );

      // 2D view only
      if ( j < 3 )
      {
        connect( m_layerCollections[keys[i]], SIGNAL(LayerAdded(Layer*)),
                 m_views[j], SLOT(UpdateAnnotation()) );
        connect( m_layerCollections[keys[i]], SIGNAL(LayerRemoved(Layer*)),
                 m_views[j], SLOT(UpdateAnnotation()) );
        connect( m_layerCollections[keys[i]], SIGNAL(LayerAdded(Layer*)),
                 m_views[j], SLOT(Update2DOverlay()) );
        connect( m_layerCollections[keys[i]], SIGNAL(LayerRemoved(Layer*)),
                 m_views[j], SLOT(Update2DOverlay()) );
        connect( m_layerCollections[keys[i]], SIGNAL(LayerMoved(Layer*)),
                 m_views[j], SLOT(Update2DOverlay()) );
        connect( m_layerCollections[keys[i]], SIGNAL(LayerVisibilityChanged()),
                 m_views[j], SLOT(Update2DOverlay()) );
      }
      // 3D view only
      else
      {
        connect( m_layerCollections[keys[i]], SIGNAL(LayerAdded(Layer*)),
                 m_views[j], SLOT(UpdateBounds()) );
        connect( m_layerCollections[keys[i]], SIGNAL(LayerRemoved(Layer*)),
                 m_views[j], SLOT(UpdateBounds()) );
      }
    }
    connect(m_layerCollections[keys[i]], SIGNAL(LayerAdded(Layer*)),
            ui->treeWidgetCursorInfo, SLOT(UpdateAll()), Qt::QueuedConnection);
    connect(m_layerCollections[keys[i]], SIGNAL(LayerRemoved(Layer*)),
            ui->treeWidgetCursorInfo, SLOT(UpdateAll()), Qt::QueuedConnection);
    connect(m_layerCollections[keys[i]], SIGNAL(LayerMoved(Layer*)),
            ui->treeWidgetCursorInfo, SLOT(UpdateAll()), Qt::QueuedConnection);
    connect(m_layerCollections[keys[i]], SIGNAL(LayerShowInfoChanged()),
            ui->treeWidgetCursorInfo, SLOT(UpdateAll()), Qt::QueuedConnection);

    connect(m_layerCollections[keys[i]], SIGNAL(LayerAdded(Layer*)),
            ui->treeWidgetMouseInfo, SLOT(UpdateAll()), Qt::QueuedConnection);
    connect(m_layerCollections[keys[i]], SIGNAL(LayerRemoved(Layer*)),
            ui->treeWidgetMouseInfo, SLOT(UpdateAll()), Qt::QueuedConnection);
    connect(m_layerCollections[keys[i]], SIGNAL(LayerMoved(Layer*)),
            ui->treeWidgetMouseInfo, SLOT(UpdateAll()), Qt::QueuedConnection);
    connect(m_layerCollections[keys[i]], SIGNAL(LayerShowInfoChanged()),
            ui->treeWidgetMouseInfo, SLOT(UpdateAll()), Qt::QueuedConnection);

    connect(m_layerCollections[keys[i]], SIGNAL(ActiveLayerChanged(Layer*)),
            this, SLOT(OnActiveLayerChanged(Layer*)), Qt::QueuedConnection);

    /*
    connect(m_layerCollections[keys[i]], SIGNAL(ActiveLayerChanged(Layer*)),
            ui->tabAllLayers, SLOT(OnActiveLayerChanged(Layer*)), Qt::QueuedConnection);
    connect(m_layerCollections[keys[i]], SIGNAL(LayerAdded(Layer*)),
            ui->tabAllLayers, SLOT(OnLayerAdded(Layer*)), Qt::QueuedConnection);
    connect(m_layerCollections[keys[i]], SIGNAL(LayerRemoved(Layer*)),
            ui->tabAllLayers, SLOT(OnLayerRemoved(Layer*)), Qt::QueuedConnection);
            */
  }
  for ( int i = 0; i < 4; i++ )
  {
    connect( this, SIGNAL(SlicePositionChanged()), m_views[i], SLOT(OnSlicePositionChanged()) );
  }

  for (int i = 0; i < 3; i++)
  {
    connect(m_layerCollections["MRI"], SIGNAL(ActiveLayerChanged(Layer*)),
            m_views[i], SLOT(UpdateAnnotation()));
  }
  connect(m_layerCollections["MRI"], SIGNAL(LayerTransformed()),
          m_views[3], SLOT(UpdateBounds()));

  QActionGroup* actionGroupMode = new QActionGroup( this );
  actionGroupMode->addAction( ui->actionNavigate );
  actionGroupMode->addAction( ui->actionMeasure );
  actionGroupMode->addAction( ui->actionVoxelEdit );
  actionGroupMode->addAction( ui->actionReconEdit );
  actionGroupMode->addAction( ui->actionROIEdit );
  actionGroupMode->addAction( ui->actionPointSetEdit );
  actionGroupMode->setExclusive( true );

  QActionGroup* actionGroupView = new QActionGroup( this );
  actionGroupView->addAction( ui->actionViewAxial );
  actionGroupView->addAction( ui->actionViewSagittal );
  actionGroupView->addAction( ui->actionViewCoronal );
  actionGroupView->addAction( ui->actionView3D );
  actionGroupView->setExclusive( true );

  QActionGroup* actionGroupLayout = new QActionGroup( this );
  actionGroupLayout->addAction( ui->actionLayout1x1 );
  actionGroupLayout->addAction( ui->actionLayout2x2 );
  actionGroupLayout->addAction( ui->actionLayout1n3 );
  actionGroupLayout->addAction( ui->actionLayout1n3h );
  actionGroupLayout->setExclusive( true );

  for ( int i = 0; i < MAX_RECENT_FILES; i++ )
  {
    QAction* act = new QAction( this );
    act->setVisible(false);
    connect( act, SIGNAL(triggered()),
             this, SLOT(OnRecentVolumeFile()) );
    m_actionRecentVolumes << act;
    act = new QAction( this );
    act->setVisible(false);
    connect( act, SIGNAL(triggered()),
             this, SLOT(OnRecentSurfaceFile()) );
    m_actionRecentSurfaces << act;
  }
  ui->menuRecentFiles->insertActions( ui->actionSurfaces, m_actionRecentVolumes );
  ui->menuRecentFiles->insertActions( 0, m_actionRecentSurfaces );
  UpdateRecentFileActions();

  connect( actionGroupView, SIGNAL( triggered(QAction*) ),
           this, SLOT(OnSetMainView(QAction*)), Qt::QueuedConnection );
  connect( actionGroupLayout, SIGNAL( triggered(QAction*) ),
           this, SLOT(OnSetViewLayout(QAction*)), Qt::QueuedConnection );

  connect( ui->actionQuickReference, SIGNAL(triggered()), this->m_wndQuickRef, SLOT(show()));
  connect( ui->actionShowSlices, SIGNAL(toggled(bool)),
           this->ui->view3D, SLOT(SetShowAllSlices(bool)));
  connect( ui->view3D, SIGNAL(VolumeTrackMouseOver(Layer*,QVariantMap)),
           ui->treeWidgetMouseInfo, SLOT(UpdateTrackVolumeAnnotation(Layer*,QVariantMap)));

  m_threadIOWorker = new ThreadIOWorker( this );
  connect( m_threadIOWorker, SIGNAL(Error(Layer*, int)), this, SLOT(OnIOError(Layer*, int)), Qt::QueuedConnection );
  connect( m_threadIOWorker, SIGNAL(Finished(Layer*, int )), this, SLOT(OnIOFinished(Layer*, int)), Qt::QueuedConnection );
  connect( m_threadIOWorker, SIGNAL(FCDLoadFinished(LayerFCD*)), this, SLOT(OnFCDLoadFinished(LayerFCD*)), Qt::QueuedConnection);
  connect( m_threadIOWorker, SIGNAL(started()), this, SLOT(SetProcessing()));

  connect( m_threadIOWorker, SIGNAL(Progress(int)), m_statusBar, SLOT(SetProgress(int)) );
  connect( m_threadIOWorker, SIGNAL(started()), m_statusBar, SLOT(ShowProgress()));
  connect( m_threadIOWorker, SIGNAL(finished()), m_statusBar, SLOT(HideProgress()));

  connect(this, SIGNAL(SlicePositionChanged()),
          ui->treeWidgetCursorInfo, SLOT(OnCursorPositionChanged()), Qt::QueuedConnection);
  connect(ui->treeWidgetCursorInfo, SIGNAL(RASChangeTriggered(double,double,double)),
          this, SLOT(SetSlicePosition(double,double,double)));
  connect(m_layerCollections["MRI"], SIGNAL(MouseRASPositionChanged()),
          ui->treeWidgetMouseInfo, SLOT(OnMousePositionChanged()), Qt::QueuedConnection);
  connect(ui->treeWidgetMouseInfo, SIGNAL(RASChangeTriggered(double,double,double)),
          m_layerCollections["MRI"], SLOT(SetMouseRASPosition(double,double,double)));

  m_threadVolumeFilter = new VolumeFilterWorkerThread(this);
  connect(m_threadVolumeFilter, SIGNAL(Finished(VolumeFilter*)),
          this, SLOT(OnVolumeFilterFinished(VolumeFilter*)));
  connect( m_threadVolumeFilter, SIGNAL(Progress(int)), m_statusBar, SLOT(SetProgress(int)));
  connect( m_threadVolumeFilter, SIGNAL(started()), m_statusBar, SLOT(ShowProgress()));
  connect( m_threadVolumeFilter, SIGNAL(finished()), m_statusBar, SLOT(HideProgress()));

  this->LoadSettings();

  // timer served as idle loop to update action status
  QTimer* timer = new QTimer( this );
  connect( timer, SIGNAL( timeout() ), this, SLOT( OnIdle() ), Qt::QueuedConnection );
  timer->start( 500 );

  ui->menuFile->setTearOffEnabled(true);
  ui->menuEdit->setTearOffEnabled(true);
  ui->menuView->setTearOffEnabled(true);
  ui->menuLayer->setTearOffEnabled(true);
  ui->menuAction->setTearOffEnabled(true);
  ui->menuTools->setTearOffEnabled(true);
  ui->menuHelp->setTearOffEnabled(true);

  addAction(ui->actionToggleAseg);
  addAction(ui->actionToggleBrainmask);
  addAction(ui->actionToggleWm);
}

MainWindow::~MainWindow()
{
  delete m_propertyBrush;
  delete m_luts;
}

MainWindow* MainWindow::GetMainWindow()
{
  foreach ( QWidget* w, QApplication::topLevelWidgets() )
  {
    if ( w->objectName() == "MainWindow" )
    {
      return (MainWindow*)w;
    }
  }
  return NULL;
}

void MainWindow::SetProgress(int n)
{
  m_statusBar->SetProgress(n);
}

void MainWindow::LoadSettings()
{
  QSettings settings;
  restoreGeometry ( settings.value("MainWindow/Geometry").toByteArray() );
  restoreState    ( settings.value("MainWindow/WindowState").toByteArray() );
  ui->splitterMain->restoreState( settings.value("MainWindow/SplitterState").toByteArray() );
  ui->splitterInfoPanel->restoreState( settings.value("InfoPanel/SplitterState").toByteArray() );
  SetViewLayout( settings.value( "MainWindow/ViewLayout", VL_2x2 ).toInt() );
  SetMainView( settings.value( "MainWindow/MainView", MV_Sagittal ).toInt() );
  m_strLastDir = settings.value( "MainWindow/LastDir").toString();
  m_strLastFsgdDir = settings.value( "MainWindow/LastFsgdDir").toString();
  m_settingsScreenshot.Magnification = settings.value("ScreenShot/Magnification", 1).toInt();
  m_settingsScreenshot.AntiAliasing = settings.value("ScreenShot/AntiAliasing", false).toBool();
  m_settingsScreenshot.HideCoords = settings.value("ScreenShot/HideAnnotation", true).toBool();
  m_settingsScreenshot.HideCursor = settings.value("ScreenShot/HideCursor", true).toBool();
  m_settings = settings.value("Settings/General").toMap();
  if (!m_settings.contains("SaveCopy"))
  {
    m_settings["SaveCopy"] = true;
  }
  if (!m_settings.contains("BackgroundColor"))
  {
    m_settings["BackgroundColor"] = Qt::black;
  }
  if (!m_settings.contains("CursorColor"))
  {
    m_settings["CursorColor"] = Qt::red;
  }
  if (!m_settings.contains("CursorStyle"))
  {
    m_settings["CursorStyle"] = 0;
  }
  if (!m_settings.contains("AnnotationColor"))
  {
    m_settings["AnnotationColor"] = Qt::white;
  }
  if (!m_settings.contains("SyncZoom"))
  {
    m_settings["SyncZoom"] = true;
  }
  if (!m_settings.contains("MacUseCommand"))
  {
    m_settings["MacUseCommand"] = false;
  }
  if (!m_settings.contains("MacUnifiedTitleBar"))
  {
    m_settings["MacUnifiedTitleBar"] = false;
  }
  if (!m_settings.contains("DarkConsole"))
  {
    m_settings["DarkConsole"] = true;
  }

  for (int i = 0; i < 4; i++)
  {
    m_views[i]->SetBackgroundColor(m_settings["BackgroundColor"].value<QColor>());
    if ( i < 3 )
    {
      ((RenderView2D*)m_views[i])->GetCursor2D()->SetColor(m_settings["CursorColor"].value<QColor>());
      ((RenderView2D*)m_views[i])->GetCursor2D()->SetStyle(m_settings["CursorStyle"].toInt());
    }
    else
    {
      ((RenderView3D*)m_views[i])->GetCursor3D()->SetColor(m_settings["CursorColor"].value<QColor>());
    }
  }
  SyncZoom(m_settings["SyncZoom"].toBool());
  m_term->SetDarkTheme(m_settings["DarkConsole"].toBool());

#ifdef Q_WS_MAC
  this->SetUnifiedTitleAndToolBar(m_settings["MacUnifiedTitleBar"].toBool());
  this->SetUseCommandControl(m_settings["MacUseCommand"].toBool());
#endif
}

void MainWindow::SaveSettings()
{
  QSettings settings;
  settings.setValue( "MainWindow/Geometry",       saveGeometry() );
  settings.setValue( "MainWindow/WindowState",    saveState() );
  settings.setValue( "MainWindow/SplitterState",  ui->splitterMain->saveState() );
  settings.setValue( "InfoPanel/SplitterState",   ui->splitterInfoPanel->saveState() );
  settings.setValue( "MainWindow/MainView",       this->m_nMainView );
  settings.setValue( "MainWindow/ViewLayout",     this->m_nViewLayout);
  settings.setValue( "MainWindow/LastDir",        m_strLastDir );
  settings.setValue( "MainWindow/LastFsgdDir",        m_strLastFsgdDir );
  if (m_dlgSaveScreenshot)
  {
    SettingsScreenshot s = m_dlgSaveScreenshot->GetSettings();
    settings.setValue("ScreenShot/Magnification", s.Magnification);
    settings.setValue("ScreenShot/AntiAliasing", s.AntiAliasing);
    settings.setValue("ScreenShot/HideAnnotation", s.HideCoords);
    settings.setValue("ScreenShot/HideCursor", s.HideCursor);
  }
  if (m_dlgPreferences)
  {
    settings.setValue("Settings/General", m_dlgPreferences->GetSettings());
  }
  /*
  QStringList tabs;
  for (int i = 0; i < ui->tabWidgetControlPanel->count(); i++)
  {
    tabs << ui->tabWidgetControlPanel->widget(i)->objectName();
  }
  settings.setValue("ControlPanel/TabOrder", tabs);
  */
}

void MainWindow::closeEvent( QCloseEvent * event )
{
  if (this->IsBusy())
  {
    if (QMessageBox::question(this, "Warning",
                              "There is on-going data processing. If you force quit, any data that is being saved can be lost. Do you really want to force quit?",
                              QMessageBox::Yes, QMessageBox::No ) == QMessageBox::No)
    {
      event->ignore();
      return;
    }
  }

  QList<LayerEditable*> layers;
  QStringList keys = m_layerCollections.keys();
  for ( int i = 0; i < keys.size(); i++ )
  {
    for ( int j = 0; j < m_layerCollections[keys[i]]->GetNumberOfLayers(); j++ )
    {
      LayerEditable* layer = qobject_cast<LayerEditable*>(m_layerCollections[keys[i]]->GetLayer(j));
      if ( layer && layer->IsModified() )
      {
        layers << layer;
      }
    }
  }
  if ( !layers.isEmpty() )
  {
    QString msg = "The following layers have been modified but not saved:\n\n";
    for (int i = 0; i < layers.size(); i++)
    {
      msg += layers[i]->GetName() + " (" + layers[i]->GetFileName() + ")\n";
    }
    msg += "\nDo you still want to quit?";
    QMessageBox msgbox(this);
    msgbox.setIcon(QMessageBox::Question);
    QAbstractButton* yesBtn = msgbox.addButton("Quit", QMessageBox::YesRole);
    msgbox.addButton("Cancel", QMessageBox::NoRole);
    msgbox.setText(msg);
    msgbox.setWindowTitle("Warning");
    msgbox.exec();
    if (msgbox.clickedButton() != yesBtn)
    {
      event->ignore();
      return;
    }
  }

  SaveSettings();
  QMainWindow::closeEvent( event );
}

void MainWindow::ReassureGeometry()
{
  // hack to reassure X11 geometry
  if (pos() == QPoint(0,0))
  {
    move(m_ptBackUpPos);
  }
}

void MainWindow::showEvent(QShowEvent *event)
{
  QMainWindow::showEvent(event);
#ifdef Q_WS_X11
  m_ptBackUpPos = this->pos();
  QTimer::singleShot(500, this, SLOT(ReassureGeometry()));
#endif
}

void MainWindow::resizeEvent(QResizeEvent *event)
{
  if (m_statusBar->isVisible())
  {
    m_statusBar->Reposition();
  }

  QMainWindow::resizeEvent(event);
}

void MainWindow::moveEvent(QMoveEvent *event)
{
  if (m_statusBar->isVisible())
  {
    m_statusBar->Reposition();
  }

  QMainWindow::moveEvent(event);
}

bool MainWindow::ParseCommand(int argc, char *argv[], bool bAutoQuit)
{
  return (m_cmdParser->Parse(argc, argv) &&
          DoParseCommand(m_cmdParser, bAutoQuit));
}

bool MainWindow::ParseCommand(const QString& cmd, bool bAutoQuit)
{
  QString strg = cmd;
  strg.replace("~", QDir::homePath());
  return (m_cmdParser->Parse(strg) &&
          DoParseCommand(m_cmdParser, bAutoQuit));
}

bool MainWindow::ParseCommand(MyCmdLineParser* parser, const QString& cmd, bool bAutoQuit)
{
  QString strg = cmd;
  strg.replace("~", QDir::homePath());
  return (parser->Parse(strg) &&
          DoParseCommand(parser, bAutoQuit));
}


bool MainWindow::DoParseCommand(MyCmdLineParser* parser, bool bAutoQuit)
{
  QStringList sa;
  QStringList floatingArgs;
  string_array tmp_ar = parser->GetFloatingArguments();
  for ( size_t i = 0; i < tmp_ar.size(); i++)
  {
    floatingArgs << tmp_ar[i].c_str();
  }

  m_bShowTransformWindow = parser->Found( "transform-volume" );

  bool bReverseOrder = parser->Found("rorder");
  if ( parser->Found("cmd", &sa))
  {
    this->AddScript( QStringList("loadcommand") << sa[0]);
  }
  if ( parser->Found("hide", &sa))
  {
    this->AddScript( QStringList("hidelayer") << sa[0]);
  }
  if ( parser->Found("unload", &sa))
  {
    this->AddScript( QStringList("unloadlayer") << sa[0]);
  }
  if ( parser->Found( "trilinear" ) )
  {
    this->SetDefaultSampleMethod( SAMPLE_TRILINEAR );
  }
  if ( parser->Found( "cubic" ) )
  {
    this->SetDefaultSampleMethod( SAMPLE_CUBIC_BSPLINE );
  }
  if ( parser->Found( "conform" ) )
  {
    this->SetDefaultConform( true );
  }
  if (parser->Found( "smoothed" ) )
  {
    m_defaultSettings["Smoothed"] = true;
  }
//  if ( parser->Found( "colormap", &sa ))
//  {
//    this->SetDefaultColorMapType(sa[0]);
//  }
  if ( parser->Found( "viewport", &sa ) )
  {
    QString strg = sa[0].toLower();
    if ( strg == "sagittal" || strg == "sag" || strg == "x" )
    {
      SetMainView( MV_Sagittal );
    }
    else if ( strg == "coronal" ||  strg == "cor" || strg == "y" )
    {
      SetMainView( MV_Coronal );
    }
    else if ( strg == "axial" || strg == "z" )
    {
      SetMainView( MV_Axial );
    }
    else if ( strg == "3d" )
    {
      SetMainView( MV_3D );
    }
    else
    {
      std::cerr << "Unrecognized viewport name '" << qPrintable(sa[0]) << "'.\n";
    //  return false;
    }
  }
  if (parser->Found("timecourse"))
  {
    ui->actionTimeCourse->setChecked(true);
  }

  if (parser->Found("nocursor"))
  {
    OnToggleCursorVisibility(false);
  }

  if ( parser->Found( "viewsize", &sa ) )
  {
    this->AddScript( QStringList("setviewsize") << sa[0] << sa[1] );
  }

  int nRepeats = parser->GetNumberOfRepeats( "recon" );
  bool bHasVolume = false;
  for ( int n = 0; n < nRepeats; n++ )
  {
    parser->Found( "recon", &sa, n );
    for ( int i = 0; i < sa.size(); i++ )
    {
      QStringList script = QStringList("loadsubject") << sa[i];
    //  AddScript(script);
    //  qDebug() << script;
      CommandLoadSubject(script);
      bHasVolume = true;
    }
  }

  QList<QStringList> cmds;
  if ( floatingArgs.size() > 0 )
  {
    for ( int i = 0; i < floatingArgs.size(); i++ )
    {
      QStringList script = QStringList("loadvolume") << floatingArgs[i];
      if ( parser->Found( "r" ) )
      {
        script << "r";
      }
      cmds << script;
      bHasVolume = true;
    }
  }

  nRepeats = parser->GetNumberOfRepeats( "v" );
  for ( int n = 0; n < nRepeats; n++ )
  {
    parser->Found( "v", &sa, n );
    for ( int i = 0; i < sa.size(); i++ )
    {
      QStringList script = QStringList("loadvolume") << sa[i];
      if ( parser->Found( "r" ) )
      {
        script << "r";
      }
      cmds << script;
      bHasVolume = true;
    }
  }

  if (bReverseOrder)
  {
    QList<QStringList> tempList;
    for (int i = cmds.size()-1; i >= 0; i--)
      tempList << cmds[i];
    cmds = tempList;
  }
  for (int i = 0; i < cmds.size(); i++)
  {
    for (int j = 0; j < cmds[i].size(); j++)
    {
      if (cmds[i][j].contains(":basis=1"))
      {
        QStringList cmd = cmds[i];
        cmd[j].replace(":basis=1", QString(":basis=%1").arg(i));
        cmds.removeAt(i);
        cmds.insert(0, cmd);
      }
    }
  }
  AddScripts(cmds);

  nRepeats = parser->GetNumberOfRepeats( "dti" );
  for ( int n = 0; n < nRepeats; n++ )
  {
    parser->Found( "dti", &sa, n );
    for ( int i = 0; i < sa.size()/2; i++ )
    {
      QStringList script("loaddti");
      script << sa[i*2] << sa[i*2+1];
      if ( parser->Found( "r" ) )
      {
        script << "r";
      }
      this->AddScript( script );
      bHasVolume = true;
    }
  }

  nRepeats = parser->GetNumberOfRepeats( "tv" );
  for ( int n = 0; n < nRepeats; n++ )
  {
    parser->Found( "tv", &sa, n );
    for ( int i = 0; i < sa.size(); i++ )
    {
      QStringList script = QStringList("loadtrackvolume") << sa[i];
      if ( parser->Found( "r" ) )
      {
        script << "r";
      }
      this->AddScript( script );
    }
  }

  cmds.clear();
  nRepeats = parser->GetNumberOfRepeats( "l" );
  if (nRepeats > 0 && !bHasVolume)
  {
    QString msg = "Can not load volume label without loading a volume first";
    ShowNonModalMessage("Warning", msg);
    std::cerr << qPrintable(msg) << std::endl;
  }
  else
  {
    for ( int n = 0; n < nRepeats; n++ )
    {
      parser->Found( "l", &sa, n );
      for (int i = 0; i < sa.size(); i++ )
      {
        cmds << (QStringList("loadroi") << sa[i]);
      }
    }
    if (bReverseOrder)
    {
      QList<QStringList> tempList;
      for (int i = cmds.size()-1; i >= 0; i--)
        tempList << cmds[i];
      cmds = tempList;
    }
    AddScripts(cmds);
  }

  if (parser->Found("fcd", &sa))
  {
    QStringList script = QStringList("loadfcd") << sa[0] << sa[1];
    this->AddScript( script );
    m_defaultSettings["Smoothed"] = true;
  }

  cmds.clear();
  nRepeats = parser->GetNumberOfRepeats( "f" );
  for ( int n = 0; n < nRepeats; n++ )
  {
    parser->Found( "f", &sa, n );
    for (int i = 0; i < sa.size(); i++ )
    {
      QStringList script("loadsurface");
      script << sa[i];
      if ( parser->Found( "r" ) )
      {
        script << "r";
      }
      cmds << script;
    }
  }
  if (bReverseOrder)
  {
    QList<QStringList> tempList;
    for (int i = cmds.size()-1; i >= 0; i--)
      tempList << cmds[i];
    cmds = tempList;
  }
  AddScripts(cmds);

  nRepeats = parser->GetNumberOfRepeats( "t" );
  for ( int n = 0; n < nRepeats; n++ )
  {
    parser->Found( "t", &sa, n );
    for (int i = 0; i < sa.size(); i++ )
    {
      this->AddScript( QStringList("loadtrack") << sa[i] );
    }
  }

  cmds.clear();
  nRepeats = parser->GetNumberOfRepeats( "w" );
  for ( int n = 0; n < nRepeats; n++ )
  {
    parser->Found( "w", &sa, n );
    for ( int i = 0; i < sa.size(); i++ )
    {
      QStringList script("loadwaypoints");
      script << sa[i];
      if ( parser->Found( "r" ) )
      {
        script << "r";
      }
      cmds << script;
    }
  }
  if (bReverseOrder)
  {
    QList<QStringList> tempList;
    for (int i = cmds.size()-1; i >= 0; i--)
      tempList << cmds[i];
    cmds = tempList;
  }
  AddScripts(cmds);

  nRepeats = parser->GetNumberOfRepeats( "c" );
  for ( int n = 0; n < nRepeats; n++ )
  {
    parser->Found( "c", &sa, n );
    for ( int i = 0; i < sa.size(); i++ )
    {
      QStringList script("loadcontrolpoints");
      script << sa[i];
      if ( parser->Found( "r" ) )
      {
        script << "r";
      }
      this->AddScript( script );
    }
  }

  nRepeats = parser->GetNumberOfRepeats( "p-labels" );
  for ( int n = 0; n < nRepeats; n++ )
  {
    parser->Found( "p-labels", &sa, n );
    QString filenames = sa.join(";");
    QStringList script("loadpvolumes");
    script << filenames;
    sa.clear();
    if ( parser->Found( "p-prefix", &sa ) )
    {
      script << sa[0];
    }
    else
    {
      script << "n/a";
    }
    sa.clear();
    if ( parser->Found( "p-lut", &sa ) )
    {
      script << sa[0];
    }
    this->AddScript( script );
  }

  if ( parser->Found( "cmat", &sa ) )
  {
    this->AddScript( QStringList("loadconnectome") << sa[0] << sa[1] );
  }

  if ( parser->Found( "ras", &sa ) )
  {
    bool bOK;
    double ras[3];
    ras[0] = sa[0].toDouble(&bOK);
    ras[1] = sa[1].toDouble(&bOK);
    ras[2] = sa[2].toDouble(&bOK);
    if ( !bOK )
    {
      std::cerr << "Invalid argument for 'ras'. Arguments must be valid float values.\n";
      return false;
    }
    this->AddScript( QStringList("ras") << sa[0] << sa[1] << sa[2] );
  }

  if ( parser->Found( "slice", &sa ) )
  {
    bool bOK;
    int slice[3];
    slice[0] = sa[0].toInt(&bOK);
    slice[1] = sa[1].toInt(&bOK);
    slice[2] = sa[2].toInt(&bOK);
    if ( !bOK )
    {
      std::cerr << "Invalid argument for 'slice'. Arguments must be valid integers.\n";
      return false;
    }

    this->AddScript( QStringList("slice") << sa[0] << sa[1] << sa[2] );
  }

  if ( parser->Found( "zoom", &sa ) )
  {
    bool bOK;
    double dValue = sa[0].toDouble(&bOK);
    if ( !bOK || dValue == 0 )
    {
      std::cerr << "Invalid argument for 'zoom'. Argument must be a valid float value.\n";
      return false;
    }
    this->AddScript( QStringList("zoom") << sa[0] );
  }

  if ( parser->Found( "camera", &sa ) )
  {
    if ( sa.size()%2 > 0)
    {
      std::cerr << "Invalid arguments for 'cam'. Arguments must be in pairs.\n";
      return false;
    }
    this->AddScript( QStringList("setcamera") << sa );
  }

  if (parser->Found("colorscale"))
  {
    this->AddScript(QStringList("showcolorscale"));
  }

  if (parser->Found("cc"))
  {
    AddScript(QStringList("center"));
  }

  if ( parser->Found( "ss", &sa ) )
  {
    QString mag_factor = "1";
    if (sa.size() > 1)
      mag_factor = sa[1];
    this->AddScript( QStringList("screencapture") << sa[0] << mag_factor );
    if (bAutoQuit && !parser->Found("noquit"))
    {
      this->AddScript( QStringList("quit") );
    }
  }

  if (parser->Found("fly", &sa))
  {

  }

  if ( parser->Found("quit"))
    AddScript(QStringList("quit") );

  return true;
}

int MainWindow::GetActiveViewId()
{
  for (int i = 0; i < 4; i++)
  {
    if (m_views[i]->hasFocus())
    {
      return i;
    }
  }
  return m_nMainView;
}

void MainWindow::SyncZoom(bool bSync)
{
  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      if ( i != j)
      {
        disconnect(m_views[i], SIGNAL(Zooming(RenderView2D*)), m_views[j], SLOT(SyncZoomTo(RenderView2D*)));
        if ( bSync )
          connect(m_views[i], SIGNAL(Zooming(RenderView2D*)),
                  m_views[j], SLOT(SyncZoomTo(RenderView2D*)), Qt::UniqueConnection);
      }
    }
  }
}

void MainWindow::SetUseCommandControl(bool b)
{
  Interactor::SetUseCommandControl(b);
}

void MainWindow::SetUnifiedTitleAndToolBar(bool b)
{
  this->setUnifiedTitleAndToolBarOnMac(b);
}

bool MainWindow::IsEmpty()
{
  QStringList keys = m_layerCollections.keys();
  keys.removeOne("Supplement");
  for ( int i = 0; i < keys.size(); i++ )
  {
    if ( !m_layerCollections[keys[i]]->IsEmpty() )
    {
      return false;
    }
  }
  return true;
}

void MainWindow::AddScript(const QStringList & command)
{
  m_scripts << command;
}

void MainWindow::AddScripts(const QList<QStringList> &cmds)
{
  m_scripts << cmds;
}

void MainWindow::OnIdle()
{
  bool bBusy = IsBusy();
//  qDebug() << "busy: " << bBusy << "  script_running: " << m_bScriptRunning
 //     << "  script empty: " << m_scripts.isEmpty();
  if ( !bBusy && !m_bScriptRunning && !m_scripts.isEmpty() )
  {
    bool last_one = (m_scripts.size() == 1);
    RunScript();
    if (last_one)
      ui->widgetAllLayers->UpdateWidgets();
  }

  ui->actionViewSagittal->setChecked( m_nMainView == this->MV_Sagittal );
  ui->actionViewCoronal->setChecked( m_nMainView == this->MV_Coronal );
  ui->actionViewAxial->setChecked( m_nMainView == this->MV_Axial );
  ui->actionView3D->setChecked( m_nMainView == this->MV_3D );

  ui->actionLayout1x1->setChecked( m_nViewLayout == this->VL_1x1 );
  ui->actionLayout2x2->setChecked( m_nViewLayout == this->VL_2x2 );
  ui->actionLayout1n3->setChecked( m_nViewLayout == this->VL_1n3 );
  ui->actionLayout1n3h->setChecked( m_nViewLayout == this->VL_1n3h );

  RenderView* view = GetMainView();
  int nMode = view->GetInteractionMode();
  ui->actionNavigate->setChecked( nMode == RenderView::IM_Navigate );
  ui->actionVoxelEdit->setChecked( nMode == RenderView::IM_VoxelEdit );
  ui->actionReconEdit->setChecked( nMode == RenderView::IM_ReconEdit );
  ui->actionROIEdit->setChecked( nMode == RenderView::IM_ROIEdit );
  ui->actionMeasure->setChecked( nMode == RenderView::IM_Measure );
  ui->actionCropVolume->setChecked( nMode == RenderView::IM_VolumeCrop );
  ui->actionPointSetEdit->setChecked( nMode == RenderView::IM_PointSetEdit );

  if ( nMode == RenderView::IM_ROIEdit )
  {
    LayerROI* roi = ( LayerROI* )GetActiveLayer("ROI");
    ui->actionUndo->setEnabled( roi && roi->IsVisible() && roi->HasUndo() );
    ui->actionRedo->setEnabled( roi && roi->IsVisible() && roi->HasRedo() );
  }
  else if ( nMode == RenderView::IM_VoxelEdit || nMode == RenderView::IM_ReconEdit || nMode == RenderView::IM_Navigate )
  {
    LayerMRI* mri = ( LayerMRI* )GetActiveLayer( "MRI");
    ui->actionUndo->setEnabled( mri && mri->IsVisible() && mri->HasUndo() );
    ui->actionRedo->setEnabled( mri && mri->IsVisible() && mri->HasRedo() );
  }
  else if ( nMode == RenderView::IM_PointSetEdit )
  {
    LayerPointSet* wp = ( LayerPointSet* )GetActiveLayer("PointSet");
    ui->actionUndo->setEnabled( wp && wp->IsVisible() && wp->HasUndo() );
    ui->actionRedo->setEnabled( wp && wp->IsVisible() && wp->HasRedo() );
  }
  else
  {
    ui->actionUndo->setEnabled( false );
    ui->actionRedo->setEnabled( false );
  }

  LayerMRI* layerVolume       = (LayerMRI*)GetActiveLayer( "MRI");
  LayerSurface* layerSurface  = (LayerSurface*)GetActiveLayer( "Surface");
  LayerROI* layerROI  = (LayerROI*)GetActiveLayer( "ROI");
  LayerPointSet* layerPointSet  = (LayerPointSet*)GetActiveLayer( "PointSet");
  LayerTrack* layerTrack  = (LayerTrack*)GetActiveLayer( "Track");
  bool bHasLabelLayer = false;
  QList<Layer*> volumes = GetLayers("MRI");
  foreach (Layer* layer, volumes)
  {
    if (((LayerMRI*)layer)->GetProperty()->GetColorMap() == LayerPropertyMRI::LUT)
    {
      bHasLabelLayer = true;
      break;
    }
  }
//  LayerCollection* lc = GetCurrentLayerCollection();
  bool bHasLayer = !IsEmpty();
  ui->actionVoxelEdit       ->setEnabled( layerVolume && layerVolume->IsEditable() );
  ui->actionReconEdit       ->setEnabled( layerVolume && layerVolume->IsEditable() );
  ui->actionROIEdit         ->setEnabled( layerROI && layerROI->IsEditable() );
  ui->actionMeasure         ->setEnabled( layerVolume );
  ui->actionPointSetEdit    ->setEnabled( layerPointSet && layerPointSet->IsEditable() );
  ui->actionCropVolume      ->setEnabled( layerVolume && layerVolume->IsEditable() );
  ui->actionIntensityProject->setEnabled( layerVolume );
  ui->actionClosePointSet   ->setEnabled( !bBusy && layerPointSet );
  ui->actionCloseROI        ->setEnabled( !bBusy && layerROI );
  ui->actionCloseSurface    ->setEnabled( !bBusy && layerSurface );
  ui->actionCloseVolume     ->setEnabled( !bBusy && layerVolume ); 
  ui->actionCloseTrack      ->setEnabled( !bBusy && layerTrack );
  ui->actionReloadVolume    ->setEnabled( !bBusy && layerVolume );
  ui->actionReloadSurface    ->setEnabled( !bBusy && layerSurface );
  ui->actionCreateOptimalCombinedVolume->setEnabled( GetLayerCollection("MRI")->GetNumberOfLayers() > 1 );
//  ui->actionCycleLayer      ->setEnabled( lc && lc->GetNumberOfLayers() > 1 );
//  ui->actionReverseCycleLayer      ->setEnabled( lc && lc->GetNumberOfLayers() > 1 );
//  ui->actionHideAllLayers   ->setEnabled( lc && !lc->IsEmpty() );
//  ui->actionShowAllLayers   ->setEnabled( lc && !lc->IsEmpty() );
  ui->actionLoadDTIVolumes  ->setEnabled( !bBusy );
  ui->actionLoadVolume      ->setEnabled( !bBusy );
  ui->actionLoadROI         ->setEnabled( !bBusy && layerVolume );
  ui->actionLoadPointSet    ->setEnabled( !bBusy && layerVolume );
  ui->actionLoadSurface     ->setEnabled( !bBusy );
  ui->actionLoadTrackVolume ->setEnabled( !bBusy );
  ui->actionLoadTrack       ->setEnabled( !bBusy );
  ui->actionNewVolume       ->setEnabled( layerVolume );
  ui->actionNewROI          ->setEnabled( layerVolume );
  ui->actionNewPointSet     ->setEnabled( layerVolume );
  ui->actionRepositionSurface->setEnabled( layerSurface );
  ui->actionSmoothSurface   ->setEnabled( layerSurface );
  ui->actionRemoveIntersectionsSurface->setEnabled(layerSurface);
  ui->actionResetView       ->setEnabled( bHasLayer );
  ui->actionResetViewNearestAxis->setEnabled( bHasLayer && ui->view3D->isVisible() );
  ui->actionSaveMovieFrames ->setEnabled( bHasLayer );
  ui->actionSaveScreenshot  ->setEnabled( bHasLayer );
  ui->actionSavePoint       ->setEnabled( bHasLayer );
  ui->actionGoToPoint       ->setEnabled( bHasLayer );
  ui->actionSaveVolume      ->setEnabled( !bBusy && layerVolume && layerVolume->IsModified() );
  ui->actionSaveVolumeAs    ->setEnabled( layerVolume );
  ui->actionSaveROI         ->setEnabled( layerROI && layerROI->IsModified() );
  ui->actionSaveROIAs       ->setEnabled( layerROI );
  ui->actionSavePointSet    ->setEnabled( layerPointSet && layerPointSet->IsModified() );
  ui->actionSavePointSetAs  ->setEnabled( layerPointSet );
  ui->actionSaveSurface     ->setEnabled( !bBusy && layerSurface && layerSurface->IsModified() );
  ui->actionSaveSurfaceAs   ->setEnabled( layerSurface );
  ui->actionShowColorScale  ->setEnabled( bHasLayer );
  ui->actionShowSliceFrames  ->setEnabled(bHasLayer && ui->view3D->GetShowSlices());
  ui->actionShowSliceFrames ->blockSignals(true);
  ui->actionShowSliceFrames  ->setChecked(ui->view3D->GetShowSliceFrames());
  ui->actionShowSliceFrames ->blockSignals(false);
  ui->actionShowSlices      ->setEnabled(bHasLayer);
  ui->actionShowSlices      ->blockSignals(true);
  ui->actionShowSlices      ->setChecked(ui->view3D->GetShowSlices());
  ui->actionShowSlices->blockSignals(false);
  ui->actionShowLabelStats  ->setEnabled( (layerVolume && bHasLabelLayer) || layerROI );
  ui->actionToggleCursorVisibility  ->setEnabled( bHasLayer );
  ui->actionTogglePointSetVisibility->setEnabled( layerPointSet );
  ui->actionToggleROIVisibility     ->setEnabled( layerROI );
  ui->actionToggleSurfaceVisibility ->setEnabled( layerSurface );
  ui->actionToggleVolumeVisibility  ->setEnabled( layerVolume );
  ui->actionToggleVoxelCoordinateDisplay->setEnabled( bHasLayer );
  ui->actionTransformVolume ->setEnabled( layerVolume );
  ui->actionThresholdVolume->setEnabled(layerVolume);
  ui->actionVolumeSegmentation->setEnabled(layerVolume);
  ui->actionVolumeFilterConvolve->setEnabled( !bBusy && layerVolume && layerVolume->IsEditable() );
  ui->actionVolumeFilterMean    ->setEnabled( !bBusy && layerVolume && layerVolume->IsEditable() );
  ui->actionVolumeFilterMedian  ->setEnabled( !bBusy && layerVolume && layerVolume->IsEditable() );
  ui->actionVolumeFilterGradient->setEnabled( !bBusy && layerVolume && layerVolume->IsEditable() );
  ui->actionVolumeFilterSobel->setEnabled( !bBusy && layerVolume && layerVolume->IsEditable() );
  ui->actionVolumeFilterErode->setEnabled( !bBusy && layerVolume && layerVolume->IsEditable() );
  ui->actionVolumeFilterDilate->setEnabled( !bBusy && layerVolume && layerVolume->IsEditable() );
  ui->actionVolumeFilterOpen->setEnabled( !bBusy && layerVolume && layerVolume->IsEditable() );
  ui->actionVolumeFilterClose->setEnabled( !bBusy && layerVolume && layerVolume->IsEditable() );
  ui->actionSetCamera->setEnabled(bHasLayer);
  ui->actionSaveCamera->setEnabled(bHasLayer && GetMainView() == ui->view3D);
  ui->actionLoadCamera->setEnabled(bHasLayer && GetMainView() == ui->view3D);

  ui->actionLoadConnectome->setEnabled( !bBusy );  
  ui->actionCloseConnectome ->setEnabled( !bBusy && GetActiveLayer( "CMAT"));

  ui->actionLoadFCD->setEnabled( !bBusy );
  ui->actionCloseFCD->setEnabled( !bBusy && GetActiveLayer( "FCD"));

  ui->actionShowCoordinateAnnotation->setChecked(ui->viewAxial->GetShowCoordinateAnnotation());
  ui->actionShowColorScale->setChecked(view->GetShowScalarBar());

  ui->actionToggleCursorVisibility->setChecked(ui->viewAxial->GetCursor2D()->IsShown());

  ui->actionCopy->setEnabled(false);
  ui->actionCopyStructure->setEnabled(false);
  ui->actionPaste->setEnabled(false);
  int nWnd = GetActiveViewId();
  if ( ui->viewAxial->GetInteractionMode() == RenderView::IM_ROIEdit )
  {
    LayerROI* roi = ( LayerROI* )GetActiveLayer("ROI");
    ui->actionCopy->setEnabled( roi && roi->IsVisible() && nWnd >= 0 && nWnd < 3 );
    ui->actionPaste->setEnabled( roi && roi->IsVisible() && roi->IsEditable() &&
                                 nWnd >= 0 && nWnd < 3 && roi->IsValidToPaste( nWnd ) );
  }
  else if ( ui->viewAxial->GetInteractionMode() == RenderView::IM_VoxelEdit ||
            ui->viewAxial->GetInteractionMode() == RenderView::IM_ReconEdit )
  {
    LayerMRI* mri = ( LayerMRI* )GetActiveLayer("MRI");
    ui->actionCopy->setEnabled( mri && mri->IsVisible() && nWnd >= 0 && nWnd < 3 );
    ui->actionCopyStructure->setEnabled(mri && mri->IsVisible() && nWnd >= 0 && nWnd < 3);
    ui->actionPaste->setEnabled( mri && mri->IsVisible() && mri->IsEditable() &&
                                 nWnd >= 0 && nWnd < 3 && mri->IsValidToPaste( nWnd ) );
  }

  for ( int i = 0; i < m_actionRecentVolumes.size(); i++ )
  {
    m_actionRecentVolumes[i]->setEnabled( !bBusy );
  }
  for ( int i = 0; i < m_actionRecentSurfaces.size(); i++ )
  {
    m_actionRecentSurfaces[i]->setEnabled( !bBusy );
  }

  bool bMeasureWindowVisible = m_toolWindowMeasure->isVisible();
  bool bEditWindowVisible = m_toolWindowEdit->isVisible();
  bool bROIEditWindowVisible = m_toolWindowROIEdit->isVisible();
  m_toolWindowMeasure->setVisible( nMode == RenderView::IM_Measure );
  m_toolWindowEdit->setVisible( nMode == RenderView::IM_VoxelEdit || nMode == RenderView::IM_ReconEdit );
  m_toolWindowROIEdit->setVisible( nMode == RenderView::IM_ROIEdit );

  if ( !m_dlgCropVolume->isVisible() && nMode == RenderView::IM_VolumeCrop )
  {
    SetMode(RenderView::IM_Navigate);
  }
  else if ( m_dlgCropVolume->isVisible() && nMode != RenderView::IM_VolumeCrop )
  {
    m_dlgCropVolume->hide();
  }

  ui->actionShowCommandConsole->setChecked(m_term->isVisible());
  ui->actionTimeCourse->setEnabled((layerVolume && layerVolume->GetNumberOfFrames() > 1 && !layerVolume->GetCorrelationSurface()) ||
                                   (layerSurface && layerSurface->GetActiveOverlay() && layerSurface->GetActiveOverlay()->GetNumberOfFrames() > 1));
  if (ui->actionTimeCourse->isEnabled())
    ui->actionTimeCourse->setChecked(m_wndTimeCourse->isVisible());

  ui->actionToggleSurfaceSpline->setEnabled(layerSurface);
  ui->actionToggleSurfaceSpline->setChecked(this->m_bSplinePicking);

  if ((!bEditWindowVisible && m_toolWindowEdit->isVisible()) ||
      (!bMeasureWindowVisible && m_toolWindowMeasure->isVisible()) ||
      (!bROIEditWindowVisible && m_toolWindowROIEdit->isVisible()))
  {
    QTimer::singleShot(50, this, SLOT(SlotActivateWindow()));
  }
}

bool MainWindow::IsBusy()
{
  return m_threadIOWorker->isRunning() || m_bProcessing || m_threadVolumeFilter->isRunning();
}

void MainWindow::RequestRedraw()
{
  for ( int i = 0; i < 4; i++ )
  {
    if ( m_views[i]->isVisible())
    {
      m_views[i]->RequestRedraw();
    }
  }
}

void MainWindow::RunScript()
{
  if (m_scripts.isEmpty())
  {
    return;
  }

  QStringList sa = m_scripts[0];
  m_scripts.removeAt(0);
  m_bScriptRunning = true;
  QString cmd = sa[0].toLower();
  if (cmd == "loadcommand")
  {
    CommandLoadCommand(sa);
  }
  else if (cmd == "loadsubject")
  {
    CommandLoadSubject(sa);
  }
  else if ( cmd == "hidelayer" )
  {
    CommandHideLayer(sa);
  }
  else if ( cmd == "unloadlayer" )
  {
    CommandUnloadLayer(sa);
  }
  else if ( cmd == "loadvolume" )
  {
    CommandLoadVolume( sa );
  }
  else if ( cmd == "loaddti" )
  {
    CommandLoadDTI( sa );
  }
  else if ( cmd == "loadtrackvolume" || cmd == "loadvolumetrack")
  {
    CommandLoadVolumeTrack( sa );
  }
  else if ( cmd == "loadsurface" )
  {
    CommandLoadSurface( sa );
  }
  else if ( cmd == "loadsurfacevector" )
  {
    CommandLoadSurfaceVector( sa );
  }
  else if ( cmd == "loadsurfacecurvature" )
  {
    CommandLoadSurfaceCurvature( sa );
  }
  else if (cmd == "setsurfacecurvaturemap")
  {
      CommandSetSurfaceCurvatureMap( sa);
  }
  else if ( cmd == "loadsurfaceoverlay" )
  {
    CommandLoadSurfaceOverlay( sa );
  }
  else if ( cmd == "loadsurfaceannotation" )
  {
    CommandLoadSurfaceAnnotation( sa );
  }
  else if ( cmd == "loadsurfacelabel" )
  {
    CommandLoadSurfaceLabel( sa );
  }
  else if ( cmd == "loadsurfacespline")
  {
    CommandLoadSurfaceSpline( sa );
  }
  else if ( cmd == "loadconnectome" )
  {
    CommandLoadConnectomeMatrix( sa );
  }
  else if ( cmd == "loadfcd")
  {
    CommandLoadFCD( sa );
  }
  else if ( cmd == "loadroi" || sa[0] == "loadlabel" )
  {
    CommandLoadROI( sa );
  }
  else if ( cmd == "loadwaypoints" )
  {
    CommandLoadWayPoints( sa );
  }
  else if ( cmd == "loadcontrolpoints" )
  {
    CommandLoadControlPoints( sa );
  }
  else if ( cmd == "loadpvolumes" )
  {
    CommandLoadPVolumes( sa );
  }
  else if ( cmd == "loadtrack")
  {
    CommandLoadTrack(sa);
  }
  else if ( cmd == "screencapture" )
  {
    CommandScreenCapture( sa );
  }
  else if ( cmd == "quit" || cmd == "exit" )
  {
    close();
  }
  else if (cmd == "center")
  {
    for (int i = 0; i < 3; i++)
      ((RenderView2D*)m_views[i])->CenterAtCursor();
  }
  else if ( cmd == "setviewport" )
  {
    CommandSetViewport( sa );
  }
  else if ( cmd == "setviewsize" )
  {
    CommandSetViewSize( sa );
  }
  else if ( cmd == "zoom" )
  {
    CommandZoom( sa );
  }
  else if ( cmd == "setcamera")
  {
    CommandSetCamera(sa);
  }
  else if ( cmd == "ras" )
  {
    CommandSetRAS( sa );
  }
  else if ( cmd == "slice" )
  {
    CommandSetSlice( sa );
  }
  else if ( cmd == "setcolormap" )
  {
    CommandSetColorMap( sa );
  }
  else if ( cmd == "setheatscaleoptions" )
  {
    CommandSetHeadScaleOptions( sa );
  }
  else if ( cmd == "setlut" )
  {
    CommandSetLUT( sa );
  }
  else if ( cmd == "setopacity" )
  {
    CommandSetOpacity( sa );
  }
  else if ( cmd == "setsmoothed")
  {
    CommandSetSmoothed( sa );
  }
  else if ( cmd == "setdisplayoutline")
  {
    CommandSetLabelOutline(sa);
  }
  else if ( cmd == "setdisplayisosurface" )
  {
    CommandSetDisplayIsoSurface( sa );
  }
  else if ( cmd == "setisosurfacecolor" )
  {
    CommandSetIsoSurfaceColor( sa );
  }
  else if (cmd == "setisosurfaceupsample")
  {
    CommandSetIsoSurfaceUpsample( sa );
  }
  else if ( cmd == "loadisosurfaceregion" )
  {
    CommandLoadIsoSurfaceRegion( sa );
  }
  else if ( cmd == "setsurfaceoverlaymethod" )
  {
    CommandSetSurfaceOverlayMethod( sa );
  }
  else if (cmd == "setsurfaceoverlaycolormap")
  {
    CommandSetSurfaceOverlayColormap( sa );
  }
  else if ( cmd == "setsurfaceoverlayopacity" )
  {
    CommandSetSurfaceOverlayOpacity( sa );
  }
  else if ( cmd == "setsurfaceoffset" )
  {
    CommandSetSurfaceOffset( sa );
  }
  else if ( cmd == "setpointsetcolor" )
  {
    CommandSetPointSetColor( sa );
  }
  else if ( cmd == "setpointsetradius" )
  {
    CommandSetPointSetRadius( sa );
  }
  else if (cmd == "setpointsetheatmap")
  {
    CommandSetPointSetHeatmap( sa );
  }
  else if ( cmd == "setdisplayvector" )
  {
    CommandSetDisplayVector( sa );
  }
  else if ( cmd == "setdisplaytensor" )
  {
    CommandSetDisplayTensor( sa );
  }
  else if ( cmd == "setsurfacecolor" )
  {
    CommandSetSurfaceColor( sa );
  }
  else if ( cmd == "setsurfaceedgecolor" )
  {
    CommandSetSurfaceEdgeColor( sa );
  }
  else if ( cmd == "setsurfaceedgethickness" )
  {
    CommandSetSurfaceEdgeThickness( sa );
  }
  else if ( cmd == "displaysurfacevertex" )
  {
    CommandSetDisplaySurfaceVertex( sa );
  }
  else if ( cmd == "setsurfacevertexcolor" )
  {
    CommandSetSurfaceVertexColor( sa );
  }
  else if ( cmd == "setsurfacelabeloutline" )
  {
    CommandSetSurfaceLabelOutline( sa );
  }
  else if (cmd == "setsurfacelabelcolor")
  {
    CommandSetSurfaceLabelColor( sa );
  }
  else if ( cmd == "setsurfaceannotationoutline" )
  {
    CommandSetSurfaceAnnotationOutline( sa );
  }
  else if ( cmd == "setlayername" )
  {
    CommandSetLayerName( sa );
  }
  else if ( cmd == "locklayer" )
  {
    CommandLockLayer( sa );
  }
  else if ( cmd == "showlayer" )
  {
    CommandShowLayer( sa );
  }
  else if ( cmd == "gotolabel" || cmd == "gotostructure")
  {
    CommandGotoLabel( sa );
  }
  else if (cmd == "showcolorscale")
  {
    OnShowColorBar(true);
  }
  else if (cmd == "setvolumemask")
  {
    CommandSetVolumeMask( sa );
  }
  else if (cmd == "savelayer")
  {
    CommandSaveLayer(sa);
  }
  m_bScriptRunning = false;
}

void MainWindow::ClearScripts()
{
  m_scripts.clear();
  m_bScriptRunning = false;
}

void MainWindow::CommandLoadCommand(const QStringList &sa)
{
  QFile file(sa[1]);
  if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
  {
    cerr << "Can not open for read: " << qPrintable(sa[1]) << ".\n";
    return;
  }
  QStringList lines = QString(file.readAll()).trimmed().split("\n", QString::SkipEmptyParts);
  foreach (QString line, lines)
  {
    QStringList args = line.trimmed().split(QRegExp("\\s+"), QString::SkipEmptyParts);
    if (args.size() > 0 &&
        ( args[0].toLower() == "freeview" || args[0].toLower() == "fv"))
    {
      args.removeFirst();
    }
    args.prepend("freeview");
    ParseCommand(args.join(" "));
  }
}

void MainWindow::CommandLoadSubject(const QStringList &sa)
{
  QString subject_path = QProcessEnvironment::systemEnvironment().value("SUBJECTS_DIR");
  if (subject_path.isEmpty())
  {
    cerr << "SUBJECTS_DIR is not set. Can not load subject.\n";
    return;
  }
  subject_path += "/" + sa[1];
  QString args = QString("freeview -v %1/mri/wm.mgz:colormap=heat "
                         "%1/mri/brain.mgz "
                         "%1/mri/orig.mgz "
                         "%1/mri/aseg.mgz:colormap=lut "
                         "-f %1/surf/lh.white "
                         "%1/surf/rh.white "
                         "%1/surf/lh.pial:edgecolor=red "
                         "%1/surf/rh.pial:edgecolor=red "
                         "%1/surf/lh.orig:edgecolor=green "
                         "%1/surf/rh.orig:edgecolor=green "
                         "%1/surf/lh.inflated:annot=aparc:visible=0 "
                         "%1/surf/rh.inflated:annot=aparc:visible=0 ").arg(subject_path);
  MyCmdLineParser parser(m_cmdParser);
  ParseCommand(&parser, args);
}

void MainWindow::CommandHideLayer(const QStringList &sa)
{
  QString type = sa[1].toLower();
  if (type == "volume" || type == "mri")
  {
    Layer* layer = GetActiveLayer("MRI");
    if (layer)
      layer->SetVisible(false);
  }
  else if (type == "surface" || type == "surf")
  {
    Layer* layer = GetActiveLayer("Surface");
    if (layer)
      layer->SetVisible(false);
  }
  else if (type == "label" || type == "roi")
  {
    Layer* layer = GetActiveLayer("ROI");
    if (layer)
      layer->SetVisible(false);
  }
}

void MainWindow::CommandUnloadLayer(const QStringList &sa)
{
  QString type = sa[1].toLower();
  if (type == "volume" || type == "mri")
  {
    OnCloseVolume();
  }
  else if (type == "surface" || type == "surf")
  {
    OnCloseSurface();
  }
  else if (type == "label" || type == "roi")
  {
    OnCloseROI();
  }
}

void MainWindow::CommandLoadVolume( const QStringList& sa )
{
  QStringList sa_vol =sa[1].split(":");
  QString fn = sa_vol[0];
  QString reg_fn;
  QStringList scales;
  QString colormap = m_strDefaultColorMapType;
  QString colormap_scale = "grayscale";
  QString lut_name;
  QString vector_display = "no",
          vector_inversion = "none",
          vector_render = "line",
          tensor_display = "no",
          tensor_render = "boxoid";
  int nSampleMethod = m_nDefaultSampleMethod;
  bool bConform = m_bDefaultConform;
  QString gotoLabelName;
  QVariantMap sup_data;
  for ( int i = 1; i < sa_vol.size(); i++ )
  {
    QString strg = sa_vol[i];
    int n = strg.indexOf( "=" );
    if ( n != -1 )
    {
      QString subOption = strg.left(n).toLower();
      QString subArgu = strg.mid( n + 1 );
      if ( subOption == "colormap" )
      {
        colormap = subArgu.toLower();
      }
      else if ( subOption == "grayscale" ||
                subOption == "heatscale" ||
                subOption == "colorscale" )
      {
        colormap_scale = subOption;    // colormap scale might be different from colormap!
        scales = subArgu.split(",");
      }
      else if ( subOption == "heatscaleoption" ||
                subOption == "heatscaleoptions" )
      {
        QStringList script("setheatscaleoptions");
        script << subArgu.split(",");
        m_scripts.insert( 0, script );
      }
      else if ( subOption == "lut" )
      {
        lut_name = subArgu;
        if ( lut_name.isEmpty() )
        {
          cerr << "Missing lut name.\n";
        }
      }
      else if ( subOption == "vector" )
      {
        vector_display = subArgu.toLower();
        if ( vector_display.isEmpty() )
        {
          cerr << "Missing vector display argument.\n";
        }
      }
      else if ( subOption == "tensor" )
      {
        tensor_display = subArgu.toLower();
        if ( tensor_display.isEmpty() )
        {
          cerr << "Missing tensor display argument.\n";
        }
      }
      else if ( subOption == "inversion" ||
                subOption == "invert")
      {
        vector_inversion = subArgu.toLower();
        if ( vector_inversion.isEmpty() )
        {
          cerr << "Missing inversion argument.\n";
          vector_inversion = "none";
        }
      }
      else if ( subOption == "render" )
      {
        vector_render = subArgu.toLower();
        tensor_render = vector_render;
        if ( vector_render.isEmpty() )
        {
          cerr << "Missing render argument.\n";
          vector_render = "line";
          tensor_render = "boxoid";
        }
      }
      else if ( subOption == "reg" )
      {
        reg_fn = subArgu;
      }
      else if ( subOption == "sample" )
      {
        if ( subArgu.toLower() == "nearest" )
        {
          nSampleMethod = SAMPLE_NEAREST;
        }
        else if ( subArgu.toLower() == "trilinear" )
        {
          nSampleMethod = SAMPLE_TRILINEAR;
        }
        else if ( subArgu.toLower() == "cubic" )
        {
          nSampleMethod = SAMPLE_CUBIC_BSPLINE;
        }
      }
      else if ( subOption == "opacity" )
      {
        m_scripts.insert( 0, QStringList() << "setopacity" << subArgu );
      }
      else if ( subOption == "outline")
      {
        m_scripts.insert( 0, QStringList() << "setdisplayoutline" << subArgu );
      }
      else if ( subOption == "isosurface" )
      {
        QStringList script("setdisplayisosurface");
        QStringList args = subArgu.split( ",");
        if ( args.size() > 0 && args[0].size() > 0 )
        {
          script << args[0];
        }
        if ( args.size() > 1 && args[1].size() > 0 )
        {
          script << args[1];
        }
        m_scripts.insert( 0, script );
      }
      else if ( subOption == "upsample_isosurface")
      {
        m_scripts.insert( 0,  (QStringList("setisosurfaceupsample") << subArgu) );
      }
      else if (subOption == "color")
      {
        m_scripts.insert( 0,  (QStringList("setisosurfacecolor") << subArgu) );
      }
      else if ( subOption == "surface_region" || subOption == "surface_regions" )
      {
        m_scripts.insert( 0, (QStringList("loadisosurfaceregion") << QFileInfo(subArgu).absoluteFilePath()) );
      }
      else if ( subOption == "name" )
      {
        m_scripts.insert( 0, QStringList("setlayername") << "MRI" << subArgu );
      }
      else if ( subOption == "lock" || subOption == "locked" )
      {
        m_scripts.insert( 0, QStringList("locklayer") << "MRI" << subArgu );
      }
      else if ( subOption == "visible" )
      {
        m_scripts.insert( 0, QStringList("showlayer") << "MRI" << subArgu );
      }
      else if ( subOption == "gotolabel" || subOption == "structure")
      {
        m_scripts.insert(0, QStringList("gotolabel") << subArgu);
        gotoLabelName = subArgu;
      }
      else if (subOption == "mask")
      {
        m_scripts.insert(0, QStringList("setvolumemask") << subArgu);
      }
      else if (subOption == "basis")
      {
        sup_data["Basis"] = subArgu.toInt();
      }
      else if (subOption == "percentile")
      {
        sup_data["Percentile"] = true;
      }
      else if (subOption == "smoothed" || subOption == "smooth")
      {
        m_scripts.insert(0, QStringList("setsmoothed") << subArgu);
      }
      else if (!subOption.isEmpty())
      {
        cerr << "Unrecognized sub-option flag '" << strg.toAscii().constData() << "'.\n";
        return;
      }
    }
    else
    {
      cerr << "Unrecognized sub-option flag '" << strg.toAscii().constData() << "'.\n";
      return;
    }
  }
  bool bResample = false;
  if ( sa[sa.size()-1] == "r")
  {
    bResample = true;
  }

  if ( scales.size() > 0 || (colormap != "grayscale" && colormap != "greyscale") )
  {
    QStringList script("setcolormap");
    script << colormap << colormap_scale
              << scales;
    m_scripts.insert( 0, script );
  }

  if ( !lut_name.isEmpty() )
  {
    m_scripts.insert( 0, QStringList("setlut") << lut_name );
  }

  if ( !tensor_display.isEmpty() && tensor_display != "no" )
  {
    QStringList script = QStringList("setdisplaytensor") <<
                     tensor_display <<
                     tensor_render <<
                     vector_inversion;
    m_scripts.insert( 0, script );
  }
  else if ( !vector_display.isEmpty() && vector_display != "no" )
  {
    QStringList script = QStringList("setdisplayvector") <<
                     vector_display <<
                     vector_render <<
                     vector_inversion;
    m_scripts.insert( 0, script );
  }

  int nView = this->GetMainViewId();
  if (nView > 2)
  {
    nView = 0;
  }
  int orientation = nView;
  if (orientation == 0 )
  {
    orientation = 1;
  }
  else if (orientation == 1)
  {
    orientation = 0;
  }
  LoadVolumeFile( fn, reg_fn, bResample, nSampleMethod, bConform, orientation, gotoLabelName, sup_data );
}

void MainWindow::CommandSetColorMap( const QStringList& sa )
{
  int nColorMap = LayerPropertyMRI::Grayscale;;
  QString strg = sa[1];
  if ( strg == "heat" || strg == "heatscale" )
  {
    nColorMap = LayerPropertyMRI::Heat;
  }
  else if ( strg == "jet" || strg == "jetscale" )
  {
    nColorMap = LayerPropertyMRI::Jet;
  }
  else if ( strg == "lut" )
  {
    nColorMap = LayerPropertyMRI::LUT;
  }
  else if ( strg == "gecolor" || strg == "ge_color" )
  {
    nColorMap = LayerPropertyMRI::GEColor;
  }
  else if ( strg == "nih" )
  {
    nColorMap = LayerPropertyMRI::NIH;
  }
  else if ( strg == "pet" )
  {
    nColorMap = LayerPropertyMRI::PET;
  }
  else if ( strg != "grayscale" )
  {
    cerr << "Unrecognized colormap name '" << strg.toAscii().constData() << "'.\n";
  }

  int nColorMapScale = LayerPropertyMRI::Grayscale;
  strg = sa[2];
  if ( strg == "heatscale" )
  {
    nColorMapScale = LayerPropertyMRI::Heat;
  }
  else if ( strg == "colorscale" )
  {
    nColorMapScale = LayerPropertyMRI::Jet;
  }
  else if ( strg == "lut" )
  {
    nColorMapScale = LayerPropertyMRI::LUT;
  }

  QList<double> pars;
  for ( int i = 3; i < sa.size(); i++ )
  {
    bool bOK;
    double dValue = sa[i].toDouble(&bOK);
    if ( !bOK )
    {
      cerr << "Invalid color scale value(s). \n";
      break;
    }
    else
    {
      pars << dValue;
    }
  }

  SetVolumeColorMap( nColorMap, nColorMapScale, pars );
}

void MainWindow::CommandSetHeadScaleOptions( const QStringList& sa )
{
  if ( GetLayerCollection( "MRI" )->GetActiveLayer() )
  {
    LayerPropertyMRI* p = ( (LayerMRI*)GetLayerCollection( "MRI" )->GetActiveLayer() )->GetProperty();
    for ( int i = 1; i < sa.size(); i++ )
    {
      if ( sa[i] == "invert" )
      {
        p->SetHeatScaleInvert( true );
      }
      else if ( sa[i] == "truncate" )
      {
        p->SetHeatScaleTruncate( true );
      }
    }
  }
}

void MainWindow::CommandSetLayerName( const QStringList& cmd )
{
  if ( cmd.size() > 2 )
  {
    LayerCollection* lc = GetLayerCollection( cmd[1] );
    if ( lc && !lc->IsEmpty() )
    {
      lc->GetActiveLayer()->SetName( cmd[2] );
    }
  }
}

void MainWindow::CommandLockLayer( const QStringList& cmd )
{
  if ( cmd.size() > 2 && ( cmd[2] == "1" || cmd[2].toLower() == "true" ) )
  {
    LayerCollection* lc = GetLayerCollection( cmd[1] );
    if ( lc && !lc->IsEmpty() )
    {
      lc->GetActiveLayer()->Lock( true );
    }
  }
}

void MainWindow::CommandShowLayer( const QStringList& cmd )
{
  if ( cmd.size() > 2 && ( cmd[2] == "0" || cmd[2].toLower() == "false" ) )
  {
    LayerCollection* lc = GetLayerCollection( cmd[1] );
    if ( lc && !lc->IsEmpty() )
    {
      lc->GetActiveLayer()->SetVisible( false );
    }
  }
}

void MainWindow::CommandSetLabelOutline(const QStringList &cmd)
{
  QString stemp = cmd[1].toLower();
  if ( stemp == "yes"|| stemp == "true" || stemp == "1" || stemp == "on")
  {
    LayerMRI* mri = (LayerMRI*)GetLayerCollection( "MRI" )->GetActiveLayer();
    if ( mri )
    {
      mri->GetProperty()->SetShowLabelOutline(true);
    }
  }
}

void MainWindow::CommandSetSmoothed(const QStringList &cmd)
{
  QString stemp = cmd[1].toLower();
  if ( stemp == "yes"|| stemp == "true" || stemp == "1" || stemp == "on")
  {
    LayerMRI* mri = (LayerMRI*)GetLayerCollection( "MRI" )->GetActiveLayer();
    if ( mri )
    {
      mri->GetProperty()->SetTextureSmoothing(1);
    }
  }
}

void MainWindow::CommandSetDisplayVector( const QStringList& cmd )
{
  if ( cmd[1].toLower() == "yes" || cmd[1].toLower() == "true" || cmd[1] == "1" )
  {
    LayerMRI* mri = (LayerMRI*)GetLayerCollection( "MRI" )->GetActiveLayer();
    if ( mri )
    {
      if ( !mri->IsTypeOf( "DTI" ) && mri->GetNumberOfFrames() < 3 )
      {
        cerr << "Volume has less than 3 frames. Can not display as vectors.\n";
      }
      else
      {
        mri->GetProperty()->SetDisplayVector( true );

        if ( cmd[2].toLower() == "line" )
        {
          mri->GetProperty()->SetVectorRepresentation( LayerPropertyMRI::VR_Line );
        }
        else if ( cmd[2].toLower() == "bar" )
        {
          mri->GetProperty()->SetVectorRepresentation( LayerPropertyMRI::VR_Bar );
        }
        else
        {
          cerr << "Unrecognized argument '" << cmd[2].toAscii().constData() << "' for vector rendering.\n";
        }

        if ( cmd[3].toLower() != "none" )
        {
          if ( cmd[3].toLower() == "x" )
          {
            mri->GetProperty()->SetVectorInversion( LayerPropertyMRI::VI_X );
          }
          else if ( cmd[3].toLower() == "y" )
          {
            mri->GetProperty()->SetVectorInversion( LayerPropertyMRI::VI_Y );
          }
          else if ( cmd[3].toLower() == "z" )
          {
            mri->GetProperty()->SetVectorInversion( LayerPropertyMRI::VI_Z );
          }
          else
          {
            cerr << "Unknown inversion flag '" << cmd[2].toAscii().constData() << "'.\n";
          }
        }
      }
    }
  }
}


void MainWindow::CommandSetDisplayTensor( const QStringList& cmd )
{
  if ( cmd[1].toLower() == "yes" || cmd[1].toLower() == "true" || cmd[1].toLower() == "1" )
  {
    LayerMRI* mri = (LayerMRI*)GetLayerCollection( "MRI" )->GetActiveLayer();
    if ( mri )
    {
      if ( mri->GetNumberOfFrames() < 9 )
      {
        cerr << "Volume has less than 9 frames. Can not display as tensor.\n";
      }
      else
      {
        mri->GetProperty()->SetDisplayTensor( true );

        if ( cmd[2].toLower().indexOf( "box" ) != -1 )
        {
          mri->GetProperty()->SetVectorRepresentation( LayerPropertyMRI::TR_Boxoid );
        }
        else if ( cmd[2].toLower().indexOf( "ellips" ) != -1 )
        {
          mri->GetProperty()->SetVectorRepresentation( LayerPropertyMRI::TR_Ellipsoid );
        }
        else
        {
          cerr << "Unrecognized argument '" << cmd[2].toAscii().constData() << "' for tensor rendering.\n";
        }

        if ( cmd[3].toLower() != "none" )
        {
          if ( cmd[3].toLower() == "x" )
          {
            mri->GetProperty()->SetVectorInversion( LayerPropertyMRI::VI_X );
          }
          else if ( cmd[3].toLower() == "y" )
          {
            mri->GetProperty()->SetVectorInversion( LayerPropertyMRI::VI_Y );
          }
          else if ( cmd[3].toLower() == "z" )
          {
            mri->GetProperty()->SetVectorInversion( LayerPropertyMRI::VI_Z );
          }
          else
          {
            cerr << "Unknown inversion flag '" << cmd[2].toAscii().constData() << "'.\n";
          }
        }
      }
    }
  }
}


void MainWindow::CommandSetLUT( const QStringList& sa )
{
  LayerMRI* mri = (LayerMRI*)GetLayerCollection( "MRI" )->GetActiveLayer();
  if ( mri )
  {
    COLOR_TABLE* ct = m_luts->LoadColorTable( sa[1] );
    if ( ct )
    {
      mri->GetProperty()->SetLUTCTAB( ct );
    }
  }
}

void MainWindow::CommandSetOpacity( const QStringList& sa )
{
  LayerMRI* mri = (LayerMRI*)GetLayerCollection( "MRI" )->GetActiveLayer();
  if ( mri )
  {
    bool bOK;
    double dValue = sa[1].toDouble(&bOK);
    if ( bOK )
    {
      mri->GetProperty()->SetOpacity( dValue );
    }
    else
    {
      cerr << "Opacity value is not valid.\n";
    }
  }
}

void MainWindow::CommandSetDisplayIsoSurface( const QStringList& sa )
{
  LayerMRI* mri = (LayerMRI*)GetLayerCollection( "MRI" )->GetActiveLayer();
  if ( mri )
  {
    bool bOK;
    double dValue;
    if ( sa.size() > 1 )
    {
      dValue = sa[1].toDouble(&bOK);
      if ( bOK )
      {
        mri->GetProperty()->SetContourMinThreshold( dValue );
      }
      else if ( sa[1].toLower() != "on" )
      {
        cerr << "Isosurface threshold value is not valid.\n";
      }
    }
    if ( sa.size() > 2 )
    {
      dValue = sa[2].toDouble(&bOK);
      if ( bOK )
      {
        mri->GetProperty()->SetContourMaxThreshold( dValue );
      }
      else
      {
        cerr << "Isosurface threshold value is not valid.\n";
      }
    }
    mri->GetProperty()->SetShowAsContour( true );
  }
}

void MainWindow::CommandSetIsoSurfaceColor(const QStringList &cmd)
{
  LayerMRI* mri = (LayerMRI*)GetLayerCollection( "MRI" )->GetActiveLayer();
  if ( mri )
  {
    QColor color = ParseColorInput( cmd[1] );
    if ( color.isValid() )
    {
      mri->GetProperty()->SetContourColor( color.redF(), color.greenF(), color.blueF() );
    }
    else
    {
      cerr << "Invalid color name or value " << cmd[1].toAscii().constData() << ".\n";
    }
  }
}

void MainWindow::CommandSetIsoSurfaceUpsample(const QStringList &cmd)
{
  LayerMRI* mri = (LayerMRI*)GetLayerCollection( "MRI" )->GetActiveLayer();
  if ( mri )
  {
    if (cmd[1].toLower() == "on" || cmd[1].toLower() == "true" || cmd[1].toLower() == "1")
    {
      mri->GetProperty()->SetContourUpsample(true);
    }
  }
}

void MainWindow::CommandLoadIsoSurfaceRegion( const QStringList& sa )
{
  LayerMRI* mri = (LayerMRI*)GetLayerCollection( "MRI" )->GetActiveLayer();
  if ( mri )
  {
    if ( sa.size() > 1 )
    {
      if ( !mri->LoadSurfaceRegions( sa[1] ) )
      {
        cerr << "Can not load surfacer region(s) from " << sa[1].toAscii().constData() << ".\n";
      }
    }
  }
}

void MainWindow::CommandLoadDTI( const QStringList& sa )
{
  bool bResample = false;
  if ( sa.size() > 3 && sa[3] == "r" )
  {
    bResample = true;
  }

  if ( sa.size() > 2 )
  {
    QStringList sa_vol = sa[1].split(":");
    QString fn = sa_vol[0];
    QString strg, reg_fn;
    QString vector_display = "no",
            vector_inversion = "none",
            vector_render = "line";

    for ( int i = 1; i < sa_vol.size(); i++ )
    {
      QString strg = sa_vol[i];
      int n = strg.indexOf( "=" );
      if ( n != -1 )
      {
        if ( strg.left( n ).toLower() == "vector" )
        {
          vector_display = strg.mid( n + 1 ).toLower();
          if ( vector_display.isEmpty() )
          {
            cerr << "Missing vector display argument.\n";
          }
        }
        else if ( strg.left( n ).toLower() == "inversion" ||
                  strg.left( n ).toLower() == "invert" )
        {
          vector_inversion = strg.mid( n + 1 ).toLower();
          if ( vector_inversion.isEmpty() )
          {
            cerr << "Missing vector inversion argument.\n";
            vector_inversion = "none";
          }
        }
        else if ( strg.left( n ).toLower() == "render" )
        {
          vector_render = strg.mid( n + 1 ).toLower();
          {
            if ( vector_render.isEmpty() )
            {
              cerr << "Missing vector render argument.\n";
              vector_render = "line";
            }
          }
        }
        else if ( strg.left( n ).toLower() == "reg" )
        {
          reg_fn = strg.mid( n + 1 );
        }
      }
    }

    if ( !vector_display.isEmpty() && vector_display != "no" )
    {
      QStringList script("setdisplayvector");
      script << vector_display << vector_render << vector_inversion;

      m_scripts.insert( 0, script );
    }

    this->LoadDTIFile( fn, sa[2], sa.size() > 3 ? sa[3] : "", reg_fn, bResample );
  }
}

void MainWindow::CommandLoadVolumeTrack( const QStringList& sa )
{
  bool bResample = false;
  if ( sa.last() == "r" )
  {
    bResample = true;
  }
  this->LoadVolumeTrackFile(sa[1], bResample);
}

void MainWindow::CommandLoadPVolumes( const QStringList& cmd )
{
  QStringList files = cmd[1].split(";");
  QString lut = "";
  if ( cmd.size() > 3 )
  {
    lut = cmd[3];
    COLOR_TABLE* ct = m_luts->LoadColorTable( lut );
    if ( !ct )
    {
      cerr << "Can not load look up table " << lut.toAscii().constData() << ".\n";
      return;
    }
  }
  QString prefix = cmd[2];
  if (prefix == "n/a")
  {
    prefix = "";
  }
  this->LoadPVolumeFiles( files, prefix, lut );
}

void MainWindow::CommandLoadConnectomeMatrix(const QStringList& cmd )
{
  if (cmd.size() < 3)
    return;

  QStringList options = cmd[1].split(":");
  QString fn = options[0];
  QString lut;
  for ( int i = 1; i < options.size(); i++ )
  {
    QString strg = options[i];
    int n = strg.indexOf( "=" );
    if ( n != -1 )
    {
      QString option = strg.left( n ).toLower();
      QString argu = strg.mid( n+1 );
      if ( option == "lut" )
      {
        lut = argu;
      }
      else
      {
        cerr << "Unrecognized sub-option flag '" << strg.toAscii().constData() << "'.\n";
      }
    }
  }

  this->LoadConnectomeMatrixFile( fn, cmd[2], lut );
}

void MainWindow::LoadConnectomeMatrixFile(const QString &fn_cmat, const QString &fn_parcel, const QString &fn_ctab)
{
  LayerConnectomeMatrix* layer = new LayerConnectomeMatrix(m_layerVolumeRef);
  layer->SetFileName(fn_cmat);
  layer->SetParcelFilename(fn_parcel);
  layer->SetName(QFileInfo(fn_cmat).completeBaseName());  
  COLOR_TABLE* ct = NULL;
  if (!fn_ctab.isEmpty())
    m_luts->LoadColorTable( fn_ctab );
  if (ct)
    layer->SetColorTable(ct);
  else
    layer->SetColorTable(m_luts->GetColorTable(0));
  m_threadIOWorker->LoadConnectomeMatrix( layer );
//  m_statusBar->StartTimer();
}

void MainWindow::OnCloseConnectomeMatrix()
{
  LayerConnectomeMatrix* layer = (LayerConnectomeMatrix*)GetActiveLayer( "CMAT" );
  if ( !layer )
  {
    return;
  }

  GetLayerCollection( "CMAT" )->RemoveLayer( layer );
}


void MainWindow::CommandLoadROI( const QStringList& cmd )
{
  QStringList options = cmd[1].split(":");
  QString fn = options[0];
  QString ref;
  QColor color = Qt::yellow;
  double opacity = 1;
  double threshold = 0;
  for ( int i = 1; i < options.size(); i++ )
  {
    QString strg = options[i];
    int n = strg.indexOf( "=" );
    if ( n != -1 )
    {
      QString option = strg.left( n ).toLower();
      QString argu = strg.mid( n+1 );
      if ( option == "ref" || option == "template" )
      {
        ref = argu;
      }
      else if (option == "color")
      {
        color = ParseColorInput(argu);
        if (!color.isValid())
          color = Qt::yellow;
      }
      else if (option == "opacity")
      {
        opacity = argu.toDouble();
      }
      else if (option == "threshold")
      {
        threshold = argu.toDouble();
      }
      else
      {
        cerr << "Unrecognized sub-option flag '" << strg.toAscii().constData() << "'.\n";
      }
    }
  }

  LoadROIFile( fn, ref, color, opacity, threshold );
}

void MainWindow::CommandLoadTrack(const QStringList &cmd)
{
  QString fn = cmd[1];
  LoadTrackFile( fn );
}

void MainWindow::CommandLoadSurface( const QStringList& cmd )
{
  QStringList rawoverlay_list = cmd[1].split("overlay=", QString::SkipEmptyParts, Qt::CaseInsensitive);
  QStringList overlay_list;
  for (int i = 1; i < rawoverlay_list.size(); i++)
  {
    QStringList sublist = rawoverlay_list[i].split("correlation=", QString::SkipEmptyParts, Qt::CaseInsensitive);
    overlay_list << QString("overlay=") + sublist[0];
    for (int j = 1; j < sublist.size(); j++)
      overlay_list << QString("correlation=") + sublist[i];
  }

  overlay_list.insert(0, rawoverlay_list[0]);

  QString surface_fn;
  QString fn_patch = "";
  QString fn_target = "";
  QStringList sup_files;
  QStringList valid_overlay_options;
  valid_overlay_options << "overlay_reg" << "overlay_method" << "overlay_threshold"
                        << "overlay_rh" << "overlay_opacity" << "overlay_colormap";
  for (int nOverlay = 0; nOverlay < overlay_list.size(); nOverlay++)
  {
    QStringList sa_fn = overlay_list[nOverlay].split(":");
    if (nOverlay == 0)    // first one is not overlay file but actually surface file
      surface_fn = sa_fn[0];
    bool bLoadAll = false;
    bool bLabelOutline = false;
    QString labelColor;
    QString overlay_reg;
    QString overlay_opacity;
    QString overlay_method = "linearopaque";
    QStringList overlay_colormap;
    QStringList overlay_thresholds;
    bool bSecondHalfData = false;
    for ( int k = sa_fn.size()-1; k >= 0; k-- )
    {
      int n = sa_fn[k].indexOf( "=" );
      if ( n != -1  )
      {
        QString subOption = sa_fn[k].left( n ).toLower();
        QString subArgu = sa_fn[k].mid( n+1 );
        if ( subOption == "overlay_reg" )
          overlay_reg = subArgu;
        else if (subOption == "overlay_method")
          overlay_method = subArgu;
        else if (subOption == "overlay_threshold")
          overlay_thresholds = subArgu.split(",", QString::SkipEmptyParts);
        else if (subOption == "overlay_rh" && (subArgu == "1" || subArgu == "true"))
          bSecondHalfData = true;
        else if (subOption == "overlay_opacity")
          overlay_opacity = subArgu;
        else if (subOption == "overlay_colormap")
          overlay_colormap = subArgu.split(",", QString::SkipEmptyParts);
      }
    }
    if (overlay_reg.isEmpty())
      overlay_reg = "n/a";

  //  for ( int k = sa_fn.size()-1; k >= 0; k-- )
    for (int k = 0; k < sa_fn.size(); k++)
    {
      int n = sa_fn[k].indexOf( "=" );
      if ( n != -1  )
      {
        QString subOption = sa_fn[k].left( n ).toLower();
        QString subArgu = sa_fn[k].mid( n+1 );
        if ( subOption == "color" )
        {
          m_scripts.insert( 0, QStringList("setsurfacecolor") << subArgu );
        }
        else if ( subOption == "edgecolor" || subOption == "edge_color")
        {
          m_scripts.insert( 0, QStringList("setsurfaceedgecolor") << subArgu );
        }
        else if ( subOption == "edgethickness"|| subOption == "edge_thickness" )
        {
          m_scripts.insert( 0, QStringList("setsurfaceedgethickness") << subArgu );
        }
        else if ( subOption == "vertex" )
        {
          m_scripts.insert( 0, QStringList("displaysurfacevertex") << subArgu);
        }
        else if ( subOption == "vertexcolor" || subOption == "vertex_color" )
        {
          m_scripts.insert( 0, QStringList("setsurfacevertexcolor") << subArgu );
        }
        else if ( subOption == "curv" || subOption == "curvature" )
        {
          m_scripts.insert( 0, QStringList("loadsurfacecurvature") << subArgu );
        }
        else if ( subOption == "curvature_method" || subOption == "curvature_map")
        {
            m_scripts.insert(0, QStringList("setsurfacecurvaturemap") << subArgu);
        }
        else if ( subOption == "overlay" || subOption == "correlation" )
        {
          // add script to load surface overlay files
          QStringList script("loadsurfaceoverlay");
          script << subArgu;

          script << overlay_reg;
          if (subOption == "correlation")
            script << "correlation";
          else
            script << "n/a";

          if (bSecondHalfData)
              script << "rh";
          m_scripts.insert( 0, script );

          if (overlay_method != "linearopaque" || !overlay_thresholds.isEmpty())
          {
            script = QStringList("setsurfaceoverlaymethod") << overlay_method <<
                     overlay_thresholds;
            // insert right AFTER loadsurfaceoverlay command
            m_scripts.insert( 1, script );
          }

          if (!overlay_opacity.isEmpty())
            m_scripts.insert(1, QStringList("setsurfaceoverlayopacity") << overlay_opacity);

          if (!overlay_colormap.isEmpty())
            m_scripts.insert(1, QStringList("setsurfaceoverlaycolormap") << overlay_colormap);
        }
        else if ( subOption == "annot" || subOption == "annotation" )
        {
          // add script to load surface annotation files
          QStringList annot_fns =subArgu.split(",");
          for ( int i = annot_fns.size()-1; i >= 0; i-- )
          {
            m_scripts.insert( 0, QStringList("loadsurfaceannotation") << annot_fns[i] );
          }
        }
        else if ( subOption == "annot_outline" || subOption == "annotation_outline")
        {
          if ( subArgu.toLower() == "true" || subArgu.toLower() == "yes" || subArgu == "1")
          {
            for (int i = 0; i < m_scripts.size(); i++)
            {
              if (m_scripts[i][0] == "loadsurfaceannotation")
              {
                m_scripts.insert(i+1, QStringList("setsurfaceannotationoutline") << "1");
                break;
              }
            }
          }
        }
        else if ( subOption == "label" )
        {
          // add script to load surface label files
          QStringList fns = subArgu.split(",");
          for ( int i = fns.size()-1; i >= 0; i-- )
          {
            m_scripts.insert(0, QStringList("loadsurfacelabel") << fns[i]);
          }
        }
        else if ( subOption == "label_outline" || subOption == "labeloutline")
        {
          if ( subArgu.toLower() == "true" || subArgu.toLower() == "yes" || subArgu == "1")
          {
            for (int i = 0; i < m_scripts.size(); i++)
            {
              if (m_scripts[i][0] == "loadsurfacelabel")
              {
                m_scripts.insert(i+1, QStringList("setsurfacelabeloutline") << "1");
                break;
              }
            }
          }
        }
        else if (subOption == "label_color" || subOption == "labelcolor")
        {
          if (!subArgu.isEmpty())
          {
            for (int i = 0; i < m_scripts.size(); i++)
            {
              if (m_scripts[i][0] == "loadsurfacelabel")
              {
                m_scripts.insert(i+1, QStringList("setsurfacelabelcolor") << subArgu);
                break;
              }
            }
          }
        }
        else if ( subOption == "vector" )
        {
          // add script to load surface vector files
          QStringList vector_fns = subArgu.split(",");
          for ( int i = vector_fns.size() - 1 ; i >= 0 ; i-- )
          {
            m_scripts.insert(0, QStringList("loadsurfacevector") << vector_fns[i]);
          }
        }
        else if (subOption == "spline")
        {
          m_scripts.insert(0, QStringList("loadsurfacespline") << subArgu);
        }
        else if ( subOption == "patch" )
        {
          if ( subArgu.contains( "/" ) )
          {
            subArgu = QFileInfo( subArgu ).absoluteFilePath();
          }
          fn_patch = subArgu;
        }
        else if ( subOption == "target" || subOption == "target_surf")
        {
          if ( subArgu.contains( "/" ) )
          {
            subArgu = QFileInfo( subArgu ).absoluteFilePath();
          }
          fn_target = subArgu;
        }
        else if ( subOption == "name" )
        {
          m_scripts.insert( 0, QStringList("setlayername") << "Surface" << subArgu );
        }
        else if ( subOption == "lock" || subOption == "locked")
        {
          m_scripts.insert( 0, QStringList("locklayer") << "Surface" << subArgu );
        }
        else if ( subOption == "visible" )
        {
          m_scripts.insert( 0, QStringList("showlayer") << "Surface" << subArgu );
        }
        else if ( subOption == "offset" )
        {
          m_scripts.insert( 0, QStringList("setsurfaceoffset") << subArgu.split(",") );
        }
        else if ( subOption == "all")
        {
          if ( subArgu.toLower() == "true" || subArgu.toLower() == "yes" || subArgu == "1")
            bLoadAll = true;
        }
        else if (subOption == "sup_files")
        {
          sup_files = subArgu.split(",",  QString::SkipEmptyParts);
        }
        else if ( !valid_overlay_options.contains(subOption) )
        {
          cerr << "Unrecognized sub-option flag '" << subOption.toAscii().constData() << "'.\n";
          return;
        }
      }
    }
  }
  LoadSurfaceFile( surface_fn, fn_patch, fn_target, sup_files );
}

/*
void MainWindow::CommandLoadSurface( const QStringList& cmd )
{
  QString fullfn = cmd[1];
  int nIgnoreStart = fullfn.indexOf( "#" );
  int nIgnoreEnd = fullfn.indexOf( "#", nIgnoreStart+1 );
  QStringList sa_fn = MyUtils::SplitString( fullfn, ":", nIgnoreStart, nIgnoreEnd - nIgnoreStart + 1 );
  QString fn = sa_fn[0];
  QString fn_patch = "";
  QString fn_target = "";
  bool bLoadAll = false;
  QStringList sup_files;
  bool bLabelOutline = false;
  QString labelColor;
  QString overlay_reg;
  QString overlay_method = "linearopaque";
  QStringList overlay_thresholds;
  bool bSecondHalfData = false;
  for ( int k = sa_fn.size()-1; k >= 1; k-- )
  {
    int n = sa_fn[k].indexOf( "=" );
    if ( n != -1  )
    {
      QString subOption = sa_fn[k].left( n ).toLower();
      QString subArgu = sa_fn[k].mid( n+1 );
      if ( subOption == "overlay_reg" )
        overlay_reg = subArgu;
      else if (subOption == "overlay_method")
        overlay_method = subArgu;
      else if (subOption == "overlay_threshold")
        overlay_thresholds = subArgu.split(",", QString::SkipEmptyParts);
      else if (subOption == "overlay_rh" && (subArgu == "1" || subArgu == "true"))
        bSecondHalfData = true;
    }
  }
  if (overlay_reg.isEmpty())
    overlay_reg = "n/a";


  for ( int k = sa_fn.size()-1; k >= 1; k-- )
  {
    int n = sa_fn[k].indexOf( "=" );
    if ( n != -1  )
    {
      QString subOption = sa_fn[k].left( n ).toLower();
      QString subArgu = sa_fn[k].mid( n+1 );
      if ( subOption == "color" )
      {
        m_scripts.insert( 0, QString("setsurfacecolor ") + subArgu );
      }
      else if ( subOption == "edgecolor" || subOption == "edge_color")
      {
        m_scripts.insert( 0, QString("setsurfaceedgecolor ") + subArgu );
      }
      else if ( subOption == "edgethickness"|| subOption == "edge_thickness" )
      {
        m_scripts.insert( 0, QString("setsurfaceedgethickness ") + subArgu );
      }
      else if ( subOption == "vertex" )
      {
        m_scripts.insert( 0, QString("displaysurfacevertex ") + subArgu);
      }
      else if ( subOption == "vertexcolor" || subOption == "vertex_color" )
      {
        m_scripts.insert( 0, QString("setsurfacevertexcolor ") + subArgu );
      }
      else if ( subOption == "curv" || subOption == "curvature" )
      {
        m_scripts.insert( 0, QString("loadsurfacecurvature ") + subArgu );
      }
      else if ( subOption == "overlay" || subOption == "correlation" )
      {
        // add script to load surface overlay files
        nIgnoreStart = subArgu.indexOf( "#" );
        nIgnoreEnd = subArgu.indexOf( "#", nIgnoreStart+1 );
        QStringList overlay_fns = MyUtils::SplitString( subArgu, ",", nIgnoreStart, nIgnoreEnd - nIgnoreStart + 1 );
        for ( int i = overlay_fns.size() - 1 ; i >= 0 ; i-- )
        {
          QString script = "loadsurfaceoverlay ";
          int nSubStart = overlay_fns[i].indexOf( "#" );
          int nSubEnd = overlay_fns[i].indexOf( "#", nSubStart+1 );
          if ( nSubEnd == -1 )
          {
            nSubEnd = overlay_fns[i].length() - 1;
          }
          if ( nSubStart != -1 )
          {
            script += overlay_fns[i].left( nSubStart );
          }
          else
          {
            script += overlay_fns[i];
          }

          script += " " + overlay_reg;
          if (subOption == "correlation")
            script += QString(" correlation");
          else
            script += QString(" n/a");

          if (bSecondHalfData)
            script += " rh";
          m_scripts.insert( 0, script );

          // if there are sub-options attached with overlay file, parse them
          QString opt_strg;
          if ( nSubStart != -1 )
          {
            opt_strg = overlay_fns[i].mid( nSubStart+1, nSubEnd - nSubStart - 1 );
          }

          QStringList overlay_opts = opt_strg.split(":");

          QString method;
          QStringList thresholds;
          for ( int j = 0; j < overlay_opts.size(); j++ )
          {
            QString strg = overlay_opts[j];
            if ( ( n = strg.indexOf( "=" ) ) != -1 && strg.left( n ).toLower() == "method" )
            {
              method = strg.mid( n+1 ).toLower();
            }
            else if ( ( n = strg.indexOf( "=" ) ) != -1 && strg.left( n ).toLower() == "threshold" )
            {
              thresholds = strg.mid( n+1 ).split(",");
            }
          }

          if ( !method.isEmpty() || !thresholds.isEmpty() )
          {
            if (method.isEmpty())
              method = overlay_method;
            script = QString("setsurfaceoverlaymethod ") + method + " " +
                     thresholds.join(" ");
            // insert right AFTER loadsurfaceoverlay command
            m_scripts.insert( 1, script );
          }
          else if (overlay_method != "linearopaque" || !overlay_thresholds.isEmpty())
          {
            script = QString("setsurfaceoverlaymethod ") + overlay_method + " " +
                     overlay_thresholds.join(" ");
            // insert right AFTER loadsurfaceoverlay command
            m_scripts.insert( 1, script );
          }
        }
      }
      else if ( subOption == "annot" || subOption == "annotation" )
      {
        // add script to load surface annotation files
        QStringList annot_fns =subArgu.split(",");
        for ( int i = annot_fns.size()-1; i >= 0; i-- )
        {
          m_scripts.insert( 0, QString("loadsurfaceannotation ") + annot_fns[i] );
        }
      }
      else if ( subOption == "label" )
      {
        // add script to load surface label files
        QStringList fns = subArgu.split(",");
        for ( int i = fns.size()-1; i >= 0; i-- )
        {
          m_scripts.insert(0, QString("loadsurfacelabel ")+fns[i]);
        }
      }
      else if ( subOption == "vector" )
      {
        // add script to load surface vector files
        QStringList vector_fns = subArgu.split(",");
        for ( int i = vector_fns.size() - 1 ; i >= 0 ; i-- )
        {
          m_scripts.insert(0, QString("loadsurfacevector ")+vector_fns[i]);
        }
      }
      else if (subOption == "spline")
      {
        m_scripts.insert(0, QString("loadsurfacespline ")+subArgu);
      }
      else if ( subOption == "patch" )
      {
        if ( subArgu.contains( "/" ) )
        {
          subArgu = QFileInfo( subArgu ).absoluteFilePath();
        }
        fn_patch = subArgu;
      }
      else if ( subOption == "target" || subOption == "target_surf")
      {
        if ( subArgu.contains( "/" ) )
        {
          subArgu = QFileInfo( subArgu ).absoluteFilePath();
        }
        fn_target = subArgu;
      }
      else if ( subOption == "name" )
      {
        m_scripts.insert( 0, QString("setlayername Surface ")+subArgu );
      }
      else if ( subOption == "lock" )
      {
        m_scripts.insert( 0, QString("locklayer Surface ")+subArgu );
      }
      else if ( subOption == "visible" )
      {
        m_scripts.insert( 0, QString("showlayer Surface ")+subArgu );
      }
      else if ( subOption == "offset" )
      {
        QString script = "setsurfaceoffset ";
        script += subArgu.replace(",", " ");
        m_scripts.insert( 0, script );
      }
      else if ( subOption == "all")
      {
        if ( subArgu.toLower() == "true" || subArgu.toLower() == "yes" || subArgu == "1")
          bLoadAll = true;
      }
      else if ( subOption == "label_outline" || subOption == "labeloutline")
      {
        if ( subArgu.toLower() == "true" || subArgu.toLower() == "yes" || subArgu == "1")
          bLabelOutline = true;
      }
      else if (subOption == "label_color" || subOption == "labelcolor")
      {
        labelColor = subArgu;
      }
      else if (subOption == "sup_files")
      {
        sup_files = subArgu.split(",",  QString::SkipEmptyParts);
      }
      else if (subOption != "overlay_reg" && subOption != "overlay_method" && subOption != "overlay_threshold" &&
               subOption != "overlay_rh")
      {
        cerr << "Unrecognized sub-option flag '" << subOption.toAscii().constData() << "'.\n";
        return;
      }
    }
  }
  if (bLabelOutline)
  {
    for (int i = 0; i < m_scripts.size(); i++)
    {
      if (m_scripts[i].indexOf("loadsurfacelabel") == 0)
        m_scripts.insert(i+1, "setsurfacelabeloutline 1");
    }
  }
  if (!labelColor.isEmpty())
  {
    for (int i = 0; i < m_scripts.size(); i++)
    {
      if (m_scripts[i].indexOf("loadsurfacelabel") == 0)
        m_scripts.insert(i+1, QString("setsurfacelabelcolor ") + labelColor);
    }
  }
  LoadSurfaceFile( fn, fn_patch, fn_target, sup_files );
}
*/

void MainWindow::CommandSetSurfaceLabelOutline(const QStringList &cmd)
{
  LayerSurface* surf = (LayerSurface*)GetLayerCollection( "Surface" )->GetActiveLayer();
  if ( surf )
  {
    if (cmd[1] == "1")
    {
      surf->SetActiveLabelOutline(true);
    }
  }
}

void MainWindow::CommandSetSurfaceAnnotationOutline(const QStringList &cmd)
{
  LayerSurface* surf = (LayerSurface*)GetLayerCollection( "Surface" )->GetActiveLayer();
  if ( surf )
  {
    if (cmd[1] == "1")
    {
      surf->SetActiveAnnotationOutline(true);
    }
  }
}

void MainWindow::CommandSetSurfaceOverlayOpacity(const QStringList &cmd)
{
  LayerSurface* surf = (LayerSurface*)GetLayerCollection( "Surface" )->GetActiveLayer();
  if ( surf )
  {
    SurfaceOverlay* overlay = surf->GetActiveOverlay();
    if ( overlay )
    {
      bool ok;
      double opacity = cmd[1].toDouble(&ok);
      if (ok)
      {
        overlay->GetProperty()->SetOpacity(opacity);
        surf->UpdateOverlay(true);
        overlay->EmitDataUpdated();
      }
      else
      {
        cerr << "Invalid input for overlay opacity.\n";
      }
    }
  }
}

void MainWindow::CommandSetSurfaceOverlayMethod( const QStringList& cmd_in )
{
  QStringList cmd = cmd_in;
  LayerSurface* surf = (LayerSurface*)GetLayerCollection( "Surface" )->GetActiveLayer();
  if ( surf )
  {
    SurfaceOverlay* overlay = surf->GetActiveOverlay();
    if ( overlay )
    {
      int nMethod = SurfaceOverlayProperty::CM_LinearOpaque;
      if ( cmd[1] == "linear" )
      {
        nMethod = SurfaceOverlayProperty::CM_Linear;
      }
      else if ( cmd[1] == "piecewise" )
      {
        nMethod = SurfaceOverlayProperty::CM_Piecewise;
      }
      else if ( cmd[1] != "linearopaque" )
      {
        cerr << "Unrecognized overlay method name '" << cmd[1].toAscii().constData() << "'.\n";
        return;
      }

      overlay->GetProperty()->SetColorMethod( nMethod );

      bool bPercentile = false;
      if (cmd.last() == "percentile")
      {
        cmd.removeLast();
        bPercentile = true;
      }
      double values[3];
      if ( cmd.size() - 2 >= 3 )   // 3 values
      {
        bool bOK;
        values[0] = cmd[2].toDouble(&bOK);
        values[1] = cmd[3].toDouble(&bOK);
        values[2] = cmd[4].toDouble(&bOK);
        if (bPercentile)
        {
          overlay->GetProperty()->SetUsePercentile(bPercentile);
          for (int i = 0; i < 3; i++)
            values[i] = overlay->PercentileToPosition(values[i]);
        }
        if ( bOK )
        {
          overlay->GetProperty()->SetMinPoint( values[0] );
          overlay->GetProperty()->SetMidPoint( values[1] );
          overlay->GetProperty()->SetMaxPoint( values[2] );
        }
        else
        {
          cerr << "Invalid input for overlay threshold.\n";
        }
      }
      else if ( cmd.size() - 2 == 2 )   // 2 values
      {
        bool bOK;
        values[0] = cmd[2].toDouble(&bOK);
        values[1] = cmd[3].toDouble(&bOK);
        if ( bOK )
        {
          if (bPercentile)
          {
            overlay->GetProperty()->SetUsePercentile(bPercentile);
            for (int i = 0; i < 2; i++)
              values[i] = overlay->PercentileToPosition(values[i]);
          }
          overlay->GetProperty()->SetMinPoint( values[0] );
          overlay->GetProperty()->SetMaxPoint( values[1] );
          overlay->GetProperty()->SetMidPoint( ( values[0] + values[1] ) / 2 );
        }
        else
        {
          cerr << "Invalid input for overlay threshold.\n";
        }
      }
      surf->UpdateOverlay(true);
      overlay->EmitDataUpdated();
    }
  }
}

void MainWindow::CommandSetSurfaceOverlayColormap(const QStringList &cmd)
{
  LayerSurface* surf = (LayerSurface*)GetLayerCollection( "Surface" )->GetActiveLayer();
  if ( surf )
  {
    SurfaceOverlay* overlay = surf->GetActiveOverlay();
    if ( overlay )
    {
      if (cmd[1] == "colorwheel")
        overlay->GetProperty()->SetColorScale(SurfaceOverlayProperty::CS_ColorWheel);
      if (cmd.size() > 2)
      {
        if (cmd[2] == "inverse")
          overlay->GetProperty()->SetColorInverse(true);
        else if (cmd[2] == "truncate")
          overlay->GetProperty()->SetColorTruncate(true);
      }
      surf->UpdateOverlay(true);
    }
  }
}

QColor MainWindow::ParseColorInput(const QString &strg)
{
  QColor color;
  if (!strg.contains(','))
  {
    if (QColor::isValidColor(strg))
      color = QColor(strg);
    else if (strg.indexOf("light", 0, Qt::CaseInsensitive) == 0)
    {
      QString name = strg.mid(5);
      if (QColor::isValidColor(name))
        color = QColor(name).lighter();
    }
    else if (strg.indexOf("dark", 0, Qt::CaseInsensitive) == 0)
    {
      QString name = strg.mid(4);
      if (QColor::isValidColor(name))
        color = QColor(name).darker();
    }
  }
  if ( !color.isValid() )
  {
    int rgb[3];
    QStringList rgb_strs = strg.split(",");
    rgb_strs << "n/a" << "n/a";
    bool bOK;
    rgb[0] = rgb_strs[0].toInt(&bOK);
    if (bOK)
    {
      rgb[1] = rgb_strs[1].toInt(&bOK);
    }
    if (bOK)
    {
      rgb[2] = rgb_strs[2].toInt(&bOK);
    }
    if ( bOK )
    {
      color.setRgb(rgb[0], rgb[1], rgb[2]);
    }
  }
  return color;
}

void MainWindow::CommandSetSurfaceColor( const QStringList& cmd )
{
  LayerSurface* surf = (LayerSurface*)GetLayerCollection( "Surface" )->GetActiveLayer();
  if ( surf && cmd[1] != "null" )
  {
    QColor color = ParseColorInput( cmd[1] );
    if ( color.isValid() )
    {
      surf->GetProperty()->SetBinaryColor( color.redF(), color.greenF(), color.blueF() );
    }
    else
    {
      cerr << "Invalid color name or value " << cmd[1].toAscii().constData() << ".\n";
    }
  }
}

void MainWindow::CommandSetSurfaceEdgeColor( const QStringList& cmd )
{
  LayerSurface* surf = (LayerSurface*)GetLayerCollection( "Surface" )->GetActiveLayer();
  if ( surf && cmd[1] != "null" )
  {
    QColor color = ParseColorInput( cmd[1] );
    if ( color.isValid() )
    {
      surf->GetProperty()->SetEdgeColor( color.redF(), color.greenF(), color.blueF() );
    }
    else
    {
      cerr << "Invalid color name or value " << cmd[1].toAscii().constData() << ".\n";
    }
  }
}

void MainWindow::CommandSetSurfaceEdgeThickness( const QStringList& cmd )
{
  LayerSurface* surf = (LayerSurface*)GetLayerCollection( "Surface" )->GetActiveLayer();
  if ( surf )
  {
    bool bOK;
    int thickness = cmd[1].toInt(&bOK);
    if ( !bOK )
    {
      cerr << "Invalid edge thickness value. Must be a integer.\n";
    }
    else
    {
      surf->GetProperty()->SetEdgeThickness( thickness );
    }
  }
}

void MainWindow::CommandSetDisplaySurfaceVertex(const QStringList &cmd)
{
  LayerSurface* surf = (LayerSurface*)GetLayerCollection( "Surface" )->GetActiveLayer();
  if ( surf && (cmd[1].toLower() == "on" ||
                cmd[1].toLower() == "true" ||
                cmd[1].toLower() == "yes" ||
                cmd[1] == "1") )
  {
    surf->GetProperty()->ShowVertices(true);
  }
}

void MainWindow::CommandSetSurfaceVertexColor(const QStringList &cmd)
{
  LayerSurface* surf = (LayerSurface*)GetLayerCollection( "Surface" )->GetActiveLayer();
  if ( surf && cmd[1] != "null" )
  {
    QColor color = ParseColorInput( cmd[1] );
    if ( color.isValid() )
    {
      surf->GetProperty()->SetVertexColor( color.redF(), color.greenF(), color.blueF() );
    }
    else
    {
      cerr << "Invalid color name or value " << cmd[1].toAscii().constData() << ".\n";
    }
  }
}

void MainWindow::CommandSetSurfaceLabelColor(const QStringList &cmd)
{
  LayerSurface* surf = (LayerSurface*)GetLayerCollection( "Surface" )->GetActiveLayer();
  if ( surf && cmd[1] != "null" )
  {
    QColor color = ParseColorInput( cmd[1] );
    if ( color.isValid() )
    {
      surf->SetActiveLabelColor(color);
    }
    else
    {
      cerr << "Invalid color name or value " << cmd[1].toAscii().constData() << ".\n";
    }
  }
}

void MainWindow::CommandSetSurfaceOffset( const QStringList& cmd )
{
  LayerSurface* surf = (LayerSurface*)GetLayerCollection( "Surface" )->GetActiveLayer();
  if ( surf )
  {
    bool bOK;
    double pos[3];
    if (cmd.size() < 4 )
    {
      cerr << "Invalid surface offset inputs. Need 3 numbers.\n";
    }
    else
    {
      pos[0] = cmd[1].toDouble(&bOK);
      pos[1] = cmd[2].toDouble(&bOK);
      pos[2] = cmd[3].toDouble(&bOK);
      if ( !bOK )
      {
        cerr << "Invalid surface offset inputs. Need 3 numbers.\n";
      }
      else
      {
        surf->GetProperty()->SetPosition( pos );
      }
    }
  }
}

void MainWindow::CommandGotoLabel(const QStringList &cmd)
{
  LayerMRI* mri = (LayerMRI*)GetActiveLayer("MRI");
  if ( mri )
  {
    int nSlice = mri->GetGotoLabelSlice();
    if (nSlice >= 0)
    {
      double pos[3];
      int n[3];
      mri->GetSlicePosition(pos);
      mri->TargetToRAS(pos, pos);
      mri->RASToOriginalIndex(pos, n);
      QString ostr = mri->GetOrientationString();
      int nView = this->GetMainViewId();
      if (nView > 2)
        nView = 0;
      int nOrigPlane = nView;
      char ch[3][3] = {"RL", "AP", "IS"};
      for (int i = 0; i < 3; i++)
      {
        if (ostr[i] == ch[nView][0] || ostr[i] == ch[nView][1])
        {
          nOrigPlane = i;
          break;
        }
      }
      n[nOrigPlane] = nSlice;
      mri->OriginalIndexToRAS(n, pos);
      mri->RASToTarget(pos, pos);
      GetLayerCollection("MRI")->SetCursorRASPosition( pos );
      SetSlicePosition(pos);
    }
    else
    {
      cerr << "Did not find slice with most structure of " << qPrintable(cmd[1]);
    }
  }
}

void MainWindow::CommandLoadSurfaceVector( const QStringList& cmd )
{
  LoadSurfaceVectorFile( cmd[1] );
}

void MainWindow::CommandLoadSurfaceCurvature( const QStringList& cmd )
{
  LoadSurfaceCurvatureFile( cmd[1] );
}

void MainWindow::CommandSetSurfaceCurvatureMap(const QStringList &cmd)
{
    LayerSurface* layer = qobject_cast<LayerSurface*>(GetActiveLayer("Surface"));
    if ( layer )
    {
        int nMap = LayerPropertySurface::CM_Threshold;
        if (cmd[1].toLower() == "off")
            nMap = LayerPropertySurface::CM_Off;
        else if (cmd[1].toLower() == "binary")
            nMap = LayerPropertySurface::CM_Binary;
        layer->GetProperty()->SetCurvatureMap(nMap);
    }
}

void MainWindow::CommandLoadSurfaceOverlay( const QStringList& cmd )
{
  QString reg_file = cmd[2];
  if (reg_file == "n/a")
    reg_file = "";
  LoadSurfaceOverlayFile( cmd[1], reg_file, cmd.size() > 3 && cmd[3] == "correlation", cmd.size() > 4 && cmd[4] == "rh" );
}

void MainWindow::CommandLoadSurfaceAnnotation( const QStringList& cmd )
{
  LoadSurfaceAnnotationFile( cmd[1] );
}

void MainWindow::CommandLoadSurfaceLabel( const QStringList& cmd )
{
  if (!LoadSurfaceLabelFile( cmd[1] ))
  {
    if (!m_scripts.isEmpty())
    {
      if (m_scripts[0].at(0).toLower() == "loadsurfacelabel")
        m_scripts.removeAt(0);
    }
  }
}

void MainWindow::CommandLoadSurfaceSpline(const QStringList &cmd)
{
  LoadSurfaceSplineFile( cmd[1]);
}

void MainWindow::CommandLoadWayPoints( const QStringList& cmd )
{
  QStringList options = cmd[1].split(":");
  QString fn = options[0];
  QString color = "null";
  QString spline_color = "null";
  QString radius = "0";
  QString spline_radius = "0";
  QString spline_heatmap;
  for ( int i = 1; i < options.size(); i++ )
  {
    QString strg = options[i];
    int n = strg.indexOf( "=" );
    if ( n != -1 )
    {
      QString option = strg.left( n ).toLower();
      QString argu = strg.mid( n+1 );
      if ( option == "color" )
      {
        color = argu;
      }
      else if ( option == "splinecolor" )
      {
        spline_color = argu;
      }
      else if ( option == "radius" )
      {
        radius = argu;
      }
      else if ( option == "splineradius" )
      {
        spline_radius = argu;
      }
      else if (option == "splineheatmap" || option == "heatmap")
      {
        spline_heatmap = argu;
      }
      else if ( option == "name" )
      {
        m_scripts.insert( 0, QStringList("setlayername") << "PointSet" << argu );
      }
      else if ( option == "visible" )
      {
        m_scripts.insert( 0, QStringList("showlayer") << "PointSet" << argu );
      }
      else
      {
        cerr << "Unrecognized sub-option flag '" << strg.toAscii().constData() << "'.\n";
      }
    }
  }

  if ( color != "null" || spline_color != "null" )
  {
    m_scripts.insert( 0, QStringList("setpointsetcolor") << color << spline_color );
  }

  if ( radius != "0" || spline_radius != "0" )
  {
    m_scripts.insert( 0, QStringList("setpointsetradius") << radius << spline_radius );
  }
  if (!spline_heatmap.isEmpty())
  {
    m_scripts.insert( 0, QStringList("setpointsetheatmap") << spline_heatmap.split(",", QString::SkipEmptyParts));
  }

  LoadWayPointsFile( fn );
}

void MainWindow::CommandLoadControlPoints( const QStringList& cmd )
{
  QStringList options = cmd[1].split(":");
  QString fn = options[0];
  QString color = "null";
  QString radius = "0";
  for ( int i = 1; i < options.size(); i++ )
  {
    QString strg = options[i];
    int n = strg.indexOf( "=" );
    if ( n != -1 )
    {
      QString option = strg.left( n ).toLower();
      QString argu = strg.mid( n+1 );
      if ( option == "color" )
      {
        color = argu;
      }
      else if ( option == "radius" )
      {
        radius = argu;
      }
      else if ( option == "name" )
      {
        m_scripts.insert( 0, QStringList("setlayername") << "PointSet" << argu );
      }
      else if ( option == "visible" )
      {
        m_scripts.insert( 0, QStringList("showlayer") << "PointSet" << argu );
      }
      else
      {
        cerr << "Unrecognized sub-option flag '" << strg.toAscii().constData() << "'.\n";
      }
    }
  }

  if ( color != "null" )
  {
    m_scripts.insert( 0, QStringList("setpointsetcolor") << color );
  }

  if ( radius != "0" )
  {
    m_scripts.insert( 0, QStringList("setpointsetradius") << radius);
  }
  LoadControlPointsFile( fn );
}

void MainWindow::CommandSetPointSetColor( const QStringList& cmd )
{
  LayerPointSet* wp = (LayerPointSet*)GetLayerCollection( "PointSet" )->GetActiveLayer();
  if ( wp )
  {
    if ( cmd[1] != "null" )
    {
      QColor color = ParseColorInput( cmd[1] );

      if ( color.isValid() )
      {
        wp->GetProperty()->SetColor( color.redF(), color.greenF(), color.blueF() );
      }
      else
      {
        cerr << "Invalid color name or value " << cmd[1].toAscii().constData() << ".\n";
      }
    }

    if ( cmd.size() > 2 && cmd[2] != "null" )
    {
      QColor color = ParseColorInput( cmd[2] );
      if ( color.isValid() )
      {
        wp->GetProperty()->SetSplineColor( color.redF(), color.greenF(), color.blueF() );
      }
      else
      {
        cerr << "Invalid color name or value " << cmd[2].toAscii().constData() << ".\n";
      }
    }
  }
}

void MainWindow::CommandSetPointSetRadius( const QStringList& cmd )
{
  LayerPointSet* wp = (LayerPointSet*)GetLayerCollection( "PointSet" )->GetActiveLayer();
  if ( wp )
  {
    if ( cmd[1] != "0" )
    {
      bool bOK;
      double dvalue = cmd[1].toDouble(&bOK);
      if ( bOK)
      {
        wp->GetProperty()->SetRadius( dvalue );
      }
      else
      {
        cerr << "Invalid way points radius.\n";
      }
    }

    if ( cmd.size() > 2 && cmd[2] != "0" )
    {
      bool bOK;
      double dvalue = cmd[2].toDouble(&bOK);
      if ( bOK )
      {
        wp->GetProperty()->SetSplineRadius( dvalue );
      }
      else
      {
        cerr << "Invalid spline radius.\n";
      }
    }
  }
}

void MainWindow::CommandSetPointSetHeatmap(const QStringList &cmd)
{
  LayerPointSet* wp = (LayerPointSet*)GetLayerCollection( "PointSet" )->GetActiveLayer();
  if ( wp )
  {
    if (wp->GetProperty()->LoadScalarsFromFile(cmd[1]))
    {
      if (cmd.size() >= 3)
      {
        wp->GetProperty()->SetHeatScaleMin(cmd[2].toDouble());
      }
      if (cmd.size() == 4)
      {
        wp->GetProperty()->SetHeatScaleMax(cmd[3].toDouble());
        wp->GetProperty()->SetHeatScaleMid((cmd[2].toDouble()+cmd[3].toDouble())/2.0);
      }
      else if (cmd.size() > 4)
      {
        wp->GetProperty()->SetHeatScaleMid(cmd[3].toDouble());
        wp->GetProperty()->SetHeatScaleMax(cmd[4].toDouble());
      }
    }
    else
    {
      cerr << "Could not load scalar map from file " << qPrintable(cmd[1]) << "\n";
    }
  }
}

void MainWindow::CommandScreenCapture( const QStringList& cmd )
{
  double mag_factor = 1.0;
  bool bOK;
  mag_factor = cmd[2].toDouble(&bOK);
  if (bOK && mag_factor < 1)
    mag_factor = 1;

  if (!m_views[m_nMainView]->SaveScreenShot( cmd[1],
      m_settingsScreenshot.AntiAliasing,
      (int)mag_factor ))
  {
    cerr << "Failed to save screen shot to " << cmd[1].toAscii().constData() << ".\n";
  }
}

void MainWindow::CommandFlyThrough(const QStringList &cmd)
{
    if (GetMainViewId() > 2)
    {
      cerr << "Can not fly through. Please set main viewport to 2D slice view.\n";
      return;
    }
}

void MainWindow::CommandSetViewport( const QStringList& cmd )
{
  if ( cmd[1].toLower() == "x" )
  {
    SetMainView( MV_Sagittal );
  }
  else if ( cmd[1].toLower() == "y" )
  {
    SetMainView( MV_Coronal );
  }
  else if ( cmd[1].toLower() == "z" )
  {
    SetMainView( MV_Axial );
  }
  else if ( cmd[1].toLower() == "3d" )
  {
    SetMainView( MV_3D );
  }
}

void MainWindow::CommandSetViewSize( const QStringList& cmd )
{
  bool bOK;
  int x = cmd[1].toInt(&bOK);
  int y = cmd[2].toInt(&bOK);
  if ( !bOK )
  {
    cerr << "Invalid view size.\n";
    return;
  }

  QSize sz = m_views[m_nMainView]->size();
  int offsetx = x - sz.width(), offsety = y - sz.height();
  switch( m_nViewLayout )
  {
  case VL_2x2:
    offsetx *= 2;
    offsety *= 2;
    break;
  case VL_1n3:
    offsety += (int)(offsety/2.0+0.5);
    break;
  case VL_1n3h:
    offsetx += (int)(offsetx/2.0+0.5);
    break;
  default:
    break;
  }
  this->resize( this->size() + QSize(offsetx,offsety) );

  // now fine adjust the size again.
  QSize soffset = QSize(x,y)-m_views[m_nMainView]->size();
  this->resize( this->size() + soffset);
}

void MainWindow::CommandZoom( const QStringList& cmd )
{
  bool bOK;
  double dValue = cmd[1].toDouble(&bOK);
  if ( bOK && m_nMainView >= 0 )
  {
    m_views[m_nMainView]->Zoom( dValue );
  }
}

void MainWindow::CommandSetCamera(const QStringList &cmd)
{
  if (cmd[1].toLower() == "load")
  {
    OnToolLoadCamera(cmd[2]);
  }
  else
  {
    bool bOK;
    CameraOperations ops;
    for (int i = 1; i < cmd.size(); i+=2)
    {
      double dValue = cmd[i+1].toDouble(&bOK);
      if (!bOK)
      {
        cerr << "Invalid input value for " << cmd[i].toAscii().constData() << ".\n";
        return;
      }
      else
      {
        ops << CameraOperation(cmd[i], dValue);
      }
    }
    m_views[3]->SetCameraOperations(ops);
    m_views[3]->ResetCameraClippingRange();
  }
  m_views[3]->RequestRedraw();
}

void MainWindow::CommandSetRAS( const QStringList& cmd )
{
  bool bOK;
  double ras[3];
  ras[0] = cmd[1].toDouble(&bOK);
  if (bOK)
  {
    ras[1] = cmd[2].toDouble(&bOK);
  }
  if (bOK)
  {
    ras[2] = cmd[3].toDouble(&bOK);
  }
  if ( bOK )
  {
    LayerCollection* lc = GetLayerCollection( "MRI" );
    LayerMRI* layer = (LayerMRI*)lc->GetLayer( 0 );
    if ( layer )
    {
      if ( cmd.size() > 4 && cmd[4] == "tkreg" )
      {
        layer->TkRegToNativeRAS( ras, ras );
      }
      layer->RASToTarget( ras, ras );
    }
    this->GetMainView()->CenterAtWorldPosition(ras);
    lc->SetCursorRASPosition( ras );
    SetSlicePosition( ras );
  }
  else
  {
    cerr << "Invalid input values for RAS coordinates. \n";
  }
}


void MainWindow::CommandSetSlice( const QStringList& cmd )
{
  LayerCollection* lc_mri = GetLayerCollection( "MRI" );
  if ( !lc_mri->IsEmpty() )
  {
    LayerMRI* mri = (LayerMRI*)lc_mri->GetLayer( lc_mri->GetNumberOfLayers()-1 );
    int x, y, z;
    bool bOK;
    x = cmd[1].toInt(&bOK);
    y = cmd[2].toInt(&bOK);
    z = cmd[3].toInt(&bOK);
    if ( bOK )
    {
      int slice[3] = { x, y, z };
      double ras[3];
      mri->OriginalIndexToRAS( slice, ras );
      mri->RASToTarget( ras, ras );

      lc_mri->SetCursorRASPosition( ras );
      SetSlicePosition( ras );
    }
    else
    {
      cerr << "Invalide slice number(s). \n";
    }
  }
  else
  {
    cerr << "No volume was loaded. Set slice failed.\n";
  }
}

void MainWindow::SetCurrentFile( const QString &fileName, int type )
{
  QString key = "MainWindow/RecentVolumeFiles";
  if ( type == 1 )
  {
    key = "MainWindow/RecentSurfaceFiles";
  }
  QSettings settings;
  QStringList files = settings.value( key ).toStringList();
  files.removeAll(fileName);
  files.prepend(fileName);
  while ( files.size() > MAX_RECENT_FILES )
  {
    files.removeLast();
  }

  settings.setValue( key, files);

  foreach (QWidget *widget, QApplication::topLevelWidgets())
  {
    MainWindow *mainWin = qobject_cast<MainWindow *>(widget);
    if (mainWin)
    {
      mainWin->UpdateRecentFileActions();
    }
  }
}

void MainWindow::SetAction( int nAction )
{
  for ( int i = 0; i < 4; i++ )
  {
    m_views[i]->releaseMouse();
    m_views[i]->SetAction( nAction );
  }

  if ( (m_views[0]->GetInteractionMode() == RenderView::IM_VoxelEdit ||
        m_views[0]->GetInteractionMode() == RenderView::IM_ReconEdit )
      && nAction == Interactor::EM_Contour )
  {
    BrushProperty* bp =GetBrushProperty();
    LayerMRI* layer = (LayerMRI*)GetActiveLayer( "MRI" );
    if ( layer && layer->GetProperty()->GetColorMap() == LayerPropertyMRI::LUT )
    {
      layer->GetProperty()->SetShowLabelOutline( true );
    }

    if ( !bp->GetReferenceLayer() )
    {
      LayerCollection* lc_mri = GetLayerCollection( "MRI" );
      for ( int i = 0; i < lc_mri->GetNumberOfLayers(); i++ )
      {
        LayerMRI* mri = (LayerMRI*)lc_mri->GetLayer( i );
        if ( mri->GetProperty()->GetColorMap() != LayerPropertyMRI::LUT && mri->IsVisible() && mri != layer )
        {
          bp->SetReferenceLayer( mri );
          break;
        }
      }
    }
  }
}

void MainWindow::UpdateRecentFileActions()
{
  QSettings settings;
  QStringList files = settings.value("MainWindow/RecentVolumeFiles").toStringList();

  int numRecentFiles = qMin(files.size(), (int)MAX_RECENT_FILES);
  for (int i = 0; i < numRecentFiles; ++i)
  {
    QString text = tr("&%1 %2").arg(i + 1).arg(MyUtils::Win32PathProof(files[i]));
    m_actionRecentVolumes[i]->setText(text);
    m_actionRecentVolumes[i]->setData(files[i]);
    m_actionRecentVolumes[i]->setVisible(true);
  }
  for (int j = numRecentFiles; j < MAX_RECENT_FILES; ++j)
  {
    m_actionRecentVolumes[j]->setVisible(false);
  }

  files = settings.value("MainWindow/RecentSurfaceFiles").toStringList();
  numRecentFiles = qMin(files.size(), (int)MAX_RECENT_FILES);
  for (int i = 0; i < numRecentFiles; ++i)
  {
    QString text = tr("&%1 %2").arg(i + 1).arg(MyUtils::Win32PathProof(files[i]));
    m_actionRecentSurfaces[i]->setText(text);
    m_actionRecentSurfaces[i]->setData(files[i]);
    m_actionRecentSurfaces[i]->setVisible(true);
  }
  for (int j = numRecentFiles; j < MAX_RECENT_FILES; ++j)
  {
    m_actionRecentSurfaces[j]->setVisible(false);
  }

  bool bHasVolumes = m_actionRecentVolumes[0]->isVisible();
  ui->actionVolumes->setVisible( bHasVolumes );
  bool bHasSurfaces = m_actionRecentSurfaces[0]->isVisible();
  ui->actionSurfaces->setVisible( bHasSurfaces );
}

void MainWindow::OnRecentVolumeFile()
{
  QAction* action = qobject_cast<QAction*>(sender());
  if ( !action )
  {
    return;
  }
  QString filename = action->data().toString();
  if ( !filename.isEmpty() )
  {
    LoadVolumeFile( filename );
  }
}

void MainWindow::OnRecentSurfaceFile()
{
  QAction* action = qobject_cast<QAction*>(sender());
  if ( !action )
  {
    return;
  }
  QString filename = action->data().toString();
  if ( !filename.isEmpty() )
  {
    LoadSurfaceFile( filename );
  }
}

LayerCollection* MainWindow::GetLayerCollection( const QString& strType )
{
  if ( m_layerCollections.contains( strType ) )
  {
    return m_layerCollections[strType];
  }
  else
  {
    return NULL;
  }
}

Layer* MainWindow::GetActiveLayer( const QString& strType )
{
  return GetLayerCollection( strType )->GetActiveLayer();
}

Layer* MainWindow::GetTopVisibleLayer(const QString &strType)
{
  QList<Layer*> layers = GetLayers(strType);
  foreach (Layer* layer, layers)
  {
    if (layer->IsVisible())
      return layer;
  }
  return NULL;
}

QList<Layer*> MainWindow::GetLayers(const QString &strType)
{
  return GetLayerCollection( strType )->GetLayers();
}

void MainWindow::OnSetViewLayout( QAction* action )
{
  if ( action == ui->actionLayout1x1 )
  {
    SetViewLayout( VL_1x1 );
  }
  else if ( action == ui->actionLayout2x2 )
  {
    SetViewLayout( VL_2x2 );
  }
  else if ( action == ui->actionLayout1n3 )
  {
    SetViewLayout( VL_1n3 );
  }
  else if ( action == ui->actionLayout1n3h )
  {
    SetViewLayout( VL_1n3h );
  }
}

void MainWindow::OnSetMainView( QAction* action )
{
  if ( action == ui->actionViewSagittal )
  {
    SetMainView( MV_Sagittal );
  }
  else if ( action == ui->actionViewCoronal )
  {
    SetMainView( MV_Coronal );
  }
  else if ( action == ui->actionViewAxial )
  {
    SetMainView( MV_Axial );
  }
  else if ( action == ui->actionView3D )
  {
    SetMainView( MV_3D );
  }
}

void MainWindow::SetMainView( int nView )
{
  if ( nView == m_nMainView )
  {
    return;
  }
  m_nMainView = nView;
  emit MainViewChanged( nView );
  SetViewLayout( m_nViewLayout );
}

void MainWindow::SetViewLayout( int nLayout )
{
  RenderView* view[4] = { 0 };
  switch ( m_nMainView )
  {
  case MV_Coronal:
    view[0] = ui->viewCoronal;
    view[1] = ui->viewSagittal;
    view[2] = ui->viewAxial;
    view[3] = ui->view3D;
    break;
  case MV_Axial:
    view[0] = ui->viewAxial;
    view[1] = ui->viewSagittal;
    view[2] = ui->viewCoronal;
    view[3] = ui->view3D;
    break;
  case MV_3D:
    view[0] = ui->view3D;
    view[1] = ui->viewSagittal;
    view[2] = ui->viewCoronal;
    view[3] = ui->viewAxial;
    break;
  default:
    view[0] = ui->viewSagittal;
    view[1] = ui->viewCoronal;
    view[2] = ui->viewAxial;
    view[3] = ui->view3D;
    break;
  }

  ui->gridLayoutViewport->setRowStretch( 0, 0 );
  ui->gridLayoutViewport->setRowStretch( 1, 0 );
  ui->gridLayoutViewport->setColumnStretch( 0, 0 );
  ui->gridLayoutViewport->setColumnStretch( 1, 0 );

  for ( int i = 0; i < 4; i++ )
  {
    ui->gridLayoutViewport->removeWidget( view[i] );
  }

  switch ( nLayout )
  {
  case VL_2x2:
    for ( int i = 0; i < 4; i++ )
    {
      ui->gridLayoutViewport->addWidget( view[i], i/2, i%2 );
    }
    break;
  case VL_1n3:
    ui->gridLayoutViewport->addWidget( view[0], 0, 0, 1, 3 );
    ui->gridLayoutViewport->addWidget( view[1], 1, 0 );
    ui->gridLayoutViewport->addWidget( view[2], 1, 1 );
    ui->gridLayoutViewport->addWidget( view[3], 1, 2 );
    ui->gridLayoutViewport->setRowStretch( 0, 2 );
    ui->gridLayoutViewport->setRowStretch( 1, 1 );
    break;
  case VL_1n3h:
    ui->gridLayoutViewport->addWidget( view[0], 0, 0, 3, 1 );
    ui->gridLayoutViewport->addWidget( view[1], 0, 1 );
    ui->gridLayoutViewport->addWidget( view[2], 1, 1 );
    ui->gridLayoutViewport->addWidget( view[3], 2, 1 );
    ui->gridLayoutViewport->setColumnStretch( 0, 2 );
    ui->gridLayoutViewport->setColumnStretch( 1, 1 );
    break;
  default:
    ui->gridLayoutViewport->addWidget( view[0], 0, 0 );
    break;
  }

  view[0]->show();
  for ( int i = 1; i < 4; i++ )
  {
    view[i]->setVisible( nLayout != VL_1x1 );
  }

  if ( m_nViewLayout != nLayout )
  {
    m_nViewLayout = nLayout;
    emit ViewLayoutChanged( nLayout );
  }
}

LayerCollection* MainWindow::GetCurrentLayerCollection()
{
  LayerCollection* lc = NULL;
  /*
  QString name = ui->tabWidgetControlPanel->tabText( ui->tabWidgetControlPanel->currentIndex() );
  if ( name == "Volumes" )
  {
    lc = GetLayerCollection( "MRI" );
  }
  else if ( name == "ROIs" )
  {
    lc = GetLayerCollection( "ROI" );
  }
  else if ( name == "Surfaces" )
  {
    lc = GetLayerCollection( "Surface" );
  }
  else if ( name == "Point Sets" )
  {
    lc = GetLayerCollection( "PointSet" );
  }
  else if ( name == "Tracks" )
  {
    lc = GetLayerCollection( "Track");
  }
  else if ( name == "All")
  {
    lc = GetLayerCollection(ui->tabAllLayers->GetCurrentLayerType());
  }
  */

  return lc;
}

QString MainWindow::GetCurrentLayerType()
{
  return ui->widgetAllLayers->GetCurrentLayerType();
}

bool MainWindow::SetSlicePosition( int nPlane, double dPos, bool bRoundToGrid )
{
  bool bRet = false;
  QStringList keys = m_layerCollections.keys();
  for ( int i = 0; i < keys.size(); i++ )
  {
    m_layerCollections[keys[i]]->blockSignals( true );
    if ( m_layerCollections[keys[i]]->SetSlicePosition( nPlane, dPos, bRoundToGrid ) )
    {
      bRet = true;
    }
    m_layerCollections[keys[i]]->blockSignals( false );
  }
  if ( bRet )
  {
    emit SlicePositionChanged();
  }

  return bRet;
}

bool MainWindow::SetSlicePosition( double* pos )
{
  bool bRet = false;
  QStringList keys = m_layerCollections.keys();
  for ( int i = 0; i < keys.size(); i++ )
  {
    m_layerCollections[keys[i]]->blockSignals( true );
    if ( m_layerCollections[keys[i]]->SetSlicePosition( pos ) )
    {
      bRet = true;
    }
    m_layerCollections[keys[i]]->blockSignals( false );
  }
  if ( bRet )
  {
    emit SlicePositionChanged();
  }

  return bRet;
}

bool MainWindow::OffsetSlicePosition( int nPlane, double dPosDiff, bool bRoundToGrid  )
{
  bool bRet = false;
  QStringList keys = m_layerCollections.keys();
  for ( int i = 0; i < keys.size(); i++ )
  {
    m_layerCollections[keys[i]]->blockSignals( true );
    if ( m_layerCollections[keys[i]]->OffsetSlicePosition( nPlane, dPosDiff, bRoundToGrid ) )
    {
      bRet = true;
    }
    m_layerCollections[keys[i]]->blockSignals( false );
  }
  if ( bRet )
  {
    emit SlicePositionChanged();
  }

  return bRet;
}

void MainWindow::OnLoadVolume()
{
  DialogLoadVolume dlg( this );
  dlg.SetLastDir( m_strLastDir );
  QStringList recentFiles;
  for ( int i = 0; i < m_actionRecentVolumes.size(); i++ )
  {
    if ( m_actionRecentVolumes[i]->isVisible() )
    {
      recentFiles << m_actionRecentVolumes[i]->data().toString();
    }
  }
  dlg.SetRecentFiles( recentFiles );
  if ( dlg.exec() == QDialog::Accepted )
  {
    QStringList filenames = dlg.GetVolumeFileNames();
    QString reg_fn =  dlg.GetRegFileName();
    bool bHasVolume = !GetLayerCollection( "MRI" )->IsEmpty();
    bool bHasSurface = !GetLayerCollection( "Surface" )->IsEmpty();
    for (int i = 0; i < filenames.size(); i++)
    {
      QStringList script("loadvolume");
      QString fn = filenames[i];
      if ( !reg_fn.isEmpty() )
      {
        fn += ":reg=" + reg_fn;
      }

      if ( dlg.GetSampleMethod() == SAMPLE_TRILINEAR )
        fn += ":sample=trilinear";
      else if (dlg.GetSampleMethod() == SAMPLE_CUBIC_BSPLINE)
        fn += ":sample=cubic";

      fn += ":colormap=" + dlg.GetColorMap();
      if ( dlg.GetColorMap() == "lut" )
      {
        fn += ":lut=" + dlg.GetLUT();
      }

      script << fn;

      if ( (!bHasVolume && bHasSurface) || (!bHasVolume && dlg.IsToResample()) || m_bResampleToRAS )
      {
        script << "r";
      }

      AddScript( script );
    }
  }
}

void MainWindow::LoadVolumeFile( const QString& filename,
                                 const QString& reg_filename,
                                 bool bResample_in, int nSampleMethod,
                                 bool bConform,
                                 int nGotoLabelOrientation,
                                 const QString& strGotoLabelName,
                                 const QVariantMap& sup_data)
{
  QFileInfo fi(filename);
  bool bResample = bResample_in;
  if ( GetLayerCollection( "MRI")->IsEmpty())
  {
    if (!GetLayerCollection( "Surface" )->IsEmpty() || !GetLayerCollection("Track")->IsEmpty())
    {
      bResample = true;
    }
  }

  m_bResampleToRAS = bResample;
  LayerMRI* layer = new LayerMRI( m_layerVolumeRef );
  layer->SetResampleToRAS( bResample );
  layer->SetSampleMethod( nSampleMethod );
  layer->SetConform( bConform );
  layer->GetProperty()->SetLUTCTAB( m_luts->GetColorTable( 0 ) );
  layer->SetName( fi.completeBaseName() );
  layer->SetGotoLabel(nGotoLabelOrientation, strGotoLabelName);
  QString fullpath = fi.absoluteFilePath();
  if ( fullpath.isEmpty() )
  {
    fullpath = filename;
  }
  layer->SetFileName( fullpath );
  if ( !reg_filename.isEmpty() )
  {
    layer->SetRegFileName( QFileInfo( reg_filename ).absoluteFilePath() );
  }

  if ( !bResample && reg_filename.isEmpty())
  {
    LayerMRI* mri = (LayerMRI* )GetLayerCollection( "MRI" )->GetLayer( 0 );
    if ( mri )
    {
      layer->SetRefVolume( mri->GetSourceVolume() );
    }
  }

  if (sup_data.contains("Basis"))
    layer->SetLayerIndex(sup_data["Basis"].toInt());
  if (sup_data.contains("Percentile"))
    layer->GetProperty()->SetUsePercentile(sup_data["Percentile"].toBool());

  m_threadIOWorker->LoadVolume( layer );
}

bool MainWindow::OnCloseVolume()
{
  LayerMRI* layer = (LayerMRI*)GetActiveLayer( "MRI" );
  if ( !layer )
  {
    return false;
  }
  if ( layer->IsModified() )
  {
    if ( QMessageBox::question( this, "Volume Not Saved",
                                "Volume has been modifed and not been saved. Do you still want to continue?",
                                QMessageBox::Yes, QMessageBox::No ) == QMessageBox::No )
    {
      return false;
    }
  }
  GetLayerCollection( "MRI" )->RemoveLayer( layer );
  if (layer == m_layerVolumeRef)
  {
    m_layerVolumeRef = (LayerMRI*)GetActiveLayer("MRI");
  }
  return true;
}

void MainWindow::OnNewVolume()
{
  LayerCollection* col_mri = GetLayerCollection( "MRI" );
  if ( !col_mri->GetActiveLayer() )
  {
    QMessageBox::warning(this, "Error", "Can not create new volume without any visible template volume." );
    return;
  }
  DialogNewVolume dlg( this );
  dlg.SetVolumeName( "New Volume");
  if ( dlg.exec() == QDialog::Accepted )
  {
    LayerMRI* layer_template = dlg.GetTemplate();
    if (layer_template->IsTransformed())
    {
   //   QMessageBox::information(this, "Warning",
   //                            "New volume is constructed on template that has been rotated/transformed. It is recommended that you save the transformed volume and reload it before creating new volumes.");
    }
    LayerMRI* layer_new = new LayerMRI( layer_template );

    if ( !layer_new->Create( dlg.GetTemplate(), dlg.GetCopyVoxel(), dlg.GetDataType(), dlg.GetVoxelDataOption() ) )
    {
      QMessageBox::warning( this, "Error", "Can not create new volume." );
      delete layer_new;
      return;
    }
    layer_new->GetProperty()->SetLUTCTAB( m_luts->GetColorTable( 0 ) );
    layer_new->SetName( dlg.GetVolumeName() );
    col_mri->AddLayer( layer_new );
  }
}

void MainWindow::OnSaveVolume()
{
  LayerMRI* layer = (LayerMRI*)GetActiveLayer( "MRI" );
  if ( !layer )
  {
    return;
  }
  else if ( !layer->IsVisible() )
  {
    QMessageBox::warning( this, "Error", "Current volume layer is not visible. Please turn it on before saving." );
    return;
  }

  QString fn = layer->GetFileName();
  if ( fn.isEmpty() )
  {
    QString name = layer->GetName().trimmed();
    name.replace(" ", "_");
    fn = QFileDialog::getSaveFileName( this, "Save volume file",
                                       AutoSelectLastDir( "mri" ),
                                       "Volume files (*.mgz *.mgh *.nii *.nii.gz *.img *.mnc);;All files (*)");
  }

  if ( !fn.isEmpty() )
  {
    QFileInfo fi(fn);
    if ( fi.suffix() != "nii" && fi.suffix() != "img" &&
         fi.suffix() != "mgz" && fi.suffix() != "mgh" &&
         !fi.completeSuffix().contains("nii.gz"))
    {
      fn += ".mgz";
    }
    layer->SetFileName( fn );
    m_scripts.append(QStringList("savelayer") << QString::number(layer->GetID()));
  }
  else
  {
    m_scripts.clear();
  }
}

bool MainWindow::SaveVolumeAs()
{
  LayerCollection* col_mri = GetLayerCollection( "MRI" );
  LayerMRI* layer_mri = ( LayerMRI* )col_mri->GetActiveLayer();
  if ( !layer_mri)
  {
    return false;
  }
  else if ( !layer_mri->IsVisible() )
  {
    QMessageBox::warning( this, "Error", "Current volume layer is not visible. Please turn it on before saving." );
    return false;
  }

  QString fn;
  if (layer_mri->IsTransformed())
  {
    DialogSaveVolume dlg(this, layer_mri->GetFileName());
    if (dlg.exec() == QDialog::Accepted)
    {
      fn = dlg.GetFileName();
      layer_mri->SetWriteResampled(dlg.GetResample());
      layer_mri->SetCropToOriginal(dlg.GetCrop());
    }
  }
  else
    fn = QFileDialog::getSaveFileName( this, "Save volume",
                                       QFileInfo( layer_mri->GetFileName() ).absolutePath(),
                                       "Volume files (*.mgz *.mgh *.nii *.nii.gz *.img *.mnc);;All files (*)");
  if ( !fn.isEmpty() )
  {
    layer_mri->SetFileName( fn );
    OnSaveVolume();
    ui->widgetAllLayers->UpdateWidgets();
    return true;
  }
  else
    return false;
}

void MainWindow::SaveVolumeAsAndReload()
{

}

void MainWindow::OnLoadDTI()
{
  DialogLoadDTI dlg(this);
  dlg.Initialize( m_bResampleToRAS, GetLayerCollection( "MRI" )->IsEmpty() );
  dlg.SetLastDir( m_strLastDir );;

  if ( dlg.exec() != QDialog::Accepted )
  {
    return;
  }

  this->LoadDTIFile( dlg.GetVectorFileName(), dlg.GetFAFileName(), dlg.GetEigenvalueFileName(),
                    dlg.GetRegFileName(), dlg.IsToResample() );
}

void MainWindow::OnLoadTrackVolume()
{
  QString fn = QFileDialog::getOpenFileName( this, "Open track volume",
                                     AutoSelectLastDir( "mri" ),
                                     "Volume files (*.mgz *.mgh *.nii *.nii.gz *.img *.mnc);;All files (*)");
  if (!fn.isEmpty())
    this->LoadVolumeTrackFile(fn, this->m_bResampleToRAS);
}

void MainWindow::LoadDTIFile( const QString& fn_vector,
                              const QString& fn_fa,
                              const QString& fn_scale,
                              const QString& reg_filename,
                              bool bResample )
{
  m_bResampleToRAS = bResample;

  LayerDTI* layer = new LayerDTI( m_layerVolumeRef );
  layer->SetResampleToRAS( bResample );
  QString layerName = QFileInfo( fn_vector ).completeBaseName();
  if ( QFileInfo( fn_vector ).suffix() == "gz" )
  {
    layerName = QFileInfo( layerName ).completeBaseName();
  }
  layer->SetName( layerName );
  layer->SetFileName( fn_fa );
  layer->SetVectorFileName( fn_vector );
  if ( !reg_filename.isEmpty() )
  {
    layer->SetRegFileName( QFileInfo(reg_filename).absoluteFilePath() );
  }
  if (!fn_scale.isEmpty())
  {
    layer->SetEigenvalueFileName(QFileInfo(fn_scale).absoluteFilePath() );
  }
  m_threadIOWorker->LoadVolume( layer );
}

void MainWindow::LoadVolumeTrackFile(const QString &fn, bool bResample)
{
  m_bResampleToRAS = bResample;

  LayerVolumeTrack* layer = new LayerVolumeTrack( m_layerVolumeRef );
  layer->SetResampleToRAS( bResample );
  layer->GetProperty()->SetLUTCTAB( m_luts->GetColorTable( 0 ) );
  QString layerName = QFileInfo( fn ).completeBaseName();
  if ( QFileInfo( fn ).suffix() == "gz" )
  {
    layerName = QFileInfo( layerName ).completeBaseName();
  }
  layer->SetName( layerName );
  layer->SetFileName( fn );
  m_threadIOWorker->LoadVolume( layer );
}

void MainWindow::LoadPVolumeFiles( const QStringList& filenames, const QString& prefix, const QString& lut )
{
  LayerPLabel* layer = new LayerPLabel( m_layerVolumeRef );
  layer->SetVolumeFileNames( filenames );
  layer->SetFileNamePrefix( prefix );
  layer->SetLUT( lut );
  m_threadIOWorker->LoadVolume( layer );
}

void MainWindow::OnNewROI()
{
  LayerCollection* col_mri = GetLayerCollection( "MRI" );
  LayerMRI* layer_mri = ( LayerMRI* )col_mri->GetActiveLayer();
  if ( !layer_mri)
  {
    QMessageBox::warning( this, "Error", "Can not create new ROI without volume template.");
    return;
  }

  // enter the name of the new ROI
  DialogNewROI dlg( this );
  dlg.SetROIName( "New ROI" );
  if ( dlg.exec() == QDialog::Accepted )
  {
    // finally we are about to create new ROI.
    LayerCollection* col_roi = GetLayerCollection( "ROI" );
    if ( col_roi->IsEmpty() )
    {
      col_roi->SetWorldOrigin( col_mri->GetWorldOrigin() );
      col_roi->SetWorldSize( col_mri->GetWorldSize() );
      col_roi->SetWorldVoxelSize( col_mri->GetWorldVoxelSize() );
      col_roi->SetSlicePosition( col_mri->GetSlicePosition() );
    }
    LayerROI* layer_roi = new LayerROI( dlg.GetTemplate() );
    layer_roi->SetName( dlg.GetROIName() );
    col_roi->AddLayer( layer_roi );

    SetMode( RenderView::IM_ROIEdit );
  }
}

void MainWindow::OnLoadROI()
{
  QStringList filenames = QFileDialog::getOpenFileNames( this, "Select label file",
                          AutoSelectLastDir( "label" ),
                          "Label files (*)");
  for ( int i = 0; i < filenames.size(); i++)
  {
    this->AddScript( QStringList("loadroi") << filenames[i] );
  }
}

void MainWindow::LoadROIFile( const QString& fn, const QString& ref_vol, const QColor& color, double opacity, double threshold )
{
  LayerMRI* ref = NULL;
  LayerCollection* col_mri = GetLayerCollection( "MRI" );
  if ( ref_vol.isEmpty() )
  {
 //   cout << "No template volume given, using current volume as template for ROI " << fn.toAscii().constData() << ".\n";
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
      else if ( QFileInfo( mri->GetFileName() ).fileName() == ref_vol )
      {
        ref = mri;
        break;
      }
    }
    if ( ref == NULL )
    {
      cerr << "Can not find given template volume: " << ref_vol.toAscii().constData()
           << ". Using current volume as template for ROI " << fn.toAscii().constData() << ".\n";
      ref = (LayerMRI*)col_mri->GetActiveLayer();
    }
  }
  LayerROI* roi = new LayerROI( ref );
  roi->SetName( QFileInfo(fn).completeBaseName()  );
  roi->GetProperty()->SetColor(color);
  roi->GetProperty()->SetOpacity(opacity);
  roi->GetProperty()->SetThreshold(threshold);
  if ( roi->LoadROIFromFile( fn ) )
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

    m_strLastDir = QFileInfo( fn ).canonicalPath();
  //  ui->tabWidgetControlPanel->setCurrentWidget( ui->tabROI );
  }
  else
  {
    delete roi;
    QMessageBox::warning( this, "Error", QString("Can not load ROI from %1").arg(fn) );
  }
}

void MainWindow::OnSaveROI()
{
  LayerCollection* col_roi = GetLayerCollection( "ROI" );
  LayerROI* layer_roi = ( LayerROI* )col_roi->GetActiveLayer();
  if ( !layer_roi )
  {
    return;
  }
  else if ( !layer_roi->IsVisible() )
  {
    QMessageBox::warning( this, "Error", "Current ROI layer is not visible. Please turn it on before saving.");
    return;
  }
  QString fn = layer_roi->GetFileName();
  if ( fn.isEmpty() )
  {
    QString def_fn = AutoSelectLastDir( "label" ) + "/" + layer_roi->GetName() + ".label";
    fn = QFileDialog::getSaveFileName( this, "Select label file",
                                       def_fn,
                                       "Label files (*)");
  }

  if ( !fn.isEmpty() )
  {
    if ( QFileInfo(fn).suffix() != "label" )
    {
      fn += ".label";
    }
    layer_roi->SetFileName( fn );
    layer_roi->ResetModified();
    layer_roi->SaveROI();
  }
}

void MainWindow::OnSaveROIAs()
{
  LayerCollection* col_roi = GetLayerCollection( "ROI" );
  LayerROI* layer_roi = ( LayerROI* )col_roi->GetActiveLayer();
  if ( !layer_roi )
  {
    return;
  }
  else if ( !layer_roi->IsVisible() )
  {
    QMessageBox::warning( this, "Error", "Current ROI layer is not visible. Please turn it on before saving.");
    return;
  }
  QString def_fn = AutoSelectLastDir( "label" ) + "/" + layer_roi->GetName() + ".label";
  QString fn = QFileDialog::getSaveFileName( this, "Select label file",
                                       def_fn,
                                       "Label files (*)");

  if ( !fn.isEmpty() )
  {
    if ( QFileInfo(fn).suffix() != "label" )
    {
      fn += ".label";
    }
    layer_roi->SetFileName( fn );
    layer_roi->ResetModified();
    layer_roi->SaveROI();
  }
}

void MainWindow::OnCloseROI()
{
  LayerROI* layer = (LayerROI*)GetActiveLayer( "ROI" );
  if ( !layer )
  {
    return;
  }
  if ( layer->IsModified() )
  {
    if ( QMessageBox::question( this, "ROI Not Saved",
                                "ROI has been modifed and not been saved. Do you still want to continue?",
                                QMessageBox::Yes, QMessageBox::No ) == QMessageBox::No )
    {
      return;
    }
  }
  GetLayerCollection( "ROI" )->RemoveLayer( layer );
}

void MainWindow::OnNewPointSet()
{
  LayerCollection* col_mri = GetLayerCollection( "MRI" );
  LayerMRI* layer_mri = ( LayerMRI* )col_mri->GetActiveLayer();
  if ( !layer_mri)
  {
    QMessageBox::warning( this, "Error", "Can not create new ROI without volume template.");
    return;
  }

  // enter the name of the new point set
  DialogNewPointSet dlg( this );
  dlg.SetPointSetName( tr("New Point Set %1").arg(GetLayerCollection("PointSet")->GetNumberOfLayers()));
  if ( dlg.exec() == QDialog::Accepted )
  {
    // finally we are about to create new point set.
    LayerCollection* col_wp = GetLayerCollection( "PointSet" );
    if ( col_wp->IsEmpty() )
    {
      col_wp->SetWorldOrigin( col_mri->GetWorldOrigin() );
      col_wp->SetWorldSize( col_mri->GetWorldSize() );
      col_wp->SetWorldVoxelSize( col_mri->GetWorldVoxelSize() );
      col_wp->SetSlicePosition( col_mri->GetSlicePosition() );
    }
    LayerPointSet* layer_wp = new LayerPointSet( dlg.GetTemplate(), dlg.GetType() );
    layer_wp->SetName( dlg.GetPointSetName() );
    col_wp->AddLayer( layer_wp );

    SetMode( RenderView::IM_PointSetEdit );
  }
}

void MainWindow::LoadWayPointsFile( const QString& fn )
{
  this->LoadPointSetFile( fn, LayerPropertyPointSet::WayPoint);
}

void MainWindow::LoadControlPointsFile( const QString& fn )
{
  this->LoadPointSetFile( fn, LayerPropertyPointSet::ControlPoint );
}

void MainWindow::LoadPointSetFile( const QString& fn, int type )
{
  LayerCollection* col_mri = GetLayerCollection( "MRI" );
  LayerMRI* mri = ( LayerMRI* )col_mri->GetActiveLayer();
  LayerPointSet * wp = new LayerPointSet( mri, type );
  wp->SetName( QFileInfo( fn ).fileName() );
  if ( wp->LoadFromFile( fn ) )
  {
    LayerCollection* col_wp = GetLayerCollection( "PointSet" );
    if ( col_wp->IsEmpty() )
    {
      col_wp->SetWorldOrigin( col_mri->GetWorldOrigin() );
      col_wp->SetWorldSize( col_mri->GetWorldSize() );
      col_wp->SetWorldVoxelSize( col_mri->GetWorldVoxelSize() );
      col_wp->SetSlicePosition( col_mri->GetSlicePosition() );
    }
    col_wp->AddLayer( wp );

    m_strLastDir = QFileInfo( fn ).canonicalPath();

  //  ui->tabWidgetControlPanel->setCurrentWidget( ui->tabPointSet );
  }
  else
  {
    delete wp;
    QMessageBox::warning( this, "Error", QString("Can not load Way Points from ") + fn);
  }
}

void MainWindow::OnLoadPointSet()
{
  DialogLoadPointSet dlg( this );
  dlg.SetLastDir(m_strLastDir);
  if ( dlg.exec() == QDialog::Accepted )
  {
    QStringList fns = dlg.GetFileNames();
    for ( int i = 0; i < fns.size(); i++ )
    {
      int nType = dlg.GetPointSetType();
      if ( nType == -1 )  // auto
      {
        if ( FSPointSet::IsLabelFormat( fns[i] ) )
        {
          nType = LayerPropertyPointSet::WayPoint;
        }
        else
        {
          nType = LayerPropertyPointSet::ControlPoint;
        }
      }
      if ( nType == LayerPropertyPointSet::WayPoint )
      {
        AddScript(QStringList("loadwaypoints")<<fns[i]);
      }
      else if ( nType == LayerPropertyPointSet::ControlPoint )
      {
        AddScript(QStringList("loadcontrolpoints")<<fns[i]);
      }
    }
  }
}

void MainWindow::OnSavePointSet()
{
  LayerPointSet* layer = (LayerPointSet*)GetActiveLayer("PointSet");
  if ( !layer )
  {
    return;
  }
  else if ( !layer->IsVisible() )
  {
    QMessageBox::information( this, "Error", "Current Way Points layer is not visible. Please turn it on before saving.");
    return;
  }

  QString fn = layer->GetFileName();
  if ( fn.isEmpty() )
  {
    OnSavePointSetAs();
  }
  else
  {
    if ( layer->GetProperty()->GetType() == LayerPropertyPointSet::WayPoint &&
         QFileInfo(fn).suffix() != "label")
    {
      fn += ".label";
    }
    if ( layer->GetProperty()->GetType() == LayerPropertyPointSet::ControlPoint &&
         QFileInfo(fn).suffix() != "dat" )
    {
      fn += ".dat";
    }
    layer->SetFileName( fn );
    layer->ResetModified();
    if (layer->Save())
    {
      layer->ResetModified();
    }
    else
    {
      QMessageBox::warning( this, "Error", "Could not save point set. Check console for more information.");
      return;
    }
  }
}

void MainWindow::OnSavePointSetAs()
{
  LayerPointSet* layer = (LayerPointSet*)GetActiveLayer("PointSet");
  if ( !layer )
  {
    return;
  }
  else if ( !layer->IsVisible() )
  {
    QMessageBox::information( this, "Error", "Current Way Points layer is not visible. Please turn it on before saving.");
    return;
  }

  DialogSavePointSet dlg( this );
  dlg.SetType( layer->GetProperty()->GetType() );
  QString fn = layer->GetFileName();
  if (fn.isEmpty())
    fn = layer->GetName();
  dlg.SetFileName(fn);
  dlg.SetLastDir(m_strLastDir);
  if ( dlg.exec() == QDialog::Accepted )
  {
    layer->SetFileName( dlg.GetFileName() );
    layer->GetProperty()->SetType( dlg.GetType() );
    OnSavePointSet();
    ui->widgetAllLayers->UpdateWidgets();
  }
}

void MainWindow::OnClosePointSet()
{
  LayerPointSet* layer = (LayerPointSet*)GetActiveLayer( "PointSet" );
  if ( !layer )
  {
    return;
  }
  if ( layer->IsModified() )
  {
    if ( QMessageBox::question( this, "Point Set Not Saved",
                                "Point set has been modifed and not been saved. Do you still want to continue?",
                                QMessageBox::Yes, QMessageBox::No ) == QMessageBox::No )
    {
      return;
    }
  }
  GetLayerCollection( "PointSet" )->RemoveLayer( layer );
}

void MainWindow::OnLoadTrack()
{
  QStringList filenames = QFileDialog::getOpenFileNames( this, "Select track file",
                          m_strLastDir,
                          "Track files (*.trk);;All files (*)");
  if ( !filenames.isEmpty() )
  {
    for ( int i = 0; i < filenames.size(); i++ )
    {
      AddScript(QStringList("loadtrack") << filenames[i]);
    }
  }
}

void MainWindow::LoadTrackFile(const QString &fn)
{
  LayerTrack* layer = new LayerTrack( m_layerVolumeRef );
  layer->SetFileName( QFileInfo(fn).absoluteFilePath() );
  m_threadIOWorker->LoadTrack( layer );
}

void MainWindow::OnCloseTrack()
{
  LayerTrack* layer = (LayerTrack*)GetActiveLayer( "Track" );
  if ( !layer )
  {
    return;
  }

  GetLayerCollection( "Track" )->RemoveLayer( layer );
}

void MainWindow::OnLoadSurface()
{
  // user getSaveFilename as a hack to allow adding options next to filename
  QStringList filenames = QFileDialog::getOpenFileNames( this, "Select surface file",
                          AutoSelectLastDir( "surf" ),
                          "Surface files (*)", 0, QFileDialog::DontConfirmOverwrite);
  if ( !filenames.isEmpty() )
  {
    for ( int i = 0; i < filenames.size(); i++ )
    {
      if (!filenames[i].trimmed().isEmpty())
        AddScript(QStringList("loadsurface") << filenames[i]);
    }
  }
}

void MainWindow::LoadSurfaceFile( const QString& filename, const QString& fn_patch, const QString& fn_target,
                                  const QStringList& sup_files)
{
  QFileInfo fi( filename );
  m_strLastDir = fi.absolutePath();
  LayerSurface* layer = new LayerSurface( m_layerVolumeRef );
  connect(layer, SIGNAL(CurrentVertexChanged(int)), m_wndGroupPlot, SLOT(SetCurrentVertex(int)), Qt::UniqueConnection);
  connect(ui->treeWidgetCursorInfo, SIGNAL(VertexChangeTriggered(int)), m_wndGroupPlot, SLOT(SetCurrentVertex(int)), Qt::UniqueConnection);
  layer->SetName( fi.fileName() );
  QString fullpath = fi.absoluteFilePath();
  if ( fullpath.isEmpty() )
  {
    fullpath = filename;
  }
  layer->SetFileName( fullpath );
  layer->SetPatchFileName( fn_patch );
  layer->SetTargetFileName( fn_target );
  layer->SetLoadSupSurfaces(sup_files);

  m_threadIOWorker->LoadSurface( layer );
  m_statusBar->StartTimer();
}

void MainWindow::OnCloseSurface()
{
  Layer* layer = GetActiveLayer( "Surface" );
  if ( !layer )
  {
    return;
  }

  GetLayerCollection( "Surface" )->RemoveLayer( layer );
}

void MainWindow::OnIOError( Layer* layer, int jobtype )
{
  ClearScripts();
  QString msg = QString("Failed to load %1 ").arg(layer->GetEndType());
  if (jobtype != ThreadIOWorker::JT_LoadSurfaceOverlay)
  {
    if ( jobtype == ThreadIOWorker::JT_SaveVolume )
    {
      msg = "Failed to save volume to ";
    }
    else if ( jobtype == ThreadIOWorker::JT_SaveSurface )
    {
      msg = "Failed to save surface to ";
    }
    QMessageBox::warning( this, "Error", msg + layer->GetFileName() );
    if ( jobtype != ThreadIOWorker::JT_SaveVolume && jobtype != ThreadIOWorker::JT_SaveSurface )
      delete layer;
  }
  else
  {
    QMessageBox::warning( this, "Error", msg + "overlay." );
  }
  m_bProcessing = false;
  m_volumeSettings.clear();
  m_surfaceSettings.clear();
}

void MainWindow::OnIOFinished( Layer* layer, int jobtype )
{
  m_statusBar->StopTimer();
  LayerCollection* lc_mri = GetLayerCollection( "MRI" );
  LayerCollection* lc_surface = GetLayerCollection( "Surface" );
  LayerCollection* lc_track = GetLayerCollection( "Track" );
  if ( jobtype == ThreadIOWorker::JT_LoadVolume && layer->IsTypeOf( "MRI" ) )
  {
    LayerMRI* mri = qobject_cast<LayerMRI*>( layer );
    if ( lc_mri->IsEmpty() )
    {
      double worigin[3], wsize[3];
      mri->GetWorldOrigin( worigin );
      mri->GetWorldSize( wsize );
      for ( int i = 0; i < 4; i++ )
      {
        m_views[i]->SetWorldCoordinateInfo( worigin, wsize, lc_surface->IsEmpty() );
      }
      if ( lc_surface->IsEmpty() )
      {
        mri->SetSlicePositionToWorldCenter();
      }
      else
      {
        mri->SetSlicePosition( lc_surface->GetSlicePosition() );
        lc_surface->SetWorldVoxelSize( mri->GetWorldVoxelSize() );
        lc_surface->SetWorldOrigin( mri->GetWorldOrigin() );
        lc_surface->SetWorldSize( mri->GetWorldSize() );
      }

      lc_mri->AddLayer( layer, true );
      lc_mri->SetCursorRASPosition( lc_mri->GetSlicePosition() );
      SetSlicePosition( lc_mri->GetSlicePosition() );
      m_layerVolumeRef = mri;

      if ( lc_surface->IsEmpty() )
      {
        for ( int i = 0; i < 4; i++ )
        {
          this->m_views[i]->SetWorldCoordinateInfo( worigin, wsize );
        }
      }
      m_views[3]->ResetCameraClippingRange();

      if (m_bShowTransformWindow)
      {
        this->OnTransformVolume();
        m_bShowTransformWindow = false;
      }
    }
    else
    {
      lc_mri->AddLayer( layer );
    }

    m_strLastDir = QFileInfo( layer->GetFileName() ).canonicalPath();
    SetCurrentFile( layer->GetFileName(), 0 );
//    ui->tabWidgetControlPanel->setCurrentWidget( ui->tabVolume );
    if (!m_volumeSettings.isEmpty())
    {
      mri->GetProperty()->RestoreFullSettings(m_volumeSettings);
      m_volumeSettings.clear();
    }
  }
  else if ( jobtype == ThreadIOWorker::JT_SaveVolume && layer->IsTypeOf( "MRI" ) )
  {
    SetCurrentFile( layer->GetFileName(), 0 );
  }
  else if ( jobtype == ThreadIOWorker::JT_LoadSurface && layer->IsTypeOf( "Surface" ) )
  {
    LayerSurface* sf = qobject_cast<LayerSurface*>( layer );
    if ( lc_surface->IsEmpty() )
    {
      double worigin[3], wsize[3];
      sf->GetWorldOrigin( worigin );
      sf->GetWorldSize( wsize );
      sf->SetSlicePositionToWorldCenter();
      if ( lc_mri->IsEmpty() )
      {
        for ( int i = 0; i < 4; i++ )
        {
          this->m_views[i]->SetWorldCoordinateInfo( worigin, wsize );
        }
        lc_surface->AddLayer( sf, true );
        SetSlicePosition( sf->GetSlicePosition() );
      }
      else
      {
        double mri_origin[3], mri_size[3], vs[3];
        lc_mri->GetWorldOrigin(mri_origin);
        lc_mri->GetWorldSize(mri_size);
        lc_mri->GetWorldVoxelSize(vs);

        for (int i = 0; i < 3; i++)
        {
          double upper = worigin[i] + wsize[i];
          if (worigin[i] >= mri_origin[i])
            worigin[i] = mri_origin[i];
          else
            worigin[i] = mri_origin[i] - ((int)((mri_origin[i]-worigin[i])/vs[i]+1))*vs[i];
          if (upper <= mri_origin[i]+mri_size[i])
            wsize[i] = mri_origin[i]+mri_size[i] - worigin[i];
          else
            wsize[i] = mri_origin[i]+mri_size[i] + ((int)((upper-mri_origin[i]-mri_size[i])/vs[i]+1))*vs[i];
        }
        lc_surface->SetWorldOrigin( worigin );
        lc_surface->SetWorldSize( wsize );
        lc_surface->SetWorldVoxelSize( vs );
        lc_surface->SetSlicePosition( lc_mri->GetSlicePosition() );
        lc_surface->AddLayer( sf );
      }
    }
    else
    {
      lc_surface->AddLayer( layer );
    }

    if ( !sf->HasValidVolumeGeometry() )
    {
    //  ShowNonModalMessage("Warning",
    //                      "Either this surface does not contain valid volume geometry information, or freeview failed to read the information. This surface may not align with volumes and other surfaces.");
      cerr << "Did not find any volume geometry information in the surface" << endl;
    }

    m_strLastDir = QFileInfo( layer->GetFileName() ).canonicalPath();
    SetCurrentFile( layer->GetFileName(), 1 );
//    ui->tabWidgetControlPanel->setCurrentWidget( ui->tabSurface );
    if (!m_surfaceSettings.isEmpty())
    {
      sf->GetProperty()->RestoreFullSettings(m_surfaceSettings);
      m_surfaceSettings.clear();
    }
    if (UpdateSurfaceCorrelation((LayerSurface*)layer) )
    {
      emit SlicePositionChanged();
    }
  }
  else if ( jobtype == ThreadIOWorker::JT_LoadSurfaceOverlay && layer->IsTypeOf("Surface") )
  {
    UpdateSurfaceCorrelation((LayerSurface*)layer);
    lc_surface->SetActiveLayer(layer);
    m_strLastDir = QFileInfo(((LayerSurface*)layer)->GetActiveOverlay()->GetFileName()).absolutePath();
    emit SlicePositionChanged();
  }
  else if ( jobtype == ThreadIOWorker::JT_LoadTrack && layer->IsTypeOf("Track"))
  {
    LayerTrack* track = qobject_cast<LayerTrack*>( layer );
    lc_track->AddLayer( track );
    m_strLastDir = QFileInfo( layer->GetFileName() ).canonicalPath();
//    ui->tabWidgetControlPanel->setCurrentWidget( ui->tabTrack );
    if (lc_surface->IsEmpty() && lc_mri->IsEmpty())
    {
      double worigin[3], wsize[3];
      track->GetWorldOrigin( worigin );
      track->GetWorldSize( wsize );
      for ( int i = 0; i < 4; i++ )
      {
        m_views[i]->SetWorldCoordinateInfo( worigin, wsize, true );
      }
      m_views[3]->ResetCameraClippingRange();
    }
  }
  else if (jobtype == ThreadIOWorker::JT_LoadConnectome && layer->IsTypeOf("CMAT"))
  {
    LayerConnectomeMatrix* cmat = qobject_cast<LayerConnectomeMatrix*>( layer );
    LayerCollection* lc_cmat = GetLayerCollection( "CMAT" );
    lc_cmat->AddLayer(cmat);
    double worigin[3], wsize[3];
    cmat->GetWorldOrigin(worigin);
    cmat->GetWorldSize(wsize);
    if (lc_mri->IsEmpty() && lc_surface->IsEmpty())
    {
      for ( int i = 0; i < 4; i++ )
      {
        m_views[i]->SetWorldCoordinateInfo( worigin, wsize, true );
      }
    }
  }
  else if (jobtype == ThreadIOWorker::JT_LoadFCD && layer->IsTypeOf("FCD"))
  {
    LayerFCD* fcd = qobject_cast<LayerFCD*>( layer );
    LayerCollection* lc = GetLayerCollection( "FCD" );
    LayerCollection* col_mri = GetLayerCollection("MRI");
    if ( lc->IsEmpty() )
    {
      lc->SetWorldOrigin( col_mri->GetWorldOrigin() );
      lc->SetWorldSize( col_mri->GetWorldSize() );
      lc->SetWorldVoxelSize( col_mri->GetWorldVoxelSize() );
      lc->SetSlicePosition( col_mri->GetSlicePosition() );
    }
    lc->AddLayer( fcd );
  }
  m_bProcessing = false;

  if ( jobtype == ThreadIOWorker::JT_SaveVolume)
  {
    std::cout << qPrintable(qobject_cast<LayerMRI*>(layer)->GetFileName()) << " saved successfully.\n";
  }
  else if ( jobtype == ThreadIOWorker::JT_SaveSurface)
  {
    std::cout << qPrintable(qobject_cast<LayerSurface*>(layer)->GetFileName()) << " saved successfully.\n";
    m_dlgRepositionSurface->UpdateUI();
  }

  if (m_scripts.isEmpty())
  {
    GetLayerCollection("MRI")->ClearLayerIndices();
  }
}

void MainWindow::OnFCDLoadFinished(LayerFCD *fcd)
{
  QList<LayerMRI*> mri_layers = fcd->GetMRILayers();
  foreach (LayerMRI* mri, mri_layers)
  {
    OnIOFinished(mri, ThreadIOWorker::JT_LoadVolume);
  }
  QList<LayerSurface*> surf_layers = fcd->GetSurfaceLayers();
  foreach (LayerSurface* surf, surf_layers)
  {
    OnIOFinished(surf, ThreadIOWorker::JT_LoadSurface);
  }
  OnIOFinished(fcd, ThreadIOWorker::JT_LoadFCD);
}

bool MainWindow::UpdateSurfaceCorrelation(LayerSurface *layer)
{
  QList<Layer*> layers = GetLayerCollection("Surface")->GetLayers();
  SurfaceOverlay* overlay = layer->GetActiveOverlay();
  LayerSurface* src = 0;
  if (overlay && overlay->HasCorrelationData())
  {
    src = layer;
  }
  else
  {
    for (int i = 0; i < layers.size(); i++)
    {
      overlay = ((LayerSurface*)layers[i])->GetActiveOverlay();
      if (overlay && overlay->HasCorrelationData())
      {
        src = (LayerSurface*)layers[i];
        break;
      }
    }
  }
  if (!src)
  {
    return false;
  }
  overlay = src->GetActiveOverlay();
  if ( overlay->HasSharedCorrelationData() )
  {
    return false;
  }

  for (int i = 0; i < layers.size(); i++)
  {
    LayerSurface* pair = (LayerSurface*)layers[i];
    if (src != pair && src->GetHemisphere() != pair->GetHemisphere() &&
        src->GetNumberOfVertices() == pair->GetNumberOfVertices())
    {
      pair->CopyCorrelationOverlay(src);
    }
  }

  return true;
}

void MainWindow::OnCycleLayer()
{
  LayerCollection* lc = GetLayerCollection(ui->widgetAllLayers->GetCurrentLayerType());
  if ( lc )
  {
    lc->CycleLayer(true, ui->viewAxial->GetInteractionMode() == RenderView2D::IM_ReconEdit);
  }
}

void MainWindow::OnReverseCycleLayer()
{
  LayerCollection* lc = GetLayerCollection(ui->widgetAllLayers->GetCurrentLayerType());
  if ( lc )
  {
    lc->CycleLayer( false, ui->viewAxial->GetInteractionMode() == RenderView2D::IM_ReconEdit );
  }
}

void MainWindow::OnEditUndo()
{
  if ( ui->viewAxial->GetInteractionMode() == RenderView2D::IM_ROIEdit )
  {
    LayerROI* roi = ( LayerROI* )GetLayerCollection( "ROI" )->GetActiveLayer();
    if ( roi )
    {
      roi->Undo();
    }
  }
  else if ( ui->viewAxial->GetInteractionMode() == RenderView2D::IM_VoxelEdit ||
            ui->viewAxial->GetInteractionMode() == RenderView2D::IM_ReconEdit ||
            ui->viewAxial->GetInteractionMode() == RenderView2D::IM_Navigate )
  {
    LayerMRI* mri = ( LayerMRI* )GetLayerCollection( "MRI" )->GetActiveLayer();
    if ( mri )
    {
      mri->Undo();
    }
  }
  else if ( ui->viewAxial->GetInteractionMode() == RenderView2D::IM_PointSetEdit )
  {
    LayerPointSet* wp = ( LayerPointSet* )GetLayerCollection( "PointSet" )->GetActiveLayer();
    if ( wp )
    {
      wp->Undo();
    }
  }
}

void MainWindow::OnEditRedo()
{
  if ( ui->viewAxial->GetInteractionMode() == RenderView2D::IM_ROIEdit )
  {
    LayerROI* roi = ( LayerROI* )GetLayerCollection( "ROI" )->GetActiveLayer();
    if ( roi )
    {
      roi->Redo();
    }
  }
  else if ( ui->viewAxial->GetInteractionMode() == RenderView2D::IM_VoxelEdit ||
            ui->viewAxial->GetInteractionMode() == RenderView2D::IM_ReconEdit ||
            ui->viewAxial->GetInteractionMode() == RenderView2D::IM_Navigate)
  {
    LayerMRI* mri = ( LayerMRI* )GetLayerCollection( "MRI" )->GetActiveLayer();
    if ( mri )
    {
      mri->Redo();
    }
  }
  else if ( ui->viewAxial->GetInteractionMode() == RenderView2D::IM_PointSetEdit )
  {
    LayerPointSet* wp = ( LayerPointSet* )GetLayerCollection( "PointSet" )->GetActiveLayer();
    if ( wp )
    {
      wp->Redo();
    }
  }
}

QString MainWindow::AutoSelectLastDir( const QString& subdirectory )
{
  return AutoSelectLastDir( m_strLastDir, subdirectory );
}

QString MainWindow::AutoSelectLastDir( const QString& lastDir, const QString& subdir )
{
  QDir dir( lastDir );
  QStringList stockdirs;
  stockdirs << "mri" << "label" << "scripts" << "surf" << "stats";
  if ( !stockdirs.contains( dir.dirName() ) )
  {
    while ( dir.cdUp() )
    {
      if ( stockdirs.contains( dir.dirName() ) )
      {
        break;
      }
    }
  }

  if ( stockdirs.contains( dir.dirName()) )
  {
    dir.cdUp();
    if ( !dir.cd( subdir ) )
    {
      dir.setPath( lastDir );
    }
  }
  else
  {
    dir.setPath( lastDir );
    dir.cd( subdir );
  }

//  qDebug() << dir.absolutePath();
  return dir.absolutePath();
}

void MainWindow::LoadSurfaceCurvature()
{
  QString filename = QFileDialog::getOpenFileName( this, "Select curvature file",
                     AutoSelectLastDir( "surf" ),
                     "Curvature files (*)");
  if ( !filename.isEmpty() )
  {
    this->LoadSurfaceCurvatureFile( filename );
  }
}

void MainWindow::LoadSurfaceCurvatureFile( const QString& filename )
{
  LayerSurface* layer = ( LayerSurface* )GetLayerCollection( "Surface" )->GetActiveLayer();
  if ( layer )
  {
    layer->LoadCurvatureFromFile( filename );
  //  m_strLastDir = fi.absoluteFilePath();
  }
}

void MainWindow::LoadSurfaceOverlay(bool bCorrelation)
{
  /*
  QString filename = QFileDialog::getOpenFileName( this, "Select overlay file",
                     AutoSelectLastDir( "surf" ),
                     "Overlay files (*)");
*/
  DialogLoadSurfaceOverlay dlg;
  dlg.SetLastDir(AutoSelectLastDir("surf"));
  if (dlg.exec() == QDialog::Accepted)
  {
    QString filename = dlg.GetFileName();
    QString reg_file = dlg.GetRegistration();
    this->LoadSurfaceOverlayFile( filename, reg_file, bCorrelation );
  }
}

void MainWindow::LoadSurfaceOverlayFile( const QString& filename, const QString& reg_file, bool bCorrelation, bool bSecondHalfData )
{
  LayerSurface* layer = ( LayerSurface* )GetLayerCollection( "Surface" )->GetActiveLayer();
  if ( layer )
  {
    QVariantMap args;
    args["FileName"] = filename;
    args["Correlation"] = bCorrelation;
    args["Registration"] = reg_file;
    args["SecondHalfData"] = bSecondHalfData;
    this->m_threadIOWorker->LoadSurfaceOverlay(layer, args);
//   m_strLastDir = QFileInfo(filename).absoluteFilePath();
  }
}

void MainWindow::LoadSurfaceAnnotation()
{
  QString filename = QFileDialog::getOpenFileName( this, "Select annotation file",
                     AutoSelectLastDir( "label" ),
                     "Annotation files (*)");
  if ( !filename.isEmpty() )
  {
    this->LoadSurfaceAnnotationFile( filename );
  }
}


void MainWindow::LoadSurfaceAnnotationFile( const QString& filename )
{
  LayerSurface* layer = ( LayerSurface* )GetLayerCollection( "Surface" )->GetActiveLayer();
  if ( layer && !layer->LoadAnnotationFromFile( filename ))
  {
    QMessageBox::warning(this, "Error", QString("Could not load annotation from %1").arg(filename));
  }
}

void MainWindow::LoadSurfaceLabel()
{
  QStringList filenames = QFileDialog::getOpenFileNames( this, "Select label files",
                     AutoSelectLastDir( "label" ),
                     "Label files (*)");
  if ( !filenames.isEmpty())
  {
    for (int i = 0; i < filenames.size(); i++)
     this->LoadSurfaceLabelFile( filenames[i] );
  }
}

bool MainWindow::LoadSurfaceLabelFile( const QString& filename )
{
  LayerSurface* layer = ( LayerSurface* )GetLayerCollection( "Surface" )->GetActiveLayer();
  if ( layer && !layer->LoadLabelFromFile( filename ))
  {
    QMessageBox::warning(this, "Error", QString("Could not load label from %1").arg(filename));
    return false;
  }
  return true;
}

void MainWindow::LoadSurfaceSpline()
{
  QString filename = QFileDialog::getOpenFileName( this, "Select surface spline file",
                     AutoSelectLastDir( "surf" ),
                     "All files (*)");
  if ( !filename.isEmpty() )
  {
    this->LoadSurfaceSplineFile( filename );
  }
}

void MainWindow::LoadSurfaceSplineFile(const QString &filename)
{
  LayerSurface* layer = ( LayerSurface* )GetLayerCollection( "Surface" )->GetActiveLayer();
  if ( layer && ! layer->LoadSplineFromFile( filename ))
  {
    QMessageBox::warning(this, "Error", QString("Could not load spline data from %1").arg(filename));;
  }
}

void MainWindow::LoadSurfaceVector()
{
  QString filename = QFileDialog::getOpenFileName( this, "Select surface vector file",
                     AutoSelectLastDir( "surf" ),
                     "All files (*)");
  if ( !filename.isEmpty() )
  {
    this->LoadSurfaceVectorFile( filename );
  }
}

void MainWindow::LoadSurfaceVectorFile( const QString& filename )
{
  LayerSurface* layer = ( LayerSurface* )GetLayerCollection( "Surface" )->GetActiveLayer();
  if ( layer )
  {
    layer->SetVectorFileName( filename );

    if ( !layer->LoadVectorFromFile() )
    {
      ShowNonModalMessage("Error", "Can not load vector file.");
    }
  }
}

void MainWindow::LoadLUT()
{
  QString filename = QFileDialog::getOpenFileName( this, "Select lookup table file",
                     m_strLastDir,
                     "LUT files (*)");
  if ( !filename.isEmpty() )
  {
    LayerMRI* mri = (LayerMRI*)GetLayerCollection( "MRI" )->GetActiveLayer();
    if ( mri )
    {
      COLOR_TABLE* ct = m_luts->LoadColorTable( filename );
      if ( ct )
      {
        mri->GetProperty()->SetLUTCTAB( ct );
      }
    }
  }
}

void MainWindow::SetMode( int nMode )
{
  int nOldMode = m_views[0]->GetInteractionMode();
  for ( int i = 0; i < 4; i++ )
  {
    m_views[i]->releaseMouse();
    m_views[i]->SetInteractionMode( nMode );
  }
  if (nOldMode == RenderView::IM_VolumeCrop && nMode != nOldMode)
  {
    m_volumeCropper->SetEnabled(false);
    RequestRedraw();
  }
  m_toolWindowEdit->UpdateReconMode();
}

int MainWindow::GetMode()
{
  return m_views[0]->GetInteractionMode();
}

void MainWindow::OnSetModeNavigate()
{
  SetMode( RenderView::IM_Navigate );
}

void MainWindow::OnSetModeMeasure()
{
  SetMode( RenderView::IM_Measure );
}

void MainWindow::OnSetModeVoxelEdit()
{
  SetMode( RenderView::IM_VoxelEdit );
}

void MainWindow::OnSetModeReconEdit()
{
  SetMode( RenderView::IM_ReconEdit );
}

void MainWindow::OnSetModeROIEdit()
{
  SetMode( RenderView::IM_ROIEdit );
}

void MainWindow::OnSetModePointSetEdit()
{
  SetMode( RenderView::IM_PointSetEdit );
}

void MainWindow::OnPreferences()
{
  if (!m_dlgPreferences)
  {
    m_dlgPreferences = new DialogPreferences( this );
    m_dlgPreferences->SetSettings(m_settings);
  }
  m_dlgPreferences->show();
}

void MainWindow::SetVolumeColorMap( int nColorMap, int nColorMapScale, const QList<double>& scales_in )
{
  LayerMRI* layer = qobject_cast<LayerMRI*>(GetActiveLayer("MRI"));
  if ( layer )
  {
    LayerPropertyMRI* p = layer->GetProperty();
    p->SetColorMap( (LayerPropertyMRI::ColorMapType) nColorMap );
    QList<double> scales = scales_in;
    if (p->GetUsePercentile())
    {
      for (int i = 0; i < scales.size(); i++)
        scales[i] = layer->GetHistoValueFromPercentile(scales[i]/100.0);
    }
    switch ( nColorMapScale )
    {
    case LayerPropertyMRI::Grayscale:
      if ( scales.size() >= 2 )
      {
        p->SetMinMaxGrayscaleWindow( scales[0], scales[1] );
      }
      else if ( !scales.empty() )
      {
        cerr << "Need 2 values for grayscale.\n";
      }
      break;
    case LayerPropertyMRI::Heat:
      if ( scales.size() >= 3 )
      {
        p->SetHeatScaleAutoMid(false);
        p->SetHeatScaleMinThreshold( scales[0] );
        p->SetHeatScaleMidThreshold( scales[1] );
        p->SetHeatScaleMaxThreshold( scales[2] );
      }
      else if ( scales.size() == 2 )
      {
        p->SetHeatScaleAutoMid(true);
        p->SetHeatScaleMinThreshold( scales[0] );
        p->SetHeatScaleMaxThreshold( scales[1] );
      }
      else if ( !scales.empty() )
      {
        cerr << "Need 2 or 3 values for heatscale.\n";
      }
      break;
    case LayerPropertyMRI::LUT:
      if ( scales.size() >= 1 )
      {
      }
      else if ( !scales.empty() )
      {
        cerr << "Need a value for lut.\n";
      }
      break;
    default:
      if ( scales.size() >= 2 )
      {
        p->SetMinMaxGenericThreshold( scales[0], scales[1] );
      }
      else if ( !scales.empty() )
      {
        cerr << "Need 2 values for colorscale.\n";
      }
      break;
    }
  }
}

void MainWindow::OnTransformVolume()
{
  if ( !m_dlgTransformVolume->isVisible() )
  {
    cout << "Warning: Transformation can only apply to volumes for now. If your data includes ROI/Surface/Way Points, please do not use this feature yet.\n";
    m_dlgTransformVolume->show();
    m_dlgTransformVolume->UpdateUI();
  }
}

void MainWindow::OnCropVolume()
{
  LayerMRI* mri = (LayerMRI*)GetActiveLayer( "MRI" );
  m_dlgCropVolume->SetVolume( mri );
  m_dlgCropVolume->show();
  m_volumeCropper->SetEnabled( true );
  m_volumeCropper->SetVolume( mri );
  m_volumeCropper->Show();
  SetMode( RenderView::IM_VolumeCrop );
}

void MainWindow::OnThresholdVolume()
{
  m_dlgThresholdVolume->show();
}

void MainWindow::OnSegmentVolume()
{
  m_dlgVolumeSegmentation->show();
}

void MainWindow::RotateVolume( std::vector<RotationElement>& rotations, bool bAllVolumes )
{
  // first update ROI and waypoints before their reference volume is rotated
  QList<Layer*> layers = GetLayerCollection("ROI")->GetLayers();
  bool bSuccess = true;
  for ( int i = 0; i < layers.size(); i++ )
  {
    ( (LayerROI*)layers[i] )->UpdateLabelData();
  }
  layers = GetLayerCollection("PointSet")->GetLayers();
  for ( int i = 0; i < layers.size(); i++ )
  {
    ( (LayerPointSet*)layers[i] )->UpdateLabelData();
  }

  layers = GetLayerCollection("MRI")->GetLayers();
  if ( bAllVolumes )
  {
    // then rotate MRI volumes
    for ( int i = 0; i < layers.size(); i++ )
    {
      if ( !layers[i]->Rotate( rotations ) )
      {
        bSuccess = false;
        break;
      }
    }
    // at last rotate others
    layers = GetLayerCollection("Surface")->GetLayers();
    for ( int i = 0; i < layers.size() && bSuccess; i++ )
    {
      if ( !layers[i]->Rotate( rotations) )
      {
        bSuccess = false;
        break;
      }
    }
  }
  else
  {
    LayerMRI* layer = (LayerMRI*) GetActiveLayer( "MRI" );
    if ( !layer->Rotate( rotations ) )
    {
      bSuccess = false;
    }
  }
  if ( !bSuccess )
  {
    ShowNonModalMessage("Error", "Error occured while rotating volumes.");
  }
}

void MainWindow::TransformVolume(double *mat, int sample_method)
{
  // first update ROI and waypoints before their reference volume is rotated
  QList<Layer*> layers = GetLayerCollection("ROI")->GetLayers();
  for ( int i = 0; i < layers.size(); i++ )
  {
    ( (LayerROI*)layers[i] )->UpdateLabelData();
  }
  layers = GetLayerCollection("PointSet")->GetLayers();
  for ( int i = 0; i < layers.size(); i++ )
  {
    ( (LayerPointSet*)layers[i] )->UpdateLabelData();
  }

  LayerMRI* layer = (LayerMRI*) GetActiveLayer( "MRI" );
  if (layer)
    layer->Transform(mat, sample_method);
  LayerLandmarks* lm = (LayerLandmarks*)GetSupplementLayer("Landmarks");
  lm->Transform(mat, sample_method);
}

void MainWindow::OnSaveScreenshot()
{
  if (!m_dlgSaveScreenshot)
  {
    m_dlgSaveScreenshot = new DialogSaveScreenshot(this);
    m_dlgSaveScreenshot->SetLastDir(m_strLastDir);
    m_dlgSaveScreenshot->SetSettings(m_settingsScreenshot);
  }
  m_dlgSaveScreenshot->show();
}

void MainWindow::OnVolumeFilterMean()
{
  LayerMRI* mri = (LayerMRI*)GetActiveLayer( "MRI" );
  if ( mri )
  {
    VolumeFilterMean* filter = new VolumeFilterMean( mri, mri );
    DialogVolumeFilter dlg( this );
    dlg.SetFilter( filter );
    dlg.ShowSigma( false );
    if ( dlg.exec() == QDialog::Accepted )
    {
      filter->SetKernelSize( dlg.GetKernelSize() );
      m_threadVolumeFilter->ExecuteFilter(filter);
    }
  }
}

void MainWindow::OnVolumeFilterMedian()
{
  LayerMRI* mri = (LayerMRI*)GetActiveLayer( "MRI" );
  if ( mri )
  {
    VolumeFilterMedian* filter = new VolumeFilterMedian( mri, mri );
    DialogVolumeFilter dlg( this );
    dlg.SetFilter( filter );
    dlg.ShowSigma( false );
    if ( dlg.exec() == QDialog::Accepted )
    {
      filter->SetKernelSize( dlg.GetKernelSize() );
      m_threadVolumeFilter->ExecuteFilter(filter);
    }
  }
}

void MainWindow::OnVolumeFilterConvolve()
{
  LayerMRI* mri = (LayerMRI*)GetActiveLayer( "MRI" );
  if ( mri )
  {
    VolumeFilterConvolve* filter = new VolumeFilterConvolve( mri, mri );
    DialogVolumeFilter dlg( this );
    dlg.SetFilter( filter );
    dlg.SetSigma( filter->GetSigma() );
    if ( dlg.exec() == QDialog::Accepted )
    {
      filter->SetKernelSize( dlg.GetKernelSize() );
      filter->SetSigma( dlg.GetSigma() );
      m_threadVolumeFilter->ExecuteFilter(filter);
    }
  }
}

void MainWindow::OnVolumeFilterGradient()
{
  LayerMRI* mri = (LayerMRI*)GetActiveLayer( "MRI" );
  if ( mri )
  {
    VolumeFilterGradient* filter = new VolumeFilterGradient( mri, mri );
    DialogGradientFilter dlg(this);
    dlg.SetSmoothing(filter->GetSmoothing());
    dlg.SetSD(filter->GetStandardDeviation());
    if ( dlg.exec() == QDialog::Accepted )
    {
      filter->SetSmoothing(dlg.GetSmoothing());
      filter->SetStandardDeviation(dlg.GetSD());
      m_threadVolumeFilter->ExecuteFilter(filter);
      mri->ResetWindowLevel();
    }
  }
}

void MainWindow::OnVolumeFilterThreshold()
{
//  LayerMRI* mri = (LayerMRI*)GetActiveLayer( "MRI" );
//  if ( mri )
//  {
//    VolumeFilterThreshold* filter = new VolumeFilterThreshold( mri, mri );
//    DialogThresholdFilter dlg(this);
//    if ( dlg.exec() == QDialog::Accepted )
//    {
//      double th[2];
//      dlg.GetThreshold(th);
//      filter->SetThreshold(th);
//      filter->SetReplaceIn(dlg.GetReplaceIn());
//      filter->SetReplaceOut(dlg.GetReplaceOut());
//      filter->SetInValue(dlg.GetInValue());
//      filter->SetOutValue(dlg.GetOutValue());
//      m_threadVolumeFilter->ExecuteFilter(filter);
//      if (dlg.GetReplaceIn() || dlg.GetReplaceOut())
//        mri->ResetWindowLevel();
//    }
//  }
}

void MainWindow::OnVolumeFilterSobel()
{
  LayerMRI* mri = (LayerMRI*)GetActiveLayer( "MRI" );
  if ( mri )
  {
    VolumeFilterSobel* filter = new VolumeFilterSobel( mri, mri );
    m_threadVolumeFilter->ExecuteFilter(filter);
    mri->ResetWindowLevel();
  }
}

void MainWindow::OnVolumeFilterErode()
{
  LayerMRI* mri = (LayerMRI*)GetActiveLayer( "MRI" );
  if ( mri )
  {
    VolumeFilterErode* filter = new VolumeFilterErode( mri, mri );
    m_threadVolumeFilter->ExecuteFilter(filter);
    mri->ResetWindowLevel();
  }
}


void MainWindow::OnVolumeFilterDilate()
{
  LayerMRI* mri = (LayerMRI*)GetActiveLayer( "MRI" );
  if ( mri )
  {
    VolumeFilterDilate* filter = new VolumeFilterDilate( mri, mri );
    m_threadVolumeFilter->ExecuteFilter(filter);
    mri->ResetWindowLevel();
  }
}


void MainWindow::OnVolumeFilterOpen()
{
  LayerMRI* mri = (LayerMRI*)GetActiveLayer( "MRI" );
  if ( mri )
  {
    VolumeFilterOpen* filter = new VolumeFilterOpen( mri, mri );
    m_threadVolumeFilter->ExecuteFilter(filter);
    mri->ResetWindowLevel();
  }
}


void MainWindow::OnVolumeFilterClose()
{
  LayerMRI* mri = (LayerMRI*)GetActiveLayer( "MRI" );
  if ( mri )
  {
    VolumeFilterClose* filter = new VolumeFilterClose( mri, mri );
    m_threadVolumeFilter->ExecuteFilter(filter);
    mri->ResetWindowLevel();
  }
}


void MainWindow::OnVolumeFilterFinished(VolumeFilter *filter)
{
  if (filter)filter->deleteLater();
}

void MainWindow::OnResetView()
{
  for (int i = 0; i < 4; i++)
  {
    m_views[i]->Reset();
  }
}

void MainWindow::OnSavePoint()
{
  QString fn;
  LayerCollection* lc = GetLayerCollection( "MRI" );
  Layer* layer = NULL;
  for ( int i = 0; i < lc->GetNumberOfLayers(); i++ )
  {
    fn = ( (LayerMRI*)lc->GetLayer( i ) )->GetFileName();
    if ( !fn.isEmpty() )
    {
      layer = lc->GetLayer(i);
      break;
    }
  }
  if ( fn.isEmpty() )
  {
    lc = GetLayerCollection( "Surface" );
    for ( int i = 0; i < lc->GetNumberOfLayers(); i++ )
    {
      fn = ( (LayerSurface*)lc->GetLayer( i ) )->GetFileName();
      if ( !fn.isEmpty() )
      {
        layer = lc->GetLayer(i);
        break;
      }
    }
  }

  bool bError = false;
  QString msg;
  if ( !fn.isEmpty() )
  {
    QString dir = AutoSelectLastDir( QFileInfo(fn).absolutePath(), "tmp" );
    fn = dir + "/edit.dat";
    QString path = getenv( "FS_SAVE_GOTO_POINT" );
    QString subjectName = layer->GetSubjectName();
    if (!path.isEmpty() && !subjectName.isEmpty())
    {
      fn = path + "-" + subjectName;
      dir = QFileInfo(fn).absolutePath();
    }
    if ( QDir(dir).exists() )
    {
      double ras[3];
      GetCursorRAS( ras, true );  // in tkReg coordinate
      QFile file(fn);
      if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
      {
        bError = true;
        msg = QString("Can not write to file ") + file.fileName();
      }
      else
      {
        QTextStream out(&file);
        out << ras[0] << " " << ras[1] << " " << ras[2] << "\n";
        file.close();
      }
    }
    else
    {
      bError = true;
      msg = QString("Directory ") + dir + " does not exist.";
    }
  }
  else
  {
    bError = true;
    msg = "Layer file name is empty. Can not decide where to save.";
  }
  if ( bError )
  {
    QMessageBox::warning(this, "Error", msg);
  }
}

void MainWindow::OnGoToPoint()
{
  LayerCollection* lc = GetLayerCollection( "MRI" );
  QString fn;
  QString path = getenv( "FS_SAVE_GOTO_POINT" );
  for ( int i = 0; i < lc->GetNumberOfLayers(); i++ )
  {
    fn = ( (LayerMRI*)lc->GetLayer( i ) )->GetFileName();
    QString subjectName = lc->GetLayer(i)->GetSubjectName();
    if (!path.isEmpty() && !subjectName.isEmpty())
      fn = path + "-" + subjectName;
    else
      fn = AutoSelectLastDir( QFileInfo(fn).absolutePath(), "tmp") + "/edit.dat";
    if ( !QFile::exists( fn ) )
    {
      fn = "";
    }
  }

  if ( !fn.isEmpty() )
  {
    QFile file( fn );
    file.open(QIODevice::ReadOnly | QIODevice::Text);
    QString strg = file.readAll();
    QStringList args = strg.trimmed().split(QRegExp("\\s+"));
    while (args.size() > 3)
      args.removeLast();
    args.insert( args.begin(), "setras" );
    args << "tkreg";
    CommandSetRAS( args );
  }
  else
  {
    QMessageBox::warning( this, "Error", "Could not find saved point.");
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
    {
      mri->NativeRASToTkReg( ras_out, ras_out );
    }
    return true;
  }
  else if ( !lc_surf->IsEmpty() )
  {
    lc_surf->GetCurrentRASPosition( ras_out );
    return true;
  }
  else
  {
    return false;
  }
}

void MainWindow::OnShowAnnotation(bool bShow)
{
  for ( int i = 0; i < 3; i++ )
  {
    ( (RenderView2D*)m_views[i] )->ShowCoordinateAnnotation(bShow);
  }

  RequestRedraw();
}

void MainWindow::OnShowColorBar(bool bShow)
{
  for ( int i = 0; i < 4; i++ )
  {
    if ( i != m_nMainView )
    {
      m_views[i]->ShowScalarBar( false );
    }
  }
  m_views[m_nMainView]->ShowScalarBar( bShow );
  m_views[m_nMainView]->UpdateScalarBar();

  RequestRedraw();
}

void MainWindow::OnCopy()
{
  int nWnd = GetActiveViewId();
  if ( ui->viewAxial->GetInteractionMode() == RenderView::IM_ROIEdit )
  {
    LayerROI* roi = ( LayerROI* )GetActiveLayer( "ROI" );
    if ( roi && nWnd >= 0 && nWnd < 3 )
    {
      roi->Copy( nWnd );
    }
  }
  else if ( ui->viewAxial->GetInteractionMode() == RenderView::IM_VoxelEdit ||
            ui->viewAxial->GetInteractionMode() == RenderView::IM_ReconEdit )
  {
    LayerMRI* mri = ( LayerMRI* )GetActiveLayer( "MRI" );
    if ( mri && nWnd >= 0 && nWnd < 3 )
    {
      mri->Copy( nWnd );
    }
  }
}

void MainWindow::OnCopyStructure()
{
  int nWnd = GetActiveViewId();
  if ( ui->viewAxial->GetInteractionMode() == RenderView::IM_VoxelEdit ||
       ui->viewAxial->GetInteractionMode() == RenderView::IM_ReconEdit )
  {
    LayerMRI* mri = ( LayerMRI* )GetActiveLayer( "MRI" );
    if ( mri && nWnd >= 0 && nWnd < 3 )
    {
      double* pos = mri->GetSlicePosition();
      if ( !mri->CopyStructure( nWnd, pos ) )
      {
        QMessageBox::warning( this, "Copy Structure Failed",
                              "Please move the cursor to the structure you want to copy and try again." );
      }
    }
  }
}

void MainWindow::OnPaste()
{
  int nWnd = GetActiveViewId();
  if ( ui->viewAxial->GetInteractionMode() == RenderView::IM_ROIEdit )
  {
    LayerROI* roi = ( LayerROI* )GetActiveLayer( "ROI" );
    if ( roi && nWnd >= 0 && nWnd < 3 )
    {
      roi->Paste( nWnd );
    }
  }
  else if ( ui->viewAxial->GetInteractionMode() == RenderView::IM_VoxelEdit ||
            ui->viewAxial->GetInteractionMode() == RenderView::IM_ReconEdit )
  {
    LayerMRI* mri = ( LayerMRI* )GetActiveLayer( "MRI" );
    if ( mri && nWnd >= 0 && nWnd < 3 )
    {
      mri->Paste( nWnd );
    }
  }
}

void MainWindow::ToggleSpecialVolume(const QString &name)
{
  QList<Layer*> layers = GetLayers("MRI");
  foreach (Layer* layer, layers)
  {
    if (layer->GetFileName().contains(name, Qt::CaseInsensitive))
    {
      layer->SetVisible(!layer->IsVisible());
    }
  }
}

void MainWindow::OnToggleWm()
{
  ToggleSpecialVolume("wm.mgz");
}

void MainWindow::OnToggleAseg()
{
  ToggleSpecialVolume("aseg.mgz");
}

void MainWindow::OnToggleBrainmask()
{
  ToggleSpecialVolume("brainmask.mgz");
}

void MainWindow::ToggleShowLayer(const QString& type )
{
  Layer* layer = GetActiveLayer(type);
  if (layer)
  {
    layer->SetVisible(!layer->IsVisible());
  }
}

void MainWindow::OnToggleShowVolume()
{
  ToggleShowLayer("MRI");
}

void MainWindow::OnToggleShowSurface()
{
  ToggleShowLayer("Surface");
}

void MainWindow::OnToggleShowROI()
{
  ToggleShowLayer("ROI");
}

void MainWindow::OnToggleShowPointSet()
{
  ToggleShowLayer("PointSet");
}

void MainWindow::OnToggleAllSurfaces()
{
  QList<Layer*> layers = GetLayers("Surface");
  bool bVisible = false;
  foreach (Layer* layer, layers)
  {
    if (layer->IsVisible())
    {
      bVisible = true;
      break;
    }
  }
  foreach (Layer* layer, layers)
    layer->SetVisible(!bVisible);
}

void MainWindow::OnAbout()
{
  DialogAbout dlg(this);
  dlg.exec();
}

void MainWindow::OnActiveLayerChanged(Layer* layer)
{
  if (!layer)
  {
    this->setWindowTitle("FreeView");
    m_wndTimeCourse->hide();
  }
  else
  {
    this->setWindowTitle(QString("FreeView (%1)")
                         .arg(MyUtils::Win32PathProof(layer->GetFileName())));
    if (layer->IsTypeOf("MRI") && !layer->IsTypeOf("DTI") && !layer->IsTypeOf("PLabel"))
    {
      if (((LayerMRI*)layer)->GetNumberOfFrames() > 1 && !((LayerMRI*)layer)->GetCorrelationSurface())
      {
        connect(layer, SIGNAL(ActiveFrameChanged(int)), m_wndTimeCourse, SLOT(SetCurrentFrame(int)), Qt::UniqueConnection);
        connect(layer, SIGNAL(ActiveFrameChanged(int)),
                ui->treeWidgetCursorInfo, SLOT(OnCursorPositionChanged()), Qt::UniqueConnection);
        connect(layer, SIGNAL(ActiveFrameChanged(int)),
                ui->treeWidgetMouseInfo, SLOT(OnMousePositionChanged()), Qt::UniqueConnection);
        connect(layer, SIGNAL(ActiveFrameChanged(int)), m_dlgLabelStats, SLOT(UpdateStats()), Qt::UniqueConnection);
        m_wndTimeCourse->UpdateData();
        if (ui->actionTimeCourse->isChecked() && !layer->IsTypeOf("VolumeTrack"))
        {
            m_wndTimeCourse->show();
        }
      }
    //  else
    //    m_wndTimeCourse->hide();
    }
    else if (layer->IsTypeOf("Surface"))
    {
      LayerSurface* surf = (LayerSurface*)layer;
      if (surf->GetActiveOverlay() && surf->GetActiveOverlay()->GetNumberOfFrames() > 1)
      {
        m_wndTimeCourse->UpdateData();
        if (ui->actionTimeCourse->isChecked())
        {
            m_wndTimeCourse->show();
        }
      }
  //    else
  //      m_wndTimeCourse->hide();
    }
  }
}

void MainWindow::OnLoadCommand()
{
  QString filename = QFileDialog::getOpenFileName( this, "Select command file",
                     m_strLastDir,
                     "Command files (*)");
  if (!filename.isEmpty())
  {
    AddScript(QStringList("loadcommand") << filename);
  }
}

void MainWindow::OnWriteMovieFrames()
{
  m_dlgWriteMovieFrames->show();
}

Layer* MainWindow::GetSupplementLayer(const QString &type)
{
  return m_layerCollections["Supplement"]->GetLayer(type);
}

void MainWindow::OnIncreaseOpacity()
{
  LayerMRI* mri = (LayerMRI*)this->GetActiveLayer( "MRI" );
  if ( mri )
  {
    double dOpacity = mri->GetProperty()->GetOpacity();
    mri->GetProperty()->SetOpacity( qMin(1., dOpacity+0.1) );
  }
}

void MainWindow::OnDecreaseOpacity()
{
  LayerMRI* mri = (LayerMRI*)this->GetActiveLayer( "MRI" );
  if ( mri )
  {
    double dOpacity = mri->GetProperty()->GetOpacity();
    mri->GetProperty()->SetOpacity( qMax(0., dOpacity-0.1) );
  }
}

void MainWindow::OnToggleCursorVisibility(bool bShow)
{
  ui->viewAxial->GetCursor2D()->Show( bShow );
  ui->viewSagittal->GetCursor2D()->Show( bShow );
  ui->viewCoronal->GetCursor2D()->Show( bShow );
  ui->view3D->GetCursor3D()->Show( bShow );
  this->RequestRedraw();
}

void MainWindow::ShowNonModalMessage(const QString &title, const QString &msg)
{
  m_dlgMessage->setWindowTitle(title);
  m_dlgMessage->setText(msg);
  m_dlgMessage->show();
}

void MainWindow::OnRepositionSurface()
{
  m_dlgRepositionSurface->show();
}

void MainWindow::OnSmoothSurface()
{
  m_dlgSmoothSurface->show();
}

void MainWindow::OnRemoveIntersectionsFromSurface()
{
  LayerSurface* surf = ( LayerSurface* )GetActiveLayer( "Surface" );
  if (surf)
  {
    surf->RemoveIntersections();
    emit SlicePositionChanged();
  }
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
    QMessageBox::warning( this, "Error", "Current surface layer is not visible. Please turn it on before saving.");
    return;
  }

  QString fn = layer_surf->GetFileName();
  if ( fn.isEmpty() )
  {
    QString name = layer_surf->GetName().trimmed();
    name.replace( " ", "_" );
    fn = QFileDialog::getSaveFileName( this, "Save surface",
                                       AutoSelectLastDir("surf"),
                                      "Surface files (*)");
  }

  if ( !fn.isEmpty() )
  {
    layer_surf->SetFileName( fn );
    m_threadIOWorker->SaveSurface( layer_surf );
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
    QMessageBox::warning( this, "Error", "Current surface layer is not visible. Please turn it on before saving.");
    return;
  }

  QString fn = QFileDialog::getSaveFileName( this, "Save surface as",
                                    layer_surf->GetFileName(),
                                    "Surface files (*)");
  if ( !fn.isEmpty() )
  {
    layer_surf->SetFileName(fn );
    SaveSurface();
    ui->widgetAllLayers->UpdateWidgets();
  }
}

void MainWindow::OnShowLabelStats()
{
  m_dlgLabelStats->show();
}

void MainWindow::OnLineProfile()
{
  m_dlgLineProfile->show();
}

void MainWindow::OnSaveIsoSurface()
{
  QString fn = QFileDialog::getSaveFileName(this, "Save IsoSurface As",
                                            m_strLastDir, "VTK files (*.vtk)");
  if (fn.isEmpty())
    return;

  LayerMRI* mri = qobject_cast<LayerMRI*>(GetActiveLayer("MRI"));
  if (mri && mri->GetProperty()->GetShowAsContour())
  {
    if (!mri->SaveIsoSurface(fn))
    {
      QMessageBox::warning(this, "Error", QString("Could not save iso surface to %1.").arg(fn));
    }
  }
}

bool MainWindow::IsRepositioningSurface()
{
  return this->m_dlgRepositionSurface->isVisible();
}

void MainWindow::OnPlot()
{
  QString fn = QFileDialog::getOpenFileName(this, "Select FSGD File", m_strLastFsgdDir,
                                            "FSGD files (*.fsgd);;All files (*)");
  if (fn.isEmpty())
    return;

  FSGroupDescriptor* fsgd = new FSGroupDescriptor(this);
  if (!fsgd->Read(fn))
  {
    QMessageBox::warning(this, "Error", "Failed to load FSGD file.");
    fsgd->deleteLater();
    return;
  }

  this->m_wndGroupPlot->SetFsgdData(fsgd);
  this->m_wndGroupPlot->show();
  this->m_wndGroupPlot->SetCurrentVertex(0);
  m_strLastFsgdDir = QFileInfo(fn).absolutePath();
}

void MainWindow::ToggleSplinePicking()
{
  m_bSplinePicking = !m_bSplinePicking;
  qDebug() << QString("Surface spline picking %1").arg(m_bSplinePicking?"enabled":"disabled");
}

void MainWindow::SetSplinePicking(bool b)
{
  if (b != m_bSplinePicking)
    ToggleSplinePicking();
}

void MainWindow::OnReloadVolume()
{
  LayerMRI* mri = qobject_cast<LayerMRI*>(this->GetActiveLayer("MRI"));
  if (mri)
  {
    DialogReloadLayer dlg;
    QString name = mri->GetName();
    QString filename = mri->GetFileName();
    QString reg_fn = mri->GetRegFileName();
    if (dlg.Execute(name, "Volume", filename) == QDialog::Accepted)
    {
      m_volumeSettings = mri->GetProperty()->GetFullSettings();
      if (dlg.GetCloseLayerFirst())
      {
        if (!OnCloseVolume())
        {
          m_volumeSettings.clear();
          return;
        }
      }
      this->LoadVolumeFile(filename, reg_fn);
    }
  }
}

void MainWindow::OnReloadSurface()
{
  LayerSurface* surf = qobject_cast<LayerSurface*>(this->GetActiveLayer("Surface"));
  if (surf)
  {
    DialogReloadLayer dlg;
    QString name = surf->GetName();
    QString filename = surf->GetFileName();
    if (dlg.Execute(name, "Surface", filename) == QDialog::Accepted)
    {
      m_surfaceSettings = surf->GetProperty()->GetFullSettings();
      if (dlg.GetCloseLayerFirst())
      {
        OnCloseSurface();
      }
      this->LoadSurfaceFile(filename);
    }
  }
}

void MainWindow::UpdateInfoPanel()
{
  QTimer::singleShot(0, ui->treeWidgetCursorInfo, SLOT(OnCursorPositionChanged()));
  QTimer::singleShot(0, ui->treeWidgetMouseInfo, SLOT(OnMousePositionChanged()));
}

void MainWindow::OnLoadConnectomeMatrix()
{
  DialogLoadConnectome dlg(this);
  if (dlg.exec() == QDialog::Accepted)
  {
    AddScript(QStringList("loadconnectome")
              << QString("%1:lut=%2").arg(dlg.GetCMATFilename()).arg(dlg.GetCTABFilename()) << dlg.GetParcelFilename());
  }
}

void MainWindow::CommandSetVolumeMask(const QStringList &cmd)
{
  LayerMRI* mri = qobject_cast<LayerMRI*>(this->GetActiveLayer("MRI"));
  LayerMRI* mask = qobject_cast<LayerMRI*>(GetLayerCollection("MRI")->GetLayerByName(cmd[1]));
  if (!mask)
  {
    cerr << "Can not find volume by name of " << qPrintable(cmd[1]) << endl;
    return;
  }
  if (mri && mask)
  {
    mri->SetMaskLayer(mask);
  }
}

void MainWindow::OnCycleSurfaceLabel()
{
  LayerSurface* surf = qobject_cast<LayerSurface*>(GetActiveLayer("Surface"));
  if (surf && surf->GetNumberOfLabels() > 0 )
  {
    int n = surf->GetActiveLabelIndex() + 1;
    if (n >= surf->GetNumberOfLabels())
      n = 0;
    surf->SetActiveLabel(n);
  }
}

QList<Layer*> MainWindow::GetSelectedLayers(const QString &layerType)
{
  return ui->widgetAllLayers->GetSelectedLayers(layerType);
}

void MainWindow::OnGoToROI()
{
  LayerROI* roi = (LayerROI*)GetActiveLayer("ROI");
  double pos[3];
  if (roi && roi->GetCentroidPosition(pos))
    SetSlicePosition(pos);
}

void MainWindow::OnGoToSurfaceLabel()
{
  LayerSurface* surf = (LayerSurface*)GetActiveLayer("Surface");
  double pos[3];
  if (surf && surf->GetActiveLabelCentroidPosition(pos))
  {
    SetSlicePosition(pos);
  }
}

void MainWindow::CommandLoadFCD(const QStringList& cmd )
{
  if (cmd.size() < 3)
    return;

  LoadFCD(cmd[1], cmd[2]);
}

void MainWindow::LoadFCD(const QString &subdir, const QString &subject)
{
  LayerFCD* layer = new LayerFCD(m_layerVolumeRef);
  connect( layer->GetWorkerThread(), SIGNAL(Progress(int)), m_statusBar, SLOT(SetProgress(int)));
  connect( layer->GetWorkerThread(), SIGNAL(started()), m_statusBar, SLOT(ShowProgress()));
  connect( layer->GetWorkerThread(), SIGNAL(finished()), m_statusBar, SLOT(HideProgress()));
  layer->SetName(subject);
  layer->SetMRILayerCTAB(m_luts->GetColorTable(0));
  QVariantMap map;
  map["SubjectDir"] = subdir;
  map["Subject"] = subject;
  m_threadIOWorker->LoadFCD( layer, map );
}

void MainWindow::OnLoadFCD()
{
  QString subject_dir = QFileDialog::getExistingDirectory(this, "Select Subject", m_strLastDir);
  if (!subject_dir.isEmpty())
  {
    m_strLastDir = subject_dir;
    QDir dir(subject_dir);
    QString subject = dir.dirName();
    dir.cdUp();
    subject_dir = dir.absolutePath();
    AddScript(QStringList("loadfcd") << subject_dir << subject);
  }
}

void MainWindow::OnCloseFCD()
{
  LayerFCD* layer = (LayerFCD*)GetActiveLayer( "FCD" );
  if ( !layer )
  {
    return;
  }

  GetLayerCollection( "FCD" )->RemoveLayer( layer );
}

QVariant MainWindow::GetSetting(const QString &key)
{
  if (m_settings.contains(key))
    return m_settings[key];
  else
    return QVariant();
}

void MainWindow::SetSetting(const QString &key, const QVariant &value)
{
  m_settings[key] = value;
}

void MainWindow::UpdateSettings()
{
  if (m_dlgPreferences)
  {
    QVariantMap map = m_dlgPreferences->GetSettings();
    QStringList keys = map.keys();
    foreach (QString key, keys)
      m_settings[key] = map[key];
  }
}

void MainWindow::CommandSaveLayer(const QStringList &cmd)
{
  if (cmd.size() < 2)
    return;

  QList<Layer*> layers = GetLayers("MRI");
  bool bOK;
  int nID = cmd[1].toInt(&bOK);
  if (!bOK)
    return;
  foreach (Layer* layer, layers)
  {
    if (layer->GetID() == nID)
    {
      m_threadIOWorker->SaveVolume(layer);
      return;
    }
  }
}

void MainWindow::SaveLayers(const QList<Layer *> &layers)
{
  foreach (Layer* layer, layers)
  {
    m_scripts.append(QStringList("savelayer") << QString::number(layer->GetID()));
  }
}

void MainWindow::OnViewSetCamera()
{
  m_dlgSetCamera->show();
  m_dlgSetCamera->raise();
}

void MainWindow::OnToolSaveCamera()
{
  QString fn = QFileDialog::getSaveFileName(this, "Save Camera", m_strLastDir, "All files (*)");
  if (!fn.isEmpty())
  {
    Json json;
    QVariantMap cam = ui->view3D->GetCamera();
    QFile file(fn);
    if (file.open(QIODevice::WriteOnly))
    {
      file.write(json.encode(cam).toUtf8());
      file.close();
    }
  }
}

void MainWindow::OnToolLoadCamera(const QString& fn_in)
{
  QString fn = fn_in;
  if (fn.isEmpty())
    fn = QFileDialog::getOpenFileName(this, "Load Camera", m_strLastDir, "All files (*)");
  if (!fn.isEmpty())
  {
    QFile file(fn);
    if (file.open(QIODevice::ReadOnly))
    {
      Json json;
      QVariantMap cam = json.decode(file.readAll());
      file.close();
      ui->view3D->SetCamera(cam);
    }
    else
    {
      qWarning() << "Can not open camera file " << fn;
    }
  }
}


