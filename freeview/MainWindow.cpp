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
#include "MainWindow.h"
#include "ui_MainWindow.h"
#include <cstddef>
#include <QtCore>
#include <QtGui>
#include <QFileInfo>
#include "LayerMRI.h"
#include "LayerPropertyMRI.h"
#include "LayerSurface.h"
#include "LayerROI.h"
#include "LayerTrack.h"
#include "LayerDTI.h"
#include "LayerVolumeTrack.h"
#include "LayerODF.h"
#include "LayerCollection.h"
#include "BrushProperty.h"
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
#include "SurfaceLabel.h"
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
#include "VolumeFilterBoundary.h"
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

#if !defined(ARM64)
#include "DialogLineProfile.h"
#include "LayerLineProfile.h"
#endif

#include "DialogLoadConnectome.h"
#include "LayerConnectomeMatrix.h"
#include "LayerFCD.h"
#include "LayerPropertyFCD.h"
#include "DialogSetCamera.h"
#include "DialogThresholdVolume.h"
#include "DialogVolumeSegmentation.h"
#include <QProcessEnvironment>
#include <QJsonDocument>
#include "DialogThresholdFilter.h"
#include "DialogLoadTransform.h"
#include "LayerPropertyTrack.h"
#include "BinaryTreeView.h"
#include "SurfaceAnnotation.h"
#include "Annotation2D.h"
#include "PanelLayer.h"
#include "WindowLayerInfo.h"
#include "DialogTransformSurface.h"
#include <QFileSystemWatcher>
#include <QClipboard>
#include <QDebug>
#ifdef Q_OS_MAC
#include "MacHelper.h"
#endif
#include "DialogMovePoint.h"
#include "VolumeFilterOptimal.h"
#include <QRegularExpression>
#include "MigrationDefs.h"
#include <QDragEnterEvent>
#include <QMimeData>

#if (QT_VERSION >= QT_VERSION_CHECK(5, 0, 0))
#include <QtWidgets>
#endif

#define LAYER_ID_OFFSET 10000
#define SETTING_VERSION 1

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
  m_cmdParser(cmdParser),
  m_bHadError(false)
{
  m_dlgSaveScreenshot = NULL;
  m_dlgPreferences = NULL;
  m_syncFileWatcher = new QFileSystemWatcher(this);
  m_sSyncFilePath = "/tmp/freeview_coord_sync.json";

  m_defaultSettings["no_autoload"] = true;  // default no autoload

  // must create layer collections first before setupui()
  m_layerCollections["MRI"] = new LayerCollection( "MRI", this );
  m_layerCollections["ROI"] = new LayerCollection( "ROI", this );
  m_layerCollections["Surface"] = new LayerCollection( "Surface", this );
  m_layerCollections["PointSet"] = new LayerCollection( "PointSet", this );
  m_layerCollections["Tract"] = new LayerCollection( "Tract", this );
  m_layerCollections["CMAT"] = new LayerCollection("CMAT", this);
  m_layerCollections["FCD"] = new LayerCollection("FCD", this);
  m_layerCollections["ODF"] = new LayerCollection("ODF", this);

  // supplemental layers will not show on control panel
  m_layerCollections["Supplement"] = new LayerCollection( "Supplement", this);
  LayerLandmarks* landmarks = new LayerLandmarks(this);
  m_layerCollections["Supplement"]->AddLayer(landmarks);

  // hidden surface layers
  m_layerCollections["HiddenSurface"] = new LayerCollection("Surface", this);

  m_luts = new LUTDataHolder();
  m_propertyBrush = new BrushProperty();
  m_volumeCropper = new VolumeCropper( this );

  connect(m_volumeCropper, SIGNAL(CropBoundChanged(LayerMRI*)), this, SLOT(RequestRedraw()));
  connect(m_layerCollections["MRI"], SIGNAL(LayerRemoved(Layer*)),
      m_propertyBrush, SLOT(OnLayerRemoved(Layer*)));

  ui->setupUi(this);
  setAcceptDrops(true);

  addAction(ui->actionReloadROI);
  addAction(ui->actionReloadPointSet);
  addAction(ui->actionLockOthers);

  ui->treeWidgetCursorInfo->SetForCursor(true);

  addAction(ui->actionIncreaseOpacity);
  addAction(ui->actionDecreaseOpacity);
  addAction(ui->actionCycleSurfaceLabel);

  addAction(ui->actionResetViewAnterior);
  addAction(ui->actionResetViewPosterior);
  addAction(ui->actionResetViewLeft);
  addAction(ui->actionResetViewRight);
  addAction(ui->actionResetViewSuperior);
  addAction(ui->actionResetViewInferior);
  addAction(ui->actionCopyView);
  addAction(ui->actionDeleteLayer);

  addAction(ui->actionNextLabelPoint);
  addAction(ui->actionShowLabelOutline);

#ifdef DISABLE_LINEPROF
  ui->actionLineProfile->setVisible(false);
#endif

  QWidget *spacer = new QWidget();
  spacer->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Preferred);
  ui->mainToolBar->addWidget(spacer);

  ui->mainToolBar->addAction(ui->actionNeurologicalView);

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
    connect(m_views[i], SIGNAL(CursorLocationClicked()), this, SLOT(On2DCursorClicked()));
  }

  m_dlgTransformSurface = new DialogTransformSurface(this);
  m_dlgTransformSurface->hide();

  connect(m_layerCollections["MRI"], SIGNAL(LayerModified()),this, SLOT(RequestRedraw()));
  connect(m_layerCollections["Supplement"], SIGNAL(LayerModified()),this, SLOT(RequestRedraw()));

  connect(ui->widgetAllLayers, SIGNAL(ToReorderLayers(QList<Layer*>)), this, SLOT(ReorderLayers(QList<Layer*>)));
  for (int i = 0; i < 4; i++)
    connect(ui->widgetAllLayers, SIGNAL(CurrentLayerSelected(Layer*)), m_views[i], SLOT(SetScalarBarLayer(Layer*)));

  connect(m_layerCollections["MRI"], SIGNAL(LayerAdded(Layer*)), m_views[3], SLOT(OnLayerVisibilityChanged()));
  connect(m_layerCollections["MRI"], SIGNAL(LayerRemoved(Layer*)), m_views[3], SLOT(OnLayerVisibilityChanged()));
  connect(m_layerCollections["MRI"], SIGNAL(LayerVisibilityChanged()), m_views[3], SLOT(OnLayerVisibilityChanged()));
  connect(m_layerCollections["Surface"], SIGNAL(LayerAdded(Layer*)), m_views[3], SLOT(OnLayerVisibilityChanged()));
  connect(m_layerCollections["Surface"], SIGNAL(LayerRemoved(Layer*)), m_views[3], SLOT(OnLayerVisibilityChanged()));
  connect(m_layerCollections["Surface"], SIGNAL(LayerVisibilityChanged()), m_views[3], SLOT(OnLayerVisibilityChanged()));

  m_dlgCropVolume = new DialogCropVolume(this);
  m_dlgCropVolume->hide();
  connect(m_volumeCropper, SIGNAL(CropBoundChanged(LayerMRI*)),
          m_dlgCropVolume, SLOT(OnCropBoundChanged(LayerMRI*)));
  connect(m_layerCollections["MRI"], SIGNAL(LayerRemoved(Layer*)),
      m_dlgCropVolume, SLOT(OnLayerRemoved(Layer*)));

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

  m_wndLayerInfo = new WindowLayerInfo(this);
  m_wndLayerInfo->hide();

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
  connect(ui->view3D, SIGNAL(SurfaceVertexClicked(LayerSurface*)),
          this, SLOT(OnSurfaceVertexClicked(LayerSurface*)));
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
  connect(m_layerCollections["MRI"], SIGNAL(LayerAdded(Layer*)), m_wndTimeCourse, SLOT(UpdateUI()));
  connect(m_layerCollections["MRI"], SIGNAL(LayerRemoved(Layer*)), m_wndTimeCourse, SLOT(UpdateUI()));
  connect(m_layerCollections["MRI"], SIGNAL(ActiveLayerChanged(Layer*)), m_wndTimeCourse, SLOT(UpdateUI()));
  connect(m_wndTimeCourse, SIGNAL(OverlayFrameChanged(int)), ui->widgetAllLayers->GetPanel("Surface"), SLOT(SetOverlayFrame(int)));

  m_wndGroupPlot = new WindowGroupPlot(this);
  m_wndGroupPlot->hide();

  m_dlgLabelStats = new DialogLabelStats(this);
  m_dlgLabelStats->hide();
  connect(this, SIGNAL(SlicePositionChanged()), m_dlgLabelStats, SLOT(UpdateStats()), Qt::QueuedConnection);
  connect(m_layerCollections["MRI"], SIGNAL(ActiveLayerChanged(Layer*)), m_dlgLabelStats, SLOT(UpdateStats()), Qt::QueuedConnection);
  connect(m_layerCollections["ROI"], SIGNAL(ActiveLayerChanged(Layer*)), m_dlgLabelStats, SLOT(UpdateStats()), Qt::QueuedConnection);

#if !defined(ARM64)
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
#endif

  m_dlgSetCamera = new DialogSetCamera(this);
  m_dlgSetCamera->hide();

  m_wndTractCluster = new BinaryTreeView(this);
  m_wndTractCluster->setWindowFlags(Qt::Window);
  m_wndTractCluster->setWindowTitle("Tract Cluster");
  m_wndTractCluster->hide();
  connect(m_wndTractCluster, SIGNAL(TreeDataLoaded(QVariantMap)), SLOT(OnTractClusterLoaded(QVariantMap)));

  m_dlgMovePoint = new DialogMovePoint(this);
  m_dlgMovePoint->hide();
  for (int i = 0; i < 3; i++)
    connect(m_views[i], SIGNAL(PointSetPicked(LayerPointSet*, int)), m_dlgMovePoint, SLOT(OnPointSetPicked(LayerPointSet*,int)));

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
      connect( m_layerCollections[keys[i]], SIGNAL(LayersReordered()),
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
        connect( m_layerCollections[keys[i]], SIGNAL(LayersReordered()),
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
    connect(m_layerCollections[keys[i]], SIGNAL(LayersReordered()),
        ui->treeWidgetCursorInfo, SLOT(UpdateAll()), Qt::QueuedConnection);

    connect(m_layerCollections[keys[i]], SIGNAL(LayerAdded(Layer*)),
        ui->treeWidgetMouseInfo, SLOT(UpdateAll()), Qt::QueuedConnection);
    connect(m_layerCollections[keys[i]], SIGNAL(LayerRemoved(Layer*)),
        ui->treeWidgetMouseInfo, SLOT(UpdateAll()), Qt::QueuedConnection);
    connect(m_layerCollections[keys[i]], SIGNAL(LayerMoved(Layer*)),
        ui->treeWidgetMouseInfo, SLOT(UpdateAll()), Qt::QueuedConnection);
    connect(m_layerCollections[keys[i]], SIGNAL(LayerShowInfoChanged()),
        ui->treeWidgetMouseInfo, SLOT(UpdateAll()), Qt::QueuedConnection);
    connect(m_layerCollections[keys[i]], SIGNAL(LayersReordered()),
        ui->treeWidgetMouseInfo, SLOT(UpdateAll()), Qt::QueuedConnection);

    if (keys[i] != "Supplement")
      connect(m_layerCollections[keys[i]], SIGNAL(ActiveLayerChanged(Layer*)),
          this, SLOT(OnActiveLayerChanged(Layer*)), Qt::QueuedConnection);
  }
  connect(m_views[3], SIGNAL(MouseIn()), ui->treeWidgetMouseInfo, SLOT(ShowHeaderText()));
  connect(m_views[3], SIGNAL(MouseOut()), ui->treeWidgetMouseInfo, SLOT(ClearHeaderText()));

  for ( int i = 0; i < 4; i++ )
  {
    connect( this, SIGNAL(SlicePositionChanged(bool)), m_views[i], SLOT(OnSlicePositionChanged(bool)) );
  }

  for (int i = 0; i < 3; i++)
  {
    connect(m_layerCollections["MRI"], SIGNAL(ActiveLayerChanged(Layer*)),
        m_views[i], SLOT(UpdateAnnotation()));
  }
  connect(m_layerCollections["MRI"], SIGNAL(LayerTransformed()),
      m_views[3], SLOT(UpdateBounds()));
  connect(m_layerCollections["Surface"], SIGNAL(ActiveLayerChanged(Layer*)),
      m_views[3], SLOT(UpdateAxesActor()));
  connect(m_layerCollections["MRI"], SIGNAL(ActiveLayerChanged(Layer*)),
      this, SLOT(UpdateLayerInfo(Layer*)), Qt::QueuedConnection);
  connect(m_layerCollections["Surface"], SIGNAL(ActiveLayerChanged(Layer*)),
      this, SLOT(UpdateLayerInfo(Layer*)), Qt::QueuedConnection);

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

  connect(m_layerCollections["Surface"], SIGNAL(LayerAdded(Layer*)), SLOT(UpdateSurfaceContralateralInfo()), Qt::QueuedConnection);
  connect(m_layerCollections["Surface"], SIGNAL(LayerRemoved(Layer*)), SLOT(UpdateSurfaceContralateralInfo()), Qt::QueuedConnection);
  connect(m_layerCollections["HiddenSurface"], SIGNAL(LayerAdded(Layer*)), SLOT(UpdateSurfaceContralateralInfo()), Qt::QueuedConnection);
  connect(m_layerCollections["HiddenSurface"], SIGNAL(LayerRemoved(Layer*)), SLOT(UpdateSurfaceContralateralInfo()), Qt::QueuedConnection);

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

  connect(ui->actionNextLabelPoint, SIGNAL(triggered()), ui->widgetAllLayers->GetPanel("MRI"), SLOT(OnGoToNextPoint()));
  addAction(ui->actionCycleOverlay);
  connect(ui->actionCycleOverlay, SIGNAL(triggered()), SIGNAL(CycleOverlayRequested()));

  addAction(ui->actionCycleAnnotation);
  connect(ui->actionCycleAnnotation, SIGNAL(triggered()), SIGNAL(CycleAnnotationRequested()));

  addAction(ui->actionViewLayerInfo);
  connect(ui->actionViewLayerInfo, SIGNAL(triggered(bool)), SLOT(OnViewLayerInfo()));

  m_widgetFloatControlPanel = new QWidget(this, Qt::Tool | Qt::WindowTitleHint | Qt::CustomizeWindowHint );
  QVBoxLayout* layout = new QVBoxLayout;
  layout->setContentsMargins(0,0,0,0);
  m_widgetFloatControlPanel->setLayout(layout);
  m_widgetFloatControlPanel->hide();
  m_widgetFloatControlPanel->setWindowTitle("Layers");

  m_widgetFloatInfoPanel = new QWidget(this, Qt::Tool | Qt::WindowTitleHint | Qt::CustomizeWindowHint );
  layout = new QVBoxLayout;
  layout->setContentsMargins(0,0,0,0);
  m_widgetFloatInfoPanel->setLayout(layout);
  m_widgetFloatInfoPanel->hide();
  m_widgetFloatInfoPanel->setWindowTitle("Info");

#ifdef Q_OS_MAC
  if (MacHelper::IsDarkMode())
  {
    ui->actionShowCoordinateAnnotation->setIcon(MacHelper::InvertIcon(ui->actionShowCoordinateAnnotation->icon(), QSize(), true));
    ui->actionNeurologicalView->setIcon(MacHelper::InvertIcon(ui->actionNeurologicalView->icon(), QSize(), true));
  }
#endif

  ui->actionTransformSurface->setVisible(false);

  qRegisterMetaType<SurfaceOverlay*>("SurfaceOverlay");
}

MainWindow::~MainWindow()
{
  UpdateSyncIds(false);

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
  m_settingsScreenshot.HideScaleBar = settings.value("ScreenShot/HideScaleBar", true).toBool();
  m_settingsScreenshot.HideCursor = settings.value("ScreenShot/HideCursor", true).toBool();
  m_settingsScreenshot.AutoTrim = settings.value("ScreenShot/AutoTrim", false).toBool();
  m_settings = settings.value("Settings/General").toMap();
  if (!m_settings.contains("SaveCopy"))
  {
    m_settings["SaveCopy"] = true;
  }
  if (!m_settings.contains("BackgroundColor"))
  {
    m_settings["BackgroundColor"] = QColor(Qt::black);
  }
  if (!m_settings.contains("CursorColor"))
  {
    m_settings["CursorColor"] = QColor(Qt::red);
  }
  if (!m_settings.contains("CursorSize"))
  {
    m_settings["CursorSize"] = 5;
    m_settings["CursorSize3D"] = 5;
  }
  if (!m_settings.contains("CursorThickness"))
  {
    m_settings["CursorThickness"] = 1;
    m_settings["CursorThickness3D"] = 1;
  }
  if (!m_settings.contains("AnnotationColor"))
  {
    m_settings["AnnotationColor"] = QColor(Qt::white);
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
  if (!m_settings.contains("AutoReorientView"))
  {
    m_settings["AutoReorientView"] = false;
  }
  if (!m_settings.contains("TextSize"))
  {
    m_settings["TextSize"] = 12;
  }
  if (!m_settings.contains("Precision"))
  {
    m_settings["Precision"] = 2;
  }
  if (m_settings["Version"].toInt() < 1 )
  {
    m_settings["3DAxesFlyMode"] = 4;
  }

  ui->actionDeleteLayer->setVisible(m_settings["AllowDeleteKey"].toBool());

  m_settings["Version"] = SETTING_VERSION;
  if (!m_settings.contains("UseComma"))
    m_settings["UseComma"] = true;

  //  OnPreferences();
  //  m_dlgPreferences->hide();

  for (int i = 0; i < 4; i++)
  {
    m_views[i]->SetBackgroundColor(m_settings["BackgroundColor"].value<QColor>());
    if ( i < 3 )
    {
      ((RenderView2D*)m_views[i])->GetCursor2D()->SetColor(m_settings["CursorColor"].value<QColor>());
      ((RenderView2D*)m_views[i])->GetCursor2D()->SetSize(m_settings["CursorSize"].toInt());
      ((RenderView2D*)m_views[i])->GetCursor2D()->SetThickness(m_settings["CursorThickness"].toInt());
      ((RenderView2D*)m_views[i])->SetAutoScaleText(m_settings["AutoScaleText"].toBool());
      ((RenderView2D*)m_views[i])->SetTextSize(m_settings["TextSize"].toInt());
      ((RenderView2D*)m_views[i])->GetAnnotation2D()->SetColor(m_settings["AnnotationColor"].value<QColor>());
    }
    else
    {
      ((RenderView3D*)m_views[i])->GetCursor3D()->SetColor(m_settings["CursorColor"].value<QColor>());
      ((RenderView3D*)m_views[i])->GetInflatedSurfCursor()->SetColor(m_settings["CursorColor"].value<QColor>());
      ((RenderView3D*)m_views[i])->GetCursor3D()->SetSize(m_settings["CursorSize3D"].toInt());
      ((RenderView3D*)m_views[i])->GetInflatedSurfCursor()->SetSize(m_settings["CursorSize3D"].toInt());
      ((RenderView3D*)m_views[i])->GetCursor3D()->SetThickness(m_settings["CursorThickness3D"].toInt());
      ((RenderView3D*)m_views[i])->GetInflatedSurfCursor()->SetThickness(m_settings["CursorThickness3D"].toInt());
      ((RenderView3D*)m_views[i])->SetAxesFlyMode(m_settings["3DAxesFlyMode"].toInt());
    }
  }
  SyncZoom(m_settings["SyncZoom"].toBool());
  m_term->SetDarkTheme(m_settings["DarkConsole"].toBool());

  QString val = m_settings.value("ShortcutCycleLayer").toString();
  if (!val.isEmpty() && val != "Default")
  {
    DialogPreferences::SetActionShortcut(ui->actionCycleLayer, val);
  }
  val = m_settings.value("ShortcutToggleVolume").toString();
  if (!val.isEmpty() && val != "Default")
  {
    DialogPreferences::SetActionShortcut(ui->actionToggleVolumeVisibility, val);
  }
  val = m_settings.value("ShortcutToggleSurface").toString();
  if (!val.isEmpty() && val != "Default")
  {
    DialogPreferences::SetActionShortcut(ui->actionToggleSurfaceVisibility, val);
  }

#ifdef Q_OS_MAC
  //  this->SetUnifiedTitleAndToolBar(m_settings["MacUnifiedTitleBar"].toBool());
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
    settings.setValue("ScreenShot/HideScaleBar", s.HideScaleBar);
    settings.setValue("ScreenShot/HideCursor", s.HideCursor);
  }
  if (m_dlgPreferences)
  {
    settings.setValue("Settings/General", m_dlgPreferences->GetSettings());
  }
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
  keys.removeOne("Supplement");
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
    QAbstractButton* yesBtn = msgbox.addButton("Quit without Saving", QMessageBox::YesRole);
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

  QList<Layer*> ps_layers = GetLayers("PointSet");
  foreach (Layer* layer, ps_layers)
  {
    if (layer->property("remind_edit").toBool() && !((LayerPointSet*)layer)->IsEdited())
    {
      QMessageBox box(QMessageBox::Question, tr("Unedited Point Set"),
                      tr("Point set %1 has not been edited. Do you still want to close it?").arg(layer->GetName()),
                      QMessageBox::Yes | QMessageBox::Cancel);
      box.setButtonText(QMessageBox::Yes, tr("Close It"));
      box.setDefaultButton(QMessageBox::Cancel);
      if (box.exec() != QMessageBox::Yes)
      {
        event->ignore();
        return;
      }
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
#ifdef Q_OS_LINUX
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
  if (parser->Found( "auto-load-surf"))
  {
    m_defaultSettings["no_autoload"] = false;
  }
  if ( parser->Found( "viewport", &sa ) )
  {
    this->AddScript( QStringList("setviewport") << sa[0]);
  }
  if (parser->Found("layout", &sa))
  {
    int n = sa[0].toInt()-1;
    if (n < 0)
      n = 0;
    else if (n > 3)
      n = 3;
    SetViewLayout(n);
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
      QStringList fn_list;
      if (floatingArgs[i].contains("*"))
      {
        QStringList sublist = floatingArgs[i].split(":");
        QFileInfoList fi_list = QDir().entryInfoList(QStringList(sublist[0]));
        //  qDebug() << fi_list;
      }
      else
      {
        fn_list << floatingArgs[i];
      }
      for (int j = 0; j < fn_list.size(); j++)
      {
        QStringList script = QStringList("loadvolume") << fn_list[j];
        if ( parser->Found( "r" ) )
        {
          script << "r";
        }
        cmds << script;
      }
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

  if (parser->Found("fcd", &sa))
  {
    QStringList script = QStringList("loadfcd") << sa[0] << sa[1];
    if (sa.size() > 2)
      script << sa[2];
    this->AddScript( script );
    bHasVolume = true;
    m_defaultSettings["Smoothed"] = true;
  }

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

  nRepeats = parser->GetNumberOfRepeats( "odf" );
  if (nRepeats > 0 && !bHasVolume)
  {
    QString msg = "Cannot load ODF without loading a matching volume first";
    ShowNonModalMessage("Warning", msg);
    std::cerr << qPrintable(msg) << std::endl;
  }
  else
  {
    for ( int n = 0; n < nRepeats; n++ )
    {
      parser->Found( "odf", &sa, n );
      QStringList script("loadodf");
      script << sa;
      this->AddScript( script );
    }
  }

  if ( parser->Found( "cmat", &sa ) )
  {
    this->AddScript( QStringList("loadconnectome") << sa[0] << sa[1] );
  }

  if (parser->Found("tc", &sa))
  {
    this->AddScript(QStringList("loadtractcluster") << sa[0]);
  }

  if ( parser->Found( "ras", &sa ) )
  {
    bool bOK = true;
    for (int i = 0; i < 3 && bOK; i++)
      sa[i].toDouble(&bOK);
    if ( !bOK )
    {
      std::cerr << "Invalid argument for 'ras'. Arguments must be valid float values.\n";
      return false;
    }
    QStringList script("ras");
    script << sa[0] << sa[1] << sa[2];
    if (sa.size() > 3)
      script << sa[3];
    this->AddScript( script );
  }

  if ( parser->Found( "slice", &sa ) )
  {
    bool bOK = true;
    for (int i = 0; i < 3 && bOK; i++)
      sa[i].toInt(&bOK);
    if ( !bOK )
    {
      std::cerr << "Invalid argument for 'slice'. Arguments must be valid integers.\n";
      return false;
    }

    this->AddScript( QStringList("slice") << sa[0] << sa[1] << sa[2] );
  }

  if ( parser->Found("write-slice-intersection", &sa))
  {
    int start = 0, end = 0;
    bool bOK;
    start = sa[2].toInt(&bOK);
    end = sa[3].toInt(&bOK);
    for (int i = start; i <= end; i++)
    {
      //    slice[n] = i;
      //    this->AddScript(QStringList("slice") << QString::number(slice[0]) << QString::number(slice[1]) << QString::number(slice[2]));
      this->AddScript(QStringList("writesurfaceintersection") << sa[0] << sa[1].replace("%d", "%1").arg(i) << QString::number(i));
    }
  }

  if (parser->Found("view", &sa))
  {
    this->AddScript( QStringList("view") << sa[0]);
  }

  if (parser->Found("neuro-view"))
    SetNeurologicalView(true);

  if ( parser->Found("orthographic"))
    AddScript(QStringList("setorthographic"));

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

  if (parser->Found("fly", &sa))
  {

  }

  if (parser->Found("hide-3d-slices", &sa) )
  {
    ((RenderView3D*)m_views[3])->HideSlices();
  }

  if (parser->Found("hide-x-slice", &sa))
    ((RenderView3D*)m_views[3])->ShowSlice(0, false);
  if (parser->Found("hide-y-slice", &sa))
    ((RenderView3D*)m_views[3])->ShowSlice(1, false);
  if (parser->Found("hide-z-slice", &sa))
    ((RenderView3D*)m_views[3])->ShowSlice(2, false);

  if (parser->Found("hide-3d-frames", &sa) )
  {
    ((RenderView3D*)m_views[3])->SetShowSliceFrames(false);
  }

  if (parser->Found("lineprofile", &sa))
  {
    this->AddScript(QStringList("exportlineprofile") << sa[0]);
  }

  m_bVerbose = parser->Found("verbose");
  m_bContinue = parser->Found("continue");

  if (parser->Found("stdin"))
    m_term->EnableListeningStdin();

  if (parser->Found("subtitle", &sa))
  {
    m_sTitle = sa[0];
    setWindowTitle("FreeView: " + m_sTitle);
  }

  if ( parser->Found( "ss", &sa ) )
  {
    QString mag_factor = "1", auto_trim = "0";
    if (sa.size() > 1)
      mag_factor = sa[1];
    if (sa.size() > 2)
      auto_trim = sa[2];

    this->AddScript( QStringList("screenshot") << sa[0] << mag_factor << auto_trim);

    if (bAutoQuit && !parser->Found("noquit"))
    {
      this->AddScript( QStringList("quit") );
    }
  }

  nRepeats = parser->GetNumberOfRepeats( "prefix" );
  for ( int n = 0; n < nRepeats; n++ )
  {
    parser->Found( "prefix", &sa, n );
    QStringList script("setnameprefix");
    script << sa;
    this->AddScript( script );
  }

  if ( parser->Found("quit"))
    AddScript(QStringList("quit") );

  QFileInfo fi("/tmp");
  if (!fi.isWritable())
  {
    m_sSyncFilePath = QFileInfo(QStandardPaths::locate(QStandardPaths::HomeLocation, "", QStandardPaths::LocateDirectory),
                                ".freeview_coord_sync").absoluteFilePath();
  }

  if (parser->Found("sync", &sa))
  {
    if (sa.size() > 0)
      m_sSyncFilePath = sa[0];
    ui->actionSyncInstances->setChecked(true);
  }

  if (parser->Found("rotate-around-cursor"))
    ((RenderView3D*)m_views[3])->SetFocalPointAtCursor(true);

  if (QFile::exists(m_sSyncFilePath) && QFileInfo(m_sSyncFilePath).lastModified().addDays(1) < QDateTime::currentDateTime())
  {
    QFile file(m_sSyncFilePath);
    file.remove();
  }
  if(parser->Found("sphere-ignore-vg"))    setenv("FV_SPHERE_IGNORE_VG","1",1);
  if(parser->Found("no-sphere-ignore-vg")) setenv("FV_SPHERE_IGNORE_VG","0",1);

  if (!QFile::exists(m_sSyncFilePath))
  {
    QFile file(m_sSyncFilePath);
    if (file.open(QIODevice::WriteOnly))
    {
      file.write(QJsonDocument::fromVariant(QVariantMap()).toJson());
      file.flush();
      file.close();
    }
  }
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
  if (!bBusy && m_scripts.isEmpty() && property("from_cmd").toBool())
  {
    setProperty("from_cmd", false);
    cout << "Commands finished" << endl;
  }
  if ( !bBusy && !m_bScriptRunning && !m_scripts.isEmpty() )
  {
    bool last_one = (m_scripts.size() == 1);
    RunScript();
    if (last_one)
    {
      ui->widgetAllLayers->UpdateWidgets();
    }
  }

  ui->actionShowToolbar->setChecked(ui->mainToolBar->isVisible());

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

  QString type = GetCurrentLayerType();
  if ( nMode == RenderView::IM_ROIEdit || (nMode == RenderView::IM_Navigate && type == "ROI"))
  {
    LayerROI* roi = ( LayerROI* )GetActiveLayer("ROI");
    ui->actionUndo->setEnabled( roi && roi->IsVisible() && roi->HasUndo() );
    ui->actionRedo->setEnabled( roi && roi->IsVisible() && roi->HasRedo() );
  }
  else if ( nMode == RenderView::IM_VoxelEdit || nMode == RenderView::IM_ReconEdit || (nMode == RenderView::IM_Navigate && type == "MRI"))
  {
    LayerMRI* mri = ( LayerMRI* )GetActiveLayer( "MRI");
    ui->actionUndo->setEnabled( mri && mri->IsVisible() && mri->HasUndo() );
    ui->actionRedo->setEnabled( mri && mri->IsVisible() && mri->HasRedo() );
  }
  else if ( nMode == RenderView::IM_PointSetEdit || (nMode == RenderView::IM_Navigate && type == "PointSet"))
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
  LayerTrack* layerTrack  = (LayerTrack*)GetActiveLayer( "Tract");
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
  ui->actionLoadPointSet    ->setEnabled( !bBusy && (layerVolume  || layerSurface));
  ui->actionLoadSurface     ->setEnabled( !bBusy );
  ui->actionLoadTrackVolume ->setEnabled( !bBusy );
  ui->actionLoadTrack       ->setEnabled( !bBusy );
  ui->actionLoadTractCluster->setEnabled( !bBusy );
  ui->actionNewVolume       ->setEnabled( layerVolume );
  ui->actionNewROI          ->setEnabled( layerVolume );
  ui->actionNewPointSet     ->setEnabled( layerVolume || layerSurface );
  ui->actionRepositionSurface->setEnabled( layerSurface );
  ui->actionSmoothSurface   ->setEnabled( layerSurface );
  ui->actionRemoveIntersectionsSurface->setEnabled(layerSurface);
  ui->actionResetView       ->setEnabled( bHasLayer );
  ui->actionResetViewNearestAxis->setEnabled( bHasLayer && ui->view3D->isVisible() );
  ui->actionRotateView90    ->setEnabled( bHasLayer && ui->view3D->isVisible() );
  ui->actionRotateView180   ->setEnabled( bHasLayer && ui->view3D->isVisible() );
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
  ui->actionLoadPatch       ->setEnabled( layerSurface );
  ui->actionSavePatchAs     ->setEnabled( layerSurface );
  ui->actionLoadParameterization->setEnabled( layerSurface );
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
  ui->actionTransformSurface->setEnabled(layerSurface);
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
  ui->actionVolumeFilterBoundary->setEnabled( !bBusy && layerVolume && layerVolume->IsEditable() );
  ui->actionSetCamera->setEnabled(bHasLayer);
  ui->actionSaveCamera->setEnabled(bHasLayer && GetMainView() == ui->view3D);
  ui->actionLoadCamera->setEnabled(bHasLayer && GetMainView() == ui->view3D);

  ui->actionLoadConnectome->setEnabled( !bBusy );
  ui->actionCloseConnectome ->setEnabled( !bBusy && GetActiveLayer( "CMAT"));

  ui->actionLoadFCD->setEnabled( !bBusy );
  ui->actionCloseFCD->setEnabled( !bBusy && GetActiveLayer( "FCD"));

  ui->actionLoadODF->setEnabled( !bBusy && layerVolume );
  ui->actionCloseODF->setEnabled( !bBusy && GetActiveLayer( "ODF"));

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

  if (!layerSurface)
    m_dlgTransformSurface->hide();

  if (!layerVolume)
    m_dlgTransformVolume->hide();

  ui->actionNeurologicalView->setChecked(((RenderView2D*)m_views[0])->GetNeurologicalView());
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
  else if ( cmd == "settrackvolumeframe" )
  {
    CommandSetVolumeTrackFrame(sa);
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
  else if (cmd == "loadtractcluster")
  {
    CommandLoadTractCluster(sa);
  }
  else if ( cmd == "loadodf")
  {
    CommandLoadODF( sa );
  }
  else if ( cmd == "loadfcd")
  {
    CommandLoadFCD( sa );
  }
  else if ( cmd == "loadroi" || sa[0] == "loadlabel" )
  {
    CommandLoadROI( sa );
  }
  else if ( cmd == "gotoroi")
  {
    OnGoToROI(true);
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
  else if ( cmd == "settrackcolor")
  {
    CommandSetTrackColor( sa );
  }
  else if ( cmd == "settrackrender")
  {
    CommandSetTrackRender( sa );
  }
  else if (cmd == "setorthographic")
  {
    m_views[3]->SetParallelProjection(true);
  }
  else if ( cmd == "screenshot" )
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
    double pos[3];
    GetLayerCollection("MRI")->GetSlicePosition(pos);
    m_views[3]->CenterAtWorldPosition(pos);
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
  else if (cmd == "writesurfaceintersection")
  {
    CommandWriteSurfaceIntersection( sa );
  }
  else if ( cmd == "setcolormap" )
  {
    CommandSetColorMap( sa );
  }
  else if ( cmd == "setheatscaleoptions" )
  {
    CommandSetHeatScaleOptions( sa );
  }
  else if ( cmd == "setheatscaleoffset")
  {
    CommandSetHeatScaleOffset( sa );
  }
  else if ( cmd == "setlut" )
  {
    CommandSetLUT( sa );
  }
  else if ( cmd == "setselectedlabels")
  {
    CommandSetSelectedLabels( sa );
  }
  else if ( cmd == "setopacity" )
  {
    CommandSetOpacity( sa );
  }
  else if ( cmd == "setsmoothed")
  {
    CommandSetSmoothed( sa );
  }
  else if ( cmd == "setrgb" )
  {
    CommandSetRgb( sa );
  }
  else if ( cmd == "setdisplayoutline")
  {
    CommandSetLabelOutline(sa);
  }
  else if ( cmd == "setdisplayisosurface" )
  {
    CommandSetDisplayIsoSurface( sa );
  }
  else if ( cmd == "saveisosurface")
  {
    OnSaveIsoSurface(sa.last());
  }
  else if ( cmd == "setisosurfacecolor" )
  {
    CommandSetIsoSurfaceColor( sa );
  }
  else if (cmd == "setisosurfacesmooth")
  {
    CommandSetIsoSurfaceSmooth( sa );
  }
  else if (cmd == "setisosurfaceupsample")
  {
    CommandSetIsoSurfaceUpsample( sa );
  }
  else if (cmd == "setextractallregions")
  {
    CommandSetExtractAllRegions( sa );
  }
  else if ( cmd == "loadisosurfaceregion" )
  {
    CommandLoadIsoSurfaceRegion( sa );
  }
  else if ( cmd == "setsurfaceoverlaymethod" )
  {
    CommandSetSurfaceOverlayMethod( sa );
  }
  else if ( cmd == "setsurfaceoverlaycustom" )
  {
    CommandSetSurfaceOverlayCustom( sa );
  }
  else if (cmd == "setsurfaceoverlaycolormap")
  {
    CommandSetSurfaceOverlayColormap( sa );
  }
  else if ( cmd == "setsurfaceoverlayopacity" )
  {
    CommandSetSurfaceOverlayOpacity( sa );
  }
  else if ( cmd == "setsurfaceoverlayoffset")
  {
    CommandSetSurfaceOverlayOffset(sa);
  }
  else if ( cmd == "setsurfaceoverlayframe")
  {
    CommandSetSurfaceOverlayFrame( sa );
  }
  else if (cmd == "setsurfaceoverlaysmooth")
  {
    CommandSetSurfaceOverlaySmooth( sa );
  }
  else if (cmd == "setsurfaceoverlaymask")
  {
    CommandSetSurfaceOverlayMask( sa );
  }
  else if ( cmd == "setsurfaceoffset" )
  {
    CommandSetSurfaceOffset( sa );
  }
  else if ( cmd == "gotosurfacevertex")
  {
    CommandGoToSurfaceVertex( sa );
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
  else if ( cmd == "setsurfaceopacity" )
  {
    CommandSetSurfaceOpacity( sa );
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
  else if ( cmd == "hidesurfacein3d")
  {
    CommandHideSurfaceIn3D( sa );
  }
  else if ( cmd == "setsurfacevertexcolor" )
  {
    CommandSetSurfaceVertexColor( sa );
  }
  else if ( cmd == "setsurfacelabeloutline" )
  {
    CommandSetSurfaceLabelOutline( sa );
  }
  else if ( cmd == "setsurfacelabelopacity" )
  {
    CommandSetSurfaceLabelOpacity( sa );
  }
  else if (cmd == "setsurfacelabelcolor")
  {
    CommandSetSurfaceLabelColor( sa );
  }
  else if (cmd == "setsurfacelabelthreshold")
  {
    CommandSetSurfaceLabelThreshold( sa );
  }
  else if (cmd == "gotosurfacelabel")
  {
    OnGoToSurfaceLabel(true);
  }
  else if (cmd == "hidesurfacelabel")
  {
    CommandHideSurfaceLabel(sa);
  }
  else if ( cmd == "setsurfaceannotationoutline" )
  {
    CommandSetSurfaceAnnotationOutline( sa );
  }
  else if ( cmd == "loadsurfaceparameterization")
  {
    CommandLoadSurfaceCoordsFromParameterization( sa );
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
  else if (cmd == "linkmri")
  {
    CommandLinkVolume( sa );
  }
  else if ( cmd == "gotolabel" || cmd == "gotostructure")
  {
    CommandGoToLabel( sa );
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
  else if (cmd == "gotocontralateralsurface")
  {
    LayerSurface* surf = reinterpret_cast<LayerSurface*>(sa[1].toULongLong());
    if (surf)
    {
      GetLayerCollection("Surface")->SetActiveLayer(surf);
      GoToContralateralPoint(surf);
    }
  }
  else if (cmd == "setcurrentvertex")
  {
    LayerSurface* surf = qobject_cast<LayerSurface*>(GetActiveLayer("Surface"));
    if (surf)
    {
      bool bOk;
      int n = sa[1].toInt(&bOk);
      if (bOk && n >= 0)
        surf->SetCurrentVertex(n);
    }
  }
  else if (cmd == "reorderlayers")
  {
    CommandReorderLayers(sa);
  }
  else if (cmd == "unloadlayers")
  {
    CommandUnloadLayers(sa);
  }
  else if (cmd == "setactiveframe")
  {
    CommandSetActiveFrame(sa);
  }
  else if (cmd == "setautoadjustframecontrast")
  {
    CommandSetAutoAdjustFrameContrast(sa);
  }
  else if (cmd == "setactivelayer")
  {
    CommandSetActiveLayer(sa);
  }
  else if (cmd == "view")
  {
    if (sa[1] == "left")
      ui->view3D->ResetViewLeft();
    else if (sa[1] == "right")
      ui->view3D->ResetViewRight();
    else if (sa[1] == "anterior")
      ui->view3D->ResetViewAnterior();
    else if (sa[1] == "posterior")
      ui->view3D->ResetViewPosterior();
    else if (sa[1] == "inferior")
      ui->view3D->ResetViewInferior();
    else if (sa[1] == "posterior")
      ui->view3D->ResetViewPosterior();
    else if (sa[1] == "lateral")
      ui->view3D->ResetViewLateral();
    else if (sa[1] == "medial")
      ui->view3D->ResetViewMedial();
  }
  else if (cmd == "resetview")
  {
    OnResetView();
  }
  else if (cmd == "exportlineprofile")
  {
    CommandExportLineProfileThickness(sa);
  }
  else if (cmd == "setnameprefix")
  {
    if (sa.size() > 2)
    {
      QList<Layer*> layers = GetLayers("MRI");
      QString prefix = sa[1];
      for (int j = 2; j < sa.size(); j++)
      {
        foreach (Layer* layer, layers)
        {
          if (layer->GetFileName() == QFileInfo(sa[j]).absoluteFilePath())
          {
            layer->SetName(prefix + "/" + layer->GetName());
            layers.removeOne(layer);
            sa.removeAt(j);
            j--;
          }
        }
      }
    }
  }
  else
  {
    cerr << "Command '" << qPrintable(cmd) << "' was not recognized.\n";
  }
  m_bScriptRunning = false;
}

void MainWindow::ClearScripts()
{
  m_scripts.clear();
  m_bScriptRunning = false;
  setProperty("from_cmd", false);
}

void MainWindow::CommandLoadCommand(const QStringList &sa)
{
  if (sa.size() < 2)
  {
    cerr << "No filename specified.\n";
    return;
  }
  QFile file(sa[1]);
  if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
  {
    cerr << "Can not open for read: " << qPrintable(sa[1]) << ".\n";
    return;
  }
  QStringList lines = QString(file.readAll()).trimmed().split("\n", MD_SkipEmptyParts);
  foreach (QString line, lines)
  {
    if (line.trimmed().indexOf("#") == 0)
      continue;

    QStringList args = line.trimmed().split(QRegularExpression("\\s+"), MD_SkipEmptyParts);
    if (args.size() > 0 &&
        ( args[0].toLower() == "freeview" || args[0].toLower() == "fv"))
    {
      args.removeFirst();
    }
    if (!args.isEmpty() && args.first().at(0) == '-')
    {
      args.prepend("freeview");
      ParseCommand(args.join(" "));
    }
    else
    {
      AddScript(args);
    }
  }

  cout << "Executing commands from " << qPrintable(sa[1]) << endl;
  setProperty("from_cmd", true);
}

void MainWindow::CommandLoadSubject(const QStringList &sa)
{
  QString subject_path = QProcessEnvironment::systemEnvironment().value("SUBJECTS_DIR");
  if (subject_path.isEmpty() || !QDir(subject_path).isReadable())
  {
    cerr << "SUBJECTS_DIR is not set or not accessible. Can not load subject.\n";
    return;
  }
  subject_path += "/" + sa[1];
  QString args = QString("freeview -v "
                         "%1/mri/norm.mgz "
                         "%1/mri/T1.mgz "
                         "%1/mri/brainmask.mgz "
                         "%1/mri/wm.mgz:colormap=heat:visible=0:opacity=0.4 "
                         "%1/mri/aseg.mgz:colormap=lut:opacity=0.22 "
                         "-f %1/surf/lh.white "
                         "%1/surf/rh.white "
                         "%1/surf/lh.pial:edgecolor=red "
                         "%1/surf/rh.pial:edgecolor=red "
                         //                         "%1/surf/lh.orig:edgecolor=green:visible=0 "
                         //                         "%1/surf/rh.orig:edgecolor=green:visible=0 "
                         "%1/surf/lh.inflated:annot=aparc:visible=0 "
                         "%1/surf/rh.inflated:annot=aparc:visible=0 ").arg(subject_path);
  if (QFile::exists(QString("%1/surf/lh.orig.nofix").arg(subject_path)))
      args += QString("%1/surf/lh.orig.nofix:overlay=%1/surf/lh.defect_labels:edgecolor=overlay:overlay_threshold=0.01,100,percentile:visible=0 ").arg(subject_path);
  if (QFile::exists(QString("%1/surf/rh.orig.nofix").arg(subject_path)))
      args += QString("%1/surf/rh.orig.nofix:overlay=%1/surf/rh.defect_labels:edgecolor=overlay:overlay_threshold=0.01,100,percentile:visible=0 ").arg(subject_path);
  args +=  "-viewport coronal ";
  QString control_pt_file = QString("%1/tmp/control.dat").arg(subject_path);
  if (QFile::exists(control_pt_file))
    args += "-c " + control_pt_file;
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
  if (sa.size() < 2)
    return;

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
      tensor_render = "boxoid",
      vector_width = "1",
      vector_norm_th = "0";
  int nSampleMethod = m_nDefaultSampleMethod;
  bool bConform = m_bDefaultConform;
  QString gotoLabelName;
  QVariantMap sup_data;
  QString selected_labels;
  bool bLinked = false;
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
                subOption == "heatscaleoptions" ||
                subOption == "heatscale_option" ||
                subOption == "heatscale_options")
      {
        QStringList script("setheatscaleoptions");
        script << subArgu.split(",");
        m_scripts.insert( 0, script );
      }
      else if ( subOption == "heatscale_offset" )
      {
        m_scripts.insert( 0, QStringList() << "setheatscaleoffset" << subArgu );
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
      else if ( subOption == "vector_width")
      {
        vector_width = subArgu;
        if ( vector_width.isEmpty() )
        {
          cerr << "Missing vector width argument.\n";
        }
      }
      else if ( subOption == "vector_norm_threshold")
      {
        vector_norm_th = subArgu;
        if ( vector_norm_th.isEmpty() )
        {
          cerr << "Missing vector norm threshold argument.\n";
        }
      }
      else if ( subOption == "vector_skip" )
      {
        if ( subArgu.isEmpty() )
        {
          cerr << "Missing vector_skip argument.\n";
        }
        else
          sup_data["VectorSkip"] = subArgu;
      }
      else if ( subOption == "vector_normalize" )
      {
        if ( subArgu.isEmpty() )
        {
          cerr << "Missing vector_normalize argument.\n";
        }
        else
          sup_data["VectorNormalize"] = (subArgu.toLower() == "true" || subArgu.toLower() == "yes" || subArgu == "1");
      }
      else if ( subOption == "vector_scale" )
      {
        if ( subArgu.isEmpty() )
        {
          cerr << "Missing vector_scale argument.\n";
        }
        else
          sup_data["VectorLengthScale"] = subArgu;
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
      else if ( subOption == "reg" || subOption == "transform")
      {
        reg_fn = subArgu;
      }
      else if ( subOption == "sample" || subOption == "resample" || subOption == "interpolation" )
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
        QStringList args = subArgu.split(",", MD_SkipEmptyParts);
        script << args;
        m_scripts.insert( 0, script );
      }
      else if ( subOption == "isosurface_output")
      {
        m_scripts.insert(m_scripts.size()-1, (QStringList("saveisosurface") << subArgu));
      }
      else if ( subOption == "upsample_isosurface")
      {
        m_scripts.insert( 0,  (QStringList("setisosurfaceupsample") << subArgu) );
      }
      else if (subOption == "isosurface_color")
      {
        m_scripts.insert( 0,  (QStringList("setisosurfacecolor") << subArgu) );
      }
      else if (subOption == "isosurface_smooth")
      {
        m_scripts.insert(0, (QStringList("setisosurfacesmooth") << subArgu) );
      }
      else if (subOption == "extract_all_regions")
      {
        m_scripts.insert( 0,  (QStringList("setextractallregions") << subArgu) );
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
      else if ( subOption == "link" || subOption == "linked")
      {
        bLinked = true;
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
      else if (subOption == "rgb" || subOption == "RGB")
      {
        m_scripts.insert(0, QStringList("setrgb") << subArgu);
      }
      else if (subOption == "id")
      {
        sup_data["ID"] = subArgu.toInt();
      }
      else if ( subOption == "frame")
      {
        m_scripts.insert( 0, QStringList() << "setactiveframe" << subArgu );
      }
      else if ( subOption == "auto_adjust_frame_contrast" )
      {
        m_scripts.insert( 0, QStringList() << "setautoadjustframecontrast" << subArgu );
      }
      else if (subOption == "ignore_header")
      {
        sup_data["IgnoreHeader"] = true;
      }
      else if (subOption == "binary_color")
      {
        QColor color = ParseColorInput( subArgu );
        if ( color.isValid() )
          sup_data["BinaryColor"] = color;
        else
          cerr << "Unrecognized color input for :binary_color.\n";
      }
      else if (subOption == "select_label")
      {
        selected_labels = subArgu;
      }
      else if (!subOption.isEmpty())
      {
        cerr << "Unrecognized sub-option flag '" << strg.toLatin1().constData() << "'.\n";
        return;
      }
    }
    else
    {
      cerr << "Unrecognized sub-option flag '" << strg.toLatin1().constData() << "'.\n";
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

    if (colormap == "lut" && !selected_labels.isEmpty())
    {
      m_scripts.insert(1, QStringList("setselectedlabels") << selected_labels);
    }
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
                                                            vector_inversion <<
                                                            vector_width;
    m_scripts.insert( 0, script );
  }
  else if ( !vector_display.isEmpty() && vector_display != "no" )
  {
    QStringList script = QStringList("setdisplayvector") <<
                                                            vector_display <<
                                                            vector_render <<
                                                            vector_inversion <<
                                                            vector_width <<
                                                            vector_norm_th << "new";
    m_scripts.insert( 0, script );
  }

  if (bLinked)
    m_scripts.insert(0, QStringList("linkmri") << "1" );

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
  int nColorMap = LayerPropertyMRI::Grayscale;
  int nColorMapScale = LayerPropertyMRI::Grayscale;
  QString strg = sa[1];
  if ( strg == "heat" || strg == "heatscale" )
  {
    nColorMap = LayerPropertyMRI::Heat;
  }
  else if ( strg == "jet" || strg == "jetscale" )
  {
    nColorMap = LayerPropertyMRI::Jet;
  }
  else if ( strg == "turbo" || strg == "turboscale" )
  {
    nColorMap = LayerPropertyMRI::Turbo;
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
  else if ( strg == "binary" )
  {
    nColorMap = LayerPropertyMRI::Binary;
  }
  else if ( strg != "grayscale" )
  {
    cerr << "Unrecognized colormap name '" << strg.toLatin1().constData() << "'.\n";
  }

  QList<double> pars;
  if (sa.size() > 2)
  {
    strg = sa[2];
    bool bOK;
    strg.toDouble((&bOK));
    int nStart = 3;
    if (bOK)
    {
      nColorMapScale = nColorMap;
      nStart = 2;
    }
    else if ( strg == "heatscale" )
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

    for ( int i = nStart; i < sa.size(); i++ )
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
  }

  SetVolumeColorMap( nColorMap, nColorMapScale, pars );
}

void MainWindow::CommandSetSelectedLabels(const QStringList &cmd)
{
  if ( GetLayerCollection( "MRI" )->GetActiveLayer() )
  {
    LayerPropertyMRI* p = ( (LayerMRI*)GetLayerCollection( "MRI" )->GetActiveLayer() )->GetProperty();
    QStringList list = cmd[1].split(",");
    p->SetUnselectAllLabels();
    foreach (QString str, list)
    {
      p->SetSelectLabel(str.toInt(), true);
    }
    emit RefreshLookUpTableRequested();
  }
}

void MainWindow::CommandSetHeatScaleOptions( const QStringList& sa )
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

void MainWindow::CommandSetHeatScaleOffset(const QStringList &sa)
{
  if ( GetLayerCollection( "MRI" )->GetActiveLayer() )
  {
    LayerPropertyMRI* p = ( (LayerMRI*)GetLayerCollection( "MRI" )->GetActiveLayer() )->GetProperty();
    double dValue;
    if (sa[1].toLower() == "mean")
    {
      dValue = p->GetFrameMeanValue();
    }
    else
    {
      bool bOK;
      dValue = sa[1].toDouble(&bOK);
      if ( !bOK )
      {
        cerr << "Heatscale offset value is not valid.\n";
        return;
      }
    }
    p->SetHeatScaleOffset(dValue);
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
  if ( cmd.size() > 2 )
  {
    LayerCollection* lc = GetLayerCollection( cmd[1] );
    if ( lc && !lc->IsEmpty() )
    {
      if ( cmd[2] == "1" || cmd[2].toLower() == "true" )
        lc->GetActiveLayer()->SetVisible(true);
      else if (cmd[2] == "0" || cmd[2].toLower() == "false" )
        lc->GetActiveLayer()->SetVisible(false);
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

void MainWindow::CommandSetRgb(const QStringList &cmd)
{
  QString stemp = cmd[1].toLower();
  if ( stemp == "yes"|| stemp == "true" || stemp == "1" || stemp == "on")
  {
    LayerMRI* mri = (LayerMRI*)GetLayerCollection( "MRI" )->GetActiveLayer();
    if ( mri && mri->GetNumberOfFrames() == 3)
    {
      mri->GetProperty()->SetDisplayRGB(true);
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
        else if ( cmd[2].toLower() == "direction" || cmd[2].toLower() == "directional" )
        {
          mri->GetProperty()->SetVectorRepresentation( LayerPropertyMRI::VR_Direction_Line );
        }
        else if ( cmd[2].toLower() == "bar" )
        {
          mri->GetProperty()->SetVectorRepresentation( LayerPropertyMRI::VR_Bar );
        }
        else
        {
          cerr << "Unrecognized argument '" << cmd[2].toLatin1().constData() << "' for vector rendering.\n";
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
            cerr << "Unknown inversion flag '" << cmd[2].toLatin1().constData() << "'.\n";
          }
        }

        bool ok;
        double val = cmd[4].toDouble(&ok);
        if (ok)
        {
          mri->GetProperty()->SetVectorLineWidth(val);
        }
        else
        {
          cerr << "Unknown vector width value '" << cmd[4].toLatin1().constData() << "'.\n";
        }

        val = cmd[5].toDouble(&ok);
        if (ok)
        {
          mri->GetProperty()->SetVectorNormThreshold(val);
        }
        else
        {
          cerr << "Unknown vector norm threshold value '" << cmd[5].toLatin1().constData() << "'.\n";
        }

        if (val == 1 && cmd.size() > 6)
        {
          QList<Layer*> list = GetLayers("MRI");
          foreach (Layer* layer, list)
          {
            LayerMRI* mlayer = qobject_cast<LayerMRI*>(layer);
            if (mlayer != mri && mlayer->GetProperty()->GetDisplayVector())
            {
              mri->GetProperty()->SetVectorDisplayScale(mlayer->GetProperty()->GetVectorDisplayScale());
              mri->GetProperty()->SetVectorLineWidth(mlayer->GetProperty()->GetVectorLineWidth());
              mri->GetProperty()->SetNormalizeVector(mlayer->GetProperty()->GetNormalizeVector());
              mri->GetProperty()->SetVectorRepresentation(mlayer->GetProperty()->GetVectorRepresentation());
              break;
            }
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
          cerr << "Unrecognized argument '" << cmd[2].toLatin1().constData() << "' for tensor rendering.\n";
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
            cerr << "Unknown inversion flag '" << cmd[3].toLatin1().constData() << "'.\n";
          }
        }

        bool ok;
        double val = cmd[4].toDouble(&ok);
        if (ok)
        {
          mri->GetProperty()->SetVectorLineWidth(val);
        }
        else
        {
          cerr << "Unknown vector width value '" << cmd[4].toLatin1().constData() << "'.\n";
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
  QString val_strg = (sa.size()>2?sa[2]:sa[1]);
  bool bOK;
  double dValue = val_strg.toDouble(&bOK);
  if ( !bOK )
  {
    cerr << "Opacity value is not valid.\n";
    return;
  }
  if (sa.size() > 2)
  {
    LayerCollection* lc = NULL;
    QString type = sa[1].toLower();
    if (type == "mri")
    {
      lc = GetLayerCollection("MRI");
      if (lc && lc->GetActiveLayer())
        ((LayerMRI*)lc->GetActiveLayer())->GetProperty()->SetOpacity(dValue);
    }
    else if (type == "surface")
    {
      lc = GetLayerCollection("Surface");
      if (lc && lc->GetActiveLayer())
        ((LayerSurface*)lc->GetActiveLayer())->GetProperty()->SetOpacity(dValue);
    }
    else if (type == "roi")
    {
      lc = GetLayerCollection("ROI");
      if (lc && lc->GetActiveLayer())
        ((LayerROI*)lc->GetActiveLayer())->GetProperty()->SetOpacity(dValue);
    }
    else if (type == "pointset")
    {
      lc = GetLayerCollection("PointSet");
      if (lc && lc->GetActiveLayer())
        ((LayerPointSet*)lc->GetActiveLayer())->GetProperty()->SetOpacity(dValue);
    }
    else if (type == "tract")
    {
      lc = GetLayerCollection("Tract");
      if (lc && lc->GetActiveLayer())
        ((LayerTrack*)lc->GetActiveLayer())->GetProperty()->SetOpacity(dValue);
    }
  }
  else
  {
    LayerMRI* mri = (LayerMRI*)GetLayerCollection( "MRI" )->GetActiveLayer();
    if ( mri )
    {
      mri->GetProperty()->SetOpacity( dValue );
    }
  }
}

void MainWindow::CommandSetActiveFrame( const QStringList& sa )
{
  LayerMRI* mri = (LayerMRI*)GetLayerCollection( "MRI" )->GetActiveLayer();
  if ( mri )
  {
    bool bOK;
    int val = sa[1].toInt(&bOK);
    if ( bOK && val >= 0 && val < mri->GetNumberOfFrames())
    {
      mri->SetActiveFrame(val);
    }
    else
    {
      cerr << "Frame value is not valid.\n";
    }
  }
}

void MainWindow::CommandSetAutoAdjustFrameContrast( const QStringList& sa )
{
  LayerMRI* mri = (LayerMRI*)GetLayerCollection( "MRI" )->GetActiveLayer();
  if ( mri )
  {
    mri->GetProperty()->SetAutoAdjustFrameLevel(sa[1].toLower() == "true" || sa[1] == "1");
  }
}


void MainWindow::CommandSetDisplayIsoSurface( const QStringList& sa_in )
{
  LayerMRI* mri = (LayerMRI*)GetLayerCollection( "MRI" )->GetActiveLayer();
  if ( mri )
  {
    bool bOK;
    double dValue;
    QStringList sa = sa_in;
    sa.removeFirst();
    if ( sa_in.size() > 1 )
    {
      dValue = sa_in[1].toDouble(&bOK);
      if ( bOK )
      {
        mri->GetProperty()->SetContourMinThreshold( dValue );
        sa.removeFirst();
      }
    }
    if ( sa_in.size() > 2 )
    {
      dValue = sa_in[2].toDouble(&bOK);
      if ( bOK )
      {
        mri->GetProperty()->SetContourMaxThreshold( dValue );
        sa.removeFirst();
      }
    }
    for (int i = 0; i < sa.size(); i++)
    {
      if (sa[i] == "voxelize")
        mri->GetProperty()->SetShowVoxelizedContour(true);
      else if (sa[i] != "on" || sa[i] != "1")
        cerr << "Unrecognized option(s) for isosurface";
    }
    connect(mri, SIGNAL(IsoSurfaceUpdating()), SLOT(SetProcessing()));
    connect(mri, SIGNAL(IsoSurfaceUpdated()), SLOT(SetProcessingFinished()));
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
      cerr << "Invalid color name or value " << cmd[1].toLatin1().constData() << ".\n";
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

void MainWindow::CommandSetIsoSurfaceSmooth(const QStringList &cmd)
{
  LayerMRI* mri = (LayerMRI*)GetLayerCollection( "MRI" )->GetActiveLayer();
  if ( mri )
  {
    bool bOk;
    int nIterations = cmd[1].toInt(&bOk);
    if (bOk && nIterations > 0)
    {
      mri->GetProperty()->SetContourSmoothIterations(nIterations);
    }
  }
}

void MainWindow::CommandSetExtractAllRegions(const QStringList &cmd)
{
  LayerMRI* mri = (LayerMRI*)GetLayerCollection( "MRI" )->GetActiveLayer();
  if ( mri )
  {
    if (cmd[1].toLower() == "off" || cmd[1].toLower() == "false" || cmd[1].toLower() == "0")
    {
      mri->GetProperty()->SetContourExtractAllRegions(false);
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
        cerr << "Can not load surfacer region(s) from " << sa[1].toLatin1().constData() << ".\n";
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
        else if ( strg.left( n ).toLower() == "reg" ||
                  strg.left( n ).toLower() == "transform")
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

  QStringList list = sa[1].split(":");
  if (list.size() > 1)
  {
    QStringList sublist = list[1].split("=");
    if (sublist.size() > 1 && sublist[0] == "frame")
      m_scripts.insert(0, QStringList("settrackvolumeframe") << sublist[1]);
  }
  this->LoadVolumeTrackFile(list[0], bResample);
}

void MainWindow::CommandSetVolumeTrackFrame(const QStringList &cmd)
{
  Layer* layer = GetActiveLayer("MRI");
  if (layer && layer->IsTypeOf("VolumeTrack"))
  {
    LayerVolumeTrack* vt = (LayerVolumeTrack*)layer;
    QStringList frames = cmd[1].split(",");
    vt->ShowAllLabels(false);
    for (int i = 0; i < frames.size(); i++)
    {
      int nFrame = frames[i].toInt();
      if (nFrame >= 0)
        vt->SetFrameVisible(nFrame, true);
    }
    emit RefreshLookUpTableRequested();
  }
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
      cerr << "Can not load look up table " << lut.toLatin1().constData() << ".\n";
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
        cerr << "Unrecognized sub-option flag '" << strg.toLatin1().constData() << "'.\n";
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
  QVariantMap args;
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
        QColor color = ParseColorInput(argu);
        if (color.isValid())
          args["color"] = color;
      }
      else if (option == "opacity")
      {
        args["opacity"] = argu.toDouble();
      }
      else if (option == "threshold")
      {
        args["threshold"] = argu.toDouble();
      }
      else if (option == "id")
      {
        args["id"] = argu.toInt();
      }
      else if (option == "centroid")
      {
        AddScript(QStringList("gotoroi"));
        AddScript(QStringList("center"));
      }
      else if (option == "name")
      {
        args["name"] = argu;
      }
      else
      {
        cerr << "Unrecognized sub-option flag '" << strg.toLatin1().constData() << "'.\n";
      }
    }
  }

  LoadROIFile( fn, ref, args );
}

void MainWindow::CommandLoadTrack(const QStringList &cmd)
{
  QStringList list = cmd[1].split(":");
  QString fn = list[0];
  LoadTrackFile( fn );
  if (list.size() > 1)
  {
    for (int i = 1; i < list.size(); i++)
    {
      QStringList sublist = list[i].split("=");
      if (sublist.size() > 1)
      {
        if (sublist[0] == "color")
          m_scripts.insert(0, QStringList("settrackcolor") << sublist[1]);
        else if (sublist[0] == "render")
          m_scripts.insert(0, QStringList("settrackrender") << sublist[1]);
        else
          cerr << "Unrecognized sub-option flag '" << sublist[0].toLatin1().constData() << "'.\n";
      }
    }
  }
}

void MainWindow::CommandSetTrackColor(const QStringList &cmd)
{
  LayerTrack* layer = (LayerTrack*)GetActiveLayer("Tract");
  if (layer && cmd.size() > 1)
  {
    QColor color = ParseColorInput(cmd[1]);
    if (color.isValid())
    {
      layer->GetProperty()->SetColorCode(LayerPropertyTrack::SolidColor);
      layer->GetProperty()->SetSolidColor(color);
    }
  }
}

void MainWindow::CommandSetTrackRender(const QStringList &cmd)
{
  LayerTrack* layer = (LayerTrack*)GetActiveLayer("Tract");
  if (layer && cmd.size() > 1)
  {
    if (cmd[1] == "tube" || cmd[1] == "tubes")
      layer->GetProperty()->SetRenderRep(LayerPropertyTrack::Tube);
  }
}

void MainWindow::CommandLoadSurface( const QStringList& cmd )
{
  QStringList rawoverlay_list = cmd[1].split("overlay=", MD_SkipEmptyParts, Qt::CaseInsensitive);
  QStringList overlay_list;
  for (int i = 1; i < rawoverlay_list.size(); i++)
  {
    QStringList sublist = rawoverlay_list[i].split("correlation=", MD_SkipEmptyParts, Qt::CaseInsensitive);
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
  QVariantMap sup_options;
  valid_overlay_options << "overlay_reg" << "overlay_method" << "overlay_threshold" << "overlay_color"
                        << "overlay_rh" << "overlay_opacity" << "overlay_frame" << "overlay_smooth" << "overlay_custom"
                        << "overlay_mask" << "overlay_offset";
  bool bNoAutoLoad = m_defaultSettings["no_autoload"].toBool();
  for (int nOverlay = 0; nOverlay < overlay_list.size(); nOverlay++)
  {
    QStringList sa_fn = overlay_list[nOverlay].split(":");
    if (nOverlay == 0)    // first one is not overlay file but actually surface file
      surface_fn = sa_fn[0];
    bool bLoadAll = false;
    //    bool bLabelOutline = false;
    //    QString labelColor;
    QString overlay_reg;
    QString overlay_opacity;
    QString overlay_frame;
    QString overlay_smooth_steps;
    QString overlay_method;
    QStringList overlay_color;
    QStringList overlay_thresholds;
    QStringList overlay_custom;
    QStringList overlay_mask;
    QString overlay_offset;
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
          overlay_thresholds = subArgu.split(",", MD_SkipEmptyParts);
        else if (subOption == "overlay_rh" && (subArgu == "1" || subArgu == "true"))
          bSecondHalfData = true;
        else if (subOption == "overlay_opacity")
          overlay_opacity = subArgu;
        else if (subOption == "overlay_color")
          overlay_color = subArgu.split(",", MD_SkipEmptyParts);
        else if (subOption == "overlay_frame")
          overlay_frame = subArgu;
        else if (subOption == "overlay_smooth")
          overlay_smooth_steps = subArgu;
        else if (subOption == "overlay_custom")
          overlay_custom = subArgu.split(",", MD_SkipEmptyParts);
        else if (subOption == "overlay_mask")
          overlay_mask = subArgu.split(",", MD_SkipEmptyParts);
        else if (subOption == "overlay_offset")
          overlay_offset = subArgu;
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
        else if (subOption == "opacity")
        {
          m_scripts.insert( 0, QStringList("setsurfaceopacity") << subArgu );
        }
        else if ( subOption == "id")
        {
          bool ok;
          subArgu.toInt(&ok);
          if (ok)
            sup_options["ID"] = subArgu.toInt();
        }
        else if ( subOption == "edgecolor" || subOption == "edge_color")
        {
          m_scripts.insert( 0, QStringList("setsurfaceedgecolor") << subArgu );
        }
        else if( subOption == "affinexfm" || subOption == "reg")
        {
          // The LTA can point in either direction as MRISltaMultiply()
          // will determine the right direction if it can
          sup_options["affinexform_filename"] = subArgu;
        }
        else if ( subOption == "edgethickness"|| subOption == "edge_thickness" )
        {
          m_scripts.insert( 0, QStringList("setsurfaceedgethickness") << subArgu );
        }
        else if ( subOption == "vertex" )
        {
          m_scripts.insert( 0, QStringList("displaysurfacevertex") << subArgu);
        }
        else if ( subOption == "current_vertex")
        {
          m_scripts.insert(0, QStringList("setcurrentvertex") << subArgu);
        }
        else if ( subOption == "hide_in_3d")
        {
          m_scripts.insert( 0, QStringList("hidesurfacein3d") << subArgu);
        }
        else if ( subOption == "no_auto_load" || subOption == "no_autoload")
        {
          bNoAutoLoad = true;
        }
        else if ( subOption == "vertexcolor" || subOption == "vertex_color" )
        {
          m_scripts.insert( 0, QStringList("setsurfacevertexcolor") << subArgu );
        }
        else if ( subOption == "curv" || subOption == "curvature" )
        {
          m_scripts.insert( 0, QStringList("loadsurfacecurvature") << subArgu );
        }
        else if ( subOption == "curvature_method" || subOption == "curvature_map" || subOption == "curvature_setting")
        {
          m_scripts.insert(0, QStringList("setsurfacecurvaturemap") << subArgu);
        }
        else if ( subOption == "overlay" || subOption == "correlation" ||
                  subOption == "mrisp" || subOption == "parameterization_overlay")
        {
          // add script to load surface overlay files
          QStringList script("loadsurfaceoverlay");
          script << subArgu;

          script << overlay_reg;
          if (subOption == "correlation")
            script << "correlation";
          else if (subOption == "mrisp" || subOption == "parameterization_overlay")
            script << "mrisp";
          else
            script << "n/a";

          if (bSecondHalfData)
            script << "rh";
          m_scripts.insert( 0, script );

          if (overlay_method.isEmpty())
            overlay_method = "linearopaque";
          if (overlay_method != "linearopaque" || !overlay_thresholds.isEmpty())
          {
            script = QStringList("setsurfaceoverlaymethod") << overlay_method;
            if (!overlay_thresholds.isEmpty())
                 script << overlay_thresholds.join(",");
            // insert right AFTER loadsurfaceoverlay command
            m_scripts.insert( 1, script );
          }

          if (!overlay_custom.isEmpty())
            m_scripts.insert(1, QStringList("setsurfaceoverlaycustom") << overlay_custom);

          if (!overlay_opacity.isEmpty())
            m_scripts.insert(1, QStringList("setsurfaceoverlayopacity") << overlay_opacity);

          if (!overlay_color.isEmpty())
            m_scripts.insert(1, QStringList("setsurfaceoverlaycolormap") << overlay_color);

          if (!overlay_frame.isEmpty())
            m_scripts.insert(1, QStringList("setsurfaceoverlayframe") << overlay_frame);

          if (!overlay_smooth_steps.isEmpty())
            m_scripts.insert(1, QStringList("setsurfaceoverlaysmooth") << overlay_smooth_steps);

          if (!overlay_mask.isEmpty())
            m_scripts.insert(1, QStringList("setsurfaceoverlaymask") << overlay_mask);

          if (!overlay_offset.isEmpty())
            m_scripts.insert(1, QStringList("setsurfaceoverlayoffset") << overlay_offset);
        }
        else if ( subOption == "mrisps" )
        {
          m_scripts.insert( 0, QStringList("loadsurfaceparameterization") << subArgu );
        }
        else if ( subOption == "annot" || subOption == "annotation" || subOption == "aparc" )
        {
          // add script to load surface annotation files
          QStringList annot_fns = subArgu.split(",");
          for ( int i = annot_fns.size()-1; i >= 0; i-- )
          {
            m_scripts.insert( 0, QStringList("loadsurfaceannotation") << annot_fns[i] );
          }
        }
        else if ( subOption == "annot_outline" || subOption == "annotation_outline" || subOption == "aparc_outline")
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
        else if ( subOption == "label_opacity" || subOption == "labelopacity")
        {
          if (!subArgu.isEmpty())
          {
            for (int i = 0; i < m_scripts.size(); i++)
            {
              if (m_scripts[i][0] == "loadsurfacelabel")
              {
                m_scripts.insert(i+1, QStringList("setsurfacelabelopacity") << subArgu);
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
        else if (subOption == "label_threshold" || subOption == "labelthreshold")
        {
          if (!subArgu.isEmpty())
          {
            for (int i = 0; i < m_scripts.size(); i++)
            {
              if (m_scripts[i][0] == "loadsurfacelabel")
              {
                m_scripts.insert(i+1, QStringList("setsurfacelabelthreshold") << subArgu);
                break;
              }
            }
          }
        }
        else if (subOption == "label_centroid" || subOption == "labelcentroid")
        {
          if (!subArgu.isEmpty())
          {
            for (int i = 0; i < m_scripts.size(); i++)
            {
              if (m_scripts[i][0] == "loadsurfacelabel")
              {
                m_scripts.insert(i+1, QStringList("gotosurfacelabel"));
                break;
              }
            }
          }
        }
        else if (subOption == "label_visible")
        {
          if (subArgu == "0" || subArgu == "false" || subArgu == "no")
          {
            for (int i = 0; i < m_scripts.size(); i++)
            {
              if (m_scripts[i][0] == "loadsurfacelabel")
              {
                m_scripts.insert(i+1, QStringList("hidesurfacelabel"));
                break;
              }
            }
          }
        }
        else if (subOption == "annot_zorder" || subOption == "annotation_zorder" )
        {
          sup_options["ZOrderAnnotation"] = subArgu;
        }
        else if (subOption == "label_zorder")
        {
          sup_options["ZOrderLabel"] = subArgu;
        }
        else if (subOption == "overlay_zorder" )
        {
          sup_options["ZOrderOverlay"] = subArgu;
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
          if ( !subArgu.contains( "/" ) )
          {
            subArgu = "./" + subArgu;
          }
          fn_patch = QFileInfo( subArgu ).absoluteFilePath();
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
          sup_files = subArgu.split(",",  MD_SkipEmptyParts);
        }
        else if (subOption == "goto")
        {
          m_scripts.insert(0, QStringList("gotosurfacevertex") << subArgu);
        }
        else if (subOption == "sphere")
        {
          sup_options["sphere"] = subArgu;
        }
        else if ( subOption == "ignore_vg" || subOption == "ignore_volume_geometry")
        {
          if ( subArgu.toLower() == "true" || subArgu.toLower() == "yes" || subArgu == "1")
          {
            sup_options["ignore_vg"] = true;
          }
        }
        else if (subOption == "no_shading")
        {
          if ( subArgu.toLower() == "true" || subArgu.toLower() == "yes" || subArgu == "1")
          {
            sup_options["no_shading"] = true;
          }
        }
        else if ( !valid_overlay_options.contains(subOption) )
        {
          cerr << "Unrecognized sub-option flag '" << subOption.toLatin1().constData() << "'.\n";
          return;
        }
      }
    }
    if (bLoadAll)
    {
      sup_files << "white" << "inflated" << "pial" << "orig";
    }
  }
  if (bNoAutoLoad)
    sup_options["no_autoload"] = true;

  LoadSurfaceFile( surface_fn, fn_patch, fn_target, sup_files, sup_options );
}

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

void MainWindow::CommandSetSurfaceLabelOpacity(const QStringList &cmd)
{
  LayerSurface* surf = (LayerSurface*)GetLayerCollection( "Surface" )->GetActiveLayer();
  if ( surf )
  {
    bool ok;
    cmd[1].toDouble(&ok);
    if (ok && surf->GetActiveLabel())
    {
      surf->GetActiveLabel()->SetOpacity(cmd[1].toDouble());
    }
  }
}

void MainWindow::CommandSetSurfaceLabelThreshold(const QStringList &cmd)
{
  LayerSurface* surf = (LayerSurface*)GetLayerCollection( "Surface" )->GetActiveLayer();
  if ( surf )
  {
    bool ok;
    cmd[1].toDouble(&ok);
    if (ok && surf->GetActiveLabel())
    {
      surf->GetActiveLabel()->SetThreshold(cmd[1].toDouble());
    }
  }
}

void MainWindow::CommandHideSurfaceLabel(const QStringList &cmd)
{
  Q_UNUSED(cmd);
  LayerSurface* surf = (LayerSurface*)GetLayerCollection( "Surface" )->GetActiveLayer();
  if ( surf && surf->GetActiveLabel())
  {
    surf->GetActiveLabel()->SetVisible(false);
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

void MainWindow::CommandSetSurfaceOverlayOffset(const QStringList &cmd)
{
  LayerSurface* surf = (LayerSurface*)GetLayerCollection( "Surface" )->GetActiveLayer();
  if ( surf )
  {
    SurfaceOverlay* overlay = surf->GetActiveOverlay();
    if ( overlay )
    {
      bool ok;
      double val = cmd[1].toDouble(&ok);
      if (ok)
      {
        overlay->GetProperty()->SetOffset(val);
        surf->UpdateOverlay(true);
        overlay->EmitDataUpdated();
      }
      else
      {
        cerr << "Invalid input for overlay offset.\n";
      }
    }
  }
}

void MainWindow::CommandSetSurfaceOverlayFrame(const QStringList &cmd)
{
  LayerSurface* surf = (LayerSurface*)GetLayerCollection( "Surface" )->GetActiveLayer();
  if ( surf )
  {
    SurfaceOverlay* overlay = surf->GetActiveOverlay();
    if ( overlay )
    {
      bool ok;
      int frame = cmd[1].toInt(&ok);
      if (ok)
      {
        overlay->SetActiveFrame(frame);
        surf->UpdateOverlay(true);
        overlay->EmitDataUpdated();
      }
      else
      {
        cerr << "Invalid input for overlay frame.\n";
      }
    }
  }
}

void MainWindow::CommandSetSurfaceOverlaySmooth(const QStringList &cmd)
{
  LayerSurface* surf = (LayerSurface*)GetLayerCollection( "Surface" )->GetActiveLayer();
  if ( surf )
  {
    SurfaceOverlay* overlay = surf->GetActiveOverlay();
    if ( overlay )
    {
      bool ok;
      int steps = cmd[1].toInt(&ok);
      if (ok && steps > 0)
      {
        overlay->GetProperty()->SetSmooth(true);
        overlay->GetProperty()->SetSmoothSteps(steps);
        overlay->UpdateSmooth();
      }
      else
      {
        if (!ok)
          cerr << "Invalid input for overlay smoothing.\n";
      }
    }
  }
}

void MainWindow::CommandSetSurfaceOverlayMask(const QStringList &cmd)
{
  LayerSurface* surf = (LayerSurface*)GetLayerCollection( "Surface" )->GetActiveLayer();
  if ( surf )
  {
    SurfaceOverlay* overlay = surf->GetActiveOverlay();
    if ( overlay )
    {
      if (cmd.size() > 2 && (cmd[2].toLower() == "invert" || cmd[2].toLower() == "inverse"))
        overlay->GetProperty()->SetMaskInverse(true);
      emit OverlayMaskRequested(cmd[1]);
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
      QStringList methods = cmd[1].split(",");
      if (methods.contains("linear", Qt::CaseInsensitive))
      {
        nMethod = SurfaceOverlayProperty::CM_Linear;
      }
      else if (methods.contains("piecewise", Qt::CaseInsensitive))
      {
        nMethod = SurfaceOverlayProperty::CM_Piecewise;
      }
      else if (!methods.contains("linearopaque", Qt::CaseInsensitive) && !methods.contains("mid_to_min", Qt::CaseInsensitive) &&
               !methods.contains("midtomin", Qt::CaseInsensitive))
      {
        cerr << "Unrecognized overlay method name '" << cmd[1].toLatin1().constData() << "'.\n";
        return;
      }

      overlay->GetProperty()->SetColorMethod( nMethod );

      if (cmd.size() > 2)
      {
        cmd = cmd[2].split(",", MD_SkipEmptyParts);
        bool bPercentile = false, bIgnoreZeros = false;
        while (cmd.last() == "percentile" || cmd.last() == "ignore_zeros")
        {
          if (cmd.last() == "percentile")
            bPercentile = true;
          else
            bIgnoreZeros = true;
          cmd.removeLast();
        }

        double values[3];
        if (bIgnoreZeros)
          overlay->GetProperty()->SetIgnoreZeros(bIgnoreZeros);
        if (bPercentile)
          overlay->GetProperty()->SetUsePercentile(bPercentile);

        bool bOK;
        if ( cmd.size() >= 3 )   // 3 values
        {
          values[0] = cmd[0].toDouble(&bOK);
          values[1] = cmd[1].toDouble(&bOK);
          values[2] = cmd[2].toDouble(&bOK);
        }
        else if (cmd.size() == 2)
        {
          values[0] = cmd[0].toDouble(&bOK);
          values[2] = cmd[1].toDouble(&bOK);
          values[1] = (values[0]+values[2])/2;
        }
        if (bPercentile)
        {
          for (int i = 0; i < 3; i++)
            values[i] = overlay->PercentileToPosition(values[i], bIgnoreZeros);
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

      if (methods.contains("mid_to_min", Qt::CaseInsensitive) || methods.contains("midtomin", Qt::CaseInsensitive))
      {
        overlay->GetProperty()->SetAutoMidToMin(true);
      }

      surf->UpdateOverlay(true);
      overlay->EmitDataUpdated();
    }
  }
}

void MainWindow::CommandSetSurfaceOverlayCustom( const QStringList& cmd_in )
{
  QStringList cmd = cmd_in;
  LayerSurface* surf = (LayerSurface*)GetLayerCollection( "Surface" )->GetActiveLayer();
  if ( surf )
  {
    SurfaceOverlay* overlay = surf->GetActiveOverlay();
    if ( overlay )
    {
      if (cmd.size() < 2)
      {
        cerr << "Insufficient overlay_custom argments\n";
        return;
      }
      else if (cmd.size() == 2)
      {
        overlay->GetProperty()->LoadCustomColorScale(cmd[1]);
      }
      else
      {
        QGradientStops stops;
        QColor c;
        bool bOK;
        for (int i = 1; i < cmd.size(); i++)
        {
          double dval = cmd[i].toDouble(&bOK);
          if (!bOK)
            break;

          c = QColor(cmd[i+1]);
          if (c.isValid())
          {
            i++;
          }
          else
          {
            int r,g,b;
            r = cmd[i+1].toInt(&bOK);
            if (!bOK)
              break;
            g = cmd[i+2].toInt(&bOK);
            if (!bOK)
              break;
            b = cmd[i+3].toInt(&bOK);
            if (!bOK)
              break;
            c = QColor(r,g,b);
            if (!c.isValid())
              break;
            else
              i+=3;
          }
          stops << QGradientStop(dval, c);
        }

        if (!bOK || !c.isValid())
        {
          cerr << "Invalid input for customized overlay color.\n";
          return;
        }
        overlay->GetProperty()->SetColorScale(SurfaceOverlayProperty::CS_Custom);
        overlay->GetProperty()->SetCustomColorScale(stops);
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
      for (int i = 1; i < cmd.size(); i++)
      {
        if (cmd[i] == "colorwheel")
          overlay->GetProperty()->SetColorScale(SurfaceOverlayProperty::CS_ColorWheel);
        else if (cmd[i] == "inverse")
          overlay->GetProperty()->SetColorInverse(true);
        else if (cmd[i] == "truncate")
          overlay->GetProperty()->SetColorTruncate(true);
        else if (cmd[i] == "clearlower")
          overlay->GetProperty()->SetClearLower(true);
        else if (cmd[i] == "clearhigher")
          overlay->GetProperty()->SetClearHigher(true);
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
      cerr << "Invalid color name or value " << cmd[1].toLatin1().constData() << ".\n";
    }
  }
}

void MainWindow::CommandSetSurfaceEdgeColor( const QStringList& cmd )
{
  LayerSurface* surf = (LayerSurface*)GetLayerCollection( "Surface" )->GetActiveLayer();
  if ( surf && cmd[1] != "null" )
  {
    if (cmd[1].toLower() == "off" || cmd[1].toLower() == "surface" || cmd[1].toLower() == "overlay")
    {
      surf->GetProperty()->SetUseSurfaceColorOn2D(true);
    }
    else
    {
      QColor color = ParseColorInput( cmd[1] );
      if ( color.isValid() )
      {
        surf->GetProperty()->SetEdgeColor( color.redF(), color.greenF(), color.blueF() );
      }
      else
      {
        cerr << "Invalid color name or value " << cmd[1].toLatin1().constData() << ".\n";
      }
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

void MainWindow::CommandSetSurfaceOpacity( const QStringList& cmd )
{
  LayerSurface* surf = (LayerSurface*)GetLayerCollection( "Surface" )->GetActiveLayer();
  if ( surf )
  {
    bool bOK;
    double opacity = cmd[1].toDouble(&bOK);
    if ( !bOK || opacity < 0 || opacity > 1)
    {
      cerr << "Invalid opacity value. Must be between 0 and 1.\n";
    }
    else
    {
      surf->GetProperty()->SetOpacity(opacity);
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

void MainWindow::CommandHideSurfaceIn3D(const QStringList &cmd)
{
  LayerSurface* surf = (LayerSurface*)GetLayerCollection( "Surface" )->GetActiveLayer();
  if ( surf && (cmd[1].toLower() == "on" ||
                cmd[1].toLower() == "true" ||
                cmd[1].toLower() == "yes" ||
                cmd[1] == "1") )
  {
    surf->SetHideIn3D(true);
  }
}

void MainWindow::CommandGoToSurfaceVertex(const QStringList &cmd)
{
  LayerSurface* surf = (LayerSurface*)GetLayerCollection( "Surface" )->GetActiveLayer();
  if ( surf )
  {
    bool bOK;
    int nVertex = cmd[1].toInt(&bOK);
    if ( !bOK )
    {
      cerr << "Invalid edge thickness value. Must be a integer.\n";
    }
    else
    {
      double pos[3];
      if (surf->GetTargetAtVertex(nVertex, pos))
      {
        this->GetMainView()->CenterAtWorldPosition(pos);
        GetLayerCollection("MRI")->SetCursorRASPosition( pos );
        SetSlicePosition(pos);
        double v[3] = {1, 0, 0};
        surf->GetSmoothedVertexNormal(nVertex, v);
        m_views[3]->AlignViewToNormal(v);
      }
    }
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
      cerr << "Invalid color name or value " << cmd[1].toLatin1().constData() << ".\n";
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
      cerr << "Invalid color name or value " << cmd[1].toLatin1().constData() << ".\n";
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

void MainWindow::CommandGoToLabel(const QStringList &cmd)
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
    bool bOK;
    QStringList list = cmd[1].split(",");
    double val = list[0].toDouble(&bOK);
    if (!bOK)
    {
      int nMap = LayerPropertySurface::CM_Threshold;
      if (cmd[1].toLower() == "off")
        nMap = LayerPropertySurface::CM_Off;
      else if (cmd[1].toLower() == "binary")
        nMap = LayerPropertySurface::CM_Binary;
      layer->GetProperty()->SetCurvatureMap(nMap);
    }
    else
    {
      layer->GetProperty()->SetThresholdMidPoint(val);
      if (list.size() > 1)
      {
        val = list[1].toDouble(&bOK);
        if (bOK)
          layer->GetProperty()->SetThresholdSlope(val);
      }
    }
  }
}

void MainWindow::CommandLoadSurfaceOverlay( const QStringList& cmd_in )
{
  QStringList cmd = cmd_in;
  while (cmd.size() < 4)
    cmd << "n/a";
  QString reg_file = cmd[2];
  if (reg_file == "n/a")
    reg_file = "";
  if (cmd[3] == "mrisp")
    LoadSurfaceParameterization(cmd[1]);
  else
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

void MainWindow::CommandLoadSurfaceCoordsFromParameterization(const QStringList &cmd)
{
  LoadSurfaceCoordsFromParameterization(cmd[1]);
}

void MainWindow::CommandLoadWayPoints( const QStringList& cmd )
{
  QStringList options = cmd[1].split(":");
  QString fn = options[0];
  QString color = "null";
  QString spline_color = "null";
  QString radius = "1";
  QString spline_radius = "0";
  QString spline_heatmap;
  QVariantMap args;
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
      else if ( option == "id")
        args["id"] = argu;
      else
      {
        cerr << "Unrecognized sub-option flag '" << strg.toLatin1().constData() << "'.\n";
      }
    }
  }

  if ( color != "null" || spline_color != "null" )
  {
    m_scripts.insert( 0, QStringList("setpointsetcolor") << color << spline_color );
  }

  if ( radius != "1" || spline_radius != "0" )
  {
    m_scripts.insert( 0, QStringList("setpointsetradius") << radius << spline_radius );
  }
  if (!spline_heatmap.isEmpty())
  {
    m_scripts.insert( 0, QStringList("setpointsetheatmap") << spline_heatmap.split(",", MD_SkipEmptyParts));
  }

  LoadWayPointsFile( fn, args );
}

void MainWindow::CommandLoadControlPoints( const QStringList& cmd )
{
  QStringList options = cmd[1].split(":");
  QString fn = options[0];
  QString color = "null";
  QString radius = "0";
  QVariantMap args;
  bool bCreateNew = false;
  bool bRemindEdit = false;
  QString name;
  if (options.contains("new", Qt::CaseInsensitive))
  {
    options.removeAll("new");
    bCreateNew = true;
    name = QFileInfo(fn).completeBaseName();
  }
  if (options.contains("remind_edit", Qt::CaseInsensitive))
  {
    options.removeAll("remind_edit");
    bRemindEdit = true;
  }
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
      else if ( option == "id")
        args["id"] = argu;
      else
      {
        cerr << "Unrecognized sub-option flag '" << strg.toLatin1().constData() << "'.\n";
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
  if (QFile::exists(fn) || !bCreateNew)
  {
    if (bRemindEdit)
      args["remind_edit"] = true;
    LoadControlPointsFile( fn, args );
  }
  else if (bCreateNew)
  {
    OnNewPointSet(true);
    SetMode(RenderView::IM_Navigate);
    LayerPointSet* ps = (LayerPointSet*)GetActiveLayer("PointSet");
    ps->SetFileName(fn);
    if (args.contains("id"))
      ps->SetID(args["id"].toInt());
    if (!name.isEmpty())
      ps->SetName(name);
    if (bRemindEdit)
      ps->setProperty("remind_edit", true);
  }
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
        cerr << "Invalid color name or value " << cmd[1].toLatin1().constData() << ".\n";
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
        cerr << "Invalid color name or value " << cmd[2].toLatin1().constData() << ".\n";
      }
    }
  }
}

void MainWindow::CommandSetPointSetRadius( const QStringList& cmd )
{
  LayerPointSet* wp = (LayerPointSet*)GetLayerCollection( "PointSet" )->GetActiveLayer();
  if ( wp )
  {
    if ( !cmd[1].isEmpty() )
    {
      bool bOK;
      double dvalue = cmd[1].toDouble(&bOK);
      if ( bOK && dvalue >= 0)
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

  bool auto_trim = false;
  if (cmd.size() > 3 && (cmd[3] == "autotrim" || cmd[3] == "true" || cmd[3] == "1"))
    auto_trim = true;

  if (cmd[1].contains("%name"))
  {
    QString type = GetCurrentLayerType();
    if (type == "MRI" || type == "Surface")
    {
      QStringList files;
      QList<Layer*> layers = GetLayers(type);
      for (int n = 0; n < layers.size(); n++)
      {
        for (int i = 0; i < layers.size(); i++)
          layers[i]->SetVisible(i == n);

        GetMainView()->RequestRedraw(true);
        QString fn = cmd[1];
        fn.replace("%name", layers[n]->GetName());
        if (!GetMainView()->SaveScreenShot( fn, m_settingsScreenshot.AntiAliasing,
                                            (int)mag_factor))
        {
          cerr << "Failed to save screen shot to " << fn.toLatin1().constData() << ".\n";
        }
        else
          files << fn;
      }
      if (auto_trim)
        GetMainView()->TrimImageFiles(files);
    }
  }
  else if (!GetMainView()->SaveScreenShot( cmd[1],
                                           m_settingsScreenshot.AntiAliasing,
                                           (int)mag_factor, auto_trim))
  {
    cerr << "Failed to save screen shot to " << cmd[1].toLatin1().constData() << ".\n";
  }
}

void MainWindow::CommandFlyThrough(const QStringList &cmd)
{
  Q_UNUSED(cmd);
  if (GetMainViewId() > 2)
  {
    cerr << "Can not fly through. Please set main viewport to 2D slice view.\n";
    return;
  }
}

void MainWindow::CommandSetViewport( const QStringList& cmd )
{
  QString strg = cmd[1].toLower();
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
    std::cerr << "Unrecognized viewport name '" << qPrintable(cmd[1]) << "'.\n";
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
  SetViewSize(x, y);
}

void MainWindow::SetViewSize(int x, int y)
{
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
    for (int i = 0; i < 4; i++)
      m_views[i]->Zoom( dValue );
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
        cerr << "Invalid input value for " << cmd[i].toLatin1().constData() << ".\n";
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
  ras[0] = cmd[1].split(",").first().toDouble(&bOK);
  if (bOK)
  {
    ras[1] = cmd[2].split(",").first().toDouble(&bOK);
  }
  if (bOK)
  {
    ras[2] = cmd[3].split(",").first().toDouble(&bOK);
  }
  if ( bOK )
  {
    LayerMRI* layer = qobject_cast<LayerMRI*>(GetActiveLayer("MRI"));
    if ( layer )
    {
      if ( cmd.size() > 4 && cmd[4] == "tkreg" )
      {
        layer->TkRegToNativeRAS( ras, ras );
      }
      layer->RASToTarget( ras, ras );
    }
    else
    {
      LayerSurface* surf = qobject_cast<LayerSurface*>(GetActiveLayer("Surface"));
      if (surf && cmd.size() > 4 && cmd[4] == "tkreg")
      {
        surf->GetTargetAtSurfaceRAS(ras, ras);
      }
    }
    this->GetMainView()->CenterAtWorldPosition(ras);
    GetLayerCollection("MRI")->SetCursorRASPosition( ras );
    SetSlicePosition( ras );
    ((RenderView3D*)m_views[3])->MapToInflatedCoords(ras);
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
    x = cmd[1].split(",").first().toInt(&bOK);
    y = cmd[2].split(",").first().toInt(&bOK);
    z = cmd[3].split(",").first().toInt(&bOK);
    if ( bOK )
    {
      int slice[3] = { x, y, z };
      double ras[3];
      mri->OriginalIndexToRAS( slice, ras );
      mri->RASToTarget( ras, ras );

      lc_mri->SetCursorRASPosition( ras );
      SetSlicePosition( ras );
      ((RenderView3D*)m_views[3])->MapToInflatedCoords(ras);
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

void MainWindow::CommandWriteSurfaceIntersection( const QStringList& cmd)
{
  LayerSurface* surf = (LayerSurface*)GetActiveLayer("Surface");
  LayerCollection* lc_mri = GetLayerCollection( "MRI" );
  LayerMRI* mri = NULL;
  if (!lc_mri->IsEmpty())
    mri = (LayerMRI*)lc_mri->GetLayer( lc_mri->GetNumberOfLayers()-1 );
  if (surf && mri)
  {
    QString ostr = mri->GetOrientationString();
    QString slice_str = "IS";
    int nPlane = 2;
    if (cmd[1].contains("sag", Qt::CaseInsensitive))
    {
      slice_str = "RL";
      nPlane = 0;
    }
    else if (cmd[1].contains("cor", Qt::CaseInsensitive))
    {
      slice_str = "AP";
      nPlane = 1;
    }
    int n = 0;
    for (int i = 0; i < 3; i++)
    {
      if (slice_str.contains(ostr[i]))
      {
        n = i;
        break;
      }
    }
    int slice[3] = { 0, 0, 0 };
    slice[n] = cmd[3].toInt();
    double ras[3];
    mri->OriginalIndexToRAS( slice, ras );
    mri->RASToTarget( ras, ras );
    lc_mri->SetCursorRASPosition( ras );
    SetSlicePosition( ras );

    surf->WriteIntersection(cmd[2], nPlane, mri);
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
  //  settings.sync();

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
  LayerCollection* lc = GetLayerCollection(strType);
  if (lc)
    return lc->GetLayers();
  else
    return QList<Layer*>();
}

QList<Layer*> MainWindow::GetVisibleLayers(const QString &strType)
{
  QList<Layer*> list = GetLayers(strType);
  QList<Layer*> list_visbile;
  foreach (Layer* l, list)
  {
    if (l->IsVisible())
      list_visbile << l;
  }
  return list_visbile;
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
    emit SlicePositionChanged(true);
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
    emit SlicePositionChanged(true);
  }

  return bRet;
}

bool MainWindow::OffsetSlicePosition( int nPlane, double dPosDiff, bool bRoundToGrid  )
{
  bool bRet = false;
  QStringList keys = m_layerCollections.keys();
  LayerCollection* lc_mri = m_layerCollections["MRI"];
  if (!lc_mri->IsEmpty())
  {
    lc_mri->blockSignals( true );
    if ( lc_mri->OffsetSlicePosition( nPlane, dPosDiff, bRoundToGrid ) )
    {
      bRet = true;
    }
    lc_mri->blockSignals( false );
    if (bRet)
    {
      keys.removeOne("MRI");
      double slicePos[3];
      lc_mri->GetSlicePosition(slicePos);
      for ( int i = 0; i < keys.size(); i++ )
      {
        m_layerCollections[keys[i]]->blockSignals( true );
        m_layerCollections[keys[i]]->SetSlicePosition(nPlane, slicePos[nPlane], false);
        m_layerCollections[keys[i]]->blockSignals( false );
      }
    }
  }
  else
  {
    for ( int i = 0; i < keys.size(); i++ )
    {
      m_layerCollections[keys[i]]->blockSignals( true );
      if ( m_layerCollections[keys[i]]->OffsetSlicePosition( nPlane, dPosDiff, bRoundToGrid ) )
      {
        bRet = true;
      }
      m_layerCollections[keys[i]]->blockSignals( false );
    }
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

      if (dlg.GetLoadAsVector())
        fn += ":vector=1";

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
    if (!GetLayerCollection( "Surface" )->IsEmpty() || !GetLayerCollection("Tract")->IsEmpty())
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

  layer->GetProperty()->blockSignals(true);
  if (sup_data.contains("Basis"))
    layer->SetLayerIndex(sup_data["Basis"].toInt());

  if (sup_data.contains("Percentile"))
    layer->GetProperty()->SetUsePercentile(sup_data["Percentile"].toBool());

  if (sup_data.contains("ID"))
    layer->SetID(sup_data["ID"].toInt());

  if (sup_data.contains("VectorSkip"))
    layer->GetProperty()->SetVectorSkip(qMax(0, sup_data["VectorSkip"].toInt()));

  if (sup_data.contains("VectorNormalize"))
    layer->GetProperty()->SetNormalizeVector(sup_data["VectorNormalize"].toBool());

  if (sup_data.contains("VectorLengthScale"))
    layer->GetProperty()->SetVectorDisplayScale(sup_data["VectorLengthScale"].toDouble());

  if (sup_data.contains("BinaryColor"))
    layer->GetProperty()->SetBinaryColor(sup_data["BinaryColor"].value<QColor>());

  layer->GetProperty()->blockSignals(false);

  if (sup_data.value("IgnoreHeader").toBool())
    layer->SetIgnoreHeader(true);

  m_threadIOWorker->LoadVolume( layer );
}

bool MainWindow::OnCloseVolume(const QList<Layer*>& layers_in)
{
  QList<Layer*> layers = layers_in;
  if (layers.isEmpty())
    layers = GetSelectedLayers( "MRI" );
  if ( layers.isEmpty() )
  {
    return false;
  }
  foreach (Layer* layer, layers)
  {
    if ( qobject_cast<LayerMRI*>(layer)->IsModified() )
    {
      QMessageBox box(QMessageBox::Question, tr("Close Volume"),
                      "Volume has been modifed and not been saved. Do you still want to continue?",
                      QMessageBox::Yes | QMessageBox::Cancel);
      box.setButtonText(QMessageBox::Yes, tr("Continue Without Saving"));
      box.setDefaultButton(QMessageBox::Cancel);
      if (box.exec() != QMessageBox::Yes)
      {
        return false;
      }
    }
  }
  foreach (Layer* layer, layers)
  {
    if (layer == m_layerVolumeRef)
    {
      m_layerVolumeRef = NULL;
    }
  }
  GetLayerCollection( "MRI" )->RemoveLayers( layers );
  if (m_layerVolumeRef == NULL)
    m_layerVolumeRef = (LayerMRI*)GetActiveLayer("MRI");

  OnSetModeNavigate();

  if (GetLayers("MRI").isEmpty())
  {
    LayerCollection* lc = GetLayerCollection("Supplement");
    QList<Layer*> layers = lc->GetLayers("MRI");
    foreach (Layer* layer, layers)
    {
      if (layer->GetName() == "GEOS_DRAW")
      {
        lc->RemoveLayer(layer);
      }
      else if (layer->GetName() == "GEOS_FILL")
      {
        lc->RemoveLayer(layer);
      }
    }
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
    col_mri->AddLayer(layer_new);
    ConnectMRILayer(layer_new);
    emit NewVolumeCreated();
  }
}

void MainWindow::ConnectMRILayer(LayerMRI *mri)
{
  for (int i = 0; i < 4; i++)
    connect(mri->GetProperty(), SIGNAL(ColorMapChanged()), m_views[i], SLOT(UpdateScalarBar()), Qt::UniqueConnection);
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

void MainWindow::LoadROIFile( const QString& fn, const QString& ref_vol, const QVariantMap& args )
{
  LayerMRI* ref = NULL;
  LayerCollection* col_mri = GetLayerCollection( "MRI" );
  if ( ref_vol.isEmpty() )
  {
    //   cout << "No template volume given, using current volume as template for ROI " << fn.toLatin1().constData() << ".\n";
    ref = (LayerMRI*)col_mri->GetActiveLayer();
  }
  else
  {
    for ( int i = 0; i < col_mri->GetNumberOfLayers(); i++ )
    {
      LayerMRI* mri = ( LayerMRI* )col_mri->GetLayer( i );
      if (ref_vol == QString::number(mri->GetID()) || ref_vol == mri->GetName()
          || QFileInfo( mri->GetFileName() ).fileName() == ref_vol )
      {
        ref = mri;
        break;
      }
    }
    if ( ref == NULL )
    {
      cerr << "Can not find given template volume: " << ref_vol.toLatin1().constData()
           << ". Using current volume as template for ROI " << fn.toLatin1().constData() << ".\n";
      ref = (LayerMRI*)col_mri->GetActiveLayer();
    }
  }
  LayerROI* roi = new LayerROI( ref );
  roi->SetName( QFileInfo(fn).completeBaseName() );
  if (args.contains("color"))
    roi->GetProperty()->SetColor(args["color"].value<QColor>());
  if (args.contains("opacity"))
    roi->GetProperty()->SetOpacity(args["opacity"].toDouble());
  if (args.contains("threshold"))
    roi->GetProperty()->SetThreshold(args["threshold"].toDouble());
  if (args.contains("id"))
    roi->SetID(args["id"].toInt());
  if (args.contains("name"))
    roi->SetName(args["name"].toString());
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
    if (!layer_roi->SaveROI())
    {
      QMessageBox::warning( this, "Error", "Could not save ROI at " + fn);
    }
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

void MainWindow::OnCloseROI(const QList<Layer*>& layers_in)
{
  QList<Layer*> layers = layers_in;
  if (layers.isEmpty())
    layers = GetSelectedLayers( "ROI" );
  if ( layers.isEmpty() )
  {
    return;
  }
  foreach (Layer* layer, layers)
  {
    if ( qobject_cast<LayerROI*>(layer)->IsModified() )
    {
      QMessageBox box(QMessageBox::Question, tr("Close ROI"),
                      "ROI has been modifed and not been saved. Do you still want to continue?",
                      QMessageBox::Yes | QMessageBox::Cancel);
      box.setButtonText(QMessageBox::Yes, tr("Continue Without Saving"));
      box.setDefaultButton(QMessageBox::Cancel);
      if (box.exec() != QMessageBox::Yes)
      {
        return;
      }
    }
  }
  GetLayerCollection( "ROI" )->RemoveLayers( layers );
  OnSetModeNavigate();
}

void MainWindow::OnNewPointSet(bool bSilent)
{
  LayerCollection* col_mri = GetLayerCollection( "MRI" );
  LayerMRI* layer_mri = ( LayerMRI* )col_mri->GetActiveLayer();
  if (!layer_mri && !GetActiveLayer("Surface"))
  {
    QMessageBox::warning( this, "Error", "Can not create new ROI without volume template.");
    return;
  }

  DialogNewPointSet dlg( this );
  dlg.SetPointSetName( tr("New Point Set %1").arg(GetLayerCollection("PointSet")->GetNumberOfLayers()));
  if (bSilent || dlg.exec() == QDialog::Accepted )
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
    LayerPointSet* layer_wp;
    if (bSilent)
    {
      if (layer_mri)
        layer_wp = new LayerPointSet( layer_mri, LayerPropertyPointSet::Enhanced);
      else
        layer_wp = new LayerPointSet( (LayerSurface*)GetActiveLayer("Surface"), LayerPropertyPointSet::Enhanced);
    }
    else
    {
      LayerMRI* mri = dlg.GetTemplate();
      if (mri)
        layer_wp = new LayerPointSet( dlg.GetTemplate(), dlg.GetType() );
      else
        layer_wp = new LayerPointSet(dlg.GetTemplateSurface(), dlg.GetType());
    }
    if (!bSilent)
      layer_wp->SetName( dlg.GetPointSetName() );
    col_wp->AddLayer( layer_wp );

    SetMode( RenderView::IM_PointSetEdit );
  }
}

void MainWindow::LoadWayPointsFile( const QString& fn, const QVariantMap& args )
{
  this->LoadPointSetFile( fn, LayerPropertyPointSet::WayPoint, args);
}

void MainWindow::LoadControlPointsFile( const QString& fn, const QVariantMap& args )
{
  this->LoadPointSetFile( fn, LayerPropertyPointSet::ControlPoint, args );
}

void MainWindow::LoadPointSetFile( const QString& fn, int type, const QVariantMap& args )
{
  LayerCollection* col_mri = GetLayerCollection( "MRI" );
  LayerMRI* mri = ( LayerMRI* )col_mri->GetActiveLayer();
  LayerSurface* surf = (LayerSurface*)GetActiveLayer("Surface");
  LayerPointSet * wp;
  if (mri)
    wp = new LayerPointSet( mri, type );
  else
    wp = new LayerPointSet(surf, type);
  wp->SetName( QFileInfo( fn ).fileName() );
  if ( wp->LoadFromFile( fn ) )
  {
    if (args.contains("id"))
      wp->SetID(args["id"].toInt());

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

    if (args["remind_edit"].toBool())
      wp->setProperty("remind_edit", true);
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

void MainWindow::OnSavePointSet(bool bForce)
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
  if ( fn.isEmpty() || (!bForce && layer->IsEnhanced() && layer->GetProperty()->GetType() != LayerPropertyPointSet::Enhanced) )
  {
    if (layer->IsEnhanced())
      QMessageBox::information(this, "Save Point Set", "To save comments or other added information, it is recommended to save point set in enhanced format in Json.");
    OnSavePointSetAs();
  }
  else
  {
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
  QString fn = layer->GetFileName();
  if (fn.isEmpty())
    fn = layer->GetName();
  int nType = layer->GetProperty()->GetType();
  if (layer->IsEnhanced() || nType == LayerPropertyPointSet::ControlPoint)
    nType = LayerPropertyPointSet::Enhanced;
  dlg.SetFileName(fn, nType);
  dlg.SetType(nType);
  dlg.SetLastDir(m_strLastDir);
  if (dlg.exec() == QDialog::Accepted)
  {
    layer->SetFileName( dlg.GetFileName() );
    layer->GetProperty()->SetType( dlg.GetType() );
    OnSavePointSet(true);
    ui->widgetAllLayers->UpdateWidgets();
  }
}

void MainWindow::OnClosePointSet(const QList<Layer*>& layers_in)
{
  QList<Layer*> layers = layers_in;
  if (layers.isEmpty())
    layers = GetSelectedLayers( "PointSet" );

  if ( layers.isEmpty() )
  {
    return;
  }
  foreach (Layer* layer, layers)
  {
    LayerPointSet* ps = qobject_cast<LayerPointSet*>(layer);
    if (ps->IsModified())
    {
      QMessageBox box(QMessageBox::Question, tr("Close Point Set"),
                      tr("Point Set %1 has been modifed and not been saved. Do you still want to close it?").arg(ps->GetName()),
                      QMessageBox::Yes | QMessageBox::Cancel);
      box.setButtonText(QMessageBox::Yes, tr("Close Without Saving"));
      box.setDefaultButton(QMessageBox::Cancel);
      if (box.exec() != QMessageBox::Yes)
      {
        return;
      }
    }
    if (ps->property("remind_edit").toBool() && !ps->IsEdited())
    {
      QMessageBox box(QMessageBox::Question, tr("Close Point Set"),
                      tr("Point Set %1 has not been edited. Do you still want to close it?").arg(ps->GetName()),
                      QMessageBox::Yes | QMessageBox::Cancel);
      box.setButtonText(QMessageBox::Yes, tr("Close It"));
      box.setDefaultButton(QMessageBox::Cancel);
      if (box.exec() != QMessageBox::Yes)
      {
        return;
      }
    }
  }

  GetLayerCollection( "PointSet" )->RemoveLayers( layers );
  OnSetModeNavigate();
}

void MainWindow::OnLoadTrack()
{
  QStringList filenames = QFileDialog::getOpenFileNames( this, "Select tract file",
                                                         m_strLastDir,
                                                         "Tract files (*.trk);;All files (*)");
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
  QList<Layer*> layers = GetSelectedLayers( "Tract" );
  if ( layers.isEmpty() )
  {
    return;
  }

  GetLayerCollection( "Tract" )->RemoveLayers( layers );
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
                                  const QStringList& sup_files_in, const QVariantMap& sup_options)
{
  QFileInfo fi( filename );
  QString ext = fi.completeSuffix(); 
  m_strLastDir = fi.absolutePath();
  LayerSurface* layer = new LayerSurface( m_layerVolumeRef );
  connect(layer, SIGNAL(CurrentVertexChanged(int)), m_wndGroupPlot, SLOT(SetCurrentVertex(int)), Qt::UniqueConnection);
  connect(ui->treeWidgetCursorInfo, SIGNAL(VertexChangeTriggered(int)), m_wndGroupPlot, SLOT(SetCurrentVertex(int)), Qt::UniqueConnection);
  connect(layer, SIGNAL(SurfaceOverlyDataUpdated()), ui->treeWidgetCursorInfo, SLOT(UpdateAll()), Qt::UniqueConnection);
  connect(layer, SIGNAL(ActiveSurfaceChanged(int)), ui->view3D, SLOT(OnLayerVisibilityChanged()), Qt::UniqueConnection);
  connect(this, SIGNAL(SlicePositionChanged(bool)), layer, SLOT(OnSlicePositionChanged3D()), Qt::UniqueConnection);
  connect(layer, SIGNAL(SurfaceOverlayAdded(SurfaceOverlay*)), m_wndTimeCourse, SLOT(UpdateUI()), Qt::UniqueConnection);
  connect(layer, SIGNAL(ActiveOverlayChanged(int)), m_wndTimeCourse, SLOT(UpdateAll()), Qt::UniqueConnection);
  connect(layer, SIGNAL(FlattenedPatchLoaded()), this, SLOT(OnFlattendSurfacePatchLoaded()), Qt::QueuedConnection);
  layer->SetName( fi.fileName() );
  QString fullpath = fi.absoluteFilePath();
  if ( fullpath.isEmpty() )
  {
    fullpath = filename;
  }
  QStringList sup_files = sup_files_in;
  if (sup_options.value("no_autoload").toBool())
  {
    layer->SetSphereFileName("");
    if (fi.fileName().contains("inflated"))
    {
      QString fn = fi.absoluteFilePath();
      fn.replace(".inflated", ".white");
      if (QFile::exists(fn) && !sup_files.contains("white"))
        sup_files << "white";
    }
  }
  else
  {
    if (fi.fileName().contains("inflated.nofix"))
    {
      if (!sup_files.contains("orig.nofix"))
        sup_files << "orig.nofix";
    }
    else if (fi.fileName().contains("inflated"))
    {
      if (!sup_files.contains("white"))
        sup_files << "white";
    }
  }
  layer->SetFileName( fullpath );
  layer->SetPatchFileName( fn_patch );
  layer->SetTargetFileName( fn_target );
  layer->SetLoadSupSurfaces(sup_files);
  if (sup_options.contains("ID"))
    layer->SetID(sup_options.value("ID").toInt());
  if (sup_options.contains("sphere"))
    layer->SetSphereFileName(sup_options["sphere"].toString());
  layer->GetProperty()->blockSignals(true);
  if (sup_options.contains("ZOrderAnnotation"))
    layer->GetProperty()->SetZOrderAnnotation(sup_options["ZOrderAnnotation"].toInt());
  if (sup_options.contains("ZOrderLabel"))
    layer->GetProperty()->SetZOrderLabel(sup_options["ZOrderLabel"].toInt());
  if (sup_options.contains("ZOrderOverlay"))
    layer->GetProperty()->SetZOrderOverlay(sup_options["ZOrderOverlay"].toInt());
  layer->GetProperty()->blockSignals(false);

  QVariantMap args;
  if(sup_options.contains("ignore_vg"))
    args["ignore_vg"] = sup_options["ignore_vg"];
  // For sphere, ignore the vol geom by default
  if(fi.fileName().contains("sphere")){
    char *a = getenv("FV_SPHERE_IGNORE_VG");
    if(a == NULL || strcmp(a,"0")!=0) args["ignore_vg"] = 1;
  }
  if(args["ignore_vg"].toInt()) {
    printf("INFO: ignoring vol geom for sphere %s\n",filename.toStdString().c_str());
    printf("  If you want the vol geom, then add --no-sphere-ignore-vg or setenv FV_SPHERE_IGNORE_VG 0\n");
    printf("  Note that freeview identifies a sphere by checking if the filename contains 'sphere'.\n");
  }

  if (sup_options.contains("affinexform_filename"))
    args["affinexform_filename"] = sup_options["affinexform_filename"];

  if (sup_options["no_shading"].toBool())
    layer->SetNoShading(true);

  m_threadIOWorker->LoadSurface( layer, args );
  m_statusBar->StartTimer();
}

void MainWindow::OnCloseSurface(const QList<Layer*>& layers_in)
{
  QList<Layer*> layers = layers_in;
  if (layers.isEmpty())
    layers = GetSelectedLayers("Surface");
  if ( layers.isEmpty() )
  {
    return;
  }

  GetLayerCollection( "Surface" )->RemoveLayers( layers );
}

void MainWindow::OnLoadPatch()
{
  LayerSurface* surf = ( LayerSurface* )GetActiveLayer("Surface");
  if ( !surf )
    return;

  QString filename = QFileDialog::getOpenFileName( this, "Select patch file",
                                                   AutoSelectLastDir( "surf" ),
                                                   "Patch files (*)", 0, QFileDialog::DontConfirmOverwrite);
  if ( !filename.isEmpty() && !surf->LoadPatch(filename) )
  {
    QMessageBox::warning(this, "Error", QString("Could not load patch from %1").arg(filename));
  }
}

void MainWindow::OnSavePatchAs()
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

  QString fn = QFileDialog::getSaveFileName( this, "Save patch as",
                                             layer_surf->GetFileName(),
                                             "Patch files (*)");
  if ( !fn.isEmpty() && !layer_surf->WritePatch(fn))
  {
    QMessageBox::warning(this, "Error", QString("Could not save patch to %1").arg(fn));
  }
}

void MainWindow::OnIOError( Layer* layer, int jobtype )
{
  bool bQuit = false;
  foreach (QStringList list, m_scripts)
  {
    if (list[0] == "quit")
    {
      bQuit = true;
      break;
    }
  }

  if (!m_bContinue)
    ClearScripts();
  QString msg = QString("Failed to load %1 ").arg(layer->GetEndType());
  if (jobtype != ThreadIOWorker::JT_LoadSurfaceOverlay)
  {
    if ( jobtype == ThreadIOWorker::JT_SaveVolume )
    {
      msg = "Failed to save volume to " + layer->property("saved_name").toString();
    }
    else if ( jobtype == ThreadIOWorker::JT_SaveSurface )
    {
      msg = "Failed to save surface to " + layer->GetFileName();
    }
    if (!bQuit)
      QMessageBox::warning( this, "Error", msg);
    if ( jobtype != ThreadIOWorker::JT_SaveVolume && jobtype != ThreadIOWorker::JT_SaveSurface )
      delete layer;
  }
  else
  {
    msg += "overlay.";
    if (!bQuit)
      QMessageBox::warning( this, "Error", msg);
  }
  m_bProcessing = false;
  m_layerSettings.clear();
  if (bQuit)
  {
    cout << qPrintable(msg) << endl;
    close();
  }
}

void MainWindow::OnIOFinished( Layer* layer, int jobtype )
{
  m_statusBar->StopTimer();
  LayerCollection* lc_mri = GetLayerCollection( "MRI" );
  LayerCollection* lc_surface = GetLayerCollection( "Surface" );
  LayerCollection* lc_track = GetLayerCollection( "Tract" );
  LayerCollection* lc_odf = GetLayerCollection( "ODF" );
  LayerCollection* lc_sup = GetLayerCollection( "Supplement");
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

      lc_sup->SetWorldVoxelSize( mri->GetWorldVoxelSize() );
      lc_sup->SetWorldOrigin( mri->GetWorldOrigin() );
      lc_sup->SetWorldSize( mri->GetWorldSize() );

      lc_mri->AddLayer( layer, true );
      ConnectMRILayer(mri);
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

      int dim[3];
      double vs[3];
      mri->GetVolumeInfo(dim, vs);
      if (dim[0] == 1)
        SetMainView(0);
      else if (dim[1] == 1)
        SetMainView(1);
      else if (dim[2] == 1)
        SetMainView(2);
    }
    else
    {
      lc_mri->AddLayer( mri );
      ConnectMRILayer(mri);
    }

    m_strLastDir = QFileInfo( layer->GetFileName() ).canonicalPath();
    SetCurrentFile( layer->GetFileName(), 0 );
    if (m_layerSettings.contains(layer->GetID()))
    {
      QVariantMap settings = m_layerSettings[layer->GetID()];
      if (settings.contains("name"))
        mri->SetName(settings["name"].toString());
      if (settings.contains("index"))
      {
        QList<Layer*> layers = lc_mri->GetLayers();
        layers.removeAt(0);
        layers.insert(settings["index"].toInt(), layer);
        lc_mri->ReorderLayers(layers);
      }
      mri->GetProperty()->RestoreFullSettings(settings);
      if (settings.contains("frame"))
        mri->SetActiveFrame(settings["frame"].toInt());
      m_layerSettings.remove(layer->GetID());
    }
  }
  else if ( jobtype == ThreadIOWorker::JT_SaveVolume && layer->IsTypeOf( "MRI" ) )
  {
    SetCurrentFile( layer->GetFileName(), 0 );
  }
  else if ( jobtype == ThreadIOWorker::JT_LoadSurface && layer->IsTypeOf( "Surface" ) )
  {
    LayerSurface* sf = qobject_cast<LayerSurface*>( layer );
    if (sf->property("hidden").toBool())
    {
      GetLayerCollection("HiddenSurface")->AddLayer(sf);
      m_bProcessing = false;
      return;
    }
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

    if ( !sf->HasValidVolumeGeometry() && !sf->property("IgnoreVG").toBool())
    {
      //  ShowNonModalMessage("Warning",
      //                      "Either this surface does not contain valid volume geometry information, or freeview failed to read the information. This surface may not align with volumes and other surfaces.");
      cout << "Did not find any volume info" << endl;
    }

    m_strLastDir = QFileInfo( layer->GetFileName() ).canonicalPath();
    SetCurrentFile( layer->GetFileName(), 1 );
    if (m_layerSettings.contains(layer->GetID()))
    {
      QVariantMap settings = m_layerSettings[layer->GetID()];
      sf->GetProperty()->RestoreFullSettings(settings);
      m_layerSettings.remove(layer->GetID());
    }
    if (UpdateSurfaceCorrelation((LayerSurface*)layer) )
    {
      emit SlicePositionChanged(true);
    }

    LoadSphereLeftRightIfNeeded(sf);

    if (sf->GetPatchFileName().contains("flat"))
      OnFlattendSurfacePatchLoaded();
  }
  else if ( jobtype == ThreadIOWorker::JT_LoadSurfaceOverlay && layer->IsTypeOf("Surface") )
  {
    UpdateSurfaceCorrelation((LayerSurface*)layer);
    lc_surface->SetActiveLayer(layer);
    m_strLastDir = QFileInfo(((LayerSurface*)layer)->GetActiveOverlay()->GetFileName()).absolutePath();
    emit SlicePositionChanged(true);
  }
  else if ( jobtype == ThreadIOWorker::JT_LoadTrack && layer->IsTypeOf("Tract"))
  {
    LayerTrack* track = qobject_cast<LayerTrack*>( layer );
    lc_track->AddLayer( track );
    m_strLastDir = QFileInfo( layer->GetFileName() ).canonicalPath();
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
  else if ( jobtype == ThreadIOWorker::JT_LoadODF && layer->IsTypeOf("ODF"))
  {
    LayerODF* odf = qobject_cast<LayerODF*>( layer );
    m_strLastDir = QFileInfo( layer->GetFileName() ).canonicalPath();
    double worigin[3], wsize[3];
    odf->GetWorldOrigin( worigin );
    odf->GetWorldSize( wsize );
    if (lc_surface->IsEmpty() && lc_mri->IsEmpty())
    {
      for ( int i = 0; i < 4; i++ )
      {
        m_views[i]->SetWorldCoordinateInfo( worigin, wsize, true );
      }
      m_views[3]->ResetCameraClippingRange();
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
      lc_odf->SetWorldOrigin( worigin );
      lc_odf->SetWorldSize( wsize );
      lc_odf->SetWorldVoxelSize( vs );
      lc_odf->AddLayer( odf );
      lc_odf->SetSlicePosition( lc_mri->GetSlicePosition() );
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
    std::cout << qPrintable(layer->property("saved_name").toString()) << " saved successfully.\n";
  }
  else if ( jobtype == ThreadIOWorker::JT_SaveSurface)
  {
    std::cout << qPrintable(qobject_cast<LayerSurface*>(layer)->GetFileName()) << " saved successfully.\n";
    m_dlgRepositionSurface->UpdateUI();
  }
  else if ( jobtype == ThreadIOWorker::JT_TransformVolume)
  {
    RequestRedraw();
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
  LayerCollection* lc = GetLayerCollection(GetCurrentLayerType());
  if ( lc )
  {
    lc->CycleLayer(true, true); //ui->viewAxial->GetInteractionMode() == RenderView2D::IM_ReconEdit);
  }
}

void MainWindow::OnReverseCycleLayer()
{
  LayerCollection* lc = GetLayerCollection(GetCurrentLayerType());
  if ( lc )
  {
    lc->CycleLayer( false, true); //ui->viewAxial->GetInteractionMode() == RenderView2D::IM_ReconEdit );
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

QString MainWindow::AutoSelectLastDir( const QString& lastDir_in, const QString& subdir )
{
  // ignore lastDir_in if there is a selected layer
  MainWindow* mainwnd = MainWindow::GetMainWindow();
  QString lastDir = lastDir_in;
  QString layerType = mainwnd->GetCurrentLayerType();
  if (!layerType.isEmpty())
  {
    Layer* layer = mainwnd->GetActiveLayer(layerType);
    if (layer && !layer->GetFileName().isEmpty())
      lastDir = QFileInfo(layer->GetFileName()).absolutePath();
  }
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
  m_dlgPreferences->raise();
}

void MainWindow::SetVolumeColorMap( int nColorMap, int nColorMapScale, const QList<double>& scales_in )
{
  LayerMRI* layer = qobject_cast<LayerMRI*>(GetActiveLayer("MRI"));
  if ( layer )
  {
    LayerPropertyMRI* p = layer->GetProperty();
    p->SetColorMap( (LayerPropertyMRI::ColorMapType) nColorMap );
    if (!scales_in.isEmpty())
    {
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
          p->SetHeatScaleSetMidToMin(true);
          p->SetHeatScaleMinThreshold( scales[0]);
          p->SetHeatScaleMaxThreshold( scales[1]);
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
}

void MainWindow::OnTransformVolume()
{
  cout << "Warning: Transformation can only apply to volumes for now. If your data includes ROI/Surface/Way Points, please do not use this feature yet.\n";
  m_dlgTransformVolume->show();
  m_dlgTransformVolume->raise();
  m_dlgTransformVolume->UpdateUI();
}

void MainWindow::OnTransformSurface()
{
  m_dlgTransformSurface->show();
  m_dlgTransformSurface->raise();
  m_dlgTransformSurface->UpdateUI();
}

void MainWindow::OnCropVolume()
{
  LayerMRI* mri = (LayerMRI*)GetActiveLayer( "MRI" );
  m_dlgCropVolume->SetVolume( mri );
  m_dlgCropVolume->show();
  m_dlgCropVolume->raise();
  m_volumeCropper->SetEnabled( true );
  m_volumeCropper->SetVolume( mri );
  m_volumeCropper->Show();
  SetMode( RenderView::IM_VolumeCrop );
  for (int i = 0; i < 4; i++)
    m_views[i]->ResetCameraClippingRange();
}

void MainWindow::OnThresholdVolume()
{
  m_dlgThresholdVolume->show();
  m_dlgThresholdVolume->raise();
}

void MainWindow::OnSegmentVolume()
{
  m_dlgVolumeSegmentation->show();
  m_dlgVolumeSegmentation->raise();
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
  m_dlgSaveScreenshot->raise();
}

void MainWindow::OnCopyView()
{
  if (m_dlgSaveScreenshot)
    SetScreenShotSettings(m_dlgSaveScreenshot->GetSettings());
  QString fn = QDir::tempPath() + "/freeview-temp-" + QString::number(QDateTime::currentMSecsSinceEpoch()) + ".png";
  GetMainView()->SaveScreenShot(fn, m_settingsScreenshot.AntiAliasing, 1.0, m_settingsScreenshot.AutoTrim);
  QClipboard *clipboard = QGuiApplication::clipboard();
  if (clipboard)
    clipboard->setImage(QImage(fn));
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

void MainWindow::OnVolumeFilterBoundary()
{
  LayerMRI* mri = (LayerMRI*)GetActiveLayer( "MRI" );
  if ( mri )
  {
    VolumeFilterBoundary* filter = new VolumeFilterBoundary( mri, mri );
    m_threadVolumeFilter->ExecuteFilter(filter);
    mri->ResetWindowLevel();
  }
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
  if (lc->IsEmpty())
    lc = GetLayerCollection("Surface");
  QString fn;
  QString path = getenv( "FS_SAVE_GOTO_POINT" );
  for ( int i = 0; i < lc->GetNumberOfLayers(); i++ )
  {
    fn = lc->GetLayer( i )->GetFileName();
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
    QStringList args = strg.trimmed().split(QRegularExpression("\\s+"));
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
    LayerSurface* surf = (LayerSurface*)lc_surf->GetActiveLayer();
    double slice_pos[3];
    lc_surf->GetSlicePosition(slice_pos);
    surf->GetSurfaceRASAtTarget( slice_pos, ras_out);
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
    if (layer->IsVisible() && !layer->IsLocked())
    {
      bVisible = true;
      break;
    }
  }
  foreach (Layer* layer, layers)
  {
    if (!layer->IsLocked())
      layer->SetVisible(!bVisible);
  }
}

void MainWindow::OnAbout()
{
  DialogAbout dlg(this);
  dlg.exec();
}

void MainWindow::OnActiveLayerChanged(Layer* layer)
{
  QString title = "FreeView";
  if (!m_sTitle.isEmpty())
    title += ": " + m_sTitle;
  if (!layer)
  {
    this->setWindowTitle(title);
    m_wndTimeCourse->hide();
  }
  else
  {
    QString fn = layer->GetFileName();
    if (layer->IsTypeOf("Tract") && ((LayerTrack*)layer)->IsCluster())
      fn = QFileInfo(fn).absolutePath() + "/*.trk";
    this->setWindowTitle(QString("%1 (%2)").arg(title)
                         .arg(fn));
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
    else if (layer->IsTypeOf("Tract"))
    {
      QList<Layer*> layers = GetLayers("Tract");
      foreach (Layer* layer, layers)
        disconnect(m_wndTractCluster, 0, layer, 0);
      LayerTrack* tract = qobject_cast<LayerTrack*>(layer);
      if (tract && tract->IsCluster())
      {
        connect(m_wndTractCluster, SIGNAL(TreeNodeActivated(QStringList)), tract,
                SLOT(LoadTrackFromFiles(QStringList)), Qt::UniqueConnection);
        m_wndTractCluster->SetData(tract->GetClusterData());
      }
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
  m_dlgWriteMovieFrames->raise();
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
  ui->view3D->ShowCursor( bShow );
  this->RequestRedraw();
}

void MainWindow::ShowNonModalMessage(const QString &title, const QString &msg)
{
  m_dlgMessage->setWindowTitle(title);
  m_dlgMessage->setText(msg);
  m_dlgMessage->show();
  m_dlgMessage->raise();
}

void MainWindow::OnRepositionSurface()
{
  m_dlgRepositionSurface->show();
  m_dlgRepositionSurface->raise();
}

void MainWindow::OnSmoothSurface()
{
  m_dlgSmoothSurface->show();
  m_dlgSmoothSurface->raise();
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
  m_dlgLabelStats->raise();
}

void MainWindow::OnLineProfile()
{
#if !defined(ARM64)
  m_dlgLineProfile->show();
  m_dlgLineProfile->raise();
#endif
}

void MainWindow::OnSaveIsoSurface(const QString& fn_in)
{
  LayerMRI* layer = qobject_cast<LayerMRI*>(GetActiveLayer("MRI"));
  if (!layer)
    return;

  QString fn = fn_in;
  QString selectedFilter;
  if (fn.isEmpty())
    fn = QFileDialog::getSaveFileName( NULL,
                                       "Save iso-surface",
                                       MainWindow::GetMainWindow()->AutoSelectLastDir("mri") + "/" + layer->GetName(),
                                       "VTK files (*.vtk);;STL files (*.stl);;All files (*)", &selectedFilter);
  else
    selectedFilter = QFileInfo(fn).suffix();
  if ( !fn.isEmpty() )
  {
    QString selected_suffix = selectedFilter.left(3).toLower();
    if (selected_suffix == "all")
      selected_suffix = "vtk";
    QFileInfo fi(fn);
    if (fi.suffix().toLower() != selected_suffix)
      fn += "." + selected_suffix;
    if ( !layer->SaveContourToFile( fn ) )
    {
      QMessageBox::warning(this, "Error", "Can not save surface to file.");
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
  this->m_wndGroupPlot->raise();
  this->m_wndGroupPlot->SetCurrentVertex(0);
  m_strLastFsgdDir = QFileInfo(fn).absolutePath();
}

void MainWindow::ToggleSplinePicking()
{
  m_bSplinePicking = !m_bSplinePicking;
  cout << qPrintable(QString("Surface spline picking %1").arg(m_bSplinePicking?"enabled":"disabled")) << endl;
}

void MainWindow::SetSplinePicking(bool b)
{
  if (b != m_bSplinePicking)
    ToggleSplinePicking();
}

void MainWindow::OnReloadVolume()
{
  QList<Layer*> sel_layers = GetSelectedLayers("MRI");
  if (!sel_layers.isEmpty())
  {
    DialogReloadLayer dlg;
    if (dlg.Execute(sel_layers) == QDialog::Accepted)
    {
      bool bCloseFirst = dlg.GetCloseLayerFirst();
      int active_layer_id = GetActiveLayer("MRI")->GetID();
      for (int i = 0; i < sel_layers.size(); i++)
      {
        LayerMRI* mri = qobject_cast<LayerMRI*>(sel_layers[i]);
        QVariantMap map = mri->GetProperty()->GetFullSettings();
        if (mri->GetActiveFrame() > 0)
          map["frame"] = mri->GetActiveFrame();
        m_layerSettings[mri->GetID()+(bCloseFirst?0:LAYER_ID_OFFSET)] = map;
      }

      QList<Layer*> all_layers = GetLayers("MRI");
      QStringList layer_order;
      foreach (Layer* layer, all_layers)
        layer_order << QString::number(layer->GetID());

      //      if (dlg.GetCloseLayerFirst())
      //      {
      //        if (!OnCloseVolume())
      //        {
      //          m_volumeSettings.clear();
      //          return;
      //        }
      //      }
      QStringList layer_ids;
      if (bCloseFirst)
      {
        for (int i = 0; i < sel_layers.size(); i++)
        {
          layer_ids << QString::number(sel_layers[i]->GetID()+LAYER_ID_OFFSET);
        }
      }
      for (int i = sel_layers.size()-1; i >= 0; i--)
      {
        LayerMRI* mri = qobject_cast<LayerMRI*>(sel_layers[i]);
        QString name = mri->GetName();
        QString filename = mri->GetFileName();
        QString reg_fn = mri->GetRegFileName();
        QString args = filename + ":name=" + name;
        if (!reg_fn.isEmpty())
          args += ":reg=" + reg_fn;
        args += QString(":id=%1").arg(mri->GetID()+(bCloseFirst?0:LAYER_ID_OFFSET));
        if (bCloseFirst)
        {
          mri->SetID(mri->GetID()+LAYER_ID_OFFSET);
        }
        AddScript(QStringList("loadvolume") << args);
        if (bCloseFirst)
        {
          mri->MarkAboutToDelete();
          AddScript(QStringList("unloadlayers") << "mri" << QString::number(mri->GetID()));
        }
      }

      if (bCloseFirst)
      {
        AddScript(QStringList("reorderlayers") << "mri" << layer_order.join(","));
        for (int i = 0; i < layer_ids.size(); i++)
          layer_ids[i] = QString::number(layer_ids[i].toInt()-LAYER_ID_OFFSET);
        AddScript(QStringList("setactivelayer") << "mri" << QString::number(active_layer_id) << layer_ids.join(","));
      }
    }
  }
  else
    QMessageBox::warning(this, "Reload Volume", "Select at least one volume to reload.");
}

void MainWindow::OnReloadROI()
{
  QList<Layer*> sel_layers = GetSelectedLayers("ROI");
  if (!sel_layers.isEmpty())
  {
    DialogReloadLayer dlg;
    if (dlg.Execute(sel_layers) == QDialog::Accepted)
    {
      bool bCloseFirst = dlg.GetCloseLayerFirst();
      int active_layer_id = GetActiveLayer("ROI")->GetID();
      //      for (int i = 0; i < sel_layers.size(); i++)
      //      {
      //        LayerROI* roi = qobject_cast<LayerROI*>(sel_layers[i]);
      //        m_layerSettings[roi->GetID()] = roi->GetProperty()->GetFullSettings();
      //      }

      QList<Layer*> all_layers = GetLayers("ROI");
      QStringList layer_order;
      foreach (Layer* layer, all_layers)
        layer_order << QString::number(layer->GetID());

      QStringList layer_ids;
      if (bCloseFirst)
      {
        for (int i = 0; i < sel_layers.size(); i++)
        {
          layer_ids << QString::number(sel_layers[i]->GetID()+LAYER_ID_OFFSET);
        }
      }
      for (int i = sel_layers.size()-1; i >= 0; i--)
      {
        LayerROI* roi = qobject_cast<LayerROI*>(sel_layers[i]);
        QString name = roi->GetName();
        QString filename = roi->GetFileName();
        QString args = filename + ":name=" + name;
        double* rgb = roi->GetProperty()->GetColor();
        args += QString(":id=%1:color=%2,%3,%4:opacity=%5:threshold=%6:ref=%7").arg(roi->GetID()+(bCloseFirst?0:LAYER_ID_OFFSET))
            .arg((int)(rgb[0]*255)).arg((int)(rgb[1]*255)).arg((int)(rgb[2]*255)).arg(roi->GetProperty()->GetOpacity())
            .arg(roi->GetProperty()->GetThreshold()).arg(roi->GetRefMRI()->GetID());

        if (bCloseFirst)
          roi->SetID(roi->GetID()+LAYER_ID_OFFSET);
        AddScript(QStringList("loadroi") << args);
        if (bCloseFirst)
          AddScript(QStringList("unloadlayers") << "roi" << QString::number(roi->GetID()));
      }
      if (bCloseFirst)
      {
        AddScript(QStringList("reorderlayers") << "roi" << layer_order.join(","));
        for (int i = 0; i < layer_ids.size(); i++)
          layer_ids[i] = QString::number(layer_ids[i].toInt()-LAYER_ID_OFFSET);
        AddScript(QStringList("setactivelayer") << "roi" << QString::number(active_layer_id) << layer_ids.join(","));
      }
    }
  }
}

void MainWindow::OnReloadPointSet()
{
  QList<Layer*> sel_layers = GetSelectedLayers("PointSet");
  if (!sel_layers.isEmpty())
  {
    DialogReloadLayer dlg;
    if (dlg.Execute(sel_layers) == QDialog::Accepted)
    {
      //      for (int i = 0; i < sel_layers.size(); i++)
      //      {
      //        LayerROI* roi = qobject_cast<LayerROI*>(sel_layers[i]);
      //        m_layerSettings[roi->GetID()] = roi->GetProperty()->GetFullSettings();
      //      }

      bool bCloseFirst = dlg.GetCloseLayerFirst();
      int active_layer_id = GetActiveLayer("PointSet")->GetID();
      QList<Layer*> all_layers = GetLayers("PointSet");
      QStringList layer_order;
      foreach (Layer* layer, all_layers)
        layer_order << QString::number(layer->GetID());

      QStringList layer_ids;
      if (bCloseFirst)
      {
        for (int i = 0; i < sel_layers.size(); i++)
        {
          layer_ids << QString::number(sel_layers[i]->GetID()+LAYER_ID_OFFSET);
        }
      }
      for (int i = sel_layers.size()-1; i >= 0; i--)
      {
        LayerPointSet* ps = qobject_cast<LayerPointSet*>(sel_layers[i]);
        QString filename = ps->GetFileName();
        double* rgb = ps->GetProperty()->GetColor();
        QString args = filename + QString(":id=%1:color=%2,%3,%4:name=%5:radius=%6:visible=%7").arg(ps->GetID()+(bCloseFirst?0:LAYER_ID_OFFSET))
            .arg((int)(rgb[0]*255)).arg((int)(rgb[1]*255)).arg((int)(rgb[2]*255)).arg(ps->GetName())
            .arg(ps->GetProperty()->GetRadius()).arg(ps->IsVisible()?1:0);

        if (bCloseFirst)
          ps->SetID(ps->GetID()+LAYER_ID_OFFSET);
        if (ps->GetProperty()->GetType() == LayerPropertyPointSet::ControlPoint )
          AddScript(QStringList("loadcontrolpoints") << args);
        else
        {
          double* rgb = ps->GetProperty()->GetSplineColor();
          args += QString(":splinecolor=%1,%2,%3:splineradius=%4").arg((int)(rgb[0]*255)).arg((int)(rgb[1]*255)).arg((int)(rgb[2]*255))
              .arg(ps->GetProperty()->GetSplineRadius());
          AddScript(QStringList("loadwaypoints") << args);
          if (bCloseFirst)
            AddScript(QStringList("unloadlayers") << "pointset" << QString::number(ps->GetID()));
        }
      }
      if (bCloseFirst)
      {
        AddScript(QStringList("reorderlayers") << "pointset" << layer_order.join(","));
        for (int i = 0; i < layer_ids.size(); i++)
          layer_ids[i] = QString::number(layer_ids[i].toInt()-LAYER_ID_OFFSET);
        AddScript(QStringList("setactivelayer") << "pointset" << QString::number(active_layer_id) << layer_ids.join(","));
      }
    }
  }
}

void MainWindow::OnReloadSurface()
{
  QList<Layer*> sel_layers = GetSelectedLayers("Surface");
  if (!sel_layers.isEmpty())
  {
    DialogReloadLayer dlg;
    if (dlg.Execute(sel_layers) == QDialog::Accepted)
    {
      bool bCloseFirst = dlg.GetCloseLayerFirst();
      int active_layer_id = GetActiveLayer("Surface")->GetID();
      for (int i = sel_layers.size()-1; i >= 0; i--)
      {
        LayerSurface* surf = qobject_cast<LayerSurface*>(sel_layers[i]);
        m_layerSettings[surf->GetID()+(bCloseFirst?0:LAYER_ID_OFFSET)] = surf->GetProperty()->GetFullSettings();
      }

      QList<Layer*> all_layers = GetLayers("Surface");
      QStringList layer_order;
      foreach (Layer* layer, all_layers)
        layer_order << QString::number(layer->GetID());

      QStringList layer_ids;
      if (bCloseFirst)
      {
        for (int i = 0; i < sel_layers.size(); i++)
        {
          layer_ids << QString::number(sel_layers[i]->GetID()+LAYER_ID_OFFSET);
        }
      }
      for (int i = sel_layers.size()-1; i >= 0; i--)
      {
        LayerSurface* surf = qobject_cast<LayerSurface*>(sel_layers[i]);
        QString args = QString("%1:name=%2:id=%3").arg(surf->GetFileName()).arg(surf->GetName()).arg(surf->GetID()+(bCloseFirst?0:LAYER_ID_OFFSET));
        if (surf->GetCurrentVertex() >= 0)
          args += QString(":current_vertex=%1").arg(surf->GetCurrentVertex());
        args += QString(":overlay_zorder=%1:label_zorder=%2:annot_zorder=%3")
            .arg(surf->GetProperty()->GetZOrderOverlay())
            .arg(surf->GetProperty()->GetZOrderLabel())
            .arg(surf->GetProperty()->GetZOrderAnnotation());

        if (bCloseFirst)
          surf->SetID(surf->GetID()+LAYER_ID_OFFSET);
        AddScript(QStringList("loadsurface") << args);

        for (int j = surf->GetNumberOfOverlays()-1; j >= 0; j--)
        {
          SurfaceOverlay* overlay = surf->GetOverlay(j);
          QStringList script("loadsurfaceoverlay");
          if (overlay)
          {
            script << overlay->GetFileName() << overlay->GetRegFileName();
            if (overlay->HasCorrelationData())
              script << "correlation";
            AddScript(script);
          }
        }
        for (int j = surf->GetNumberOfLabels()-1; j >= 0; j--)
        {
          SurfaceLabel* label = surf->GetLabel(j);
          if (label)
          {
            AddScript(QStringList("loadsurfacelabel") << label->GetFileName());
            if (label->GetShowOutline())
              AddScript(QStringList("setsurfacelabeloutline") << "1");
            if (!label->IsVisible())
              AddScript(QStringList("hidesurfacelabel"));
            if (label->GetOpacity() != 1)
              AddScript(QStringList("setsurfacelabelopacity") << QString::number(label->GetOpacity()));
            double* c = label->GetColor();
            AddScript(QStringList("setsurfacelabelcolor") << QString("%1,%2,%3").arg((int)(c[0]*255)).arg((int)(c[1]*255)).arg((int)(c[2]*255)));
            AddScript(QStringList("setsurfacelabelthreshold") << QString::number(label->GetThreshold()));
          }
        }
        for (int j = surf->GetNumberOfAnnotations()-1; j >= 0; j--)
        {
          SurfaceAnnotation* annot = surf->GetAnnotation(j);
          if (annot)
          {
            AddScript(QStringList("loadsurfaceannotation") << annot->GetFilename());
            if (annot->GetShowOutline())
              AddScript(QStringList("setsurfaceannotationoutline") << "1");
          }
        }
        surf->MarkAboutToDelete();
        if (dlg.GetCloseLayerFirst())
          AddScript(QStringList("unloadlayers") << "surface" << QString::number(surf->GetID()));
      }

      if (bCloseFirst)
      {
        AddScript(QStringList("reorderlayers") << "surface" << layer_order.join(","));
        for (int i = 0; i < layer_ids.size(); i++)
          layer_ids[i] = QString::number(layer_ids[i].toInt()-LAYER_ID_OFFSET);
        AddScript(QStringList("setactivelayer") << "surface" << QString::number(active_layer_id) << layer_ids.join(","));
      }
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

void MainWindow::OnGoToROI(bool center)
{
  LayerROI* roi = (LayerROI*)GetActiveLayer("ROI");
  double pos[3];
  if (roi && roi->GetCentroidPosition(pos))
  {
    SetSlicePosition(pos);
    if (center)
      GetMainView()->CenterAtWorldPosition(pos);
  }
}

void MainWindow::OnGoToPointSet(bool center)
{
  LayerPointSet* ps = (LayerPointSet*)GetActiveLayer("PointSet");
  double pos[3];
  if (ps && ps->GetCentroidPosition(pos))
  {
    SetSlicePosition(pos);
    if (center)
      GetMainView()->CenterAtWorldPosition(pos);
  }
}

void MainWindow::OnGoToSurfaceLabel(bool center)
{
  LayerSurface* surf = (LayerSurface*)GetActiveLayer("Surface");
  double pos[3];
  if (surf && surf->GetActiveLabelCentroidPosition(pos))
  {
    bool mappedFromInflated = ((RenderView3D*)m_views[3])->MapInflatedCoords(surf, pos, pos, true);
    SetSlicePosition(pos);
    if (center)
    {
      GetMainView()->CenterAtWorldPosition(pos);
      if (!mappedFromInflated)
      {
        int nVertex = surf->GetVertexIndexAtTarget(pos, NULL);
        if (nVertex >= 0)
        {
          double v[3];
          surf->GetSmoothedVertexNormal(nVertex, v);
          m_views[3]->AlignViewToNormal(v);
        }
      }
    }
  }
}

void MainWindow::CommandLoadFCD(const QStringList& cmd )
{
  if (cmd.size() < 3)
    return;

  LoadFCD(cmd[1], cmd[2], cmd.size() > 3? cmd[3] : "");
}

void MainWindow::LoadFCD(const QString &subdir, const QString &subject, const QString& suffix)
{
  LayerFCD* layer = new LayerFCD(m_layerVolumeRef);
  connect( layer->GetWorkerThread(), SIGNAL(Progress(int)), m_statusBar, SLOT(SetProgress(int)));
  connect( layer->GetWorkerThread(), SIGNAL(started()), m_statusBar, SLOT(ShowProgress()));
  connect( layer->GetWorkerThread(), SIGNAL(finished()), m_statusBar, SLOT(HideProgress()));
  if (suffix.isEmpty())
    layer->SetName(subject);
  else
    layer->SetName(subject + "." + suffix);
  layer->SetMRILayerCTAB(m_luts->GetColorTable(0));
  QVariantMap map;
  map["SubjectDir"] = subdir;
  map["Subject"] = subject;
  map["Suffix"] = suffix;
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
  return m_settings.value(key);
}

void MainWindow::SetSetting(const QString &key, const QVariant &value)
{
  m_settings[key] = value;
}

void MainWindow::UpdateSettings()
{
  if (m_dlgPreferences)
  {
    QVariantMap old = m_settings;
    QVariantMap map = m_dlgPreferences->GetSettings();
    QStringList keys = map.keys();
    foreach (QString key, keys)
      m_settings[key] = map[key];

    //    if (old["AutoSetMidToMin"].toBool() != m_settings["AutoSetMidToMin"].toBool())
    //    {
    //      QList<Layer*> layers = GetLayers("MRI");
    //      foreach (Layer* l, layers)
    //      {
    //        LayerMRI* mri = (LayerMRI*)l;
    //        if (mri->GetProperty()->GetHeatScaleAutoMid())
    //        {
    //          double dMin = mri->GetProperty()->GetHeatScaleMinThreshold();
    //          double dMax = mri->GetProperty()->GetHeatScaleMaxThreshold();
    //          if (m_settings["AutoSetMidToMin"].toBool())
    //            mri->GetProperty()->SetHeatScaleMidThreshold(dMin);
    //          else
    //            mri->GetProperty()->SetHeatScaleMidThreshold((dMin+dMax)/2);
    //        }
    //      }
    //    }

    ui->actionDeleteLayer->setVisible(m_settings["AllowDeleteKey"].toBool());
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

void MainWindow::CommandReorderLayers(const QStringList &cmd)
{
  if (cmd.size() < 3)
    return;

  LayerCollection* lc = NULL;
  QString type = cmd[1].toLower();
  if (type == "mri")
    lc = GetLayerCollection("MRI");
  else if (type == "surface")
    lc = GetLayerCollection("Surface");
  else if (type == "roi")
    lc = GetLayerCollection("ROI");
  else if (type == "pointset")
    lc = GetLayerCollection("PointSet");

  if (lc)
  {
    QList<int> ids;
    QStringList list = cmd[2].split(",");
    for (int i = 0; i < list.size(); i++)
    {
      ids << list[i].toInt();
    }
    lc->UpdateLayerOrder(ids);
  }
}

void MainWindow::CommandSetActiveLayer(const QStringList &cmd)
{
  if (cmd.size() < 3)
    return;

  LayerCollection* lc = NULL;
  QString type = cmd[1].toLower();
  if (type == "mri")
    lc = GetLayerCollection("MRI");
  else if (type == "surface")
    lc = GetLayerCollection("Surface");
  else if (type == "roi")
    lc = GetLayerCollection("ROI");
  else if (type == "pointset")
    lc = GetLayerCollection("PointSet");
  else if (type == "tract")
    lc = GetLayerCollection("Tract");

  if (lc)
  {
    bool bOK = false;
    int nId = cmd[2].toInt(&bOK);
    Layer* layer = NULL;
    if (bOK)
      layer = lc->GetLayerById(nId);
    else
      layer = lc->GetLayerByName(cmd[2]);
    if (layer)
      lc->SetActiveLayer(layer);
    if (cmd.size() >= 4)
    {
      QList<int> ids;
      QStringList list = cmd[3].split(",");
      for (int i = 0; i < list.size(); i++)
      {
        ids << list[i].toInt();
      }
      ui->widgetAllLayers->SetSelectedLayers(ids);
    }
  }
}

void MainWindow::CommandUnloadLayers(const QStringList &cmd)
{
  if (cmd.size() < 3)
    return;

  LayerCollection* lc = NULL;
  QString type = cmd[1].toLower();
  if (type == "mri")
    lc = GetLayerCollection("MRI");
  else if (type == "surface")
    lc = GetLayerCollection("Surface");
  else if (type == "roi")
    lc = GetLayerCollection("ROI");
  else if (type == "pointset")
    lc = GetLayerCollection("PointSet");
  else if (type == "tract")
    lc = GetLayerCollection("Tract");

  if (lc)
  {
    QList<int> ids;
    QStringList list = cmd[2].split(",");
    for (int i = 0; i < list.size(); i++)
    {
      ids << list[i].toInt();
    }
    QList<Layer*> layers;
    for (int i = 0; i < lc->GetNumberOfLayers(); i++)
    {
      if (ids.contains(lc->GetLayer(i)->GetID()))
        layers << lc->GetLayer(i);
    }
    if (!layers.isEmpty())
    {
      if (type == "mri")
        OnCloseVolume(layers);
      else if (type == "surface")
        OnCloseSurface(layers);
      else if (type == "roi")
        OnCloseROI(layers);
      else if (type == "pointset")
        OnClosePointSet(layers);
      else if (type == "tract")
        OnCloseTrack();
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
    QVariantMap cam = ui->view3D->GetCameraInfo();
    QFile file(fn);
    if (file.open(QIODevice::WriteOnly))
    {
      file.write(QJsonDocument::fromVariant(cam).toJson());
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
      QVariantMap cam = QJsonDocument::fromJson(file.readAll()).toVariant().toMap();
      file.close();
      ui->view3D->SetCamera(cam);
    }
    else
    {
      qWarning() << "Can not open camera file " << fn;
    }
  }
}

void MainWindow::GoToContralateralPoint()
{
  LayerFCD* layer = qobject_cast<LayerFCD*>(GetActiveLayer("FCD"));
  if ( layer )
  {
    double pos[3];
    layer->GetSlicePosition(pos);
    if (layer->GoToContralateralPoint(pos, pos))
    {
      SetSlicePosition(pos);
      CenterAtWorldPosition(pos);
    }
  }
  else
  {
    LayerSurface* surf = (LayerSurface*)GetActiveLayer("Surface");
    if (surf)
    {
      if (surf->IsContralateralReady())
        GoToContralateralPoint(surf);
      else
      {
        QString fn = surf->GetFileName();
        if (surf->GetHemisphere() == 0)
          fn.replace("lh.", "rh.");
        else
          fn.replace("rh.", "lh.");
        AddScript(QStringList("loadsurface") << fn);
        AddScript(QStringList("gotocontralateralsurface") << QString::number(reinterpret_cast<quintptr>(surf)));
      }
    }
  }
}

void MainWindow::GoToContralateralPoint(LayerSurface *layer_in)
{
  LayerSurface* layer = layer_in;
  double pos[3];
  layer->GetSlicePosition(pos);
  int nvo = -1;
  bool bInflated = layer->GetFileName().contains("inflated");
  if (bInflated)
    nvo = layer->GetCurrentVertex();
  else
    nvo = layer->GetVertexIndexAtTarget(pos, NULL);
  if (nvo < 0 && layer->GetContralateralSurface())
  {
    layer = layer->GetContralateralSurface();
    if (bInflated)
      nvo = layer->GetCurrentVertex();
    else
      nvo = layer->GetVertexIndexAtTarget(pos, NULL);
  }
  nvo = layer->GetContralateralVertex(nvo);
  if (nvo >= 0)
  {
    layer = layer->GetContralateralSurface();
    GetLayerCollection("Surface")->SetActiveLayer(layer);
    if (layer)
    {
      layer->SetCurrentVertex(nvo);
      if (bInflated)
      {
        layer->GetTargetAtVertex(nvo, pos);
        ((RenderView3D*)m_views[3])->MapInflatedCoords(layer, pos, pos, m_settings["AutoReorientView"].toBool(), true);
      }
      else
        layer->GetTargetAtVertex(nvo, pos);
      SetSlicePosition(pos);
      CenterAtWorldPosition(pos);
    }
  }
  else
  {
    cout << "Did not find any vertex at cursor on " << qPrintable(layer->GetName()) << endl;
  }
}

Layer* GetLayerByFilename(const QString& fn, const QList<Layer*>& layers)
{
  foreach (Layer* layer, layers)
  {
    if (layer->GetFileName() == fn)
      return layer;
  }
  return NULL;
}

LayerSurface* GetContralateralSurfaceLayer(LayerSurface* surf, const QList<Layer*>& layers)
{
  QString fn = surf->GetFileName();
  if (surf->GetHemisphere() == 0)
    fn.replace("lh.", "rh.");
  else
    fn.replace("rh.", "lh.");
  return (LayerSurface*)GetLayerByFilename(fn, layers);
}

void MainWindow::UpdateSurfaceContralateralInfo()
{
  QList<Layer*> surfs = GetLayers("Surface");
  QList<Layer*> sphere_surfs = GetLayerCollection("HiddenSurface")->GetLayers("Surface");
  foreach(Layer* layer, surfs)
    ((LayerSurface*)layer)->ResetContralateralInfo();
  if (surfs.isEmpty())
  {
    GetLayerCollection("HiddenSurface")->Clear();
    return;
  }
  else if (sphere_surfs.isEmpty())
    return;

  for (int i = 0; i < surfs.size(); i++)
  {
    LayerSurface* surf = (LayerSurface*)surfs[i];
    LayerSurface* surf2 = GetContralateralSurfaceLayer(surf, surfs);
    if (surf2)
    {
      QString path = QFileInfo(surf->GetFileName()).absolutePath();
      LayerSurface* sphere1 = (LayerSurface*)GetLayerByFilename(QFileInfo(path + "/lh.sphere.d1.left_right").absoluteFilePath(), sphere_surfs);
      LayerSurface* sphere2 = (LayerSurface*)GetLayerByFilename(QFileInfo(path + "/rh.sphere.d1.left_right").absoluteFilePath(), sphere_surfs);
      if (sphere1 && sphere2)
      {
        surf->SetContralateralLayer(surf2, sphere1, sphere2);
        surf2->SetContralateralLayer(surf, sphere1, sphere2);
      }
      surfs.removeOne(surf2);
    }
    surfs.removeAt(i);
    i--;
  }
}

void MainWindow::LoadSphereLeftRightIfNeeded(LayerSurface *sf)
{
  QList<Layer*> layers = GetLayerCollection("HiddenSurface")->GetLayers("Surface");
  QString fullpath = QFileInfo(sf->GetFileName()).absolutePath();
  if (sf->GetHemisphere() == 0)
    fullpath += "/lh.sphere.d1.left_right";
  else
    fullpath += "/rh.sphere.d1.left_right";
  fullpath = QFileInfo(fullpath).absoluteFilePath();
  if (GetLayerByFilename(fullpath, layers))
    return;

  if (QFile::exists(fullpath))
  {
    LayerSurface* layer = new LayerSurface( m_layerVolumeRef );
    layer->SetFileName( fullpath );
    QVariantMap args;
    args["hidden"] = true;
    m_threadIOWorker->LoadSurface( layer, args );
    m_statusBar->StartTimer();
  }
}

void MainWindow::OnSurfaceVertexClicked(LayerSurface *surf)
{
  if (m_bVerbose)
  {
    int nVert = surf->GetCurrentVertex();
    if (nVert >= 0)
    {
      double ras[3], tkras[3];
      //            surf->GetRASAtVertex(nVert, ras);
      //            surf->GetSurfaceRASAtVertex(nVert, tkras);
      surf->GetSlicePosition(ras);
      surf->GetRASAtTarget(ras, ras);
      surf->GetSurfaceRASAtRAS(ras, tkras);
      printf("RAS: %.4f %.4f %.4f\n", ras[0], ras[1], ras[2]);
      printf("SurfaceRAS: %.4f %.4f %.4f\n", tkras[0], tkras[1], tkras[2]);
    }
  }
  if (surf->GetFileName().contains("inflated"))
  {
    double pos[3];
    surf->GetSlicePosition(pos);
    CenterAtWorldPosition(pos);
  }
}

void MainWindow::CenterAtWorldPosition(double *pos, bool mainview_only)
{
  if (mainview_only)
    GetMainView()->CenterAtWorldPosition(pos);
  else
  {
    for (int i = 0; i < 3; i++)
      m_views[i]->CenterAtWorldPosition(pos);
  }
}

void MainWindow::On2DCursorClicked()
{
  if (m_bVerbose)
  {
    LayerMRI* mri = qobject_cast<LayerMRI*>(GetActiveLayer("MRI"));
    LayerSurface* surf = qobject_cast<LayerSurface*>(GetActiveLayer("Surface"));
    double ras[3], tkras[3];
    if (mri)
    {
      mri->GetSlicePosition(ras);
      mri->TargetToRAS(ras, ras);
      mri->NativeRASToTkReg(ras, tkras);
      printf("RAS: %.4f %.4f %.4f\n", ras[0], ras[1], ras[2]);
      printf("tkReg: %.4f %.4f %.4f\n", tkras[0], tkras[1], tkras[2]);
    }
    else if (surf)
    {
      surf->GetSlicePosition(ras);
      surf->GetRASAtTarget(ras, ras);
      surf->GetSurfaceRASAtRAS(ras, tkras);
      printf("RAS: %.4f %.4f %.4f\n", ras[0], ras[1], ras[2]);
      printf("SurfaceRAS: %.4f %.4f %.4f\n", tkras[0], tkras[1], tkras[2]);
    }
  }

  double pos[3];
  GetLayerCollection("MRI")->GetSlicePosition(pos);
  ((RenderView3D*)m_views[3])->MapToInflatedCoords(pos);
}

bool MainWindow::LoadSurfaceRGBMap(const QString& fn)
{
  QString filename = fn;
  if (filename.isEmpty())
    filename = QFileDialog::getOpenFileName( this, "Select RGB file",
                                             AutoSelectLastDir( "surf" ),
                                             "All files (*)");
  if ( !filename.isEmpty() )
  {
    LayerSurface* layer = qobject_cast<LayerSurface*>(GetActiveLayer("Surface"));
    if ( layer )
    {
      if (!layer->LoadRGBFromFile(filename))
      {
        ShowNonModalMessage("Error", "Can not load rgb file.");
        return false;
      }
    }
    return true;
  }
  else
    return false;
}

void MainWindow::ReorderLayers(const QList<Layer *> &layers)
{
  if (!layers.isEmpty())
  {
    QString type = layers[0]->GetPrimaryType();
    GetLayerCollection(type)->ReorderLayers(layers);
  }
}

void MainWindow::OnApplyVolumeTransform()
{
  LayerMRI* mri = qobject_cast<LayerMRI*>(this->GetActiveLayer("MRI"));
  if (mri)
  {
    DialogLoadTransform dlg;
    QString filename = mri->GetFileName();
    if (dlg.exec() == QDialog::Accepted)
    {
      QVariantMap map = mri->GetProperty()->GetFullSettings();
      map["name"] = mri->GetName();
      map["index"] = GetLayerCollection("MRI")->GetLayerIndex(mri);
      int layer_id = mri->GetID();
      m_layerSettings[layer_id] = map;
      if (!OnCloseVolume())
      {
        m_layerSettings.remove(layer_id);
        return;
      }
      QVariantMap sup_data;
      sup_data["ID"] = layer_id;
      this->LoadVolumeFile(filename, dlg.GetFilename(), m_bResampleToRAS, dlg.GetSampleMethod(), false,
                           -1, "", sup_data);
    }
  }
}

void MainWindow::OnLoadSurfaceLabelRequested(const QString &fn)
{
//  AddScript(QStringList("loadsurfacelabel") << fn);
//  AddScript(QStringList("hidesurfacelabel"));
  if (LoadSurfaceLabelFile(fn))
    CommandHideSurfaceLabel();
}

void MainWindow::OnLoadTractCluster()
{
  QString dirPath = QFileDialog::getExistingDirectory(this, "Select Folder");
  if (!dirPath.isEmpty())
  {
    m_wndTractCluster->Load(dirPath);
  }
}

void MainWindow::CommandLoadTractCluster(const QStringList &cmd)
{
  m_wndTractCluster->Load(cmd[1]);
}

void MainWindow::OnTractClusterLoaded(const QVariantMap& data)
{
  QStringList filenames = data.value("filenames").toStringList();
  if (filenames.isEmpty())
  {
    QMessageBox::warning(this, "Error", "Could not find any tract files in selected folder");
  }
  else
  {
    ShowTractClusterMap();
    LayerTrack* layer = new LayerTrack( m_layerVolumeRef, NULL, true );
    layer->SetClusterData(data);
    m_threadIOWorker->LoadTrack( layer );
  }
}

void MainWindow::ShowTractClusterMap()
{
  m_wndTractCluster->show();
  m_wndTractCluster->raise();
}

void MainWindow::OnLoadVolumeTransform()
{
  LayerMRI* mri = ( LayerMRI* )GetActiveLayer( "MRI");
  if (!mri)
    return;

  DialogLoadTransform dlg(this);
  if ( dlg.exec() == QDialog::Accepted )
  {
    mri->SetRegFileName(dlg.GetFilename());
    mri->SetSampleMethod(dlg.GetSampleMethod());
    m_threadIOWorker->TransformVolume(mri);
  }
}

void MainWindow::OnUnloadVolumeTransform()
{
  LayerMRI* mri = ( LayerMRI* )GetActiveLayer( "MRI");
  if (!mri)
    return;

  QVariantMap args;
  args["unload"] = true;
  m_threadIOWorker->TransformVolume(mri, args);
}

void MainWindow::OnLoadSurfaceParameterization()
{
  QString filename = QFileDialog::getOpenFileName( this, "Select Parameterization File",
                                                   AutoSelectLastDir( "surf" ),
                                                   "Parameterization files (*)", 0, QFileDialog::DontConfirmOverwrite);
  if ( !filename.isEmpty() )
    LoadSurfaceParameterization(filename);
}

void MainWindow::LoadSurfaceParameterization(const QString &filename)
{
  LayerSurface* surf = ( LayerSurface* )GetActiveLayer("Surface");
  if (surf && !surf->LoadParameterization(filename) )
  {
    QMessageBox::warning(this, "Error", QString("Could not load parameterization from %1").arg(filename));
  }
}

void MainWindow::OnStereoRender(bool bOn)
{
  ui->view3D->SetStereoTypeToAnaglyph();
  ui->view3D->SetStereoRender(bOn);
}

Layer* MainWindow::FindSupplementLayer(const QString &name)
{
  LayerCollection* lc = GetLayerCollection("Supplement");
  QList<Layer*> layers = lc->GetLayers();
  foreach (Layer* layer, layers)
  {
    if (layer->GetName() == name)
      return layer;
  }
  return NULL;
}

void MainWindow::SetCurrentTimeCourseFrame(int nFrame)
{
  m_wndTimeCourse->SetCurrentFrame(nFrame);
}

void MainWindow::OnViewLayerInfo()
{
  QString type = GetCurrentLayerType();
  Layer* layer = GetActiveLayer(type);
  if (layer)
  {
    m_wndLayerInfo->UpdateInfo(layer);
    m_wndLayerInfo->show();
    m_wndLayerInfo->raise();
  }
}

void MainWindow::UpdateLayerInfo(Layer* layer)
{
  if (layer && m_wndLayerInfo->isVisible())
    m_wndLayerInfo->UpdateInfo(layer);
}

void MainWindow::LoadSurfaceCoordsFromParameterization( const QString& filename )
{
  LayerSurface* layer = ( LayerSurface* )GetActiveLayer("Surface");
  if ( layer && !layer->LoadCoordsFromParameterization(filename))
  {
    QMessageBox::warning(this, "Error", QString("Could not load parameterization from %1").arg(filename));
  }
}

void MainWindow::OnExportLabelStats()
{
  QString fn = QFileDialog::getSaveFileName( this, "Save Label Stats",
                                             AutoSelectLastDir( "mri" ),
                                             "CSV files (*.csv)");
  if (!fn.isEmpty() )
  {
    LayerMRI* mri = (LayerMRI*)GetActiveLayer("MRI");
    if (!mri->ExportLabelStats(fn))
      QMessageBox::warning(this, "Error", QString("Could not save label stats to %1").arg(fn));
  }
}

void MainWindow::CommandExportLineProfileThickness(const QStringList &cmd)
{
  QVariantMap opts;
  QStringList ar = cmd[1].split(":");
  for (int i = 1; i < ar.size(); i++)
  {
    QStringList list = ar[i].split("=");
    if (list.size() > 1)
    {
      if (list[0].toLower() == "spacing")
        opts["spacing"] = list[1].toDouble();
      else if (list[0].toLower() == "resolution")
        opts["resolution"] = list[1].toDouble();
      else if (list[0].toLower() == "offset")
        opts["offset"] = list[1].toDouble();
      else if (list[0].toLower() == "segments")
        opts["segments"] = list[1].toDouble();
    }
  }
  QString fn = ar[0];
  if (!ExportLineProfileThickness(fn, opts))
    cerr << "Failed to export line profile thickness to " << qPrintable(fn) << endl;
}

bool MainWindow::ExportLineProfileThickness(const QString &filename, const QVariantMap& opts)
{
  QList<Layer*> layers = GetLayers("PointSet");
  QList<LayerPointSet*> lines;
  foreach (Layer* layer, layers)
  {
    lines.insert(0, qobject_cast<LayerPointSet*>(layer));
  }

  if (lines.size() > 1)
  {
#if !defined(ARM64)
    LayerLineProfile* lp = new LayerLineProfile(GetMainViewId(), this, lines.first(), lines.last());
    lines.removeFirst();
    lines.removeLast();
#endif

    double dVoxelSize = 1.0;
    LayerMRI* mri = qobject_cast<LayerMRI*>(GetActiveLayer("MRI"));
    if (mri)
    {
      double vs[3];
      mri->GetWorldVoxelSize(vs);
      dVoxelSize = qMin(vs[0], qMin(vs[1], vs[2]));
    }

    double spacing = 1, resolution = 1, offset = 5;
    int samples = 100;
    if (opts.contains("spacing"))
      spacing = opts["spacing"].toDouble();
    if (opts.contains("resolution"))
      resolution = opts["resolution"].toDouble();
    if (opts.contains("offset"))
      offset = opts["offset"].toDouble();
    if (opts.contains("segments"))
      samples = opts["segments"].toInt();

#if !defined(ARM64)
    if (!lp->Solve(spacing, dVoxelSize, resolution, offset))
    {
      cerr << "Could not solve line profile\n";
      lp->deleteLater();
      return false;
    }
    else if (!lp->ExportThickness(filename, lines, samples))
    {
      cerr << "Could not export line profile thickness\n";
      lp->deleteLater();
      return false;
    }
    else
      lp->deleteLater();
#endif
  }
  return true;
}

void MainWindow::WriteLog(const QString &str_in, const QString &filename, bool bOverwrite)
{
  QFileInfo fi(QStandardPaths::locate(QStandardPaths::DocumentsLocation, "", QStandardPaths::LocateDirectory), filename);
  QFile file(fi.absoluteFilePath());
  file.open(bOverwrite?QFile::WriteOnly:QFile::Append);
  QString str = QString("[%1] %2\r\n").arg(QDateTime::currentDateTime().toString()).arg(str_in);
  file.write(str.toUtf8());
  file.flush();
  file.close();
}

void MainWindow::OnShowControlPanel(bool bShow)
{
  if (ui->widgetControlPanel->parentWidget() == ui->widgetControlPanelHolder)
    ui->widgetControlPanelHolder->setVisible(bShow);
  else
    m_widgetFloatControlPanel->setVisible(bShow);
}

void MainWindow::OnFloatPanels(bool bFloat)
{
  static QByteArray geometryControlPanel = ui->widgetControlPanelHolder->saveGeometry();
  static QByteArray geometryInfoPanel = ui->widgetInfoPanelHolder->saveGeometry();
  if (bFloat)
  {
    ui->widgetControlPanelHolder->hide();
    ui->layoutControlPanelHolder->removeWidget(ui->widgetControlPanel);
    m_widgetFloatControlPanel->show();
    m_widgetFloatControlPanel->layout()->addWidget(ui->widgetControlPanel);
    ui->widgetControlPanel->show();
    m_widgetFloatControlPanel->restoreGeometry(geometryControlPanel);

    ui->widgetInfoPanelHolder->hide();
    ui->layoutInfoPanelHolder->removeWidget(ui->widgetInfoPanel);
    m_widgetFloatInfoPanel->show();
    m_widgetFloatInfoPanel->layout()->addWidget(ui->widgetInfoPanel);
    ui->widgetInfoPanel->show();
    m_widgetFloatInfoPanel->restoreGeometry(geometryInfoPanel);
  }
  else
  {
    geometryControlPanel = m_widgetFloatControlPanel->saveGeometry();
    ui->widgetControlPanelHolder->show();
    m_widgetFloatControlPanel->layout()->removeWidget(ui->widgetControlPanel);
    ui->widgetControlPanel->show();
    ui->layoutControlPanelHolder->addWidget(ui->widgetControlPanel);
    m_widgetFloatControlPanel->hide();

    geometryInfoPanel = m_widgetFloatInfoPanel->saveGeometry();
    ui->widgetInfoPanelHolder->show();
    m_widgetFloatInfoPanel->layout()->removeWidget(ui->widgetInfoPanel);
    ui->widgetInfoPanel->show();
    ui->layoutInfoPanelHolder->addWidget(ui->widgetInfoPanel);
    m_widgetFloatInfoPanel->hide();
  }
}

void MainWindow::CommandLinkVolume(const QStringList &cmd)
{
  if ( cmd.size() > 1 )
  {
    LayerMRI* mri = qobject_cast<LayerMRI*>(GetActiveLayer("MRI"));
    if ( mri )
    {
      if ( cmd[1] == "1" || cmd[1].toLower() == "true" )
      {
        QList<LayerMRI*> linked_vols = ui->widgetAllLayers->GetLinkedVolumes();
        while (!linked_vols.isEmpty() && linked_vols[0]->GetProperty()->GetColorMap() == LayerPropertyMRI::LUT)
          linked_vols.removeFirst();
        if (!linked_vols.isEmpty() && linked_vols[0] != mri)
        {
          mri->GetProperty()->CopyWindowLevelSettings(linked_vols[0]->GetProperty());
        }
        emit LinkVolumeRequested(mri);
      }
    }
  }
}

void MainWindow::OnSyncInstances(bool bChecked)
{
  if (!QFileInfo(m_sSyncFilePath).isWritable())
    return;

  if (bChecked)
  {
    m_syncFileWatcher->addPath(m_sSyncFilePath);
    connect(m_syncFileWatcher, SIGNAL(fileChanged(QString)), SLOT(OnSyncFileChanged(QString)), Qt::UniqueConnection);
    connect(this, SIGNAL(SlicePositionChanged()), SLOT(UpdateSyncCoord()), Qt::ConnectionType(Qt::QueuedConnection | Qt::UniqueConnection));
    UpdateSyncIds(true);
  }
  else
  {
    m_syncFileWatcher->removePath(m_sSyncFilePath);
    disconnect(m_syncFileWatcher, SIGNAL(fileChanged(QString)), this, SLOT(OnSyncFileChanged(QString)));
    disconnect(this, SIGNAL(SlicePositionChanged()), this, SLOT(UpdateSyncCoord()));
    UpdateSyncIds(false);
  }
}

void MainWindow::UpdateSyncCoord()
{
  QFile file(m_sSyncFilePath);
  QVariantMap map;
  if (file.open(QIODevice::ReadOnly))
  {
    map = QJsonDocument::fromJson(file.readAll()).toVariant().toMap();
    file.close();
  }

  QVariantMap ras;
  double pos[3];
  GetLayerCollection("MRI")->GetSlicePosition(pos);
  ras["x"] = pos[0];
  ras["y"] = pos[1];
  ras["z"] = pos[2];
  map["ras"] = ras;
  map["instance_id"] = qApp->applicationPid();
  if (file.open(QIODevice::WriteOnly))
  {
    file.write(QJsonDocument::fromVariant(map).toJson());
    file.flush();
    file.close();
  }
  else
  {
    qWarning() << "Can not write to sync file " << m_sSyncFilePath;
  }
}

void MainWindow::UpdateSyncIds(bool bAdd)
{
  QFile file(m_sSyncFilePath);
  QVariantMap map;
  if (file.open(QIODevice::ReadOnly))
  {
    map = QJsonDocument::fromJson(file.readAll()).toVariant().toMap();
    file.close();
  }

  QStringList list = map.value("instance_list").toStringList();
  if (list.size() > 4)
    list = list.mid(0, 4);
  QString strg = QString::number(qApp->applicationPid());
  if (bAdd && !list.contains(strg))
    list.insert(list.begin(), strg);
  else if (!bAdd && list.contains(strg))
    list.removeAll(strg);
  map["instance_list"] = list;
  if (file.open(QIODevice::WriteOnly))
  {
    file.write(QJsonDocument::fromVariant(map).toJson());
    file.flush();
    file.close();
  }
}

void MainWindow::OnTileSyncedWindows()
{
  QFile file(m_sSyncFilePath);
  QVariantMap map;
  if (file.open(QIODevice::ReadOnly))
  {
    map = QJsonDocument::fromJson(file.readAll()).toVariant().toMap();
    file.close();
  }

  QStringList list = map.value("instance_list").toStringList();
  QString id_str = QString::number(qApp->applicationPid());
  TileWindow(list.indexOf(id_str));
  list.removeAll(id_str);
  map["to_be_tiled"] = list;
  if (file.open(QIODevice::WriteOnly))
  {
    file.write(QJsonDocument::fromVariant(map).toJson());
    file.flush();
    file.close();
  }
}

void MainWindow::TileWindow(int n)
{
  QRect rc = QGuiApplication::primaryScreen()->geometry();
  if (n == 0)
    rc.setWidth(rc.width()/2);
  else
    rc.setLeft(rc.left()+rc.width()/2);
  setGeometry(rc);
}

void MainWindow::OnSyncFileChanged(const QString &fn)
{
  QFile file(fn);
  if (file.open(QIODevice::ReadOnly))
  {
    QVariantMap map = QJsonDocument::fromJson(file.readAll()).toVariant().toMap();
    file.close();
    if (map["instance_id"].toLongLong() != qApp->applicationPid() && map.contains("ras"))
    {
      double pos[3];
      pos[0] = map["ras"].toMap().value("x").toDouble();
      pos[1] = map["ras"].toMap().value("y").toDouble();
      pos[2] = map["ras"].toMap().value("z").toDouble();
      disconnect(this, SIGNAL(SlicePositionChanged()), this, SLOT(UpdateSyncCoord()));
      SetSlicePosition(pos);
      connect(this, SIGNAL(SlicePositionChanged()), SLOT(UpdateSyncCoord()), Qt::ConnectionType(Qt::QueuedConnection | Qt::UniqueConnection));
    }
    QString id_str = QString::number(qApp->applicationPid());
    QStringList list = map.value("to_be_tiled").toStringList();
    if (list.contains(id_str))
    {
      list.removeAll(id_str);
      map["to_be_tiled"] = list;
      list = map.value("instance_list").toStringList();
      TileWindow(list.indexOf(id_str));
      if (file.open(QIODevice::WriteOnly))
      {
        file.write(QJsonDocument::fromVariant(map).toJson());
        file.flush();
        file.close();
      }
    }
  }
  else
  {
    qWarning() << "Can not open sync file " << fn;
  }
}

void MainWindow::OnLoadODF()
{
  QString fn = QFileDialog::getOpenFileName(this, "Load ODF", m_strLastDir);
  if (!fn.isEmpty())
  {
    AddScript(QStringList("loadodf") << fn);
  }
}

void MainWindow::OnCloseODF()
{
  LayerODF* layer = (LayerODF*)GetActiveLayer( "ODF" );
  if ( !layer )
    return;

  GetLayerCollection( "ODF" )->RemoveLayer( layer );
}

void MainWindow::CommandLoadODF(const QStringList& cmd )
{
  if (cmd.size() < 2)
    return;

  LayerODF* layer = new LayerODF(m_layerVolumeRef);
  QVariantMap map;
  QStringList list = cmd[1].split(":");
  QString fn = list.first();
  if (list.size() > 1)
  {
    list = list[1].split("=");
    if (list.size() > 1 && (list.first() == "permute" || list.first() == "permuted"))
    {
      map["permuted"] = (list[1] == "1" || list[1] == "true");
    }
  }
  map["Filename"] = QFileInfo(fn).absoluteFilePath();
  if (cmd.size() > 2)
    map["vertex_filename"] = QFileInfo(cmd[2]).absoluteFilePath();
  if (cmd.size() > 3)
    map["face_filename"] = QFileInfo(cmd[3]).absoluteFilePath();
  m_threadIOWorker->LoadODF( layer, map );
}

void MainWindow::OnSaveLabelAsVolume()
{
  LayerMRI* layer_mri = ( LayerMRI* )GetActiveLayer("MRI");
  if ( !layer_mri || !sender())
  {
    return;
  }

  int nVal = sender()->property("label_value").toInt();
  QString fn = QFileDialog::getSaveFileName( this, "Save volume",
                                             QFileInfo( layer_mri->GetFileName() ).absolutePath(),
                                             "Volume files (*.mgz *.mgh *.nii *.nii.gz *.img *.mnc);;All files (*)");
  if ( !fn.isEmpty() )
  {
    layer_mri->setProperty("label_value", nVal);
    layer_mri->setProperty("label_fn", fn);
    m_scripts.append(QStringList("savelayer") << QString::number(layer_mri->GetID()));
    ui->widgetAllLayers->UpdateWidgets();
  }
}

void MainWindow::OnPointSetToLabel()
{
  LayerPointSet* ps = qobject_cast<LayerPointSet*>(GetActiveLayer("PointSet"));
  LayerMRI* mri = qobject_cast<LayerMRI*>(GetActiveLayer("MRI"));
  RenderView2D* view = qobject_cast<RenderView2D*>(GetMainView());
  if (ps && mri && view)
  {
    mri->SaveForUndo(view->GetViewPlane());
    mri->UpdateVoxelsByPointSet(ps, view->GetViewPlane());
  }
}

void MainWindow::OnCreateOptimalVolume()
{
  QList<Layer*> mris = GetSelectedLayers("MRI");
  if (mris.size() < 2)
  {
    mris = GetVisibleLayers("MRI");
    if (mris.size() < 2)
      mris = GetLayers("MRI");
  }

  QList<Layer*> rois = GetSelectedLayers("ROI");
  if (rois.size() < 2)
  {
    rois = GetVisibleLayers("ROI");
    if (rois.size() < 2)
      rois = GetLayers("ROI");
  }
  if (rois.size() < 2)
  {
    QMessageBox::warning(this, "Optimal Combined Volume", "Need two ROIs to compute optimal combined volume");
    return;
  }

  LayerMRI* mri_template = (LayerMRI*)GetActiveLayer( "MRI" );
  LayerMRI* mri_new = new LayerMRI( mri_template );
  if ( !mri_new->Create( mri_template, false, 3))
  {
    QMessageBox::warning( this, "Error", "Can not create new volume." );
    delete mri_new;
    return;
  }
  mri_new->GetProperty()->SetLUTCTAB( m_luts->GetColorTable( 0 ) );
  mri_new->SetName( "optimal combined" );
  GetLayerCollection("MRI")->AddLayer(mri_new);
  ConnectMRILayer(mri_new);
  emit NewVolumeCreated();

  QList<LayerMRI*> input_mris;
  QList<LayerROI*> input_rois;
  for (int i = 0; i < mris.size(); i++)
    input_mris << (LayerMRI*)mris[i];
  for (int i = 0; i < 2; i++)
    input_rois << (LayerROI*)rois[i];

  VolumeFilterOptimal* filter = new VolumeFilterOptimal(input_mris, input_rois, mri_new, this);
  filter->SetResetWindowLevel();
  m_threadVolumeFilter->ExecuteFilter(filter);
}

void MainWindow::OnDeleteLayer()
{
  QString type = GetCurrentLayerType();
  Layer* layer = GetActiveLayer(type);
  if (!layer)
    return;

  if (type == "MRI")
    OnCloseVolume();
  else if (type == "Surface")
    OnCloseSurface();
  else if (type == "ROI")
    OnCloseROI();
  else if (type == "PointSet")
    OnClosePointSet();
}

void MainWindow::SetNeurologicalView(bool b)
{
  for (int i = 0; i < 3; i++)
  {
    ((RenderView2D*)m_views[i])->SetNeurologicalView(b);
  }
  QList<Layer*> layers = GetLayers("Surface");
  foreach (Layer* layer, layers)
  {
    ((LayerSurface*)layer)->SetDisplayInNeurologicalView(b);
  }
  layers = GetLayers("FCD");
  foreach (Layer* layer, layers)
  {
    ((LayerFCD*)layer)->SetDisplayInNeurologicalView(b);
  }
  layers = GetLayers("MRI");
  foreach (Layer* layer, layers)
  {
    ((LayerMRI*)layer)->SetDisplayInNeurologicalView(b);
  }
}

void MainWindow::OnShowLabelOutline(bool bShow)
{
  LayerMRI* mri = (LayerMRI*)GetLayerCollection( "MRI" )->GetActiveLayer();
  if ( mri )
  {
    mri->GetProperty()->SetShowLabelOutline(bShow);
  }
}

void MainWindow::OnFlattendSurfacePatchLoaded()
{
  if (m_nViewLayout == VL_1x1)
    SetMainView(MV_3D);

  ((RenderView3D*)m_views[3])->ResetViewSuperior();
}

void MainWindow::dragEnterEvent(QDragEnterEvent *event)
{
  QList<QUrl> urls = event->mimeData()->urls();
  for (int i = 0; i < urls.size(); i++)
  {
    if (!IsAcceptableUrl(urls[i]))
    {
      urls.removeAt(i);
      i--;
    }
  }

  if (!urls.isEmpty())
    event->acceptProposedAction();
}

bool MainWindow::IsAcceptableUrl(const QUrl &url, int url_type) // -1: all, 0: volume, 1: surface
{
  if (url.scheme().toLower() == "file")
  {
    QFileInfo fi(url.toLocalFile());
    QString suffix = fi.completeSuffix().toLower();
    if (suffix.contains("nii") || suffix.contains("mgz") || suffix.contains("mgh") ||
        suffix.contains("gii") || suffix == "dcm" || suffix == "hdr" ||
        suffix == "img" || suffix == "tif" || suffix == "tiff" || suffix == "bmp")
    {
      if (url_type <= 0)
        return true;
      else
        return false;
    }
    if (suffix == "orig" || suffix == "inflated" || suffix == "white" || suffix == "pial" ||
        suffix == "3gp" || suffix == "3g2")
    {
      if (url_type != 0)
        return true;
      else
        return false;
    }
  }
  return false;
}

void MainWindow::dropEvent(QDropEvent *event)
{
  const QMimeData* mimeData = event->mimeData();

  if (mimeData->hasUrls())
  {
    QList<QUrl> urlList = mimeData->urls();
    if (!urlList.isEmpty())
    {
      event->acceptProposedAction();
      foreach (QUrl url, urlList)
      {
        if (IsAcceptableUrl(url, 0)) // volume
        {
          AddScript(QStringList("loadvolume") << url.toLocalFile());
        }
        else if (IsAcceptableUrl(url, 1)) // surface
        {
          AddScript(QStringList("loadsurface") << url.toLocalFile());
        }
      }
    }
  }
}

bool MainWindow::GetNeurologicalView()
{
  return ((RenderView2D*)m_views[0])->GetNeurologicalView();
}
